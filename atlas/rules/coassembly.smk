


if config.get("assembler", "megahit") == "coassembly":

    assembly_params['megahit']={'default':'','meta-sensitive':'--presets meta-sensitive','meta-large':' --presets meta-large'}
    ASSEMBLY_FRACTIONS= MULTIFILE_FRACTIONS
    if PAIRED_END and config.get("merge_pairs_before_assembly", True):

        if 'se' in MULTIFILE_FRACTIONS:

            localrules: merge_se_me_for_megahit
            rule merge_se_me_for_megahit:
                input:
                    expand("{{sample}}/assembly/reads/{assembly_preprocessing_steps}_{fraction}.fastq.gz",
                    fraction=['se','me'], assembly_preprocessing_steps=assembly_preprocessing_steps)
                output:
                    temp(expand("{{sample}}/assembly/reads/{assembly_preprocessing_steps}_{fraction}.fastq.gz",
                                fraction=['co'], assembly_preprocessing_steps=assembly_preprocessing_steps))
                shell:
                    "cat {input} > {output}"

            ASSEMBLY_FRACTIONS = ['R1','R2','co']
        else:
            ASSEMBLY_FRACTIONS = ['R1','R2','me']





    def megahit_coassembly_input(wildcards):
        "returns a list of lists"

        return dict(zip(ASSEMBLY_FRACTIONS,
                [expand("{sample}/assembly/reads/{assembly_preprocessing_steps}_{fraction}.fastq.gz",
                fraction=fraction,
                assembly_preprocessing_steps=assembly_preprocessing_steps,
                sample=SAMPLES) for fraction in ASSEMBLY_FRACTIONS ]
                        ))

    def megahit_coassembly_input_parsing(input):

        out=''
        if hasattr(input,'se'):
            out+= f" --read {','.join(input.se)} "
        elif hasattr(input,'co'):
            out+= f" --read {','.join(input.co)} "
        elif hasattr(input,'me'):
            out+= f" --read {','.join(input.me)} "

        if hasattr(input,'R1'):

            out+= f" -1 {','.join(input.R1)} -2 {','.join(input.R2)} "
        return out
#    ruleorder: coassembly_megahit > run_megahit
    rule coassembly_megahit:
          input:
              unpack(megahit_coassembly_input)
          output:
              temp("co-assembly/assembly/megahit/co-assembly_prefilter.contigs.fa")
          benchmark:
              "logs/benchmarks/assembly/megahit/co-assembly.txt"
    #        shadow:
    #            "shallow" #needs to be shallow to find input files
          log:
              "co-assembly/logs/assembly/megahit.log"
          params:
              min_count = config.get("megahit_min_count", MEGAHIT_MIN_COUNT),
              k_min = config.get("megahit_k_min", MEGAHIT_K_MIN),
              k_max = config.get("megahit_k_max", MEGAHIT_K_MAX),
              k_step = config.get("megahit_k_step", MEGAHIT_K_STEP),
              merge_level = config.get("megahit_merge_level", MEGAHIT_MERGE_LEVEL),
              prune_level = config.get("megahit_prune_level", MEGAHIT_PRUNE_LEVEL),
              low_local_ratio = config.get("megahit_low_local_ratio", MEGAHIT_LOW_LOCAL_RATIO),
              min_contig_len = config.get("prefilter_minimum_contig_length", PREFILTER_MINIMUM_CONTIG_LENGTH),
              outdir = lambda wc, output: os.path.dirname(output[0]),
              inputs = lambda wc, input: megahit_coassembly_input_parsing(input),
              preset = assembly_params['megahit'][config['megahit_preset']],
              prefix = "co-assembly"
          conda:
              "%s/assembly.yaml" % CONDAENV
          threads:
              config.get("assembly_threads", ASSEMBLY_THREADS)
          resources:
              mem = config.get("assembly_memory", ASSEMBLY_MEMORY) #in GB
          shell:
              """
                  #rm -r {params.outdir} 2> {log}

                  megahit \
                  {params.inputs} \
                  --tmp-dir {TMPDIR} \
                  --num-cpu-threads {threads} \
                  --k-min {params.k_min} \
                  --k-max {params.k_max} \
                  --k-step {params.k_step} \
                  --out-dir {params.outdir} \
                  --out-prefix {params.prefix}_prefilter \
                  --min-contig-len {params.min_contig_len} \
                  --min-count {params.min_count} \
                  --merge-level {params.merge_level} \
                  --prune-level {params.prune_level} \
                  --low-local-ratio {params.low_local_ratio} \
                  --memory {resources.mem}000000000  \
                  --continue \
                  {params.preset} >> {log} 2>&1
              """

    localrules: rename_megahit_output
    rule rename_megahit_output:
        input:
            "{sample}/assembly/megahit/{sample}_prefilter.contigs.fa"
        output:
            temp("{sample}/assembly/{sample}_raw_contigs.fasta")
        shell:
            "cp {input} {output}"

    rule coassembly:
        input:
            "co-assembly/co-assembly_contigs.fasta",
            #expand("co-assembly/sequence_alignment/{sample}.bam"),
            #"co-assembly/assembly/contig_stats/postfilter_coverage_stats.txt",
            #"co-assembly/assembly/contig_stats/prefilter_contig_stats.txt",
            "co-assembly/assembly/contig_stats/final_contig_stats.txt"
        output:
            touch("co-assembly/finished_assembly")


    ruleorder: get_coassembly_metabat_depth_file > get_metabat_depth_file
    rule get_coassembly_metabat_depth_file:
        input:
            bam = lambda wc: expand("{sample}/sequence_alignment/{sample_reads}.bam",
                         sample_reads = SAMPLES,
                         sample='co-assembly')
        output:
            temp("co-assembly/binning/metabat/metabat_depth.txt")
        log:
            "co-assembly/binning/metabat/metabat.log"
        conda:
            "%s/metabat.yaml" % CONDAENV
        threads:
            config['threads']
        resources:
            mem = config["java_mem"]
        shell:
            """
            jgi_summarize_bam_contig_depths --outputDepth {output} {input.bam} \
                &> >(tee {log})
            """