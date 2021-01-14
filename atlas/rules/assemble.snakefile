import os
import re
import sys
from glob import glob
from snakemake.utils import report
import warnings
from copy import deepcopy




def get_preprocessing_steps(config):
    preprocessing_steps = ['QC']
    if config.get("normalize_reads_before_assembly", False):
        #preprocessing_steps.append("normalized")
        raise NotImplementedError("Normalization is depricated. It showed to reduce assembly performance at least for metagenomics."
                                  "remove the line 'normalize_reads_before_assembly' in your config file."
                                  " If you would like to have added the feature again write an issue.")

    if config.get("error_correction_before_assembly", True):
        preprocessing_steps.append("errorcorr")

    if config.get("merge_pairs_before_assembly", True) and PAIRED_END:
        preprocessing_steps.append("merged")

    return ".".join(preprocessing_steps)


assembly_preprocessing_steps = get_preprocessing_steps(config)


if  SKIP_QC & (len(MULTIFILE_FRACTIONS)<3):
    # I have problems with se reads

    rule init_pre_assembly_processing:
        input:
            get_quality_controlled_reads #expect SE or R1,R2 or R1,R2,SE
        output:
             temp(expand("{{sample}}/assembly/reads/QC_{fraction}.fastq.gz",fraction= MULTIFILE_FRACTIONS))
        params:
            inputs = lambda wc, input: io_params_for_tadpole(input,'in'),
            interleaved = lambda wc: "t" if (config.get('interleaved_fastqs',False) & SKIP_QC) else "f",
            outputs = lambda wc, output: io_params_for_tadpole(output,'out'),
            verifypaired = "t" if PAIRED_END else "f"
        log:
            "{sample}/logs/assembly/init.log"
        conda:
            "%s/required_packages.yaml" % CONDAENV
        threads:
            config.get("simplejob_threads", 1)
        resources:
            mem = config["simplejob_mem"],
            java_mem = int(int(config["simplejob_mem"] * JAVA_MEM_FRACTION))
        shell:
            """
            reformat.sh {params.inputs} \
                interleaved={params.interleaved} \
                {params.outputs} \
                iupacToN=t \
                touppercase=t \
                qout=33 \
                overwrite=true \
                verifypaired={params.verifypaired} \
                addslash=t \
                trimreaddescription=t \
                threads={threads} \
                -Xmx{resources.java_mem}G 2> {log}
            """

else:
    localrules: init_pre_assembly_processing
    rule init_pre_assembly_processing:
        input:
            get_quality_controlled_reads #expect SE or R1,R2 or R1,R2,SE
        output:
             temp(expand("{{sample}}/assembly/reads/QC_{fraction}.fastq.gz",fraction= MULTIFILE_FRACTIONS))
        log:
            "{sample}/logs/assembly/init.log"
        threads:
            1
        run:
            #make symlink
            assert len(input) == len(output), "Input and ouput files have not same number, can not create symlinks for all."
            for i in range(len(input)):
                os.symlink(os.path.abspath(input[i]),output[i])
#
# rule normalize_coverage_across_kmers:
#     input:
#         get_quality_controlled_reads #expect SE or R1,R2 or R1,R2,SE
#     output:
#         temp(expand("{{sample}}/assembly/reads/QC.normalized_{fraction}.fastq.gz",
#             fraction=MULTIFILE_FRACTIONS))
#     params:
#         k = config.get("normalization_kmer_length", NORMALIZATION_KMER_LENGTH),
#         t = config.get("normalization_target_depth", NORMALIZATION_TARGET_DEPTH),
#         minkmers = config.get("normalization_minimum_kmers", NORMALIZATION_MINIMUM_KMERS),
#         input_single = lambda wc, input: "in=%s" % input.se if hasattr(input, 'se') else "null",
#         extra_single = lambda wc, input: "extra={0},{1}" % (**input) if len(input)==3 else "",
#         has_paired_end_files = lambda wc, input: "t" if len(input)>1 else "f",
#         input_paired = lambda wc, input: "in=%s in2=%s" % (input.R1, input.R2) if hasattr(input, 'R1') else "null",
#         extra_paired = lambda wc, input: "extra=%s" % input.se if hasattr(input, 'se') else "",
#         output_single = lambda wc, output, input: "out=%s" % output[2] if hasattr(input, 'R1') else "out=%s" % output[0],
#         output_paired = lambda wc, output, input: "out=%s out2=%s" % (output[0], output[1]) if hasattr(input, 'R1') else "null",
#         tmpdir = "tmpdir=%s" % TMPDIR if TMPDIR else ""
#     log:
#         "{sample}/logs/assembly/pre_process/normalization.log"
#     benchmark:
#         "logs/benchmarks/assembly/pre_process/normalization/{sample}.txt"
#     conda:
#         "%s/required_packages.yaml" % CONDAENV
#     threads:
#         config.get("threads", 1)
#     resources:
#         mem = config["mem"],
#         java_mem = int(config["mem"] * JAVA_MEM_FRACTION)
#     shell:
#         """
#         if [ {params.input_single} != "null" ];
#         then
#             bbnorm.sh {params.input_single} \
#                 {params.extra_single} \
#                 {params.output_single} \
#                 {params.tmpdir} \
#                 k={params.k} target={params.t} \
#                 minkmers={params.minkmers} prefilter=t \
#                 threads={threads} \
#                 -Xmx{resources.java_mem}G 2> {log}
#         fi
#
#         if [ {params.has_paired_end_files} = "t" ];
#         then
#             bbnorm.sh {params.input_paired} \
#                 {params.extra_paired} \
#                 {params.output_paired} \
#                 {params.tmpdir} \
#                 k={params.k} target={params.t} \
#                 minkmers={params.minkmers} prefilter=t \
#                 threads={threads} \
#                 -Xmx{resources.java_mem}G 2>> {log}
#         fi
#         """



rule error_correction:
    input:
        expand("{{sample}}/assembly/reads/{{previous_steps}}_{fraction}.fastq.gz",
            fraction=MULTIFILE_FRACTIONS)
    output:
        temp(expand("{{sample}}/assembly/reads/{{previous_steps}}.errorcorr_{fraction}.fastq.gz",
            fraction=MULTIFILE_FRACTIONS))
    benchmark:
        "logs/benchmarks/assembly/pre_process/{sample}_error_correction_{previous_steps}.txt"
    log:
        "{sample}/logs/assembly/pre_process/error_correction_{previous_steps}.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    resources:
        mem = config["mem"],
        java_mem = int(config["mem"] * JAVA_MEM_FRACTION)
    params:
        inputs = lambda wc, input : io_params_for_tadpole(input),
        outputs = lambda wc, output: io_params_for_tadpole(output,key='out')
    threads:
        config.get("threads", 1)
    shell:
        """
        tadpole.sh -Xmx{resources.java_mem}G \
            prealloc=1 \
            {params.inputs} \
            {params.outputs} \
            mode=correct \
            threads={threads} \
            ecc=t ecco=t 2>> {log}
        """


rule merge_pairs:
    input:
        expand("{{sample}}/assembly/reads/{{previous_steps}}_{fraction}.fastq.gz",
            fraction=['R1','R2'])
    output:
        temp(expand("{{sample}}/assembly/reads/{{previous_steps}}.merged_{fraction}.fastq.gz",
            fraction=['R1','R2','me']))
    threads:
        config.get("threads", 1)
    resources:
        mem = config["mem"],
        java_mem = int(config["mem"] * JAVA_MEM_FRACTION)
    conda:
        "%s/required_packages.yaml" % CONDAENV
    log:
        "{sample}/logs/assembly/pre_process/merge_pairs_{previous_steps}.log"
    benchmark:
        "logs/benchmarks/assembly/pre_process/merge_pairs_{previous_steps}/{sample}.txt"
    shadow:
        "shallow"
    params:
        kmer = config.get("merging_k", MERGING_K),
        extend2 = config.get("merging_extend2", MERGING_EXTEND2),
        flags = config.get("merging_flags", MERGING_FLAGS)
    shell:
        """
        bbmerge.sh -Xmx{resources.java_mem}G threads={threads} \
            in1={input[0]} in2={input[1]} \
            outmerged={output[2]} \
            outu={output[0]} outu2={output[1]} \
            {params.flags} k={params.kmer} \
            extend2={params.extend2} 2> {log}
        """
localrules: passtrough_se_merged
rule passtrough_se_merged:
    input:
        "{sample}/assembly/reads/{previous_steps}_se.fastq.gz"
    output:
        temp("{sample}/assembly/reads/{previous_steps}.merged_se.fastq.gz")
    shell:
        "cp {input} {output}"


assembly_params={}
if config.get("assembler", "megahit") == "megahit":
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



    def megahit_input_parsing(input):
        Nfiles=len(input)

        if Nfiles==1:
            out= f"--read {input[0]}"
        else:
            out= f"-1 {input[0]} -2 {input[1]} "

            if Nfiles ==3:
                out+= f"--read {input[2]}"
        return out

    rule run_megahit:
        input:
            expand("{{sample}}/assembly/reads/{assembly_preprocessing_steps}_{fraction}.fastq.gz",
            fraction=ASSEMBLY_FRACTIONS, assembly_preprocessing_steps=assembly_preprocessing_steps)
        output:
            temp("{sample}/assembly/megahit/{sample}_prefilter.contigs.fa")
        benchmark:
            "logs/benchmarks/assembly/megahit/{sample}.txt"
#        shadow:
#            "shallow" #needs to be shallow to find input files
        log:
            "{sample}/logs/assembly/megahit.log"
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
            inputs = lambda wc, input: megahit_input_parsing(input),
            preset = assembly_params['megahit'][config['megahit_preset']],
        conda:
            "%s/assembly.yaml" % CONDAENV
        threads:
            config["assembly_threads"]
        resources:
            mem = config["assembly_memory"],
            time= config["runtime"]["assembly"]
        shell:
            """
                rm -r {params.outdir} 2> {log}

                megahit \
                {params.inputs} \
                --tmp-dir {TMPDIR} \
                --num-cpu-threads {threads} \
                --k-min {params.k_min} \
                --k-max {params.k_max} \
                --k-step {params.k_step} \
                --out-dir {params.outdir} \
                --out-prefix {wildcards.sample}_prefilter \
                --min-contig-len {params.min_contig_len} \
                --min-count {params.min_count} \
                --merge-level {params.merge_level} \
                --prune-level {params.prune_level} \
                --low-local-ratio {params.low_local_ratio} \
                --memory {resources.mem}000000000  \
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


else:


    ASSEMBLY_FRACTIONS = deepcopy(MULTIFILE_FRACTIONS)
    if config.get("merge_pairs_before_assembly", True) and PAIRED_END:
        ASSEMBLY_FRACTIONS += ['me']

    assembly_params['spades'] = {'meta':'--meta','normal':'', 'rna':'--rna'}

    def spades_parameters(wc,input):
        if not os.path.exists("{sample}/assembly/params.txt".format(sample=wc.sample)):

            params={}

            reads = dict(zip(ASSEMBLY_FRACTIONS,input))

            if not PAIRED_END:
                params['inputs']= " -s {se} ".format(**reads)
            else:
                params['inputs']= " --pe1-1 {R1} --pe1-2 {R2} ".format(**reads)

                if 'se' in ASSEMBLY_FRACTIONS:
                    params['inputs']+= "--pe1-s {se} ".format(**reads)
                if 'me' in ASSEMBLY_FRACTIONS:
                    params['inputs']+= "--pe1-m {me} ".format(**reads)

            # Long reads:

            if (config['longread_type'] is not None) & (str(config['longread_type']).lower()!='none'):

                long_read_file = get_files_from_sampleTable(wc.sample,'longreads')[0]
                params['longreads'] = " --{t} {f} ".format(t=config['longread_type'],f=long_read_file)
            else:
                params['longreads'] = ""


            params['preset'] = assembly_params['spades'][config['spades_preset']]
            params['skip_error_correction'] = "--only-assembler" if config['spades_skip_BayesHammer'] else ""
            params['extra'] = config['spades_extra']


        else:

            params = {"inputs": "--restart-from last",
                      "preset":"",
                      "skip_error_correction":"",
                      "extra":"",
                      "longreads":""}

        params['outdir']= "{sample}/assembly".format(sample=wc.sample)

        return params


    rule run_spades:
        input:
            expand("{{sample}}/assembly/reads/{assembly_preprocessing_steps}_{fraction}.fastq.gz",
                fraction=ASSEMBLY_FRACTIONS,
                assembly_preprocessing_steps=assembly_preprocessing_steps)
        output:
            "{sample}/assembly/contigs.fasta",
            "{sample}/assembly/scaffolds.fasta"
        benchmark:
            "logs/benchmarks/assembly/spades/{sample}.txt"
        params:
            p= lambda wc,input: spades_parameters(wc,input),
            k = config.get("spades_k", SPADES_K),
        log:
            "{sample}/logs/assembly/spades.log"
        conda:
            "%s/assembly.yaml" % CONDAENV
        threads:
            config["assembly_threads"]
        resources:
            mem = config["assembly_memory"],
            time= config["runtime"]["assembly"]
        shell:
            "spades.py "
            " --threads {threads} "
            " --memory {resources.mem} "
            " -o {params.p[outdir]} "
            " -k {params.k}"
            " {params.p[preset]} "
            " {params.p[extra]} "
            " {params.p[inputs]} "
            " {params.p[longreads]} "
            " {params.p[skip_error_correction]} "
            " > {log} 2>&1 "



    localrules: rename_spades_output
    rule rename_spades_output:
        input:
            "{{sample}}/assembly/{sequences}.fasta".format(sequences= 'scaffolds' if config['spades_use_scaffolds'] else 'contigs' )
        output:
            temp("{sample}/assembly/{sample}_raw_contigs.fasta")
        shell:
            "cp {input} {output}"


rule rename_contigs:
    # standardizes header labels within contig FASTAs
    input:
        "{sample}/assembly/{sample}_raw_contigs.fasta"
    output:
        "{sample}/assembly/{sample}_prefilter_contigs.fasta"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    shell:
        """rename.sh in={input} out={output} ow=t prefix={wildcards.sample}"""


rule calculate_contigs_stats:
    input:
        "{sample}/assembly/{sample}_{assembly_step}_contigs.fasta"
    output:
        "{sample}/assembly/contig_stats/{assembly_step}_contig_stats.txt"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        1
    resources:
        mem = 1
    shell:
        "stats.sh in={input} format=3 > {output}"


rule combine_sample_contig_stats:
    input:
        expand("{{sample}}/assembly/contig_stats/{assembly_step}_contig_stats.txt", assembly_step= ['prefilter','final'])
    output:
        "{sample}/assembly/contig_stats.tsv"
    run:
        import os
        import pandas as pd

        c = pd.DataFrame()
        for f in input:
            df = pd.read_csv(f,sep='\t')
            assembly_step = os.path.basename(f).replace("_contig_stats.txt", "")
            c.loc[assembly_step]

        c.to_csv(output[0], sep='\t')

if config['filter_contigs']:

    ruleorder: align_reads_to_prefilter_contigs > align_reads_to_final_contigs

    rule align_reads_to_prefilter_contigs:
        input:
            reads= get_quality_controlled_reads,
            fasta = rules.rename_contigs.output
        output:
            sam = temp("{sample}/sequence_alignment/alignment_to_prefilter_contigs.sam")
        params:
            input = lambda wc, input : input_params_for_bbwrap( input.reads),
            maxsites = config.get("maximum_counted_map_sites", MAXIMUM_COUNTED_MAP_SITES),
            max_distance_between_pairs = config.get('contig_max_distance_between_pairs', CONTIG_MAX_DISTANCE_BETWEEN_PAIRS),
            paired_only = 't' if config.get("contig_map_paired_only", CONTIG_MAP_PAIRED_ONLY) else 'f',
            min_id = config.get('contig_min_id', CONTIG_MIN_ID),
            maxindel = 100,
            ambiguous = 'all'
        shadow:
            "shallow"
        log:
            "{sample}/logs/assembly/post_process/align_reads_to_prefiltered_contigs.log"
        conda:
            "%s/required_packages.yaml" % CONDAENV
        threads:
            config.get("threads", 1)
        resources:
            mem = config["mem"],
            java_mem = int(config["mem"] * JAVA_MEM_FRACTION)
        shell:
            """
                bbwrap.sh nodisk=t \
                ref={input.fasta} \
                {params.input} \
                trimreaddescriptions=t \
                out={output.sam} \
                threads={threads} \
                pairlen={params.max_distance_between_pairs} \
                pairedonly={params.paired_only} \
                minid={params.min_id} \
                mdtag=t \
                xstag=fs \
                nmtag=t \
                sam=1.3 \
                local=t \
                ambiguous={params.ambiguous} \
                secondary=t \
                saa=f \
                append=t \
                machineout=t \
                maxsites={params.maxsites} \
                -Xmx{resources.java_mem}G 2> {log}
            """


    rule pileup_prefilter:
        input:
            fasta = "{sample}/assembly/{sample}_prefilter_contigs.fasta",
            sam = "{sample}/sequence_alignment/alignment_to_prefilter_contigs.sam"
        output:
            covstats = "{sample}/assembly/contig_stats/prefilter_coverage_stats.txt"
        params:
            pileup_secondary = 't'
        log:
            "{sample}/logs/assembly/post_process/pilup_prefilter_contigs.log"
        conda:
            "%s/required_packages.yaml" % CONDAENV
        threads:
            config.get("threads", 1)
        resources:
            mem = config["mem"],
            java_mem = int(config["mem"] * JAVA_MEM_FRACTION)
        shell:
            """
            pileup.sh ref={input.fasta} in={input.sam} \
                threads={threads} \
                -Xmx{resources.java_mem}G \
                covstats={output.covstats} \
                concise=t \
                secondary={params.pileup_secondary} 2> {log}
            """


    rule filter_by_coverage:
        input:
            fasta = "{sample}/assembly/{sample}_prefilter_contigs.fasta",
            covstats = "{sample}/assembly/contig_stats/prefilter_coverage_stats.txt"
        output:
            fasta = "{sample}/assembly/{sample}_final_contigs.fasta",
            removed_names = "{sample}/assembly/{sample}_discarded_contigs.fasta"
        params:
            minc = config.get("minimum_average_coverage", MINIMUM_AVERAGE_COVERAGE),
            minp = config.get("minimum_percent_covered_bases", MINIMUM_PERCENT_COVERED_BASES),
            minr = config.get("minimum_mapped_reads", MINIMUM_MAPPED_READS),
            minl = config.get("minimum_contig_length", MINIMUM_CONTIG_LENGTH),
            trim = config.get("contig_trim_bp", CONTIG_TRIM_BP)
        log:
            "{sample}/logs/assembly/post_process/filter_by_coverage.log"
        conda:
            "%s/required_packages.yaml" % CONDAENV
        threads:
            1
        resources:
            mem = config["mem"],
            java_mem = int(config["mem"] * JAVA_MEM_FRACTION)
        shell:
            """filterbycoverage.sh in={input.fasta} \
                   cov={input.covstats} \
                   out={output.fasta} \
                   outd={output.removed_names} \
                   minc={params.minc} \
                   minp={params.minp} \
                   minr={params.minr} \
                   minl={params.minl} \
                   trim={params.trim} \
                   -Xmx{resources.java_mem}G 2> {log}"""

# HACK: this makes two copies of the same file
else: # no filter
    localrules: do_not_filter_contigs
    rule do_not_filter_contigs:
        input:
            rules.rename_contigs.output
        output:
            "{sample}/assembly/{sample}_final_contigs.fasta"
        threads:
            1
        shell:
            "cp {input} {output}"


localrules: finalize_contigs
rule finalize_contigs:
    input:
        "{sample}/assembly/{sample}_final_contigs.fasta"
    output:
        "{sample}/{sample}_contigs.fasta"
    threads:
        1
    run:
        os.symlink(os.path.relpath(input[0],os.path.dirname(output[0])),output[0])


# generalized rule so that reads from any "sample" can be aligned to contigs from "sample_contigs"
rule align_reads_to_final_contigs:
    input:
        reads=get_quality_controlled_reads,
        fasta = "{sample_contigs}/{sample_contigs}_contigs.fasta",
    output:
        sam = temp("{sample_contigs}/sequence_alignment/{sample}.sam"),
        #unmapped = temp(expand("{{sample_contigs}}/assembly/unmapped_post_filter/{{sample}}_unmapped_{fraction}.fastq.gz",
        #                  fraction=MULTIFILE_FRACTIONS))
    params:
        input = lambda wc, input : input_params_for_bbwrap( input.reads),
        maxsites = config.get("maximum_counted_map_sites", MAXIMUM_COUNTED_MAP_SITES),
        #unmapped = lambda wc, output: "outu1={0},{2} outu2={1},null".format(*output.unmapped) if PAIRED_END else "outu={0}".format(*output.unmapped),
        max_distance_between_pairs = config.get('contig_max_distance_between_pairs', CONTIG_MAX_DISTANCE_BETWEEN_PAIRS),
        paired_only = 't' if config.get("contig_map_paired_only", CONTIG_MAP_PAIRED_ONLY) else 'f',
        ambiguous = 'all' if CONTIG_COUNT_MULTI_MAPPED_READS else 'best',
        min_id = config.get('contig_min_id', CONTIG_MIN_ID),
        maxindel = 100 # default 16000 good for genome deletions but not necessarily for alignment to contigs
    shadow:
        "shallow"
    benchmark:
        "logs/benchmarks/assembly/calculate_coverage/align_reads_to_filtered_contigs/{sample}_to_{sample_contigs}.txt"
    log:
        "{sample_contigs}/logs/assembly/calculate_coverage/align_reads_from_{sample}_to_filtered_contigs.log" # this file is udes for assembly report
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config["mem"],
        java_mem = int(config["mem"] * JAVA_MEM_FRACTION)
    shell:
        """
        bbwrap.sh nodisk=t \
            ref={input.fasta} \
            {params.input} \
            trimreaddescriptions=t \
            out={output.sam} \
            threads={threads} \
            pairlen={params.max_distance_between_pairs} \
            pairedonly={params.paired_only} \
            minid={params.min_id} \
            mdtag=t \
            xstag=fs \
            nmtag=t \
            sam=1.3 \
            local=t \
            ambiguous={params.ambiguous} \
            secondary=t \
            saa=f \
            maxsites={params.maxsites} \
            -Xmx{resources.java_mem}G \
            2> {log}
        """

ruleorder: bam_2_sam_contigs > convert_sam_to_bam

rule convert_sam_to_bam:
    input:
        "{file}.sam"
    output:
        bam="{file}.bam",
        sam=temp("{file}.forpileup.sam"),
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("simplejob_threads", 1)
    shadow:
        "shallow"
    resources:
        mem = 2 * config.get("simplejob_threads", 1)
    shell:
        """samtools view \
               -m 1G \
               -@ {threads} \
               -bSh1 {input} | samtools sort \
                                   -m 1G \
                                   -@ {threads} \
                                   -T {wildcards.file}_tmp \
                                   -o {output.bam} \
                                   -O bam -
            mv {input} {output.sam}
        """


rule bam_2_sam_contigs:
    input:
        "{sample}/sequence_alignment/{sample}.bam"
    output:
        temp("{sample}/sequence_alignment/{sample}.forpileup.sam")
    threads:
        config['simplejob_threads']
    resources:
        mem = config["simplejob_mem"],
    shadow:
        "shallow"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    shell:
        """
        reformat.sh in={input} out={output} sam=1.3
        """


rule pileup:
    input:
        fasta = "{sample}/{sample}_contigs.fasta",
        sam = "{sample}/sequence_alignment/{sample}.forpileup.sam",
    output:
        basecov = temp("{sample}/assembly/contig_stats/postfilter_base_coverage.txt.gz"),
        covhist = "{sample}/assembly/contig_stats/postfilter_coverage_histogram.txt",
        covstats = "{sample}/assembly/contig_stats/postfilter_coverage_stats.txt",
        bincov = "{sample}/assembly/contig_stats/postfilter_coverage_binned.txt"
    params:
        pileup_secondary = 't' if config.get("count_multi_mapped_reads", CONTIG_COUNT_MULTI_MAPPED_READS) else 'f',
    benchmark:
        "logs/benchmarks/assembly/calculate_coverage/pileup/{sample}.txt"
    log:
        "{sample}/logs/assembly/calculate_coverage/pilup_final_contigs.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config["mem"],
        java_mem = int(config["mem"] * JAVA_MEM_FRACTION)
    shell:
        """pileup.sh ref={input.fasta} in={input.sam} \
               threads={threads} \
               -Xmx{resources.java_mem}G \
               covstats={output.covstats} \
               hist={output.covhist} \
               basecov={output.basecov}\
               concise=t \
               secondary={params.pileup_secondary} \
               bincov={output.bincov} 2> {log}"""




rule create_bam_index:
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        1
    resources:
        mem = 2 * config.get("simplejob_threads", 1)
    shell:
        "samtools index {input}"

rule predict_genes:
    input:
        "{sample}/{sample}_contigs.fasta"
    output:
        fna = "{sample}/annotation/predicted_genes/{sample}.fna",
        faa = "{sample}/annotation/predicted_genes/{sample}.faa",
        gff = "{sample}/annotation/predicted_genes/{sample}.gff"
    conda:
        "%s/prodigal.yaml" % CONDAENV
    log:
        "{sample}/logs/gene_annotation/prodigal.txt"
    benchmark:
        "logs/benchmarks/prodigal/{sample}.txt"
    threads:
        1
    shell:
        """
        prodigal -i {input} -o {output.gff} -d {output.fna} \
            -a {output.faa} -p meta -f gff 2> {log}
        """



localrules: get_contigs_from_gene_names
rule get_contigs_from_gene_names:
    input:
        faa = "{sample}/annotation/predicted_genes/{sample}.faa",
    output:
        tsv= "{sample}/annotation/predicted_genes/{sample}.tsv"
    run:
        header = ["gene_id", "Contig", "Gene_nr", "Start", "Stop", "Strand", "Annotation"]
        with open(output.tsv, "w") as tsv:
            tsv.write("\t".join(header) + "\n")
            with open(input.faa) as fin:
                gene_idx = 0
                for line in fin:
                    if line[0] == ">":
                        text = line[1:].strip().split(" # ")
                        old_gene_name = text[0]
                        text.remove(old_gene_name)
                        old_gene_name_split = old_gene_name.split("_")
                        gene_nr = old_gene_name_split[-1]
                        contig_nr = old_gene_name_split[-2]
                        sample = "_".join(old_gene_name_split[:len(old_gene_name_split) - 2])
                        tsv.write("{gene_id}\t{sample}_{contig_nr}\t{gene_nr}\t{text}\n".format(
                            text="\t".join(text),
                            gene_id=old_gene_name,
                            i=gene_idx,
                            sample=sample,
                            gene_nr=gene_nr,
                            contig_nr=contig_nr))
                        gene_idx += 1
    #



localrules: build_assembly_report
rule build_assembly_report:
    input:
        contig_stats = expand("{sample}/assembly/contig_stats/final_contig_stats.txt", sample=SAMPLES),
        gene_tables = expand("{sample}/annotation/predicted_genes/{sample}.tsv", sample=SAMPLES),
        mapping_logs = expand("{sample}/logs/assembly/calculate_coverage/align_reads_from_{sample}_to_filtered_contigs.log", sample=SAMPLES),
        # mapping logs will be incomplete unless we wait on alignment to finish
        bams = expand("{sample}/sequence_alignment/{sample}.bam", sample=SAMPLES)
    output:
        report = "reports/assembly_report.html",
        combined_contig_stats = 'stats/combined_contig_stats.tsv'
    params:
        samples = SAMPLES
    conda:
        "%s/report.yaml" % CONDAENV
    log:
        "logs/assembly/report.log"
    script:
        "../report/assembly_report.py"
