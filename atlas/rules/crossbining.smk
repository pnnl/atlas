### VAMB




localrules: combine_contigs,vamb

rule vamb:
    input:
        "vamb/clustering"



rule combine_contigs:
    input:
        expand("{sample}/{sample}_contigs.fasta",sample=SAMPLES)
    output:
        "vamb/combined_contigs.fasta.gz"
    threads:
        1
    conda:
        "../envs/vamb.yaml"
    shell:
        "concatenate.py {output} {input} -m 2000 --keepnames"
        #"cat {input} > {output}"


rule minimap_index:
    input:
        contigs = rules.combine_contigs.output
    output:
        mmi = temp("vamb/combined_contigs.mmi")
    params:
        index_size= "12G"
    resources:
        mem= config['large_mem'],
    threads:
        1
    log:
        "log/vamb/index.log"
    benchmark:
        "log/benchmarks/vamb/mminimap_index.tsv"
    conda:
        "../envs/minimap.yaml"
    shell:
        "minimap2 -I {params.index_size} -d {output} {input} 2> {log}"

rule samtools_dict:
    input:
        contigs = rules.combine_contigs.output
    output:
        dict = "vamb/combined_contigs.dict"
    resources:
        mem = config['simplejob_mem']
    threads:
        1
    log:
        "log/vamb/samtools_dict.log"
    conda:
        "../envs/minimap.yaml"
    shell:
        "samtools dict {input} | cut -f1-3 > {output} 2> {log}"

rule minimap:
    input:
        fq = lambda wildcards: input_paired_only(get_quality_controlled_reads(wildcards)),
        mmi = "vamb/combined_contigs.mmi",
        dict = "vamb/combined_contigs.dict"
    output:
        bam = temp("vamb/mapped/{sample}.bam")
    threads:
        config['threads']
    log:
        "log/vamb/mapping/{sample}.minimap.log"
    benchmark:
        "log/benchmarks/vamb/mminimap/{sample}.tsv"
    conda:
        "../envs/minimap.yaml"
    shell:
        '''minimap2 -t {threads} -ax sr {input.mmi} {input.fq} | grep -v "^@" | cat {input.dict} - | samtools view -F 3584 -b - > {output.bam} 2>{log}'''

ruleorder: sort_bam > minimap > convert_sam_to_bam
rule sort_bam:
    input:
        "vamb/mapped/{sample}.bam"
    output:
        temp("vamb/mapped/{sample}.sort.bam")
    params:
        prefix="vamb/mapped/tmp.{sample}"
    threads:
        2
    resources:
        mem = config['simplejob_mem'],
        time = int(config["runtime"]["simple_job"])
    log:
        "log/vamb/mapping/{sample}.sortbam.log"
    conda:
        "../envs/minimap.yaml"
    shell:
        "samtools sort {input} -T {params.prefix} --threads {threads} -m 3G -o {output} 2>{log}"


rule summarize_bam_contig_depths:
    input:
        bam = expand(rules.sort_bam.output, sample=SAMPLES)
    output:
        "vamb/coverage.jgi.tsv"
    log:
        "log/vamb/combine_coverage.log"
    conda:
        "../envs/metabat.yaml"
    threads:
        config['threads']
    resources:
        mem = config["mem"]
    shell:
        "jgi_summarize_bam_contig_depths "
        " --outputDepth {output} "
        " {input.bam} &> {log} "

localrules: convert_jgi2vamb_coverage
rule convert_jgi2vamb_coverage:
    input:
        "vamb/coverage.jgi.tsv"
    output:
        "vamb/coverage.tsv"
    log:
        "log/vamb/convert_jgi2vamb_coverage.log"
    threads:
        1
    script:
        "../scripts/convert_jgi2vamb_coverage.py"

rule run_vamb:
    input:
        coverage = "vamb/coverage.tsv",
        fasta = rules.combine_contigs.output,
    output:
        directory("vamb/clustering")
    conda:
        "../envs/vamb.yaml"
    threads:
         config["threads"]
    resources:
        mem = config["mem"]
    log:
        "log/vamb/run_vamb.log"
    benchmark:
        "log/benchmarks/vamb/run_vamb.tsv"
    params:
        params= "-m 2000 --minfasta 500000"
    shell:
        "vamb --outdir {output} "
        " -o '_' "
        " --jgi {input.coverage} "
        " --fasta {input.fasta} "
        " {params.params} "
        "2> {log}"
