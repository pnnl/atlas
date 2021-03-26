wildcard_constraints:
    sra_run="[S,E,D]RR[0-9]+"

localrules: prefetch
rule prefetch:
    output:
        sra=temp(touch("SRAreads/{sra_run}_downloaded")),
        # not givins sra file as output allows for continue from the same download
    params:
        outdir= 'SRAreads' #lambda wc,output: os.path.dirname(output[0])
    log:
        "log/SRAdownload/{sra_run}.log"
    benchmark:
        "log/benchmarks/SRAdownload/prefetch/{sra_run}.tsv"
    threads:
        1
    resources:
        mem=1,
        time= int(config["runtime"]["simple_job"]),
        internet_connection=1
    conda:
        "%s/sra.yaml" % CONDAENV
    shell:
        " mkdir -p {params.outdir} 2> {log} "
        " ; "
        " prefetch "
        " --output-directory {params.outdir} "
        " -X 999999999 "
        " --progress "
        " --log-level info "
        " {wildcards.sra_run} &>> {log} "
        " ; "
        " vdb-validate {params.outdir}/{wildcards.sra_run}/{wildcards.sra_run}.sra &>> {log} "




rule extract_run:
    input:
        flag=rules.prefetch.output,
    output:
        expand("SRAreads/{{sra_run}}_{fraction}.fastq.gz",
                fraction= ['1','2']
                 )
    params:
        outdir=os.path.abspath('SRAreads'),
        sra_file= "SRAreads/{sra_run}/{sra_run}.sra",
        tmpdir= TMPDIR
    log:
        "log/SRAdownload/{sra_run}.log"
    benchmark:
        "log/benchmarks/SRAdownload/fasterqdump/{sra_run}.tsv"
    threads:
        config['simplejob_threads']
    resources:
        time= int(config["runtime"]["simple_job"]),
        mem=1 #default 100Mb
    conda:
        "%s/sra.yaml" % CONDAENV
    shell:
        " vdb-validate {params.sra_file} &>> {log} "
        " ; "
        " parallel-fastq-dump "
        " --threads {threads} "
        " --gzip --split-files "
        " --outdir {params.outdir} "
        " --tmpdir {params.tmpdir} "
        " --skip-technical --split-3 "
        " -s {params.sra_file} &> {log} "
        " ; "
        " rm -rf {params.outdir}/{wildcards.sra_run} 2>> {log} "


rule download_all_reads:
    input:
        expand("SRAreads/{sample}_{fraction}.fastq.gz",sample=SAMPLES,fraction=['1','2'])
