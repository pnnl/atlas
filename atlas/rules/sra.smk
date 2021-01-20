wildcard_constraints:
    SRR="[S,E,D]RR[0-9]+"

#localrules: prefetch
rule prefetch:
    output:
        sra=touch(temp("SRAreads/{run}/downloaded")),
        # not givins sra file as output allows for continue from the same download
    params:
        outdir= lambda wc,input: os.path.dirname(input[0])
    log:
        "log/SRAdownload/{run}.log"
    benchmark:
        "log/benchmarks/logs/SRAdownload/prefetch/{run}.tsv"
    threads:
        1
    resources:
        mem=1,
        time= int(config["runtime"]["simple_job"]),
        internet_connection=1
    conda:
        "%s/sra.yaml" % CONDAENV
    shell:
        " mkdir -p {params.outdir} 2> {log}; "
        " prefetch "
        " --output-directory {params.outdir} "
        " -X 999999999 "
        " {wildcards.run} &>> {log}"


rule extract_run:
    input:
        flag=rules.prefetch.output,
    output:
        "SRAreads/{run}_1.fastq",
        "SRAreads/{run}_2.fastq",
    params:
        outdir='SRAreads',
        tmpdir= TMPDIR
    log:
        "log/SRAdownload/{run}.log"
    threads:
        config['simplejob_threads']
    resources:
        time= lambda wildcards, attempt: attempt * 6, #h
        mem=config['simplejob_mem'],
    conda:
        "%s/sra.yaml" % CONDAENV
    shell:
        " cd {params.outdir} 2>> {log};"
        " "
        " fasterq-dump "
        " --threads {threads} "
        " --temp {params.tmpdir} "
        " --log-level info "
        " {wildcards.run} "
        " &>> ../{log} ; "
        " "
        " rm -rf {wildcards.run} 2> ../{log} "
        #" pigz -1 -p{threads} {wildcards.run}_1.fastq"
