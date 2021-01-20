
wildcard_constraints:
    run="[S,E,D]RR[0-9]+"

rule grab_seq:
    output:
        "SRAreads/{run}_1.fastq.gz",
        "SRAreads/{run}_2.fastq.gz",
    params:
        outdir='SRAreads'
    log:
        "log/SRAdownload/{run}.log"
    threads:
        config['simplejob_threads']
    resources:
        time= lambda wildcards, attempt: attempt * 12, #h
        mem=config['simplejob_mem'],
    conda:
        "%s/grabseq.yaml" % CONDAENV
    shadow:
        "shallow"
    shell:
        "grabseqs sra {wildcards.run}"
        " -t {threads} "
        " -o {params.outdir}"
        " &> {log} "
