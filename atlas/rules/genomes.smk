if 'genome_dir' in config:
    genome_dir= config['genome_dir']
    assert os.path.exists(genome_dir), f"genome_dir ({genome_dir}) doesn't exist"
else:
    genome_dir="genomes/genomes"


## dRep
localrules: get_all_bins, all_contigs2bins
rule get_all_bins:
    input:
        bins=expand("{sample}/binning/{binner}/bins",
               sample= SAMPLES, binner= config['final_binner']),
        #cluster_attribution=expand("{sample}/binning/{binner}/cluster_attribution.tsv",
        #       sample= SAMPLES, binner= config['final_binner'])
    output:
        temp(directory("genomes/all_bins"))
    run:
        os.mkdir(output[0])
        from utils.io import symlink_relative

        for bin_folder in input.bins:

            fasta_files = [f for f in os.listdir(bin_folder) if f.endswith('.fasta') ]
            symlink_relative(fasta_files,bin_folder,output[0])

rule all_contigs2bins:
        input:
            expand("{sample}/binning/{binner}/cluster_attribution.tsv",
                   sample= SAMPLES, binner= config['final_binner'])
        output:
            temp("genomes/clustering/all_contigs2bins.tsv.gz")
        run:
            from utils.io import cat_files
            cat_files(input,output[0],gzip=True)


localrules: get_quality_for_dRep_from_checkm, merge_checkm
rule get_quality_for_dRep_from_checkm:
    input:
        "reports/genomic_bins_{binner}.tsv".format(binner=config['final_binner'])
    output:
        temp("genomes/quality.csv")
    run:
        import pandas as pd

        D= pd.read_csv(input[0],index_col=0,sep='\t')

        D.index=D.index.astype(str)+".fasta"
        D.index.name="genome"
        D.columns= D.columns.str.lower()
        D.iloc[:,:3].to_csv(output[0])


rule merge_checkm:
    input:
        completeness= expand("{sample}/binning/{binner}/checkm/completeness.tsv",
               sample= SAMPLES, binner= config['final_binner']),
        taxonomy= expand("{sample}/binning/{binner}/checkm/taxonomy.tsv",
               sample= SAMPLES, binner= config['final_binner']),
        markers= expand("{sample}/binning/{binner}/checkm/storage/tree/concatenated.fasta",
               sample= SAMPLES, binner= config['final_binner'])
    output:
        checkm="genomes/checkm/checkm_all_bins.tsv",
        markers= "genomes/checkm/all_bins_markers.fasta"
    run:

        import pandas as pd
        import shutil
        D=[]

        for i in range(len(SAMPLES)):
            df= pd.read_csv(input.completeness[i],index_col=0,sep='\t')
            df= df.join(pd.read_csv(input.taxonomy[i],index_col=0,sep='\t'))
            D.append(df)

        D= pd.concat(D,axis=0)
        D.to_csv(output.checkm,sep='\t')

        with open(output.markers,'wb') as fout:
            for fasta in input.markers:
                shutil.copyfileobj(open(fasta,'rb'),fout)




rule first_dereplication:
    input:
        "genomes/all_bins",
        quality= "genomes/quality.csv"
    output:
        directory("genomes/pre_dereplication/dereplicated_genomes")
    threads:
        config['threads']
    log:
        "logs/genomes/pre_dereplication.log"
    conda:
        "%s/dRep.yaml" % CONDAENV
    shadow:
        "shallow"
    params:
        filter= " --noQualityFiltering " if config['genome_dereplication']['filter']['noFilter'] else "",
        filter_length= config['genome_dereplication']['filter']['length'],
        filter_completeness= config['genome_dereplication']['filter']['completeness'],
        filter_contamination= config['genome_dereplication']['filter']['contamination'],
        ANI= config['genome_dereplication']['ANI'],
        completeness_weight= config['genome_dereplication']['score']['completeness'] ,
        contamination_weight=config['genome_dereplication']['score']['contamination'] ,
        strain_heterogeneity_weight= config['genome_dereplication']['score']['completeness'] , #not in table
        N50_weight=config['genome_dereplication']['score']['N50'] ,
        size_weight=config['genome_dereplication']['score']['length'] ,
        opt_parameters = config['genome_dereplication']['opt_parameters'],
        work_directory= lambda wc,output: os.path.dirname(output[0]),
        sketch_size= config['genome_dereplication']['sketch_size']

    shell:
        " rm -rf {params.work_directory} ;"
        " dRep dereplicate "
        " {params.filter} "
        " --genomes {input[0]}/*.fasta "
        " --genomeInfo {input.quality} "
        " --length {params.filter_length} "
        " --completeness {params.filter_completeness} "
        " --contamination {params.filter_contamination} "
        " --SkipSecondary "
        " --P_ani {params.ANI} "
        " --completeness_weight {params.completeness_weight} "
        " --contamination_weight {params.contamination_weight} "
        " --strain_heterogeneity_weight {params.strain_heterogeneity_weight} "
        " --N50_weight {params.N50_weight} "
        " --size_weight {params.size_weight} "
        " --MASH_sketch {params.sketch_size} "
        " --processors {threads} "
        " {params.opt_parameters} "
        " {params.work_directory} "
        " &> {log} "

rule second_dereplication:
    input:
        rules.first_dereplication.output,
        quality= "genomes/quality.csv"
    output:
        directory("genomes/Dereplication/dereplicated_genomes")
    threads:
        config['threads']
    log:
        "logs/genomes/dereplication.log"
    shadow:
        "shallow"
    conda:
        "%s/dRep.yaml" % CONDAENV
    params:
        ANI= config['genome_dereplication']['ANI'],
        completeness_weight= config['genome_dereplication']['score']['completeness'] ,
        contamination_weight=config['genome_dereplication']['score']['contamination'] ,
        strain_heterogeneity_weight= config['genome_dereplication']['score']['completeness'] , #not in table
        N50_weight=config['genome_dereplication']['score']['N50'] ,
        size_weight=config['genome_dereplication']['score']['length'] ,
        opt_parameters = config['genome_dereplication']['opt_parameters'],
        work_directory= lambda wc,output: os.path.dirname(output[0]),
        sketch_size= config['genome_dereplication']['sketch_size']

    shell:
        " rm -rf {params.work_directory} ;"
        " dRep dereplicate "
        " --genomes {input[0]}/*.fasta "
        " --genomeInfo {input.quality} "
        " --noQualityFiltering "
        " --S_ani {params.ANI} "
        " --completeness_weight {params.completeness_weight} "
        " --contamination_weight {params.contamination_weight} "
        " --strain_heterogeneity_weight {params.strain_heterogeneity_weight} "
        " --N50_weight {params.N50_weight} "
        " --size_weight {params.size_weight} "
        " --MASH_sketch {params.sketch_size} "
        " --processors {threads} "
        " {params.opt_parameters} "
        " {params.work_directory} "
        " &> {log} "


localrules: rename_genomes
checkpoint rename_genomes:
    input:
        genomes="genomes/Dereplication/dereplicated_genomes",
        pre_dereplication= "genomes/pre_dereplication/dereplicated_genomes"
    output:
        dir= directory("genomes/genomes"),
        mapfile_contigs="genomes/clustering/contig2genome.tsv",
        mapfile_genomes = "genomes/clustering/old2newID.tsv",
        mapfile_bins= "genomes/clustering/allbins2genome.tsv"
    params:
        rename_contigs=config['rename_mags_contigs']
    shadow:
        "shallow"
    log:
        "logs/genomes/rename_genomes.log"
    script:
        "rename_genomes.py"




def get_genomes_(wildcards):

    if genome_dir=='genomes/genomes':
        checkpoints.rename_genomes.get() # test if checkpoint passed

    genomes= glob_wildcards(os.path.join(genome_dir, "{genome}.fasta")).genome

    if len(genomes)==0:
        logger.critical("No genomes found after dereplication. "
                        "You don't have any Metagenome assembled genomes with sufficient quality. "
                        "You may want to change the assembly, binning or filtering parameters. "
                        "Or focus on the genecatalog workflow only.")
        exit(1)

    return genomes



rule run_all_checkm_lineage_wf:
    input:
        touched_output = "logs/checkm_init.txt",
        dir = genome_dir
    output:
        "genomes/checkm/completeness.tsv",
        "genomes/checkm/storage/tree/concatenated.fasta"
    params:
        output_dir = lambda wc, output: os.path.dirname(output[0]),
    conda:
        "%s/checkm.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        time=config["runtime"]["long"],
        mem=config["large_mem"]
    log:
        "logs/genomes/checkm.log"
    benchmark:
        "logs/benchmarks/checkm_lineage_wf/all_genomes.tsv"
    shell:
        """
        rm -r {params.output_dir}
        checkm lineage_wf \
            --file {params.output_dir}/completeness.tsv \
            --tab_table \
            --quiet \
            --extension fasta \
            --threads {threads} \
            {input.dir} \
            {params.output_dir} &> {log}
        """

### Quantification




rule build_db_genomes:
    input:
        genome_dir
    output:
        index="ref/genome/3/summary.txt",
        fasta=temp("genomes/all_contigs.fasta")
    threads:
        config.get("threads", 6)
    resources:
        mem = config["mem"],
        java_mem = int(config["mem"] * JAVA_MEM_FRACTION)
    log:
        "logs/genomes/mapping/build_bbmap_index.log"
    shell:
        """
        cat {input}/*.fasta > {output.fasta} 2> {log}
        bbmap.sh build=3 -Xmx{resources.java_mem}G ref={output.fasta} threads={threads} local=f 2>> {log}

        """





# generalized rule so that reads from any "sample" can be aligned to contigs from "sample_contigs"
rule align_reads_to_MAGs:
    input:
        reads=get_quality_controlled_reads,
        ref = rules.build_db_genomes.output.index,
    output:
        sam = temp("genomes/alignments/{sample}.sam"),
        unmapped = expand("genomes/alignments/unmapped/{{sample}}_{fraction}.fastq.gz",
                          fraction=MULTIFILE_FRACTIONS)
    params:
        input = lambda wc, input : input_params_for_bbwrap(input.reads),
        unmapped = lambda wc, output: io_params_for_tadpole(output.unmapped,'outu'),
        maxsites = config.get("maximum_counted_map_sites", MAXIMUM_COUNTED_MAP_SITES),
        max_distance_between_pairs = config.get('contig_max_distance_between_pairs', CONTIG_MAX_DISTANCE_BETWEEN_PAIRS),
        paired_only = 't' if config.get("contig_map_paired_only", CONTIG_MAP_PAIRED_ONLY) else 'f',
        ambiguous = 'all' if CONTIG_COUNT_MULTI_MAPPED_READS else 'best',
        min_id = config.get('contig_min_id', CONTIG_MIN_ID),
        maxindel = 100 # default 16000 good for genome deletions but not necessarily for alignment to contigs
    shadow:
        "shallow"
    log:
        "logs/genomes/mapping/map_{sample}.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config["mem"],
        java_mem = int(config["mem"] * JAVA_MEM_FRACTION)
    shell:
        """
            bbwrap.sh \
            build=3 \
            {params.input} \
            trimreaddescriptions=t \
            outm={output.sam} \
            {params.unmapped} \
            threads={threads} \
            pairlen={params.max_distance_between_pairs} \
            pairedonly={params.paired_only} \
            minid={params.min_id} \
            mdtag=t \
            xstag=fs \
            nmtag=t \
            sam=1.3 \
            ambiguous={params.ambiguous} \
            secondary=t \
            saa=f \
            maxsites={params.maxsites} \
            -Xmx{resources.java_mem}G \
            2> {log}
        """




rule pileup_MAGs:
    input:
        sam = "genomes/alignments/{sample}.sam",
    output:
        basecov = temp("genomes/alignments/{sample}_base_coverage.txt.gz"),
        covhist = temp("genomes/alignments/{sample}_coverage_histogram.txt"),
        covstats = temp("genomes/alignments/{sample}_coverage.txt"),
        bincov = temp("genomes/alignments/{sample}_coverage_binned.txt")
    log:
        "logs/genomes/alignments/pilup_{sample}.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config["mem"],
        java_mem = int(config["mem"] * JAVA_MEM_FRACTION)
    shell:
        """pileup.sh in={input.sam} \
               threads={threads} \
               -Xmx{resources.java_mem}G \
               covstats={output.covstats} \
               hist={output.covhist} \
               basecov={output.basecov}\
               concise=t \
               bincov={output.bincov} 2> {log}"""



localrules: combine_coverages_MAGs,combine_bined_coverages_MAGs
rule combine_coverages_MAGs:
    input:
        covstats = expand("genomes/alignments/{sample}_coverage.txt",
            sample=SAMPLES),
        cluster_attribution = "genomes/clustering/contig2genome.tsv"
    output:
        "genomes/counts/median_contig_coverage.tsv",
        "genomes/counts/raw_counts_contigs.tsv",
        "genomes/counts/raw_counts_genomes.tsv",
    run:
        import pandas as pd
        from utils.parsers_bbmap import combine_coverages


        combined_cov,Counts_contigs = combine_coverages(input.covstats,SAMPLES)

        combined_cov.to_csv(output[0],sep='\t')
        Counts_contigs.to_csv(output[1],sep='\t')


        cluster_attribution = pd.read_csv(input.cluster_attribution,header=None,index_col=0,squeeze=True,sep='\t')

        Counts_genome= Counts_contigs.groupby(cluster_attribution,axis=1).sum().T
        Counts_genome.index.name='Sample'
        Counts_genome.to_csv(output[2],sep='\t')



rule combine_bined_coverages_MAGs:
    input:
        binned_coverage_files = expand("genomes/alignments/{sample}_coverage_binned.txt",
            sample=SAMPLES),
        cluster_attribution = "genomes/clustering/contig2genome.tsv"
    params:
        samples= SAMPLES
    output:
        binned_cov= "genomes/counts/binned_coverage.tsv.gz",
        median_abund = "genomes/counts/median_coverage_genomes.tsv"
    run:

        import pandas as pd
        import os

        def read_coverage_binned(covarage_binned_file):
            return pd.read_csv(covarage_binned_file,sep='\t',
                             skiprows=2,
                             index_col=[0,2],
                             usecols=[0,1,2],
                             squeeze=True)


        binCov={}
        for i, cov_file in enumerate(input.binned_coverage_files):

            sample= params.samples[i]

            binCov[sample] = read_coverage_binned(cov_file)

        binCov = pd.DataFrame(binCov)
        binCov.index.names=['Contig','Position']
        binCov.to_csv(output.binned_cov,sep='\t',compression='gzip')

        cluster_attribution = pd.read_csv(input.cluster_attribution,header=None,index_col=0,squeeze=True,sep='\t')

        Median_abund= binCov.groupby(cluster_attribution.loc[binCov.index.get_level_values(0)].values).median().T

        Median_abund.to_csv(output.median_abund,sep='\t')

# rule predict_genes_genomes:
#     input:
#         dir=genome_dir
#     output:
#         directory("genomes/annotations/genes")
#     conda:
#         "%s/prodigal.yaml" % CONDAENV
#     log:
#         "logs/genomes/prodigal.log"
#     shadow:
#         "shallow"
#     threads:
#         config.get("threads", 1)
#     script:
#         "predict_genes_of_genomes.py"



rule predict_genes_genomes:
    input:
        "genomes/genomes/{genome}.fasta"
    output:
        fna = "genomes/annotations/genes/{genome}.fna",
        faa = "genomes/annotations/genes/{genome}.faa",
        gff = "genomes/annotations/genes/{genome}.gff"
    conda:
        "%s/prodigal.yaml" % CONDAENV
    log:
        "logs/genomes/prodigal/{genome}.txt"
    threads:
        1
    resources:
        mem= config['simplejob_mem']
    shell:
        """
        prodigal -i {input} -o {output.gff} -d {output.fna} \
            -a {output.faa} -p meta -f gff 2> {log}
        """

def get_all_genes(wildcards,extension='.faa'):
    return expand("genomes/annotations/genes/{genome}{extension}",
           genome=get_genomes_(wildcards),extension=extension)

localrules: all_prodigal
rule all_prodigal:
    input:
        get_all_genes
    output:
        touch("genomes/annotations/genes/predicted")



## detached from DAG


## annotation

rule run_prokka_bins:
    input:
        "genomes/genomes/{genome}.fasta"
    output:
        expand("genomes/annotations/prokka/{{genome}}/{{genome}}.{extension}",
               extension= ["err","faa","ffn",
                           "fna","fsa","gff",
                           "tbl","tsv","txt"])
    log:
        "genomes/annotations/prokka/{genome}/{genome}.log"
    params:
        outdir = lambda wc, output: os.path.dirname(output[0]),
        kingdom = config.get("prokka_kingdom", PROKKA_KINGDOM)
    shadow:
        "shallow"
    conda:
        "%s/prokka.yaml" % CONDAENV
    threads:
        config.get("simplejob_threads", 1)
    shell:
        """prokka --outdir {params.outdir} \
               --force \
               --metagenome \
               --prefix {wildcards.genome} \
               --locustag {wildcards.genome} \
               --kingdom {params.kingdom} \
               --cpus {threads} \
               {input}
              """


def genome_all_prokka_input(wildcards):
    genome_dir = genome_dir
    return expand("genomes/annotations/prokka/{genome}/{genome}.tsv",
           genome=glob_wildcards(os.path.join(genome_dir, "{genome}.fasta")).genome)

rule all_prokka:
    input:
        genome_all_prokka_input
    output:
        "genomes/annotations/prokka.tsv"
    run:
        import pandas as pd
        out= pd.concat([pd.read_csv(file,index_col=0,sep='\t') for file in input],axis=0).sort_index()

        out.to_csv(output[0],sep='\t')
