import os



if config['genecatalog']['source']=='contigs':

    localrules: concat_genes
    rule concat_genes:
        input:
            faa= expand("{sample}/annotation/predicted_genes/{sample}.faa", sample=SAMPLES),
            fna= expand("{sample}/annotation/predicted_genes/{sample}.fna", sample=SAMPLES)
        output:
            faa=  temp("Genecatalog/all_genes_unfiltered.faa"),
            fna = temp("Genecatalog/all_genes_unfiltered.fna"),
        run:
            from utils.io import cat_files
            cat_files(input.faa,output.faa)
            cat_files(input.fna,output.fna)



else:

    localrules: concat_genes
    rule concat_genes:
        input:
            "genomes/annotations/orf2genome.tsv",
            faa= lambda wc: get_all_genes(wc,".faa"),
            fna= lambda wc: get_all_genes(wc,".fna")
        output:
            faa=  temp("Genecatalog/all_genes_unfiltered.faa"),
            fna = temp("Genecatalog/all_genes_unfiltered.fna"),
        run:
            from utils.io import cat_files
            cat_files(input.faa,output.faa)
            cat_files(input.fna,output.fna)



localrules: filter_genes
rule filter_genes:
    input:
        fna="Genecatalog/all_genes_unfiltered.fna",
        faa="Genecatalog/all_genes_unfiltered.faa"
    output:
        fna= "Genecatalog/all_genes/predicted_genes.fna",
        faa= "Genecatalog/all_genes/predicted_genes.faa",
    threads:
        1
    params:
        min_length=config['genecatalog']['minlength']
    run:
        from Bio import SeqIO
        faa = SeqIO.parse(input.faa,'fasta')
        fna = SeqIO.parse(input.fna,'fasta')

        with open(output.faa,'w') as out_faa, open(output.fna,'w') as out_fna:

            for gene in fna:
                protein = next(faa)

                if len(gene) >= params.min_length:
                    SeqIO.write(gene,out_fna,'fasta')
                    SeqIO.write(protein,out_faa,'fasta')


if (config['genecatalog']['clustermethod']=='linclust') or (config['genecatalog']['clustermethod']=='mmseqs'):

    rule cluster_genes:
        input:
            faa= "Genecatalog/all_genes/predicted_genes.faa"
        output:
            db=temp(directory("Genecatalog/all_genes/predicted_genes")),
            clusterdb = temp(directory("Genecatalog/clustering/mmseqs"))
        conda:
            "%s/mmseqs.yaml" % CONDAENV
        log:
            "logs/Genecatalog/clustering/cluster_proteins.log"
        threads:
            config.get("threads", 1)
        params:
            tmpdir= os.path.join(config['tmpdir'],"mmseqs"),
            clustermethod = 'linclust' if config['genecatalog']['clustermethod']=='linclust' else 'cluster',
            coverage=config['genecatalog']['coverage'], #0.8,
            minid=config['genecatalog']['minid'], # 0.00
            extra=config['genecatalog']['extra'],
            clusterdb= lambda wc, output: os.path.join(output.clusterdb,'clusterdb'),
            db=lambda wc, output: os.path.join(output.db,'inputdb')
        shell:
            """
                mkdir -p {params.tmpdir} {output} 2>> {log}
                mmseqs createdb {input.faa} {params.db} &> {log}

                mmseqs {params.clustermethod} -c {params.coverage} \
                --min-seq-id {params.minid} {params.extra} \
                --threads {threads} {params.db} {params.clusterdb} {params.tmpdir}  &>>  {log}

                rm -fr  {params.tmpdir} 2>> {log}
            """


    rule get_rep_proteins:
        input:
            db= rules.cluster_genes.output.db,
            clusterdb = rules.cluster_genes.output.clusterdb,
        output:
            cluster_attribution = temp("Genecatalog/orf2gene_oldnames.tsv"),
            rep_seqs_db = temp(directory("Genecatalog/protein_catalog")),
            rep_seqs = temp("Genecatalog/representatives_of_clusters.fasta")
        conda:
            "%s/mmseqs.yaml" % CONDAENV
        log:
            "logs/Genecatalog/clustering/get_rep_proteins.log"
        threads:
            config.get("threads", 1)
        params:
            clusterdb= lambda wc, input: os.path.join(input.clusterdb,'clusterdb'),
            db=lambda wc, input: os.path.join(input.db,'inputdb'),
            rep_seqs_db=lambda wc, output: os.path.join(output.rep_seqs_db,'db')
        shell:
            """
            mmseqs createtsv {params.db} {params.db} {params.clusterdb} {output.cluster_attribution}  &> {log}

            mkdir {output.rep_seqs_db} 2>> {log}

            mmseqs result2repseq {params.db} {params.clusterdb} {params.rep_seqs_db}  &>> {log}

            mmseqs result2flat {params.db} {params.db} {params.rep_seqs_db} {output.rep_seqs}  &>> {log}

            """


    localrules: rename_protein_catalog
    rule rename_protein_catalog:
        input:
            cluster_attribution = "Genecatalog/orf2gene_oldnames.tsv",
        output:
            cluster_attribution = "Genecatalog/clustering/orf2gene.tsv.gz",
        run:
            import pandas as pd
            # CLuterID    GeneID    empty third column
            orf2gene= pd.read_csv(input.cluster_attribution,index_col=1, header=None,sep='\t')

            protein_clusters_old_names= orf2gene[0].unique()

            map_names = dict(zip(protein_clusters_old_names,
                                 utils.gen_names_for_range(len(protein_clusters_old_names),'Gene')))

            orf2gene['Gene'] = orf2gene[0].map(map_names)
            orf2gene.index.name='ORF'
            orf2gene['Gene'].to_csv(output.cluster_attribution,sep='\t',header=True,compression='gzip')






elif config['genecatalog']['clustermethod']=='cd-hit-est':

# cluster genes


    rule cluster_genes:
        input:
            fna_dir="Genecatalog/all_genes/predicted_genes.fna",
        output:
            temp("Genecatalog/representatives_of_clusters.fasta"),
            temp("Genecatalog/gene_catalog_oldnames.clstr")
        conda:
            "%s/cd-hit.yaml" % CONDAENV
        log:
            "logs/Genecatalog/cluster_genes.log"
        threads:
            config.get("threads", 1)
        resources:
            mem= config["mem"]
        params:
            coverage=config['genecatalog']['coverage'],
            identity=config['genecatalog']['minid'],
            extra= config['genecatalog']['extra'],
            prefix= lambda wc,output: os.path.splitext(output[1])[0],
        shell:
            """
                cd-hit-est -i {input} -T {threads} \
                -M {resources.mem}000 -o {params.prefix} \
                -c {params.identity} -n 9  -d 0 {params.extra} \
                -aS {params.coverage} -aL {params.coverage} &> {log}

                mv {params.prefix} {output[0]} 2>> {log}
            """





    localrules:  parse_clstr_files, rename_gene_clusters

    def parse_cd_hit_file(clstr_file):
        """

        >Cluster 0
        0	342nt, >S1_83_1... *
        1	342nt, >S2_82_1... at +/100.00%
        >Cluster 1
        0	339nt, >S1_61_1... *
        1	339nt, >S2_59_1... at +/100.00%


        """
        import numpy as np
        def parse_line(line):
            _, length, name, identity = line.strip().replace('...','\t').replace(', ','\t').split('\t')

            length= int(length.replace('nt',''))
            name=name[1:]
            if '*' in identity:
                identity= np.nan
            else:
                identity= float(identity[identity.rfind('/')+1:identity.rfind('%')])

            return name,length, identity

        Clusters= []
        with open(clstr_file) as f:
            for line in f:
                if line[0]=='>': #new cluster
                    cluster= dict(elements=[],representative=None)
                    Clusters.append(cluster)
                else:
                    name,length, identity =parse_line(line)
                    cluster['elements'].append((name,length, identity))
                    if np.isnan(identity):
                        cluster['representative']= name
        return Clusters


    def write_cd_hit_clusters(Clusters,file_handle):
            for cluster in Clusters:
                for element in cluster['elements']:
                    file_handle.write(f"{element[0]}\t{element[1]}\t{element[2]}\t{cluster['representative']}\n")



    rule parse_clstr_files:
        input:
            clustered_dir= "Genecatalog/gene_catalog_oldnames.clstr"
        output:
            temp("Genecatalog/clustering/orf2gene_oldnames.tsv")
        run:
            with open(output[0],'w') as fout:
                fout.write(f"ORF\tLength\tIdentity\tGene\n")
                Clusters=parse_cd_hit_file(input[0])
                write_cd_hit_clusters(Clusters,fout)



    rule rename_gene_clusters:
        input:
            orf2gene = "Genecatalog/clustering/orf2gene_oldnames.tsv",
        output:
            orf2gene = "Genecatalog/clustering/orf2gene.tsv.gz",
        run:
            import pandas as pd
            from Bio import SeqIO

            orf2gene= pd.read_csv(input.orf2gene,index_col=0,sep='\t')

            # rename gene repr to Gene0000XX

            gene_clusters_old_names= orf2gene['Gene'].unique()

            map_names = dict(zip(gene_clusters_old_names,
                                 utils.gen_names_for_range(len(gene_clusters_old_names),'Gene')))

            orf2gene['Gene'] = orf2gene['Gene'].map(map_names)
            orf2gene.to_csv(output.orf2gene,sep='\t',header=True,compression='gzip')

else:
    raise Exception("Didn't understood the genecatalog clustermethod: {}".format(config['genecatalog']['clustermethod']))

localrules: rename_gene_catalog
rule rename_gene_catalog:
    input:
        fna = "Genecatalog/all_genes/predicted_genes.fna",
        faa= "Genecatalog/all_genes/predicted_genes.faa",
        orf2gene = "Genecatalog/clustering/orf2gene.tsv.gz",
        representatives= "Genecatalog/representatives_of_clusters.fasta"
    output:
        fna= "Genecatalog/gene_catalog.fna",
        faa= "Genecatalog/gene_catalog.faa",
    run:
        import pandas as pd
        from Bio import SeqIO

        representatives= []
        with open(input.representatives) as fasta:
            for line in fasta:
                if line[0]=='>': representatives.append(line[1:].split()[0])

        map_names= pd.read_csv(input.orf2gene,index_col=0,sep='\t').loc[representatives,'Gene']

        # rename fna
        faa_parser = SeqIO.parse(input.faa,'fasta')
        fna_parser = SeqIO.parse(input.fna,'fasta')

        with open(output.fna,'w') as fna, open(output.faa,'w') as faa :
            for gene in fna_parser:
                protein = next(faa_parser)
                if gene.name in map_names.index:
                    gene.id = map_names[gene.name]
                    protein.id = map_names[protein.name]

                    SeqIO.write(gene,fna,'fasta')
                    SeqIO.write(protein,faa,'fasta')


rule align_reads_to_Genecatalog:
    input:
        reads=get_quality_controlled_reads,
        fasta = "Genecatalog/gene_catalog.fna",
    output:
        sam = temp("Genecatalog/alignments/{sample}.sam")
    params:
        input = lambda wc, input : input_params_for_bbwrap( input.reads),
        maxsites = 4,
        ambiguous = 'all',
        minid = config['genecatalog']['minid'],
        maxindel = 1 # default 16000 good for genome deletions but not necessarily for alignment to contigs
    shadow:
        "shallow"
    log:
        "logs/Genecatalog/alignment/{sample}_map.log"
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
            local=t \
            ref={input.fasta} \
            {params.input} \
            trimreaddescriptions=t \
            out={output.sam} \
            threads={threads} \
            minid={params.minid} \
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


rule pileup_Genecatalog:
    input:
        sam = "Genecatalog/alignments/{sample}.sam",
        bam = "Genecatalog/alignments/{sample}.bam"
    output:
        covstats = temp("Genecatalog/alignments/{sample}_coverage.tsv"),
        basecov = temp("Genecatalog/alignments/{sample}_base_coverage.txt.gz"),
    params:
        pileup_secondary = 't' # a read maay map to different genes
    log:
        "logs/Genecatalog/alignment/{sample}_pileup.log"
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
               basecov={output.basecov} \
               secondary={params.pileup_secondary} \
                2> {log}
        """

localrules: combine_gene_coverages
rule combine_gene_coverages:
    input:
        covstats = expand("Genecatalog/alignments/{sample}_coverage.tsv",
            sample=SAMPLES)
    output:
        "Genecatalog/counts/median_coverage.tsv.gz",
        "Genecatalog/counts/Nmapped_reads.tsv.gz",
    run:

        import pandas as pd
        import os

        combined_cov={}
        combined_N_reads={}
        for cov_file in input:

            sample= os.path.split(cov_file)[-1].split('_')[0]
            data= pd.read_csv(cov_file,index_col=0,sep='\t')
            data.loc[data.Median_fold<0,'Median_fold']=0
            combined_cov[sample]= data.Median_fold
            combined_N_reads[sample] = data.Plus_reads+data.Minus_reads

        pd.DataFrame(combined_cov).to_csv(output[0],sep='\t',compression='gzip')
        pd.DataFrame(combined_N_reads).to_csv(output[1],sep='\t',compression='gzip')


###########
## EGG NOG
##########

# # this rule specifies the more general eggNOG rules

# output with wildcards "{folder}/{prefix}.emapper.tsv"


rule eggNOG_homology_search:
    input:
        eggnog_db_files=get_eggnog_db_file(),
        faa = "{folder}/{prefix}.faa",
    output:
        temp("{folder}/{prefix}.emapper.seed_orthologs"),
    params:
        data_dir = EGGNOG_DIR,
        prefix = "{folder}/{prefix}"
    resources:
        mem = config["mem"]
    threads:
        config["threads"]
    shadow:
        "minimal"
    conda:
        "%s/eggNOG.yaml" % CONDAENV
    log:
        "{folder}/logs/{prefix}/eggNOG_homology_search_diamond.log"
    shell:
        """
        emapper.py -m diamond --no_annot --no_file_comments \
            --data_dir {params.data_dir} --cpu {threads} -i {input.faa} \
            -o {params.prefix} --override 2> {log}
        """



rule eggNOG_annotation:
    input:
        eggnog_db_files=get_eggnog_db_file(),
        seed = rules.eggNOG_homology_search.output
    output:
        temp("{folder}/{prefix}.emapper.annotations")
    params:
        data_dir = config['virtual_disk'] if config['eggNOG_use_virtual_disk'] else EGGNOG_DIR,
        prefix = "{folder}/{prefix}",
        copyto_shm="t" if config['eggNOG_use_virtual_disk'] else 'f'
    threads:
        config.get("threads", 1)
    resources:
        mem=  2*config["simplejob_mem"] + (37 if config['eggNOG_use_virtual_disk'] else 0)
    shadow:
        "minimal"
    conda:
        "%s/eggNOG.yaml" % CONDAENV
    log:
        "{folder}/logs/{prefix}/eggNOG_annotate_hits_table.log"
    shell:
        """

            if [ {params.copyto_shm} == "t" ] ;
            then
                cp {EGGNOG_DIR}/eggnog.db {params.data_dir}/eggnog.db 2> {log}
            fi

            emapper.py --annotate_hits_table {input.seed} --no_file_comments \
              --override -o {params.prefix} --cpu {threads} --data_dir {params.data_dir} 2>> {log}
        """





#
# localrules: get_Genecatalog_annotations
# rule get_Genecatalog_annotations:
#     input:
#         Genecatalog= 'Genecatalog/gene_catalog.fna".fna',
#         eggNOG= expand('{sample}/annotation/eggNOG.tsv',sample=SAMPLES),
#         refseq= expand('{sample}/annotation/refseq/{sample}_tax_assignments.tsv',sample=SAMPLES),
#         scg= expand("Genecatalog/annotation/single_copy_genes_{domain}.tsv",domain=['bacteria','archaea'])
#     output:
#         annotations= "Genecatalog/annotations.tsv",
#     run:
#         import pandas as pd
#
#         gene_ids=[]
#         with open(input.Genecatalog) as fasta_file:
#             for line in fasta_file:
#                 if line[0]=='>':
#                     gene_ids.append(line[1:].strip().split()[0])
#
#         eggNOG=pd.DataFrame()
#         for annotation_file in input.eggNOG:
#             eggNOG=eggNOG.append(pd.read_csv(annotation_file, index_col=0,sep='\t'))
#
#         refseq=pd.DataFrame()
#         for annotation_file in input.refseq:
#             refseq=refseq.append(pd.read_csv(annotation_file, index_col=1,sep='\t'))
#
#         scg=pd.DataFrame()
#         for annotation_file in input.scg:
#             d= pd.read_csv(annotation_file, index_col=0,header=None,sep='\t')
#             d.columns = 'scg_'+ os.path.splitext(annotation_file)[0].split('_')[-1] # bacteria or archaea
#             scg=scg.append(d)
#
#
#         annotations= refseq.join(eggNOG).join(scg).loc[gene_ids]
#         annotations.to_csv(output.annotations,sep='\t')



rule predict_single_copy_genes:
    input:
        "Genecatalog/gene_catalog.faa"
    output:
        "Genecatalog/annotation/single_copy_genes_{domain}.tsv",
    params:
        script_dir = os.path.dirname(os.path.abspath(workflow.snakefile)),
        key = lambda wc: wc.domain[:3] #bac for bacteria, #arc for archaea
    conda:
        "%s/DASTool.yaml" % CONDAENV # needs pearl
    log:
        "logs/Genecatalog/annotation/predict_single_copy_genes_{domain}.log"
    shadow:
        "shallow"
    threads:
        config['threads']
    shell:
        " DIR=$(dirname $(readlink -f $(which DAS_Tool))) "
        ";"
        " ruby {params.script_dir}/rules/scg_blank_diamond.rb diamond"
        " {input} "
        " $DIR\/db/{params.key}.all.faa "
        " $DIR\/db/{params.key}.scg.faa "
        " $DIR\/db/{params.key}.scg.lookup "
        " {threads} "
        " 2> {log} "
        " ; "
        " mv {input[0]}.scg {output}"




localrules: gene_subsets,combine_egg_nogg_annotations
checkpoint gene_subsets:
    input:
        "Genecatalog/gene_catalog.faa"
    output:
        directory("Genecatalog/subsets/genes")
    params:
        subset_size=config['genecatalog']['SubsetSize'],
    run:
        from utils import fasta
        fasta.split(input[0],params.subset_size,output[0],simplify_headers=True)


def combine_genecatalog_annotations_input(wildcards):
    dir_for_subsets = checkpoints.gene_subsets.get(**wildcards).output[0]
    Subset_names,= glob_wildcards(os.path.join(dir_for_subsets, "{subset}.faa"))
    return expand("Genecatalog/subsets/genes/{subset}.emapper.annotations",
                  subset=Subset_names)

rule combine_egg_nogg_annotations:
    input:
        combine_genecatalog_annotations_input
    output:
        "Genecatalog/annotations/eggNog.tsv.gz"
    run:

        import pandas as pd

        # read input files one after the other
        for i,annotation_table in enumerate(input):
            D = pd.read_csv(annotation_table, header=None,sep='\t')
            # Add headers, to verify size
            D.columns = EGGNOG_HEADER
            # appedn to output file, header only the first time
            D.to_csv(output[0],sep="\t",index=False,header= (i==0),compression='gzip',mode='a')











rule gene2genome:
    input:
        contigs2bins= "genomes/clustering/all_contigs2bins.tsv.gz",
        contigs2mags= "genomes/clustering/contig2genome.tsv",
        old2newID= "genomes/clustering/old2newID.tsv",
        orf2gene= "Genecatalog/clustering/orf2gene.tsv.gz"
    params:
        remaned_contigs= config['rename_mags_contigs'] & (config['genecatalog']['source']=='contigs')
    output:
        "genomes/annotations/gene2genome.tsv.gz"
    run:
        import pandas as pd

        if params.remaned_contigs:

            contigs2bins= pd.read_csv(input.contigs2bins,
                                       index_col=0,squeeze=False,sep='\t',header=None)

            contigs2bins.columns=['Bin']
            old2newID = pd.read_csv(input.old2newID,
                                       index_col=0,squeeze=True,sep='\t')

            contigs2genome=contigs2bins.join(old2newID,on='Bin').dropna().drop('Bin',axis=1)
        else:
            contigs2genome= pd.read_csv(input.contigs2mags,
                                       index_col=0,squeeze=False,sep='\t',header=None)
            contigs2genome.columns=['MAG']


        orf2gene = pd.read_csv(input.orf2gene,
                                   index_col=0,squeeze=False,sep='\t',header=0)

        orf2gene['Contig']= orf2gene.index.map(lambda s: '_'.join(s.split('_')[:-1]))
        orf2gene=orf2gene.join(contigs2genome,on='Contig')
        orf2gene= orf2gene.dropna(axis=0)

        gene2genome= orf2gene.groupby(['Gene','MAG']).size()
        gene2genome.name='Ncopies'

        gene2genome.to_csv(output[0],sep='\t',header=True,compression='gzip')



# after combination need to add eggNOG headerself.
#"{folder}/{prefix}_eggNOG.tsv"
