import os
import re
import sys
from glob import glob
from snakemake.utils import report
import warnings


localrules: build_qc_report, combine_read_length_stats, \
            combine_read_counts



def get_ribosomal_rna_input(wildcards):
    data_type = config["data_type"]

    clean_reads = expand("{sample}/sequence_quality_control/{sample}_{step}_{fraction}.fastq.gz",step='clean', fraction=MULTIFILE_FRACTIONS,sample=wildcards.sample)
    rrna_reads = expand("{sample}/sequence_quality_control/contaminants/rRNA_{fraction}.fastq.gz",fraction=MULTIFILE_FRACTIONS,sample=wildcards.sample)

    if data_type == "metagenome" and 'rrna' in [c.lower() for c in config["contaminant_references"].keys()]:
        return {'clean_reads':clean_reads, 'rrna_reads':rrna_reads}
    else:
        return {'clean_reads':clean_reads}


# def get_finalize_qc_input(wildcards):
#     inputs = get_quality_controlled_reads(wildcards)
#     try:
#         # FIXME: 'decontaminated_seqs': ['{sample}/sequence_quality_control/contaminants/PhiX_R1.fastq.gz'...
#         # inputs["decontaminated_seqs"] = rules.decontamination.output.contaminants
#         inputs["decontaminated_stats"] = "{sample}/sequence_quality_control/{sample}_decontamination_reference_stats.txt".format(sample=wildcards.sample)
#     except AttributeError:
#         pass
#     return inputs


if SKIP_QC:
    PROCESSED_STEPS = ['QC']

    get_input_fastq=get_quality_controlled_reads

else:
    PROCESSED_STEPS = ['raw']



    def get_input_fastq(wildcards):

        if not config.get('interleaved_fastqs',False):
            Raw_Headers= ['Reads_raw_'+f for f in  RAW_INPUT_FRACTIONS]
        else:
            Raw_Headers= ['Reads_raw_se']

        #get file
        if not (wildcards.sample in sampleTable.index):
            return expand("Impossible/file/{sample}_{fraction}.fastq.gz",sample=wildcards.sample,fraction=RAW_INPUT_FRACTIONS)
            logger.debug(f"Searched for qc reads for inexisitng sample. wildcards: {wildcards}")
        else:

            return get_files_from_sampleTable(wildcards.sample,Raw_Headers)



rule initialize_qc:
    input:
        unpack(get_input_fastq)
    output:
        temp(expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                    fraction=RAW_INPUT_FRACTIONS, step=PROCESSED_STEPS[-1]))
    priority: 80
    params:
        inputs = lambda wc, input: io_params_for_tadpole(input,'in'),
        interleaved = lambda wc: "t" if config.get('interleaved_fastqs',False) else "f",
        outputs = lambda wc, output: io_params_for_tadpole(output,'out'),
        verifypaired = "t" if PAIRED_END else "f",
        extra = config['importqc_params']
    log:
        "{sample}/logs/QC/init.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("simplejob_threads", 1)
    resources:
        mem = config["simplejob_mem"],
        java_mem = int(int(config["simplejob_mem"] * JAVA_MEM_FRACTION))
#    group:
#        "qc"
    shell:
        """
        reformat.sh {params.inputs} \
            interleaved={params.interleaved} \
            {params.outputs} \
            {params.extra} \
            overwrite=true \
            verifypaired={params.verifypaired} \
            threads={threads} \
            -Xmx{resources.java_mem}G 2> {log}
        """


rule get_read_stats:
    # TODO: remove run block in favor of script or alternate cli
    # see http://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#snakefiles-external-scripts
    input:
        expand("{{sample}}/sequence_quality_control/{{sample}}_{{step}}_{fraction}.fastq.gz",
            fraction=RAW_INPUT_FRACTIONS)
    output:
        "{sample}/sequence_quality_control/read_stats/{step}.zip",
        read_counts = temp("{sample}/sequence_quality_control/read_stats/{step}_read_counts.tsv"),
    priority:
        30
    log:
        "{sample}/logs/QC/read_stats/{step}.log"
    # conda:
    #     "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("simplejob_threads", 1)
    resources:
        mem = config["simplejob_mem"],
        java_mem = int(config["simplejob_mem"] * JAVA_MEM_FRACTION)
    params:
        folder = lambda wc, output: os.path.splitext(output[0])[0],
        single_end_file = "{sample}/sequence_quality_control/{sample}_{step}_se.fastq.gz"
#    group:
#        "qc"
    run:
        import datetime
        import shutil
        import sys
        sys.stderr = open(log[0], 'w')

        timestamp = datetime.datetime.now().strftime('%Y-%m-%d-%X')

        def get_read_stats(fraction, params_in):

            subfolder = os.path.join(params.folder, fraction)
            tmp_file=os.path.join(subfolder,"read_stats.tmp")
            shell("""
                    mkdir -p {subfolder} 2> {log}

                    reformat.sh {params_in} \
                    bhist={subfolder}/base_hist.txt \
                    qhist={subfolder}/quality_by_pos.txt \
                    lhist={subfolder}/readlength.txt \
                    gchist={subfolder}/gc_hist.txt \
                    gcbins=auto \
                    bqhist={subfolder}/boxplot_quality.txt \
                    threads={threads} \
                    overwrite=true \
                    -Xmx{mem}G \
                    2> >(tee -a {log} {tmp_file} )
                 """.format(subfolder=subfolder, params_in=params_in, log=log,
                            threads=threads, mem=resources.java_mem,tmp_file=tmp_file))
            content = open(tmp_file).read()
            pos = content.find('Input:')
            if pos == -1:
                raise Exception("Didn't find read number in file:\n\n" + content)
            else:

                content[pos:].split()[1:4]
                        # Input:    123 reads   1234 bases
                n_reads, _, n_bases = content[pos:].split()[1:4]

                os.remove(tmp_file)
            return int(n_reads), int(n_bases)


        if PAIRED_END:
            n_reads_pe, n_bases_pe = get_read_stats('pe', "in1={0} in2={1}".format(*input))
            n_reads_pe = n_reads_pe / 2
            headers = ['Sample', 'Step', 'Total_Reads', 'Total_Bases',
                       'Reads_pe', 'Bases_pe', 'Reads_se', 'Bases_se',
                       'Timestamp']

            if os.path.exists(params.single_end_file):
                n_reads_se, n_bases_se = get_read_stats('se', "in=" + params.single_end_file)
            else:
                n_reads_se, n_bases_se = 0, 0

            values=[n_reads_pe + n_reads_se, n_bases_pe + n_bases_se,
                    n_reads_pe, n_bases_pe,
                    n_reads_se, n_bases_se]
        else:
            headers= ['Sample', 'Step', 'Total_Reads', 'Total_Bases', 'Reads',
                      'Bases', 'Timestamp']
            values = 2 * get_read_stats('', "in=" + input[0])

        with open(output.read_counts, 'w') as f:
            f.write('\t'.join(headers) + '\n')
            f.write('\t'.join([wildcards.sample, wildcards.step] + [str(v) for v in values] + [timestamp]) + '\n')

        shutil.make_archive(params.folder, 'zip', params.folder)
        shutil.rmtree(params.folder)

if not SKIP_QC:

    if config.get('deduplicate', True):
        PROCESSED_STEPS.append("deduplicated")
        rule deduplicate_reads:
            input:
                expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                    step=PROCESSED_STEPS[-2], fraction=RAW_INPUT_FRACTIONS)
            output:
                temp(expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                    fraction=RAW_INPUT_FRACTIONS, step=PROCESSED_STEPS[-1]))
            benchmark:
                "logs/benchmarks/QC/deduplicate/{sample}.txt"
            params:
                inputs = lambda wc, input: io_params_for_tadpole(input,'in'),
                outputs = lambda wc,output: io_params_for_tadpole(output,'out'),
                dupesubs = config.get('duplicates_allow_substitutions', DUPLICATES_ALLOW_SUBSTITUTIONS),
                only_optical = 't' if config.get('duplicates_only_optical', DUPLICATES_ONLY_OPTICAL) else 'f'
            log:
                "{sample}/logs/QC/deduplicate.log"
            conda:
                "%s/required_packages.yaml" % CONDAENV
            threads:
                config.get("threads", 1)
            resources:
                mem = config["mem"],
                java_mem = int(config["mem"] * JAVA_MEM_FRACTION)
        #    group:
        #        "qc"
            shell:
                """
                clumpify.sh \
                    {params.inputs} \
                    {params.outputs} \
                    overwrite=true\
                    dedupe=t \
                    dupesubs={params.dupesubs} \
                    optical={params.only_optical}\
                    threads={threads} \
                    pigz=t unpigz=t \
                    -Xmx{resources.java_mem}G 2> {log}
                """


    PROCESSED_STEPS.append("filtered")


    rule apply_quality_filter:
        input:
            reads=expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                fraction=RAW_INPUT_FRACTIONS, step=PROCESSED_STEPS[-2]),
            adapters=ancient(config["preprocess_adapters"])
        output:
            temp(expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                fraction=MULTIFILE_FRACTIONS, step=PROCESSED_STEPS[-1])),
            stats = "{sample}/logs/{sample}_quality_filtering_stats.txt"
        benchmark:
            "logs/benchmarks/QC/quality_filter/{sample}.txt"
        params:
            ref = "ref=%s" % config.get("preprocess_adapters") if config.get("preprocess_adapters") else "",
            mink = "" if not config.get("preprocess_adapters") else "mink=%d" % config.get("preprocess_adapter_min_k", PREPROCESS_ADAPTER_MIN_K),
            ktrim = "" if not config.get("preprocess_adapters") else "ktrim=%s" % config.get("preprocess_kmer_trim", PREPROCESS_KMER_TRIM),
            trimq = config.get("preprocess_minimum_base_quality", PREPROCESS_MINIMUM_BASE_QUALITY),
            hdist = "" if not config.get("preprocess_adapters") else "hdist=%d" % config.get("preprocess_allowable_kmer_mismatches", PREPROCESS_ALLOWABLE_KMER_MISMATCHES),
            k = "" if not config.get("preprocess_adapters") else "k=%d" % config.get("preprocess_reference_kmer_match_length", PREPROCESS_REFERENCE_KMER_MATCH_LENGTH),
            qtrim = config.get("qtrim", QTRIM),
            error_correction_pe = "t" if PAIRED_END and config.get('error_correction_overlapping_pairs',True) else "f",
            minlength = config.get("preprocess_minimum_passing_read_length", PREPROCESS_MINIMUM_PASSING_READ_LENGTH),
            minbasefrequency = config.get("preprocess_minimum_base_frequency", PREPROCESS_MINIMUM_BASE_FREQUENCY),
            # we require the user to reformat to R1 and R2, non-interleaved files
            interleaved = "f",
            maxns = config.get("preprocess_max_ns", PREPROCESS_MAX_NS),
            prealloc = config.get("preallocate_ram", PREALLOCATE_RAM),
            inputs = lambda wc, input: io_params_for_tadpole(input.reads),
            outputs = lambda wc, output:"out={0} out2={1} outs={2}".format(*output) if PAIRED_END else "out={0}".format(*output)
        log:
            "{sample}/logs/QC/quality_filter.log"
        conda:
            "%s/required_packages.yaml" % CONDAENV
        threads:
            config.get("threads", 1)
        resources:
            mem = config["mem"],
            java_mem = int(config["mem"] * JAVA_MEM_FRACTION)
    #    group:
    #        "qc"
        shell:
            """
            bbduk.sh {params.inputs} \
                {params.ref} \
                interleaved={params.interleaved} \
                {params.outputs} \
                stats={output.stats} \
                overwrite=true \
                qout=33 \
                trd=t \
                {params.hdist} \
                {params.k} \
                {params.ktrim} \
                {params.mink} \
                trimq={params.trimq} \
                qtrim={params.qtrim} \
                threads={threads} \
                minlength={params.minlength} \
                maxns={params.maxns} \
                minbasefrequency={params.minbasefrequency} \
                ecco={params.error_correction_pe} \
                prealloc={params.prealloc} \
                pigz=t unpigz=t \
                -Xmx{resources.java_mem}G 2> {log}
            """

    # if there are no references, decontamination will be skipped
    if len(config.get("contaminant_references", {}).keys()) > 0:
        PROCESSED_STEPS.append("clean")

        rule build_decontamination_db:
            input:
                ancient(config["contaminant_references"].values())
            output:
                "ref/genome/1/summary.txt"
            threads:
                config.get("threads", 1)
            resources:
                mem = config["mem"],
                java_mem = int(config["mem"] * JAVA_MEM_FRACTION)
            log:
                "logs/QC/build_decontamination_db.log"
            conda:
                "%s/required_packages.yaml" % CONDAENV
            params:
                k = config.get("contaminant_kmer_length", CONTAMINANT_KMER_LENGTH),
                refs_in = " ".join(["ref_%s=%s" % (n, fa) for n, fa in config["contaminant_references"].items()]),
            shell:
                """
                bbsplit.sh -Xmx{resources.java_mem}G {params.refs_in} \
                    threads={threads} k={params.k} local=t 2> {log}
                """


        rule run_decontamination:
            input:
                expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                    step=PROCESSED_STEPS[-2], fraction=MULTIFILE_FRACTIONS),
                db = "ref/genome/1/summary.txt"
            output:
                temp(expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                    fraction=MULTIFILE_FRACTIONS, step=PROCESSED_STEPS[-1])),
                contaminants = expand("{{sample}}/sequence_quality_control/contaminants/{db}_{fraction}.fastq.gz",
                        db=list(config["contaminant_references"].keys()),
                        fraction=MULTIFILE_FRACTIONS),
                stats = "{sample}/sequence_quality_control/{sample}_decontamination_reference_stats.txt"
            benchmark:
                "logs/benchmarks/QC/decontamination/{sample}.txt"
            params:
                contaminant_folder = lambda wc, output: os.path.dirname(output.contaminants[0]),
                maxindel = config.get("contaminant_max_indel", CONTAMINANT_MAX_INDEL),
                minratio = config.get("contaminant_min_ratio", CONTAMINANT_MIN_RATIO),
                minhits = config.get("contaminant_minimum_hits", CONTAMINANT_MINIMUM_HITS),
                ambiguous = config.get("contaminant_ambiguous", CONTAMINANT_AMBIGUOUS),
                k = config.get("contaminant_kmer_length", CONTAMINANT_KMER_LENGTH),
                paired= "true" if PAIRED_END else "false",
                input_single= lambda wc,input: input[2] if PAIRED_END else input[0],
                output_single= lambda wc,output: output[2] if PAIRED_END else output[0]
            log:
                "{sample}/logs/QC/decontamination.log"
            conda:
                "%s/required_packages.yaml" % CONDAENV
            threads:
                config.get("threads", 1)
            resources:
                mem = config["mem"],
                java_mem = int(config["mem"] * JAVA_MEM_FRACTION)
        #    group:
        #        "qc"
            shell:
                """
                if [ "{params.paired}" = true ] ; then
                    bbsplit.sh in1={input[0]} in2={input[1]} \
                        outu1={output[0]} outu2={output[1]} \
                        basename="{params.contaminant_folder}/%_R#.fastq.gz" \
                        maxindel={params.maxindel} minratio={params.minratio} \
                        minhits={params.minhits} ambiguous={params.ambiguous} refstats={output.stats}\
                        threads={threads} k={params.k} local=t \
                        pigz=t unpigz=t ziplevel=9 \
                        -Xmx{resources.java_mem}G 2> {log}
                fi

                bbsplit.sh in={params.input_single}  \
                    outu={params.output_single} \
                    basename="{params.contaminant_folder}/%_se.fastq.gz" \
                    maxindel={params.maxindel} minratio={params.minratio} \
                    minhits={params.minhits} ambiguous={params.ambiguous} refstats={output.stats} append \
                    interleaved=f threads={threads} k={params.k} local=t \
                    pigz=t unpigz=t ziplevel=9 \
                    -Xmx{resources.java_mem}G 2>> {log}
                """


    PROCESSED_STEPS.append("QC")

    localrules: qcreads
    rule qcreads:
        input:
            unpack(get_ribosomal_rna_input),
            sample_table_file= ancient('samples.tsv')
        output:
            expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                step=PROCESSED_STEPS[-1], fraction=MULTIFILE_FRACTIONS)
        threads:
            1
    #    group:
    #        "qc"
        run:
            import shutil
            import pandas as pd
            for i in range(len(MULTIFILE_FRACTIONS)):

                with open(output[i], 'wb') as outFile:
                    with open(input.clean_reads[i], 'rb') as infile1:
                        shutil.copyfileobj(infile1, outFile)

                        if hasattr(input, 'rrna_reads'):
                            with open(input.rrna_reads[i], 'rb') as infile2:
                                shutil.copyfileobj(infile2, outFile)

            # append to sample table
            sample_table = load_sample_table(input.sample_table_file)

            sample_table.loc[wildcards.sample,
                             [ f'Reads_QC_{fraction}' for fraction in MULTIFILE_FRACTIONS]
                             ] = output

            sample_table.to_csv(input.sample_table_file,sep='\t')




#### STATS


if PAIRED_END:
    rule calculate_insert_size:
        input:
            lambda wildcards: input_paired_only(get_quality_controlled_reads(wildcards))
        output:
            ihist = "{sample}/sequence_quality_control/read_stats/QC_insert_size_hist.txt",
            read_length= "{sample}/sequence_quality_control/read_stats/QC_read_length_hist.txt"
        threads:
            config.get("simplejob_threads", 1)
        resources:
            mem = config["mem"],
            java_mem = int(config["mem"] * JAVA_MEM_FRACTION)
        conda:
            "%s/required_packages.yaml" % CONDAENV
        log:
            "{sample}/logs/QC/stats/calculate_insert_size.log"
        params:
            kmer = config.get("merging_k", MERGING_K),
            extend2 = config.get("merging_extend2", MERGING_EXTEND2),
            flags = 'loose ecct',
            minprob = config.get("bbmerge_minprob", "0.8"),
            inputs = lambda wc,input: io_params_for_tadpole(input)
    #    group:
    #        "qc"
        shell:
            """
            bbmerge.sh -Xmx{resources.java_mem}G threads={threads} \
                {params.inputs} \
                {params.flags} k={params.kmer} \
                extend2={params.extend2} \
                ihist={output.ihist} merge=f \
                mininsert0=35 minoverlap0=8 \
                prealloc=t prefilter=t \
                minprob={params.minprob} 2> {log}

            readlength.sh {params.inputs} out={output.read_length} 2> {log}
            """


else:
    rule calculate_read_length_hist:
        input:
            get_quality_controlled_reads
        output:
            read_length= "{sample}/sequence_quality_control/read_stats/QC_read_length_hist.txt",
        params:
            kmer = config.get("merging_k", MERGING_K),
        threads:
            config.get("simplejob_threads", 1)
        resources:
            mem = config["simplejob_mem"],
            java_mem = int(config["simplejob_mem"] * JAVA_MEM_FRACTION)
        conda:
            "%s/required_packages.yaml" % CONDAENV
        log:
            "{sample}/logs/QC/stats/calculate_read_length.log"
    #    group:
    #        "qc"
        shell:
            """
            readlength.sh in={input[0]} out={output.read_length} 2> {log}
            """


rule combine_read_length_stats:
    input:
        expand("{sample}/sequence_quality_control/read_stats/QC_read_length_hist.txt", sample=SAMPLES),
    output:
        'stats/read_length_stats.tsv'
    run:
        import pandas as pd
        import os
        from utils.parsers_bbmap import parse_comments

        stats = pd.DataFrame()

        for length_file in input:
            sample = length_file.split(os.path.sep)[0]
            data = parse_comments(length_file)
            data = pd.Series(data)[['Reads', 'Bases', 'Max', 'Min', 'Avg', 'Median', 'Mode', 'Std_Dev']]
            stats[sample] = data

        stats.to_csv(output[0], sep='\t')


# rule combine_cardinality:
#     input:
#         expand("{sample}/sequence_quality_control/read_stats/QC_cardinality.txt",sample=SAMPLES),
#     output:
#         'stats/cardinality.tsv'
#     run:
#         import pandas as pd
#         import os

#         stats= pd.Series()

#         for file in input:
#             sample= file.split(os.path.sep)[0]
#             with open(file) as f:
#                 cardinality= int(f.read().strip())

#             stats.loc[sample]=cardinality

#         stats.to_csv(output[0],sep='\t')


if PAIRED_END:
    localrules: combine_insert_stats
    rule combine_insert_stats:
        input:
            expand("{sample}/sequence_quality_control/read_stats/QC_insert_size_hist.txt",sample=SAMPLES),
        output:
            'stats/insert_stats.tsv'
        run:
            import pandas as pd
            import os
            from utils.parsers_bbmap import parse_comments

            stats = pd.DataFrame()

            for insert_file in input:
                sample = insert_file.split(os.path.sep)[0]
                data = parse_comments(insert_file)
                data = pd.Series(data)[['Mean', 'Median', 'Mode', 'STDev', 'PercentOfPairs']]
                stats[sample] = data

            stats.T.to_csv(output[0], sep='\t')

localrules: combine_read_counts, write_read_counts
rule write_read_counts:
    input:
        read_count_files = expand("{{sample}}/sequence_quality_control/read_stats/{step}_read_counts.tsv", step=PROCESSED_STEPS)
    output:
        read_stats = "{sample}/sequence_quality_control/read_stats/read_counts.tsv"
    run:
        import pandas as pd

        all_read_counts = pd.DataFrame()
        for read_stats_file in input.read_count_files:
            d = pd.read_csv(read_stats_file, index_col=[0, 1],sep='\t')
            all_read_counts = all_read_counts.append(d)
        all_read_counts.to_csv(output.read_stats, sep='\t')


rule combine_read_counts:
    input:
        expand("{sample}/sequence_quality_control/read_stats/read_counts.tsv",sample=SAMPLES)
    output:
        "stats/read_counts.tsv"
    run:
        import pandas as pd

        stats = pd.DataFrame()

        for f in input:
            d = pd.read_csv(f, index_col=[0, 1],sep='\t')
            stats = stats.append(d)

        stats.to_csv(output[0], sep='\t')


rule finalize_sample_qc:
    input:
        qcreads= get_quality_controlled_reads,
        #quality_filtering_stats = "{sample}/logs/{sample}_quality_filtering_stats.txt",
        reads_stats_zip = expand("{{sample}}/sequence_quality_control/read_stats/{step}.zip", step=PROCESSED_STEPS),
        read_length_hist = "{sample}/sequence_quality_control/read_stats/QC_read_length_hist.txt",
    output:
        touch("{sample}/sequence_quality_control/finished_QC"),






def zipfile_input(wildcards):
    """
        returns dict of list with all the zipfiles of stats.
        when QC is skipped only the QC stats are returned.
    """
    output=dict(zipfiles_QC = expand('{sample}/sequence_quality_control/read_stats/QC.zip', sample=SAMPLES))

    if not SKIP_QC:
        output['zipfiles_raw'] = expand('{sample}/sequence_quality_control/read_stats/raw.zip', sample=SAMPLES)

    return output

rule build_qc_report:
    input:
        unpack(zipfile_input),
        read_counts = "stats/read_counts.tsv",
        read_length_stats = ['stats/insert_stats.tsv', 'stats/read_length_stats.tsv'] if PAIRED_END else 'stats/read_length_stats.tsv',
    output:
        report = "reports/QC_report.html"
    log:
        "logs/QC/report.log"
    params:
        min_quality = config["preprocess_minimum_base_quality"],
    conda:
        "%s/report.yaml" % CONDAENV
    script:
         "../report/qc_report.py"
