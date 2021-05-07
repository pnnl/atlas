wildcard_constraints:
    sample="[A-z0-9_-]"

import pandas as pd
ADDITIONAL_SAMPLEFILE_HEADERS=[]


def validate_sample_table(sampleTable):

    Expected_Headers =['BinGroup'] + ADDITIONAL_SAMPLEFILE_HEADERS
    for h in Expected_Headers:
        if not (h in sampleTable.columns):
         logging.error(f"expect '{h}' to be found in samples.tsv")
         exit(1)

    if not sampleTable.index.is_unique:
        duplicated_samples=', '.join(D.index.duplicated())
        logging.error( f"Expect Samples to be unique. Found {duplicated_samples} more than once")
        exit(1)


def load_sample_table(sample_table='samples.tsv'):

    sampleTable = pd.read_csv(sample_table,index_col=0,sep='\t')
    validate_sample_table(sampleTable)
    return sampleTable


sampleTable= load_sample_table(config.get('sample_table','samples.tsv'))

SAMPLES = sampleTable.index.values
SKIP_QC=False
#GROUPS = sampleTable.BinGroup.unique()
def get_alls_samples_of_group(wildcards):
    group_of_sample= sampleTable.loc[wildcards.sample,'BinGroup']

    return list(sampleTable.loc[ sampleTable.BinGroup==group_of_sample].index)



PAIRED_END = sampleTable.columns.str.contains('R2').any() or config.get('interleaved_fastqs',False)

colum_headers_QC= sampleTable.columns[sampleTable.columns.str.startswith("Reads_QC_")]
if len(colum_headers_QC)>=1:
    MULTIFILE_FRACTIONS= list(colum_headers_QC.str.replace('Reads_QC_',''))

    if (len(MULTIFILE_FRACTIONS)==1 ) and config.get('interleaved_fastqs',False):
        MULTIFILE_FRACTIONS=['R1','R2']

else:
    MULTIFILE_FRACTIONS = ['R1', 'R2', 'se'] if PAIRED_END else ['se']

colum_headers_raw= sampleTable.columns[sampleTable.columns.str.startswith("Reads_raw_")]
if len(colum_headers_raw) ==0:
    SKIP_QC=True

    logger.info("Didn't find raw reads in sampleTable - skip QC")
    RAW_INPUT_FRACTIONS = MULTIFILE_FRACTIONS
else:
    RAW_INPUT_FRACTIONS = ['R1', 'R2'] if PAIRED_END else ['se']


if (len(colum_headers_raw) ==0) and (len(colum_headers_QC) ==0):

    raise IOError("Either raw reas or QC reads need to be in the sample table. "
                  "I din't find any collums with 'Reads_raw_<fraction>' or 'Reads_QC_<fraction>'  "
                  )


class FileNotInSampleTableException(Exception):
    """
        Exception with sampleTable
    """
    def __init__(self, message):
        super(FileNotInSampleTableException, self).__init__(message)


def get_files_from_sampleTable(sample,Headers):
    """
        Function that gets some filenames form the sampleTable for a given sample and Headers.
        It checks various possibilities for errors and throws either a
        FileNotInSampleTableException or a IOError, when something went really wrong.
    """

    if not (sample in sampleTable.index):
        raise FileNotInSampleTableException(f"Sample name {sample} is not in sampleTable")


    Error_details=f"\nsample: {sample}\nFiles: {Headers}"

    if type(Headers) == str: Headers= [Headers]

    NheadersFound= sampleTable.columns.isin(Headers).sum()

    if  NheadersFound==0 :
        raise FileNotInSampleTableException(f"None of the Files ar in sampleTable, they should be added to the sampleTable later in the workflow"+Error_details)
    elif NheadersFound < len(Headers):
        raise IOError(f"Not all of the Headers are in sampleTable, found only {NheadersFound}, something went wrong."+Error_details)

    files= sampleTable.loc[sample,Headers]

    if files.isnull().all():
        raise FileNotInSampleTableException("The following files were not available for this sample in the SampleTable"+Error_details)

    elif files.isnull().any():
        raise IOError(f"Not all of the files are in sampleTable, something went wrong."+Error_details)

    return list(files)


def get_quality_controlled_reads_(wildcards,fractions):

    QC_Headers=["Reads_QC_"+f for f in fractions]

    try:
        return get_files_from_sampleTable(wildcards.sample,QC_Headers)
    except FileNotInSampleTableException:

        # return files as named by atlas pipeline
        return expand("{sample}/sequence_quality_control/{sample}_QC_{fraction}.fastq.gz",
                        fraction=MULTIFILE_FRACTIONS,sample=wildcards.sample)


def get_quality_controlled_reads(wildcards):
    """Gets quality controlled reads for two cases. When preprocessed with
    ATLAS, returns R1, R2 and se fastq files or just se. When preprocessed
    externaly and run ATLAS workflow assembly, we expect R1, R2 or se.
    """

    if config.get('interleaved_fastqs',False) and SKIP_QC:
        QC_Headers='se'
    else:
        QC_Headers=MULTIFILE_FRACTIONS
    get_quality_controlled_reads_(wildcards,QC_Headers)





def io_params_for_tadpole(io,key='in'):
    """This function generates the input flag needed for bbwrap/tadpole for all cases
    possible for get_quality_controlled_reads.

    params:
        io  input or output element from snakemake
        key 'in' or 'out'

        if io contains attributes:
            se -> in={se}
            R1,R2,se -> in1={R1},se in2={R2}
            R1,R2 -> in1={R1} in2={R2}

    """
    N= len(io)
    if N==1:
        flag = f"{key}1={io[0]}"
    elif N==2:
        flag= f"{key}1={io[0]} {key}2={io[1]}"
    elif N==3:
        flag= f"{key}1={io[0]},{io[2]} {key}2={io[1]}"
    else:
        logger.error(("File input/output expectation is one of: "
                         "1 file = single-end/ interleaved paired-end "
                         "2 files = R1,R2, or"
                         "3 files = R1,R2,se"
                         "got: {n} files:\n{}").format('\n'.join(io),
                                                       n=len(io)))
        sys.exit(1)
    return flag

def input_paired_only(files):
    """
        return only paired interleaved or not
    """

    if not PAIRED_END:
        raise IOError("Ask for paired end samples if not paired end workflow")

    return files[:min(len(files),2)]


def input_params_for_bbwrap(input):

    if len(input)==3:
        return f"in1={input[0]},{input[2]} in2={input[1]},null"
    else:
        return io_params_for_tadpole(input)
