import logging
import multiprocessing
import os
import sys
import subprocess
import click

from atlas import __version__
from atlas.conf import make_config,prepare_sample_table,load_configfile
from atlas.conf import validate_config,run_init


logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M",
    format="[%(asctime)s %(levelname)s] %(message)s",
)

def log_exception(msg):
    logging.critical(msg)
    logging.info("Documentation is available at: https://metagenome-atlas.readthedocs.io")
    logging.info("Issues can be raised at: https://github.com/metagenome-atlas/atlas/issues")
    sys.exit(1)

@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context
def cli(obj):
    """ATLAS - workflows for assembly, annotation, and genomic binning of
    metagenomic and metatranscriptomic data.

    For updates and reporting issues, see: https://github.com/metagenome-atlas/atlas
    """

cli.add_command(run_init)




def get_snakefile(file="Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf






## QC command

@cli.command(
    "run",
    context_settings=dict(ignore_unknown_options=True),
    short_help="run atlas main workflow"
)
@click.argument(
    "workflow",
    default="all",
    type=click.Choice(["qc","assembly","binning","genomes","genecatalog","None","all"]),
#    show_default=True,
#    help="Execute only subworkflow.",
)
@click.option("-w",
    "--working-dir",
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    help="location to run atlas.",
    default="."
)
@click.option("-c",
    "--config-file",
    type=click.Path(exists=True,resolve_path=True),
    help="config-file generated with 'atlas init'",
)
@click.option(
    "-j",
    "--jobs",
    type=int,
    help="use at most this many jobs in parallel (see cluster submission for mor details).",
)
@click.option(
    "--profile",
    default=None,
    help="snakemake profile e.g. for cluster execution.",
)
@click.option(
    "-n",
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="Test execution.",
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def run_workflow(workflow, working_dir,config_file, jobs, profile, dryrun, snakemake_args):
    """Runs the ATLAS pipline

    By default all steps are executed but a sub-workflow can be specified.
    Needs a config-file and expects to find a sample table in the working-directory. Both can be generated with 'atlas init'

    Most snakemake arguments can be appended to the command for more info see 'snakemake --help'

    For more details, see: https://metagenome-atlas.readthedocs.io
    """

    if config_file is None:
        config_file = os.path.join(working_dir,'config.yaml')

    if not os.path.exists(config_file):
        logging.critical(f"config-file not found: {config_file}\n"
                         "generate one with 'atlas init'")
        sys.exit(1)

    sample_file=os.path.join(working_dir,'samples.tsv')

    if not os.path.exists(sample_file):
        logging.critical(f"sample.tsv not found in the working directory. "
                         "Generate one with 'atlas init'")
        sys.exit(1)

    validate_config(config_file, workflow)

    conf = load_configfile(config_file)

    db_dir = conf['database_dir']

    cmd = (
        "snakemake --snakefile {snakefile} --directory {working_dir} "
        "{jobs} --rerun-incomplete "
        "--configfile '{config_file}' --nolock "
        " {profile} --use-conda {conda_prefix} {dryrun} --conda-frontend mamba "
        " {target_rule} "
        " {args} "
    ).format(
        snakefile=get_snakefile(),
        working_dir=working_dir,
        jobs="--jobs {}".format(jobs) if jobs is not None else "",
        config_file=config_file,
        profile="" if (profile is None) else "--profile {}".format(profile),
        dryrun="--dryrun" if dryrun else "",
        args=" ".join(snakemake_args),
        target_rule=workflow if workflow!="None" else "",
        conda_prefix= "--conda-prefix "+os.path.join(db_dir,'conda_envs')
    )
    logging.info("Executing: %s" % cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logging.critical(e)
        exit(1)




# Download
@cli.command(
    "download",
    context_settings=dict(ignore_unknown_options=True),
    short_help="download reference files (need ~50GB)",
)
@click.option("-d",
    "--db-dir",
    help="location to store databases",
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    required=True
)
@click.option(
    "-j",
    "--jobs",
    default=1,
    type=int,
    show_default=True,
    help="number of simultaneous downloads",
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def run_download(db_dir,jobs, snakemake_args):
    """Executes a snakemake workflow to download reference database files and validate based on
    their MD5 checksum.
    """

    cmd = (
        "snakemake --snakefile {snakefile} "
        "--jobs {jobs} --rerun-incomplete "
        "--nolock  --use-conda  --conda-prefix {conda_prefix} "
        "--config database_dir='{db_dir}' {add_args} "
        "{args}"
    ).format(
        snakefile=get_snakefile("rules/download.snakefile"),
        jobs=jobs,
        db_dir=db_dir,
        conda_prefix=os.path.join(db_dir,'conda_envs'),
        add_args="" if snakemake_args and snakemake_args[0].startswith("-") else "--",
        args=" ".join(snakemake_args),
    )
    logging.info("Executing: %s" % cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logging.critical(e)
        exit(1)


# from atlas.parsers import refseq_parser
# from atlas.tables import merge_tables


# @cli.command("refseq", short_help="enables tree based LCA and LCA star methods")
# @click.argument("tsv", type=click.Path(exists=True))
# @click.argument("namemap", type=click.Path(exists=True))
# @click.argument("treefile", type=click.Path(exists=True))
# @click.argument("output", type=click.File("w", atomic=True))
# @click.option(
#     "-s",
#     "--summary-method",
#     type=click.Choice(["lca", "majority", "best"]),
#     default="lca",
#     show_default=True,
#     help="summary method for annotating ORFs; when using LCA, it's recommended that one limits the number of hits using --top-fraction though function will be assigned per the best hit; 'best' is fastest",
# )
# @click.option(
#     "-a",
#     "--aggregation-method",
#     type=click.Choice(["lca", "lca-majority", "majority"]),
#     default="lca-majority",
#     show_default=True,
#     help="summary method for aggregating ORF taxonomic assignments to contig level assignment; 'lca' will result in most stringent, least specific assignments",
# )
# @click.option(
#     "--majority-threshold",
#     type=float,
#     default=0.51,
#     show_default=True,
#     help="constitutes a majority fraction at tree node for 'lca-majority' ORF aggregation method",
# )
# @click.option(
#     "--min-identity",
#     type=int,
#     default=70,
#     show_default=True,
#     help="minimum allowable percent ID of BLAST hit",
# )
# @click.option(
#     "--min-bitscore",
#     type=int,
#     default=0,
#     show_default=True,
#     help="minimum allowable bitscore of BLAST hit; 0 disables",
# )
# @click.option(
#     "--min-length",
#     type=int,
#     default=60,
#     show_default=True,
#     help="minimum allowable BLAST alignment length",
# )
# @click.option(
#     "--max-evalue",
#     type=float,
#     default=0.000001,
#     show_default=True,
#     help="maximum allowable e-value of BLAST hit",
# )
# @click.option(
#     "--max-hits",
#     type=int,
#     default=10,
#     show_default=True,
#     help="maximum number of BLAST hits to consider when summarizing ORFs; can drastically alter ORF LCA assignments if too high without further limits",
# )
# @click.option(
#     "--table-name",
#     default="refseq",
#     help="table name within namemap database; expected columns are 'name', 'function', and 'taxonomy'",
# )
# @click.option(
#     "--top-fraction",
#     type=float,
#     default=1,
#     show_default=True,
#     help="filters ORF BLAST hits by only keep hits within this fraction of the highest bitscore; this is recommended over --max-hits",
# )
# def run_refseq_parser(
#     tsv,
#     namemap,
#     treefile,
#     output,
#     summary_method,
#     aggregation_method,
#     majority_threshold,
#     min_identity,
#     min_bitscore,
#     min_length,
#     max_evalue,
#     max_hits,
#     table_name,
#     top_fraction,
# ):
#     """Parse TSV (tabular BLAST output [-outfmt 6]), grabbing taxonomy metadata from ANNOTATION to
#     compute LCAs.
#
#     The BLAST hits are assumed to be sorted by query with decreasing bitscores (best alignment first):
#
#         \b
#         sort -k1,1 -k12,12rn tsv > sorted_tsv
#
#     Annotation file should include your BLAST subject sequence ID, a function, a taxonomy name,
#     the taxonomy ID, and the parent taxonomy ID. This file is generated from `prepare-refseq`:
#
#         \b
#         gi|507051347|ref|WP_016122340.1| two-component sensor histidine kinase Bacillus cereus 1396 86661
#         gi|507052147|ref|WP_016123132.1| two-component sensor histidine kinase Bacillus cereus 1396 86661
#         gi|507053266|ref|WP_016124222.1| two-component sensor histidine kinase Bacillus cereus 1396 86661
#
#     The RefSeq function is always derived from the best BLAST hit.
#
#     The output will give contig, ORF ID, the lineage assigned to the contig based on
#     --aggregation-method, the probability of error (erfc), taxonomy assigned to the ORF, the
#     best hit's product, the best hit's evalue, and the best hit's bitscore:
#
#         \b
#         contig     orf          taxonomy erfc orf_taxonomy refseq               refseq_evalue refseq_bitscore
#         k121_52126 k121_52126_1 root     1.0  root         hypothetical protein 1.0e-41       176.8
#
#     """
#     refseq_parser(
#         tsv,
#         namemap,
#         treefile,
#         output,
#         summary_method,
#         aggregation_method,
#         majority_threshold,
#         min_identity,
#         min_bitscore,
#         min_length,
#         max_evalue,
#         max_hits,
#         table_name,
#         top_fraction,
#     )
#
#
# @cli.command(
#     "gff2tsv", short_help="writes version of Prokka TSV with contig as new first column"
# )
# @click.argument("gff", type=click.Path(exists=True))
# @click.argument("output", type=click.File("w", atomic=True))
# @click.option(
#     "--feature-type",
#     default="all",
#     show_default=True,
#     help="feature type in GFF annotation to print",
# )
# def run_gff_to_tsv(gff, output, feature_type):
#     import re, pandas as pd
#
#     res = dict(
#         gene_id=re.compile(r"ID=(.*?)(?:;|$)"),
#         EC_number=re.compile(r"eC_number=(.*?)(?:;|$)"),
#         gene=re.compile(r"gene=(.*?)(?:;|$)"),
#         product=re.compile(r"product=(.*?)(?:;|$)"),
#         partial=re.compile(r"partial=(.*?)(?:;|$)"),
#         rbs_motif=re.compile(r"rbs_motif=(.*?)(?:;|$)"),
#         gc_cont=re.compile(r"gc_cont=(.*?)(?:;|$)"),
#         confidence=re.compile(r"conf=(.*?)(?:;|$)"),
#         start_type=re.compile(r"start_type=(.*?)(?:;|$)"),
#     )
#
#     def parse_gff_annotation(toks):
#         parsed = {}
#         for key in res:
#             try:
#                 parsed[key] = res[key].findall(toks[-1])[0]
#             except IndexError:
#                 if key == "gene_id":
#                     logging.critical("Unable to locate an ID in [%s]" % toks[-1])
#                     sys.exit(1)
#                 else:
#                     pass
#                     # do not add empty values
#         return parsed
#
#     parsed_df = {}
#
#     with open(gff) as gff_fh:
#         for line in gff_fh:
#             if line.startswith("##FASTA"):
#                 break
#             if line.startswith("#"):
#                 continue
#
#             toks = line.strip().split("\t")
#             if (toks[2] != feature_type) and (feature_type != "all"):
#                 continue
#
#             parsed_line = parse_gff_annotation(toks)
#             parsed_line["feature_type"] = toks[2]
#             parsed_line["contig"] = toks[0]
#             parsed_line["start"] = toks[3]
#             parsed_line["stop"] = toks[4]
#             parsed_line["strand"] = toks[6]
#             parsed_df[parsed_line["gene_id"]] = parsed_line
#
#     parsed_df = pd.DataFrame(parsed_df).T
#     parsed_df.drop("gene_id", inplace=True, axis=1)
#     parsed_df.to_csv(output, sep="\t")
#
#
# @cli.command("munge-blast", short_help="adds contig ID to prokka annotated ORFs")
# @click.argument("tsv", type=click.Path(exists=True))
# @click.argument("gff", type=click.Path(exists=True))
# @click.argument("output", type=click.File("w", atomic=True))
# @click.option(
#     "--gene-id",
#     default="ID",
#     show_default=True,
#     help="tag in gff attributes corresponding to ORF ID",
# )
# def run_munge_blast(tsv, gff, output, gene_id):
#     """Prokka ORFs are reconnected to their origin contigs using the GFF of the Prokka output.
#     Contig output is re-inserted as column 1, altering blast hits to be tabular + an extra initial
#     column that will be used to place the ORFs into context.
#     """
#     import re
#
#     gff_map = dict()
#
#     logging.info("step 1 of 2; parsing %s" % gff)
#     # gff attrs: ID=Flavobacterium_00802;inference=ab initio prediction:Prodigal:2.60;...
#     orf_id_re = re.compile(r"%s=(.*?)\;" % gene_id)
#     with open(gff) as prokka_gff:
#         for line in prokka_gff:
#             if line.startswith("##FASTA"):
#                 break
#             if line.startswith("#"):
#                 continue
#             toks = line.strip().split("\t")
#             try:
#                 orf_id = orf_id_re.findall(toks[-1])[0]
#             except IndexError:
#                 # some, like repeat regions, will not have a locus_tag=, but they also will not
#                 # be in the .faa file that is being locally aligned
#                 logging.warning(
#                     "Unable to locate ORF ID using '%s' for line '%s'"
#                     % (gene_id, " ".join(toks))
#                 )
#                 continue
#             gff_map[orf_id] = toks[0]
#
#     logging.info("step 2 of 2; parsing %s" % tsv)
#     # example blast hit:
#     # Flavobacterium_00002	gi|500936490|ref|WP_012025625.1|	100.0	187	0	0	1	187	1	187	1.7e-99	369.8
#     with open(tsv) as blast_hits:
#         for line in blast_hits:
#             toks = line.strip().split("\t")
#             try:
#                 toks.insert(0, gff_map[toks[0]])
#             except KeyError:
#                 logging.critical("%s was not found in the GFF [%s]" % (toks[0], gff))
#                 logging.critical("processing of %s was halted" % tsv)
#                 sys.exit(1)
#             print(*toks, sep="\t", file=output)
#
#
# @cli.command(
#     "merge-tables", short_help="merge genepredictions TSV, Counts, and Taxonomy"
# )
# @click.argument("prokkatsv", type=click.Path(exists=True))
# @click.argument("refseqtsv", type=click.Path(exists=True))
# @click.argument("output")
# @click.option(
#     "--counts", type=click.Path(exists=True), help="Feature Counts result TSV"
# )
# @click.option(
#     "--eggnog", type=click.Path(exists=True), help="Gene annotations from eggNOG-mapper"
# )
# @click.option(
#     "--completeness", type=click.Path(exists=True), help="CheckM completeness TSV"
# )
# @click.option("--taxonomy", type=click.Path(exists=True), help="CheckM taxonomy TSV")
# @click.option(
#     "--fasta",
#     multiple=True,
#     type=click.Path(exists=True),
#     help="Bin fasta file path; can be specified multiple times",
# )
# def run_merge_tables(
#     prokkatsv, refseqtsv, output, counts, eggnog, completeness, taxonomy, fasta
# ):
#     """Combines gne annotations TSV, RefSeq TSV, and Counts TSV into a single table, merging on "gene_id" tag.
#     """
#     merge_tables(
#         prokkatsv, refseqtsv, output, counts, eggnog, completeness, taxonomy, fasta
#     )
#
if __name__ == "__main__":
    cli()
