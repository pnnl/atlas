#! /bin/bash
set -euo pipefail

atlas --version

databaseDir="databases"
WD='example_data/binning'
reads_dir="example_data/reads/test"
config="--config threads=2 mem=4 java_mem=4"


rm -f $WD/samples.tsv $WD/finished_assembly
touch -m $WD/finished_assembly
touch -m $WD/sample*/*/*.bam # update timestamp of sam files

#
atlas init --db-dir $databaseDir --threads 4  -w $WD --skip-qc --interleaved-fastq $reads_dir

atlas run binning $config -w $WD $@ --omit-from run_checkm_lineage_wf

echo "Copy checkm files


"

for s in sample1 sample2 ;
do
  dest_dir=$WD/$s/binning/DASTool/checkm/
  mkdir -p $dest_dir
  mv $WD/checkm_results/$s/* $dest_dir
done

atlas run binning $config -w $WD $@

# genomes need databases
atlas run None genomes/clustering/allbins2genome.tsv  $config -w $WD $@
