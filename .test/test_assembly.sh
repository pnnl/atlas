#! /bin/bash
set -euo pipefail





atlas --version



samplenames="Mycoplasma Streptococcus"
databaseDir="databases"
WD='.test/Test_assembly'
reads_dir="example_data/reads/stub"


rm -f $WD/samples.tsv
#
atlas init --db-dir $databaseDir  -w $WD $reads_dir


atlas run -w $WD qc  $@

atlas run assembly -w $WD $@

echo "copy qc reads and assemble"

WD2='.test/Test_assembly_skipQC'
reads_dir=$WD2/"reads"

rm -f $WD2/samples.tsv
mkdir -p $reads_dir
cp $WD/*/sequence_quality_control/*_QC_R?.fastq.gz $reads_dir

atlas init --db-dir $databaseDir --assembler megahit --skip-qc -w $WD2 $reads_dir

atlas run -w $WD2 assembly  $@


echo "start from interleaved QC reads"

WD3='.test/Test_assembly_interleved'
reads_dir=$WD3/"reads"

rm -f $WD3/samples.tsv
mkdir -p $reads_dir

for sample in $samplenames ;
do
reformat.sh in=$WD/$sample/sequence_quality_control/${sample}_QC_R1.fastq.gz \
  in2=$WD/$sample/sequence_quality_control/${sample}_QC_R2.fastq.gz out=$reads_dir/${sample}.fastq.gz overwrite=true
done

atlas init --db-dir $databaseDir --skip-qc --interleaved-fastq -w $WD3 $reads_dir

atlas run assembly --config  threads=2 mem=4 java_mem=4 -w $WD3 -p $@
