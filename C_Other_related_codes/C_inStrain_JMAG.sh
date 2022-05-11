#!/bin/bash
DIR=$1
ID=$2
FASTQ_DIR=$3
BT2_INDEX=$4

<<VARIABLE
DIR: Directory for the analaysis
ID: Sample ID
FASTQ_DIR: Directory for the QCed fastq files
BT2_INDEX: Directory containing the bowtie2 index file and contig2genome file.

Required files in the ${BT2_INDEX} folder:
${BT2_INDEX}/JMAG*・・・bt2 index
${BT2_INDEX}/JMAG_contig_to_genome.tsv・・・contig2genome file
VARIABLE

if [ -e ${DIR}/JMAG_inStrain/${ID} ]; then
rm -r ${DIR}/JMAG_inStrain/${ID}
fi

mkdir -p ${DIR}/JMAG_inStrain/${ID}
cd ${DIR}/JMAG_inStrain/${ID}

NUM1=`seqkit stats \
${FASTQ_DIR}/${ID}_bmt-bow_R1.fastq.gz | \
awk 'NR>1 {print $5}' | sed 's/,//g' `
NUM2=`seqkit stats \
${FASTQ_DIR}/${ID}_bmt-bow_R2.fastq.gz | \
awk 'NR>1 {print $5}' | sed 's/,//g' `
NUM_TOTAL=`expr $NUM1 + $NUM2`
echo -e ${NUM_TOTAL} > ${ID}_number_of_total_seq_bp.txt

bowtie2 \
-p 4 --no-discordant --very-sensitive \
-x ${BT2_INDEX}/JMAG \
-1 ${FASTQ_DIR}/${ID}_bmt-bow_R1.fastq.gz \
-2 ${FASTQ_DIR}/${ID}_bmt-bow_R2.fastq.gz \
-S ${ID}.sam 2> ${ID}.bt2.log.txt
samtools view -Sb -@ 4 ${ID}.sam | samtools sort -@ 4 \
> ${ID}.bam
samtools index ${ID}.bam
rm ${ID}.sam

#calculate abundance by coverM
awk 'NR > 1 {print $2"\t"$1}' \
${BT2_INDEX}/JMAG_contig_to_genome.tsv \
> genome_def_tmp.tsv

coverm genome \
--genome-definition genome_def_tmp.tsv \
--bam-files ${ID}.bam \
--threads 1 \
-m mean trimmed_mean relative_abundance covered_bases variance \
length count reads_per_base rpkm \
--min-covered-fraction 0 --trim-min 0.05 \
--trim-max 0.95 | \
gzip -f > ${ID}_coverM.tsv.gz

rm genome_def_tmp.tsv

inStrain profile -p 4 \
-s ${BT2_INDEX}/JMAG_contig_to_genome.tsv \
--skip_plot_generation --rarefied_coverage 5 \
--min_cov 5 -o inStrain --min_mapq 2 \
${ID}.bam \
${BT2_INDEX}/JMAG.fa \
--database_mode