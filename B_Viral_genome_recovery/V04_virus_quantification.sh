#!/bin/bash
DIR=$1
ID=$2
FASTQ_DIR=$3
NUM=$4
BT2_DB=$5

if [ -e ${DIR}/Viral_quant_dwn/${ID} ]; then
rm -r ${DIR}/Viral_quant_dwn/${ID}
fi

mkdir -p ${DIR}/Viral_quant_dwn/${ID}
cd ${DIR}/Viral_quant_dwn/${ID}

#calculate mapping ratio and tag count (idxstat + coverM)
#1. down sampling
seqtk \
sample -s1 ${FASTQ_DIR}/${ID}_bmt-bow_R1.fastq.gz ${NUM} | \
gzip -f > DWN_${ID}_R1.fastq.gz

seqtk \
sample -s1 ${FASTQ_DIR}/${ID}_bmt-bow_R2.fastq.gz ${NUM} | \
gzip -f > DWN_${ID}_R2.fastq.gz
#calculate total seq bp
NUM1=`seqkit stats \
DWN_${ID}_R1.fastq.gz | \
awk 'NR>1 {print $5}' | sed 's/,//g' `
NUM2=`seqkit stats \
DWN_${ID}_R2.fastq.gz | \
awk 'NR>1 {print $5}' | sed 's/,//g' `
NUM_TOTAL=`expr $NUM1 + $NUM2`
echo -e ${NUM_TOTAL} > ${ID}_number_of_total_seq_bp.txt

#2. mapping + coverM
DATABASE="Virus"

bowtie2 -p 1 \
--very-sensitive --no-discordant \
-x ${BT2_DB} \
-1 DWN_${ID}_R1.fastq.gz \
-2 DWN_${ID}_R2.fastq.gz \
-S ${ID}_${DATABASE}.sam \
2> ${ID}_${DATABASE}.bt2.txt

ALL=`grep "overall alignment rate" ${ID}_${DATABASE}.bt2.txt | \
awk '{print $1}' | sed 's/%//g'`
CONC=`grep "aligned concordantly" ${ID}_${DATABASE}.bt2.txt | \
awk 'NR>1 {print $1}' | awk '{a+=$1}END{print a*100/'$NUM'}'`
UNIQ_CONC=`grep "aligned concordantly" ${ID}_${DATABASE}.bt2.txt | \
awk 'NR==2 {print $1}' | awk '{a+=$1}END{print a*100/'$NUM'}'`

echo -e "Overall_align\tConc_align\tUniq_conc_align" > ${ID}_${DATABASE}.bt2summary.txt
echo -e "${ALL}\t${CONC}\t${UNIQ_CONC}" >> ${ID}_${DATABASE}.bt2summary.txt

samtools view -Sb -@ 1 ${ID}_${DATABASE}.sam | samtools sort -@ 1 \
> ${ID}_${DATABASE}.bam
samtools index ${ID}_${DATABASE}.bam
rm ${ID}_${DATABASE}.sam

coverm contig \
--bam-files ${ID}_${DATABASE}.bam \
--threads 1 \
-m mean trimmed_mean covered_bases variance \
length count reads_per_base rpkm \
--min-covered-fraction 0 --trim-min 0.05 \
--trim-max 0.95 | \
gzip -f > ${ID}_${DATABASE}_coverM.tsv.gz

rm ${ID}_${DATABASE}.bam
rm ${ID}_${DATABASE}.bam.bai

rm DWN_${ID}_R1.fastq.gz
rm DWN_${ID}_R2.fastq.gz