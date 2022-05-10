#!/bin/bash
DIR=$1
ID=$2
FASTQ_DIR=$3

if [ -e ${DIR}/BIN_Metrics/${ID} ]; then
rm -r ${DIR}/BIN_Metrics/${ID}
fi

mkdir -p ${DIR}/BIN_Metrics/${ID}
cd ${DIR}/BIN_Metrics/${ID}

awk -F"\t" 'NR>1 {print $1"\t"$2"\t"$3}' ${DIR}/BIN_QC/${ID}/CheckM.passed.bins | \
while read line 
do
BIN_NAME=`echo ${line} | awk '{print $1}'`
zcat ${DIR}/BIN_QC/${ID}/QCed_bin/${BIN_NAME}.fa.gz | grep ">" | sed 's/>//g' | \
awk '{print "'$BIN_NAME'\t"$0}'
done > genome.def.tsv

coverm genome \
--genome-definition genome.def.tsv \
--bam-files ${DIR}/Assembly/${ID}/${ID}.bam \
--threads 2 \
-m trimmed_mean mean \
--min-covered-fraction 0 \
--trim-min 0.05 \
--trim-max 0.95 | \
awk 'BEGIN{OFS="\t"}NR>1 {print $0}' | \
sed '1i BIN_NAME\tTrimmed_mean\tMean' \
> coverM.tsv

awk -F"\t" 'NR>1 {print $1"\t"$2"\t"$3}' ${DIR}/BIN_QC/${ID}/CheckM.passed.bins | \
while read line 
do
BIN_NAME=`echo ${line} | awk '{print $1}'`
zcat ${DIR}/BIN_QC/${ID}/QCed_bin/${BIN_NAME}.fa.gz | grep ">" | sed 's/>//g' | \
awk '{print $0"\t'$BIN_NAME'"}'
done > contig_to_genome.tsv

awk -F"\t" 'NR>1 {print $1"\t"$2"\t"$3}' ${DIR}/BIN_QC/${ID}/CheckM.passed.bins | \
while read line 
do
BIN_NAME=`echo ${line} | awk '{print $1}'`
zcat ${DIR}/BIN_QC/${ID}/QCed_bin/${BIN_NAME}.fa.gz
done > inStrain.genome.fa

bowtie2-build -f \
inStrain.genome.fa \
inStrain.genome

bowtie2 \
-p 2 --no-discordant --very-sensitive -x inStrain.genome \
-1 ${FASTQ_DIR}/${ID}_bmt-bow_R1.fastq.gz \
-2 ${FASTQ_DIR}/${ID}_bmt-bow_R2.fastq.gz \
-S ${ID}.sam 2> ${ID}.bt2.log.txt
samtools view -Sb -@ 2 ${ID}.sam | samtools sort -@ 2 \
> ${ID}.bam
samtools index ${ID}.bam
rm ${ID}.sam

prodigal -i inStrain.genome.fa -d genes.fna

inStrain profile -p 2 -s contig_to_genome.tsv \
--skip_plot_generation --rarefied_coverage 5 --min_cov 5 --output inStrain --min_mapq 2 -g genes.fna \
${ID}.bam inStrain.genome.fa

rm genes.fna
rm -r inStrain/raw_data
cd inStrain/log
gzip -f *txt
cd ../output
gzip -f *.tsv
cd ${DIR}/BIN_Metrics/${ID}
rm inStrain.genome.fa
rm inStrain.genome.*.bt2

mkdir -p tRNA
mkdir -p rRNA

echo -e "Ala\nArg\nAsn\nAsp\nCys\nGln\nGlu\nGly\nHis\nIle\nLeu\nLys\nMet\nPhe\nPro\nSer\nThr\nTrp\nTyr\nVal" \
> typical_tRNA.txt

echo -e "BIN_NAME\tTrim_mean_cov\tMean_cov\tNUM_5SrRNA\tNUM_16SrRNA\tNUM_23SrRNA\tNUM_tRNA\tCov_inStrain\tBreadth_inStrain\tnucl_div_inStrain\trarefied_nucl_div_inStrain\tSNV_number" \
> RNA_inStrain_summary.txt

awk -F"\t" 'NR>1 {print $1"\t"$2"\t"$3}' ${DIR}/BIN_QC/${ID}/CheckM.passed.bins | \
while read line 
do
BIN_NAME=`echo ${line} | awk '{print $1}'`
ORG=`echo ${line} | awk '{print $2}'`
KING=`echo ${line} | awk '{print $3}'`
#extract bins
gunzip -c ${DIR}/BIN_QC/${ID}/QCed_bin/${BIN_NAME}.fa.gz \
> ${BIN_NAME}.fa

samtools faidx ${BIN_NAME}.fa

if [ "$KING" = "Bacteria" ]; then
OBJ_BAR="bac"
OBJ_tRNA="B"
elif [ "$KING" = "Archaea" ]; then
OBJ_BAR="arc"
OBJ_tRNA="A"
fi
#rRNA
barrnap --quiet --kingdom ${OBJ_BAR} \
--reject 0.8 \
-evalue 1e-3 --outseq rRNA/${BIN_NAME}.rRNA.fa \
${BIN_NAME}.fa | grep -v "#" \
> rRNA/${BIN_NAME}.rRNA.list
NUM=`cat rRNA/${BIN_NAME}.rRNA.list | wc -l`
NUM_5S=`grep "5S_rRNA" rRNA/${BIN_NAME}.rRNA.list | wc -l`
NUM_16S=`grep "16S_rRNA" rRNA/${BIN_NAME}.rRNA.list | wc -l`
NUM_23S=`grep "23S_rRNA" rRNA/${BIN_NAME}.rRNA.list | wc -l`
if [ $NUM -eq 0 ]; then
rm rRNA/${BIN_NAME}.rRNA.list
rm rRNA/${BIN_NAME}.rRNA.fa
else
gzip -f rRNA/${BIN_NAME}.rRNA.fa
fi
#tRNA
if [ -e tRNA/${BIN_NAME}.tRNA.list ]; then
rm tRNA/${BIN_NAME}.tRNA.list
fi
tRNAscan-SE \
${BIN_NAME}.fa \
-o tRNA/${BIN_NAME}.tRNA.list -${OBJ_tRNA}
NUM_tRNA=`grep -v "pseudo" tRNA/${BIN_NAME}.tRNA.list | \
awk 'NR>3 && $5!="Undet" {print $5}'  | \
sort -k 1,1 | uniq | join -j 1 typical_tRNA.txt - | wc -l`

#coverM result
TRIM_MEAN=`awk '$1=="'$BIN_NAME'" {print $2}' coverM.tsv`
MEAN=`awk '$1=="'$BIN_NAME'" {print $3}' coverM.tsv`
# coverage	breadth	nucl_diversity
INSTRAIN=`zcat inStrain/output/inStrain_genome_info.tsv | awk -F"\t" '$1=="'$BIN_NAME'" {print $2"\t"$3"\t"$4"\t"$13"\t"$25}'`

echo -e "${BIN_NAME}\t${TRIM_MEAN}\t${MEAN}\t${NUM_5S}\t${NUM_16S}\t${NUM_23S}\t${NUM_tRNA}\t${INSTRAIN}" \
>> RNA_inStrain_summary.txt

rm ${BIN_NAME}.fa 
rm ${BIN_NAME}.fa.fai
done

cut -f 2- RNA_inStrain_summary.txt | \
paste -d "\t" ${DIR}/BIN_QC/${ID}/CheckM.passed.bins - \
> QCpassed_bin_summary_statistics.txt

rm typical_tRNA.txt