#!/bin/bash
DIR=$1
FASTQ_DIR=$2
ID=$3
THREADS=$4
LENGTH=$5

mkdir -p ${DIR}
cd ${DIR}
mkdir -p Assembly
cd Assembly

if [ -e ${ID} ]; then
rm -r ${ID}
fi

mkdir -p ${ID}
cd ${ID}

#1. assembly by METASPADES
spades.py --meta -1 ${FASTQ_DIR}/${ID}_bmt-bow_R1.fastq.gz \
-2 ${FASTQ_DIR}/${ID}_bmt-bow_R2.fastq.gz \
-m 376 \
--threads ${THREADS} \
-o metaspades_${ID}

quast.py \
-o quast \
--min-contig ${LENGTH} \
--threads ${THREADS} \
metaspades_${ID}/contigs.fasta

#extract stats
awk 'NR==16 || NR==17 || NR==18 || NR==20 {print $NF}' \
quast/report.txt  | \
tr '\n' '\t' | \
sed "1i Number_of_contigs\tMax_contig_length\tTotal_contig_length\tN50" | \
awk 'BEGIN{OFS="\t"}{print $0}' \
> ${ID}.assembly.stats.txt

cat metaspades_${ID}/contigs.fasta | \
seqkit seq \
-m ${LENGTH} | \
sed -e '/^>/s/$/</g' | tr -d "\n" | \
sed -e 's/>/\n>/g' | sed -e 's/</\n/g' | sed '/^$/d' | \
awk '{print $0}' \
> ${ID}.contigs.fa

gzip -c metaspades_${ID}/contigs.fasta \
> ${ID}.contigs.wofilter.fa.gz
gzip -c metaspades_${ID}/scaffolds.fasta \
> ${ID}.scaffolds.wofilter.fa.gz
gzip -c metaspades_${ID}/scaffolds.paths \
> ${ID}.scaffolds.wofilter.paths.gz
gzip -c metaspades_${ID}/spades.log \
> ${ID}.spades.log.gz

rm -r metaspades_${ID}

#2. mapping read to contig (made by Metaspades)
#index by bt2
bowtie2-build \
-f ${ID}.contigs.fa \
${ID}
gzip -f ${ID}.contigs.fa
#mapping (very sensitive mode)
bowtie2 \
-p 4 --no-discordant --very-sensitive -x ${ID} \
-1 ${FASTQ_DIR}/${ID}_bmt-bow_R1.fastq.gz \
-2 ${FASTQ_DIR}/${ID}_bmt-bow_R2.fastq.gz \
-S ${ID}.sam 2> ${ID}.bt2.log.txt
samtools view -Sb -@ ${THREADS} ${ID}.sam | samtools sort -@ ${THREADS} \
> ${ID}.bam
samtools index ${ID}.bam
rm ${ID}.sam

#extract stats
#1. over all alignment ratio
ALL=`grep "overall alignment rate" ${ID}.bt2.log.txt | awk '{print $1}' | sed 's/%//g'`
#2. concordantly aligned ratio
CONC=`grep "aligned concordantly" ${ID}.bt2.log.txt | awk 'NR>1 {print $2}' | tr -d '(%)' | awk '{a+=$1}END{print a}'`
#3. concordant multi map ratio
CONCM=`grep "aligned concordantly" ${ID}.bt2.log.txt | awk 'NR>2 {print $2}' | tr -d '(%)' | awk '{a+=$1}END{print a}'`
echo -e "Overall\tConc\tConc_multi\n${ALL}\t${CONC}\t${CONCM}" \
> ${ID}.mapping.stats.txt

rm *.bt2