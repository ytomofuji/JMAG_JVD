#!/bin/bash
DIR=$1
BLAST_DB=$2
GENE_NAME=$3
ID=$4
FASTQ_DIR=$5
THREADS=$6
NUM=$7
SCRIPT_DIR=$8

<<VARIABLE
DIR: Directory for the analysis
BLAST_DB: Blast database constructed from the target genes
GENE_NAME: Name of the gene set (used for the name of the output files)
ID: Sample ID
FASTQ_DIR: Directory for the QCed fastq files
THREADS: Number of the CPU cores
NUM: Number of the sequencing reads used for the analysis
SCRIPT_DIR: A directory which contains the custom scripts used in the V01_VirSorter_and_VirFinder.sh (i.e. RCODE_summarize_Gene_DMND.r).
VARIABLE

if [ -e ${DIR}/Gene_quant_dwn/${ID}/${GENE_NAME} ]; then
rm -r ${DIR}/Gene_quant_dwn/${ID}/${GENE_NAME}
fi

mkdir -p ${DIR}/Gene_quant_dwn/${ID}/${GENE_NAME}
cd ${DIR}/Gene_quant_dwn/${ID}/${GENE_NAME}

#calculate mapping ratio and tag count (idxstat + coverM)
#calculate number of sequences and total seq bp
#1. down sampling
seqtk \
sample -s1 ${FASTQ_DIR}/${ID}_bmt-bow_R1.fastq.gz ${NUM} | \
gzip -f > DWN_${ID}_R1.fastq.gz

seqtk \
sample -s1 ${FASTQ_DIR}/${ID}_bmt-bow_R2.fastq.gz ${NUM} | \
gzip -f > DWN_${ID}_R2.fastq.gz

echo -e "Number_of_reads\tTotal_seq_bp" > ${ID}_num_seq_and_total_seq_bp.txt
seqkit stats DWN_${ID}_R1.fastq.gz | \
awk 'NR>1 {print $4"\t"$5}' | sed 's/,//g' >> ${ID}_num_seq_and_total_seq_bp.txt
seqkit stats DWN_${ID}_R2.fastq.gz | \
awk 'NR>1 {print $4"\t"$5}' | sed 's/,//g' >> ${ID}_num_seq_and_total_seq_bp.txt

#2. mapping by DMND
diamond blastx --query DWN_${ID}_R2.fastq.gz \
--db ${BLAST_DB} \
--out ${ID}_${GENE_NAME}_DMND_R1.tsv \
--max-target-seqs 5 --id 90 --outfmt 6 --threads ${THREADS}

diamond blastx --query DWN_${ID}_R2.fastq.gz \
--db ${BLAST_DB} \
--out ${ID}_${GENE_NAME}_DMND_R2.tsv \
--max-target-seqs 5 --id 90 --outfmt 6 --threads ${THREADS}

gzip -f ${ID}_${GENE_NAME}_DMND_R1.tsv
gzip -f ${ID}_${GENE_NAME}_DMND_R2.tsv

rm DWN_${ID}_R2.fastq.gz
rm DWN_${ID}_R2.fastq.gz

#3. Transform to the abundance data

Rscript ${SCRIPT_DIR}/RCODE_summarize_Gene_DMND.r \
${ID}_${GENE_NAME}_DMND_R1.tsv.gz ${ID}_${GENE_NAME}_DMND_R2.tsv.gz \
${ID}_num_seq_and_total_seq_bp.txt \
${ID}_${GENE_NAME}_DMND.summary.tsv.gz