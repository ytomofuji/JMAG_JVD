#!/bin/bash
DIR=$1
ID=$2

if [ -e ${DIR}/VIRUS_PROPHAGE/${ID} ]; then
rm -r ${DIR}/VIRUS_PROPHAGE/${ID}
fi

mkdir -p ${DIR}/VIRUS_PROPHAGE/${ID}
cd ${DIR}/VIRUS_PROPHAGE/${ID}

LANG=C sort -t$'\t' -k 1,1 ${DIR}/BIN_Metrics/${ID}/contig_to_genome.tsv \
> sorted_contig_to_genome.tsv

awk -F"\t" 'NR>1 {print $1"\t"$3"\t"$4/$2"\t"$4"\t"$2"\t"$24"\t"$8"\t"$15"\t"$16"\t"$17}' \
${DIR}/VIRUS_2/${ID}/All_metrics_summary.tsv | \
LANG=C sort -t$'\t' -k 1,1 | \
LANG=C join -t$'\t' -j 1 sorted_contig_to_genome.tsv - | \
sed '1i Contig_ID\tMAG\tProphage\tViral_ratio\tproviral_length\tcontig_length\tNew_ID\tcheckv_quality\tDecontaminated_gene_count\tDecontaminated_viral_genes\tDecontaminated_host_genes' \
> ${ID}_Prophage_summary.tsv

gzip -f ${ID}_Prophage_summary.tsv
rm sorted_contig_to_genome.tsv