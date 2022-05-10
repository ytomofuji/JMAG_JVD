#!/bin/bash
DIR=$1
ID=$2
LENGTH=$3
VIR_FINDER_DIR=$4
VIR_SORTER_DB_DIR=$5

mkdir -p ${DIR}
cd ${DIR}
mkdir -p VIRUS
cd VIRUS

if [ -e ${ID} ]; then
rm -r ${ID}
fi

mkdir -p ${ID}
cd ${ID}
#extract contig >=5kb
zcat ${DIR}/Assembly/${ID}/${ID}.contigs.fa.gz | \
seqkit seq -g \
-m ${LENGTH} | \
sed -e '/^>/s/$/</g' | tr -d "\n" | \
sed -e 's/>/\n>/g' | sed -e 's/</\n/g' | sed '/^$/d' | \
awk '{print $1}' \
> ${ID}_INPUT.fa

#1. run virSorter
wrapper_phage_contigs_sorter_iPlant.pl \
-f ${ID}_INPUT.fa \
--db 2 \
--wdir ${ID}_virsorter \
--ncpu 1 \
--data-dir ${VIR_SORTER_DB_DIR}

grep -v "#" ${ID}_virsorter/VIRSorter_global-phage-signal.csv | \
awk -F"," '$1==$3 {print $3"\t"$5"\t"0}' > ${ID}_virsorter.tmp1
grep -v "#" ${ID}_virsorter/VIRSorter_global-phage-signal.csv | \
awk -F"," '$1!=$3 {print $3"\t"$5"\t"1}' > ${ID}_virsorter.tmp2
cat ${ID}_virsorter.tmp1 ${ID}_virsorter.tmp2 > ${ID}_virsorter.tmp
cat ${ID}_virsorter.tmp | \
awk 'BEGIN{OFS="\t"}{print $1}' | \
sed 's/^VIRSorter_//g' | \
sed 's/-.*$//g' | paste -d"\t" - ${ID}_virsorter.tmp | \
sed '1i Contig_ID\tViral_ID\tCategory\tProphage_in_virsorter' > ${ID}_virsorter.txt
rm ${ID}_virsorter.tmp*

#2. run VirFinder
Rscript ${VIR_FINDER_DIR}/virfinder.R \
${ID}_INPUT.fa \
${ID}_virfinder.txt

#3. extract only highly reliable phages
# Category 1, 2, 4, 5 in VirSorter
# P < 0.01 and score > 0.9 in VirFinder
Rscript ${SCRIPT_DIR}/Virome_sorter_finder_merge.R \
${ID}_virsorter.txt \
${ID}_virfinder.txt \
0.9 0.01 ${ID}_virsorter_and_finder.txt