#!/bin/bash
DIR=$1
ID=$2
SET_ID=$3
SCRIPT_DIR=$4
CHECKV_DB=$5

export CHECKVDB=${CHECKV_DB}

mkdir -p ${DIR}
cd ${DIR}
mkdir -p VIRUS_2
cd VIRUS_2

if [ -e ${ID} ]; then
rm -r ${ID}
fi

mkdir -p ${ID}
cd ${ID}

#export VirSorter viral contigs
awk -F"\t" 'NR>1 {print $1}' ${DIR}/VIRUS/${ID}/${ID}_virsorter_and_finder.txt | \
seqkit grep -f - ${DIR}/VIRUS/${ID}/${ID}_INPUT.fa | \
sed -e '/^>/s/$/</g' | tr -d "\n" | \
sed -e 's/>/\n>/g' | sed -e 's/</\n/g' | sed '/^$/d' | \
awk '{print $0}' > ${ID}_virus.fa

#4. Quality estimation by checkV
checkv end_to_end ${ID}_virus.fa checkV -t 1

#remove multiple hit case
grep -v ">1 viral region detected" checkV/quality_summary.tsv \
> quality_summary.tsv
awk 'NR>1 {print $1}' quality_summary.tsv > checkVed_viral_list.tsv 

grep ">1 viral region detected" checkV/quality_summary.tsv > More_than_one_phage.tsv

cat checkV/proviruses.fna | \
awk '{print $1}' | \
sed "s/_1$//g" | \
cat checkV/viruses.fna - | \
seqkit grep -f checkVed_viral_list.tsv - | \
sed -e '/^>/s/$/</g' | tr -d "\n" | \
sed -e 's/>/\n>/g' | sed -e 's/</\n/g' | sed '/^$/d' | \
awk '{print $0}' \
> checkVed_virus.fa

checkv end_to_end checkVed_virus.fa checkV_2nd -t 1

#summarize these result in 1 table
Rscript ${SCRIPT_DIR}/Virome_ver2_metrics_summary.R \
${DIR}/VIRUS/${ID}/${ID}_virsorter_and_finder.txt \
quality_summary.tsv \
checkV_2nd/quality_summary.tsv \
${ID} ${SET_ID}

python3 ${SCRIPT_DIR}/Virome_rename_ver2.py \
checkVed_virus.fa \
Rename_table.tsv \
Renamed_virus.fa

gzip checkVed_virus.fa
gzip ${ID}_virus.fa
gzip Renamed_virus.fa

gzip checkV/*fna
gzip checkV_2nd/*fna

rm -r checkV/tmp
rm -r checkV_2nd/tmp