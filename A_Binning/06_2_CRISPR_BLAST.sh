#!/bin/bash
DIR=$1
ID=$2
DB=$3
DB_NAME=$4

SOFT="MINCED"

if [ -e ${DIR}/BIN_CRISPR/${ID}/BLAST/${DB_NAME}_${SOFT} ]; then
rm -r ${DIR}/BIN_CRISPR/${ID}/BLAST/${DB_NAME}_${SOFT}
fi

mkdir -p ${DIR}/BIN_CRISPR/${ID}/BLAST/${DB_NAME}_${SOFT}
cd ${DIR}/BIN_CRISPR/${ID}/BLAST/${DB_NAME}_${SOFT}


#make sum fasta
awk -F"\t" 'NR>1 {print $1"\t"$2"\t"$3}' ${DIR}/BIN_QC/${ID}/CheckM.passed.bins | \
while read line 
do
BIN_NAME=`echo ${line} | awk '{print $1}'`
cat ${DIR}/BIN_CRISPR/${ID}/SPACERS/${SOFT}/${BIN_NAME}_${SOFT}_result_spacers.fa | \
sed 's/>/>'$BIN_NAME'_XXSEPXX_/g'
done > ${ID}_${SOFT}_all_spacers.fa


#RUN BLASTN
# qseqid    sseqid  pident(ANI) length  mismatch    \
# gapopen qstart  qend sstart  send    evalue  bitscore  qlen slen
if [ -s ${ID}_${SOFT}_all_spacers.fa ]; then
blastn \
-query ${ID}_${SOFT}_all_spacers.fa \
-db ${DB} \
-out ${ID}_${SOFT}.blast.tsv \
-perc_identity 90 \
-max_target_seqs 10000 \
-word_size 8 \
-outfmt "6 std qlen slen" \
-num_threads 1


gzip -f ${ID}_${SOFT}.blast.tsv
#filter in
# > 95% ANI
# end-to-end (qlen = qend - qstart +1 )
# query coverage > 95%
zcat ${ID}_${SOFT}.blast.tsv.gz | \
awk -F"\t" '$3>95 && ($5+$6)<=1 && ($8-$7+1)>$13*0.95 {print $0}' | \
sed 's/_XXSEPXX_/\t/g' > ${ID}_${SOFT}.blast.hits
awk -F"_dastool_" '{print $1}' ${ID}_${SOFT}.blast.hits > ${ID}_${SOFT}.blast.MAG_SAMPLE
awk -F"\t" '{print $3}' ${ID}_${SOFT}.blast.hits | \
sed 's/_[0-9]\+$//g' > ${ID}_${SOFT}.blast.VIR_SAMPLE

paste ${ID}_${SOFT}.blast.MAG_SAMPLE ${ID}_${SOFT}.blast.VIR_SAMPLE | \
paste - ${ID}_${SOFT}.blast.hits | \
sed '1i CRISPR_Sample\tVirus_Sample\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen' | \
gzip -c > ${ID}_${SOFT}.blast.hits.gz

rm ${ID}_${SOFT}.blast.MAG_SAMPLE
rm ${ID}_${SOFT}.blast.VIR_SAMPLE
rm ${ID}_${SOFT}.blast.hits

else

touch ${ID}_${SOFT}.blast.tsv
gzip -f ${ID}_${SOFT}.blast.tsv
touch ${ID}_${SOFT}.blast.hits
gzip -f ${ID}_${SOFT}.blast.hits

fi
grep ">" ${ID}_${SOFT}_all_spacers.fa | \
sed "s/>//g" | sed 's/_XXSEPXX_/\t/g' | gzip -c > ${ID}_${SOFT}.spacer.list.gz
rm ${ID}_${SOFT}_all_spacers.fa

