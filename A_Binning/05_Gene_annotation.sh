#!/bin/bash
DIR=$1
ID=$2
THREADS=$3
EGG_NOG_DB_DIR=$4

if [ -e ${DIR}/BIN_Gene/${ID} ]; then
rm -r ${DIR}/BIN_Gene/${ID}
fi

mkdir -p ${DIR}/BIN_Gene/${ID}
cd ${DIR}/BIN_Gene/${ID}
export EGGNOG_DATA_DIR=${EGG_NOG_DB_DIR}

mkdir -p emapper

#COG_category,KEGG_ko,KEGG_Pathway,CAZy
echo -e "BIN_NAME\tProtein_NAME\tAnnotation" > COG_category.txt
echo -e "BIN_NAME\tProtein_NAME\tAnnotation" > KEGG_ko.txt
echo -e "BIN_NAME\tProtein_NAME\tAnnotation" > KEGG_Pathway.txt
echo -e "BIN_NAME\tProtein_NAME\tAnnotation" > KEGG_Module.txt
echo -e "BIN_NAME\tProtein_NAME\tAnnotation" > CAZy.txt

awk -F"\t" 'NR>1 {print $1"\t"$2"\t"$3}' ${DIR}/BIN_QC/${ID}/CheckM.passed.bins | \
while read line 
do
BIN_NAME=`echo ${line} | awk '{print $1}'`
ORG=`echo ${line} | awk '{print $2}'`
KING=`echo ${line} | awk '{print $3}'`
#extract bins
gunzip -c ${DIR}/BIN_QC/${ID}/QCed_bin/${BIN_NAME}.fa.gz | sed 's/_length_.*$//g' \
> ${BIN_NAME}.fa

samtools faidx ${BIN_NAME}.fa

if [ "$KING" = "Bacteria" ]; then
OBJ_PROKKA="Bacteria"
elif [ "$KING" = "Archaea" ]; then
OBJ_PROKKA="Archaea"
fi

#prokka for prot and BGC
prokka ${BIN_NAME}.fa \
--kingdom ${OBJ_PROKKA} \
--outdir Prokka_${BIN_NAME} \
--prefix ${BIN_NAME} \
--force \
--addgenes \
--locustag ${BIN_NAME} \
--cpus ${THREADS}

#EGG NOG mapper
emapper.py -i Prokka_${BIN_NAME}/${BIN_NAME}.faa \
--cpu ${THREADS} \
-o emapper/${BIN_NAME} \
--dbmem \
--override 
rm -r emappertmp_*
#make tables for analysis (COG_category,GOs,KEGG_ko,KEGG_Pathway,CAZy)
grep -v "#" emapper/${BIN_NAME}.emapper.annotations | awk -F"\t" '$7!="-"{cnt=split($7,cat,"");for(i=1;i<=cnt;i++)print "'$BIN_NAME'\t"$1"\t"cat[i]}' \
>> COG_category.txt
grep -v "#" emapper/${BIN_NAME}.emapper.annotations | awk -F"\t" '$12!="-"{cnt=split($12,cat,",");for(i=1;i<=cnt;i++)print "'$BIN_NAME'\t"$1"\t"cat[i]}' | \
sed 's/ko://g' >> KEGG_ko.txt
grep -v "#" emapper/${BIN_NAME}.emapper.annotations | awk -F"\t" '$13!="-"{cnt=split($13,cat,",");for(i=1;i<=cnt;i++)print "'$BIN_NAME'\t"$1"\t"cat[i]}' | \
grep -v "map" >> KEGG_Pathway.txt
grep -v "#" emapper/${BIN_NAME}.emapper.annotations | awk -F"\t" '$14!="-"{cnt=split($14,cat,",");for(i=1;i<=cnt;i++)print "'$BIN_NAME'\t"$1"\t"cat[i]}' \
>> KEGG_Module.txt
grep -v "#" emapper/${BIN_NAME}.emapper.annotations | awk -F"\t" '$19!="-"{cnt=split($19,cat,",");for(i=1;i<=cnt;i++)print "'$BIN_NAME'\t"$1"\t"cat[i]}' \
>> CAZy.txt
#gzip raw result file
gzip -f emapper/${BIN_NAME}.emapper.annotations
gzip -f emapper/${BIN_NAME}.emapper.seed_orthologs
gzip -f emapper/${BIN_NAME}.emapper.hits

rm ${BIN_NAME}.fa 
rm ${BIN_NAME}.fa.fai
done

gzip -f COG_category.txt
gzip -f KEGG_ko.txt
gzip -f KEGG_Pathway.txt
gzip -f KEGG_Module.txt
gzip -f CAZy.txt