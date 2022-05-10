#!/bin/bash
DIR=$1
ID=$2

if [ -e ${DIR}/BIN_CRISPR/${ID}/SPACERS ]; then
rm -r ${DIR}/BIN_CRISPR/${ID}/SPACERS
fi

mkdir -p ${DIR}/BIN_CRISPR/${ID}/SPACERS
cd ${DIR}/BIN_CRISPR/${ID}/SPACERS
mkdir -p MINCED

#extract bins
awk -F"\t" 'NR>1 {print $1"\t"$2"\t"$3}' ${DIR}/BIN_QC/${ID}/CheckM.passed.bins | \
while read line 
do
BIN_NAME=`echo ${line} | awk '{print $1}'`
gunzip -c ${DIR}/BIN_QC/${ID}/QCed_bin/${BIN_NAME}.fa.gz | sed 's/_length_.*$//g' \
> ${BIN_NAME}.fa
done

cd MINCED
awk -F"\t" 'NR>1 {print $1"\t"$2"\t"$3}' ${DIR}/BIN_QC/${ID}/CheckM.passed.bins | \
while read line 
do
BIN_NAME=`echo ${line} | awk '{print $1}'`

minced -spacers ../${BIN_NAME}.fa \
${BIN_NAME}_MINCED_result.txt
done
cd ..

#remove bins
awk -F"\t" 'NR>1 {print $1"\t"$2"\t"$3}' ${DIR}/BIN_QC/${ID}/CheckM.passed.bins | \
while read line 
do
BIN_NAME=`echo ${line} | awk '{print $1}'`
rm ${BIN_NAME}.fa
done