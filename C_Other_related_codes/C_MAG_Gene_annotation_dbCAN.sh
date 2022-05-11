#!/bin/bash
DIR=$1
ID=$2
DB_CAN_REF=$3

<<VARIABLES
DIR: Directory for the analysis
ID: Sample ID
DB_CAN_REF: database for the dbCAN2 (dbCAN-HMMdb-V10.txt)
VARIABLES

if [ -e ${DIR}/BIN_Gene_dbCAN2/${ID} ]; then
rm -r ${DIR}/BIN_Gene_dbCAN2/${ID}
fi

mkdir -p ${DIR}/BIN_Gene_dbCAN2/${ID}
cd ${DIR}/BIN_Gene_dbCAN2/${ID}

awk -F"\t" 'NR>1 {print $1"\t"$2"\t"$3}' ${DIR}/BIN_QC/${ID}/CheckM.passed.bins | \
while read line 
do
BIN_NAME=`echo ${line} | awk '{print $1}'`
hmmscan --domtblout ${BIN_NAME}.out.dm \
${DB_CAN_REF} \
${DIR}/BIN_Gene/${ID}/Prokka_${BIN_NAME}/${BIN_NAME}.faa \
> ${BIN_NAME}.dbCAN2.out
done 
