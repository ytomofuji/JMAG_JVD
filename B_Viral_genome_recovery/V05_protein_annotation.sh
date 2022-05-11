#annotation of viral proteins
DIR=$1
THREADS=$2
BATCH_FILE_DIR=$3
INPUT_FASTA=$4
EGG_NOG_DB=$5
VOG_HMM_DIR=$6
SCRIPT_DIR=$7

mkdir -p ${DIR}/Batch_${SGE_TASK_ID}
cd ${DIR}/Batch_${SGE_TASK_ID}

export EGGNOG_DATA_DIR=${EGG_NOG_DB}

#0. call protein by prodigal
seqkit grep -f ${BATCH_FILE_DIR}/Virus_batch_${SGE_TASK_ID}.tsv \
${INPUT_FASTA} | \
sed -e '/^>/s/$/</g' | tr -d "\n" | \
sed -e 's/>/\n>/g' | sed -e 's/</\n/g' | sed '/^$/d' | \
awk '{print $0}' > BATCH_${SGE_TASK_ID}_virus.fa

prodigal -p meta -i BATCH_${SGE_TASK_ID}_virus.fa \
-d BATCH_${SGE_TASK_ID}_virus_genes.fna \
-a BATCH_${SGE_TASK_ID}_virus_genes.faa \
-o BATCH_${SGE_TASK_ID}_virus_genes.gff -f gff

rm BATCH_${SGE_TASK_ID}_virus.fa
gzip -f BATCH_${SGE_TASK_ID}_virus_genes.fna
gzip -f BATCH_${SGE_TASK_ID}_virus_genes.gff

awk '{print $1}' BATCH_${SGE_TASK_ID}_virus_genes.faa | \
sed 's/\*//g' | \
sed -e '/^>/s/$/</g' | tr -d "\n" | \
sed -e 's/>/\n>/g' | sed -e 's/</\n/g' | sed '/^$/d' | \
awk '{print $0}' \
> formatted_BATCH_${SGE_TASK_ID}_virus_genes.faa

gzip -f BATCH_${SGE_TASK_ID}_virus_genes.faa 

#1. Egg nog mapper
emapper.py -i formatted_BATCH_${SGE_TASK_ID}_virus_genes.faa \
--cpu ${THREADS} \
-o BATCH_${SGE_TASK_ID}_virus \
--dbmem \
--override 

#2. VOG database
hmmscan --tblout BATCH_${SGE_TASK_ID}_virus.VOG.tsv --cpu ${THREADS} \
${VOG_HMM_DIR} \
formatted_BATCH_${SGE_TASK_ID}_virus_genes.faa \
> BATCH_${SGE_TASK_ID}_virus.VOG.out

gzip -f BATCH_${SGE_TASK_ID}_virus.VOG.out

Rscript ${SCRIPT_DIR}/RCODE_VOG_processing.r \
BATCH_${SGE_TASK_ID}_virus.VOG.tsv \
${HOME}/reference/virus/VOG/vog.annotations.tsv.gz \
BATCH_${SGE_TASK_ID}_virus.VOG.tophit.tsv

gzip -f test_virus.VOG.tsv
gzip -f test_virus.VOG.tophit.tsv
gzip -f formatted_BATCH_${SGE_TASK_ID}_virus_genes.faa