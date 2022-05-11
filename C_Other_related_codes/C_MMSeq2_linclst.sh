#Priotein cluster level comparison code
#INPUT1, INPUT2 -> fasta
#THR -> threshold for AAI (95, 90, 50)
#OUTPUT -> cluster profile (+optional: clustered faa)

INPUT1=$1
INPUT2=$2
THR=$3
DIR=$4
THREADS=$5
OPTION=$6 #ON or OFF
OUT=$7

<<VARIABLES
INPUT1: Input file 1 (.gz format)
INPUT2: Input file 2 (.gz format, Null file can be specified)
THR: Threshold for the clustering (--min-seq-id parameter in the MMSeqs2)
DIR: Directory for the analysis
THREADS: Number of the CPU cores.
OPTION: ON -> output the representative sequences of the clusters as a fasta file.
OUT: A variable used for the name of the output files.
VARIABLES

mkdir -p ${DIR}
cd ${DIR}

cat ${INPUT1} ${INPUT2} > input.fasta.gz

mmseqs createdb \
input.fasta.gz PROT_MERGED --compressed 1

mmseqs linclust \
PROT_MERGED --compressed 1 \
PROT_MERGED_PC /tmp \
--cov-mode 1 -c 0.8 --kmer-per-seq 80 --min-seq-id ${THR} --threads ${THREADS}

mmseqs createtsv \
PROT_MERGED PROT_MERGED PROT_MERGED_PC \
MMSEQ2.${OUT}.clustered.at.${THR}.tsv
gzip -f MMSEQ2.${OUT}.clustered.at.${THR}.tsv

if [ "${OPTION}" == "ON" ]; then
mmseqs result2repseq \
PROT_MERGED PROT_MERGED_PC PROT_MERGED_PC_rep --compressed 1

mmseqs result2flat \
PROT_MERGED PROT_MERGED PROT_MERGED_PC_rep \
MMSEQ2.${OUT}.clustered.at.${THR}.fasta --use-fasta-header
gzip -f MMSEQ2.${OUT}.clustered.at.${THR}.fasta
fi

rm PROT_MERGED*
rm input.fasta.gz


