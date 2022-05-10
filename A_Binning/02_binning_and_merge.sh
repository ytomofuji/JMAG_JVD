#!/bin/bash
DIR=$1
ID=$2
LENGTH=$3

mkdir -p ${DIR}
cd ${DIR}
mkdir -p BIN
cd BIN

if [ -e ${ID} ]; then
rm -r ${ID}
fi

mkdir -p ${ID}
cd ${ID}
#fa.gz in ${DIR}/Assembly/${ID}/${ID}.contigs.fa.gz
#bam in ${DIR}/Assembly/${ID}/${ID}.bam
mkdir -p work_space
cd work_space
#prepare gz file
gunzip -c -f ${DIR}/Assembly/${ID}/${ID}.contigs.fa.gz \
> ${ID}.contigs.fa

#1. binning (METABAT2)
mkdir -p ../metabat2_bins
jgi_summarize_bam_contig_depths \
--outputDepth ${ID}.metabat2.depth.txt \
${DIR}/Assembly/${ID}/${ID}.bam
metabat2 \
--minContig ${LENGTH} \
--numThreads 4 \
-i ${ID}.contigs.fa \
-a ${ID}.metabat2.depth.txt \
-o ../metabat2_bins/${ID}_metabat2


#2. binning (Maxbin2)
coverm contig \
--bam-files ${DIR}/Assembly/${ID}/${ID}.bam \
--threads 4 \
-m trimmed_mean \
--trim-min 0.05 \
--trim-max 0.95 | \
awk 'NR>1 {print $1"\t"$2}' > ${ID}.abundn.maxbin2

mkdir -p ../maxbin2_bins
run_MaxBin.pl \
-min_contig_length ${LENGTH} \
-thread 4 \
-contig ${ID}.contigs.fa \
-abund ${ID}.abundn.maxbin2 \
-out ../maxbin2_bins/${ID}_maxbin2

#3. binning (CONCOCT)
#split bins
#make fa file for CONCOCT
sed 's/\./_DOT_/g' ${ID}.contigs.fa \
> ${ID}.contigs.CONCOCT.fa

cut_up_fasta.py \
${ID}.contigs.CONCOCT.fa \
-c 10000 \
-o 0 \
--merge_last \
-b ${ID}_10K.tmp.bed | \
sed 's/\./_DOT_/g' \
> ${ID}_10K.contigs.fa

cut -f 1 ${ID}_10K.tmp.bed | \
sed 's/_DOT_/\./g' \
> ${ID}.bed.tmp
cut -f 2- ${ID}_10K.tmp.bed | \
paste -d"\t" ${ID}.bed.tmp - \
> ${ID}_10K.bed
rm ${ID}.bed.tmp
rm ${ID}_10K.tmp.bed
#calculate coverage depth
concoct_coverage_table.py \
${ID}_10K.bed \
${DIR}/Assembly/${ID}/${ID}.bam \
> ${ID}.concoct.cov
#run concoct
concoct \
--composition_file \
${ID}_10K.contigs.fa \
--coverage_file \
${ID}.concoct.cov \
-b ${ID}_concoct/ \
-t 4 \
-l ${LENGTH} -s 1
#merge cutted up cluster
merge_cutup_clustering.py \
${ID}_concoct/clustering_gt${LENGTH}.csv | \
sed 's/_DOT_/\./g' \
> ${ID}_concoct/clustering.merged.csv
#write contigs in output folder
mkdir -p ../concoct_bins
extract_fasta_bins.py \
${ID}.contigs.fa \
${ID}_concoct/clustering.merged.csv \
--output_path ../concoct_bins
if [ -e ../concoct_bins/unbinned.fa ]; then
rm ../concoct_bins/unbinned.fa
fi


#4. merge binning result (DasTools)
#make contig-ID \t BIN-ID file
#A. metabat2
Fasta_to_Scaffolds2Bin.sh \
-i ../metabat2_bins \
-e fa \
> metabat2.contig2bin.tsv
#B. maxbin2
Fasta_to_Scaffolds2Bin.sh \
-i ../maxbin2_bins \
-e fasta \
> maxbin2.contig2bin.tsv
#C. concoct
awk 'NR>1 {print $0}' ${ID}_concoct/clustering.merged.csv | \
sed "s/,/\t${ID}_concoct./g" \
> concoct.contig2bin.tsv
#run dastool 
mkdir -p ../dastool_out
DAS_Tool --search_engine diamond \
-i metabat2.contig2bin.tsv,\
maxbin2.contig2bin.tsv,\
concoct.contig2bin.tsv \
-l metabat2,maxbin2,concoct \
-c ${ID}.contigs.fa \
-o ../dastool_out/${ID} \
--threads 4 \
--create_plots 1 \
--write_bins 1

#remove intermediate files
rm ${ID}.contigs.fa
rm ${ID}_10K.contigs.fa
rm ${ID}_10K.bed

#compress intermediate files
gzip -f ${ID}.metabat2.depth.txt
gzip -f ${ID}.abundn.maxbin2
gzip -f ${ID}.concoct.cov
gzip -f *contig2bin.tsv

tar -zcvf ../metabat2_bins.tar.gz ../metabat2_bins \
--remove-files
tar -zcvf ../maxbin2_bins.tar.gz ../maxbin2_bins \
--remove-files
tar -zcvf ../concoct_bins.tar.gz ../concoct_bins \
--remove-files

tar -zcvf ${ID}_concoct.tar.gz ${ID}_concoct \
--remove-files

gzip -f ../dastool_out/${ID}_proteins.faa
gzip -f ../dastool_out/${ID}.seqlength

cd ../dastool_out
tar -zcvf ${ID}_DASTool_bins.tar.gz ${ID}_DASTool_bins \
--remove-files