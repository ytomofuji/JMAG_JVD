#!/bin/bash
DIR=$1
ID=$2
SET_NAME=$3

##### requires 8 threads #####

<<DB
REFINEM_PROTEIN_DB=gtdb_r95_protein_db.2020-07-30.faa
REFINEM_TAXONOMY_DB=gtdb_r95_taxonomy.2020-07-30.tsv
DB

mkdir -p ${DIR}
cd ${DIR}
mkdir -p BIN_QC
cd BIN_QC

if [ -e ${ID} ]; then
rm -r ${ID}
fi

mkdir -p ${ID}
cd ${ID}
#fa.gz in ${DIR}/Assembly/${ID}/${ID}.contigs.fa.gz
#bam in ${DIR}/Assembly/${ID}/${ID}.bam

#prepare bin files
tar -xvzf ${DIR}/BIN/${ID}/dastool_out/${ID}_DASTool_bins.tar.gz 
tar -xvzf ${DIR}/BIN/${ID}/metabat2_bins.tar.gz
tar -xvzf ${DIR}/BIN/${ID}/maxbin2_bins.tar.gz
tar -xvzf ${DIR}/BIN/${ID}/concoct_bins.tar.gz

#0. checkM for bins before DASTools
checkm lineage_wf -t 8 -x fa \
metabat2_bins \
metabat2_checkm
checkm qa metabat2_checkm/lineage.ms \
metabat2_checkm/ -o 2 -t 8 --tab_table \
-f metabat2_checkm.tsv

checkm lineage_wf -t 8 -x fasta \
maxbin2_bins \
maxbin2_checkm
checkm qa maxbin2_checkm/lineage.ms \
maxbin2_checkm/ -o 2 -t 8 --tab_table \
-f maxbin2_checkm.tsv

checkm lineage_wf -t 8 -x fa \
concoct_bins \
concoct_checkm
checkm qa concoct_checkm/lineage.ms \
concoct_checkm/ -o 2 -t 8 --tab_table \
-f concoct_checkm.tsv

tar -zcvf metabat2_checkm.tar.gz metabat2_checkm \
--remove-files
tar -zcvf maxbin2_checkm.tar.gz maxbin2_checkm \
--remove-files
tar -zcvf concoct_checkm.tar.gz concoct_checkm \
--remove-files
rm -r metabat2_bins
rm -r maxbin2_bins
rm -r concoct_bins

#1. checkM for bins DASTool
checkm lineage_wf -t 8 -x fa \
${ID}_DASTool_bins \
DASTool_checkm
checkm qa DASTool_checkm/lineage.ms \
DASTool_checkm/ -o 2 -t 8 --tab_table \
-f DASTool_checkm.tsv

tar -zcvf DASTool_checkm.tar.gz DASTool_checkm \
--remove-files

#2. run refineM
#prepare original contig for refineM
mkdir -p refineM
gunzip -c -f ${DIR}/Assembly/${ID}/${ID}.contigs.fa.gz \
> ${ID}.contigs.fa
#refineM (contig prop module)
refinem scaffold_stats -c 8 --genome_ext fa \
${ID}.contigs.fa ${ID}_DASTool_bins refineM ${DIR}/Assembly/${ID}/${ID}.bam
refinem outliers refineM/scaffold_stats.tsv refineM
refinem filter_bins --genome_ext fa ${ID}_DASTool_bins \
refineM/outliers.tsv refined_bins_tmp
#refineM (Taxonomy)
refinem call_genes -c 8 --genome_ext fa refined_bins_tmp refineM/Gene
refinem taxon_profile -c 8 refineM/Gene \
refineM/scaffold_stats.tsv \
${REFINEM_PROTEIN_DB} \
${REFINEM_TAXONOMY_DB} \
refineM/TAX
refinem taxon_filter -c 8 refineM/TAX refineM/taxon_filter.tsv
refinem filter_bins --genome_ext fa refined_bins_tmp \
refineM/taxon_filter.tsv refined_bins
rm -r refined_bins_tmp
rm -r ${ID}_DASTool_bins
#run checkM for the result of refineM
#remove .filtered added by refineM
cd refined_bins
rename ".filtered" "" *.fa
rename ".filtered" "" *.fa
cd ..

checkm lineage_wf -t 8 -x fa \
refined_bins \
refined_checkm
checkm qa refined_checkm/lineage.ms \
refined_checkm/ -o 2 -t 8 --tab_table \
-f refined_checkm.tsv
checkm tree_qa refined_checkm/ -o 2 --tab_table \
-f refined_tree_checkm.tsv

tar -zcvf refined_checkm.tar.gz refined_checkm \
--remove-files

#3. filtering BINS based on checkM result
# completeness > 50%, contamination < 5%, QS (completeness - 5xcontamination) > 50
# these bins proceed to further QC
awk -F"\t" 'NR>1 {print $1"\t"$5}' refined_tree_checkm.tsv | \
awk -F";" '{print $1}' | \
sed 's/k__//g' | \
LANG=C sort -k 1,1 -t$'\t' \
> CheckM.tax.sort

awk -F"\t" 'NR>1 && $2!="root (UID1)" {print $1"\t"$3"\t"$6"\t"$7"\t"$6-5*$7"\t"$8"\t"$9\
"\t"$12"\t"$14"\t"$16"\t"$18"\t"$19}' \
refined_checkm.tsv | \
awk -F"\t" 'BEGIN{OFS="\t"}$3>50 && $4<5 && $5>50 {print $0}' | \
LANG=C sort -k 1,1 -t$'\t' | \
join -t$'\t' -1 1 -2 1 CheckM.tax.sort - \
> CheckM.passed.bins.pre

#make ID for QC passed bins ${SET}_${ID}_dastool_${BIN_ID}
awk -F"\t" 'BEGIN{OFS="\t"} {print "'$SET_NAME'_'$ID'_dastool_"NR,$0}' \
CheckM.passed.bins.pre | \
sed '1i BIN_ID\tOriginal_ID\tkingdom\tNumber_of_genome\tCompleteness\tContamination\tQS\tStrain_heterogeneity\tGenome_size\tNumber_of_contigs\tN50\tMean_contig_length\tLongest_contig\tGC_perc' \
> CheckM.passed.bins
rm CheckM.passed.bins.pre
mkdir -p QCed_bin

awk -F"\t" 'NR>1 {print $1"\t"$2"\t"$3}' CheckM.passed.bins | \
while read line
do
BIN_NAME=`echo ${line} | awk '{print $1}'`
ORG=`echo ${line} | awk '{print $2}' | sed 's/.filtered//g' | sed 's/.filtered//g'`
KING=`echo ${line} | awk '{print $3}'`
#extract bins
cp refined_bins/${ORG}.fa QCed_bin/${BIN_NAME}.fa
done

tar -zcvf refined_bins.tar.gz refined_bins \
--remove-files

cd QCed_bin
gzip -f *.fa
