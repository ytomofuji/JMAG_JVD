# JMAG_JVD
This is a repository of the codes used in the Tomofuji et al (prokaryotic and viral genomes recovered from 787 Japanese gut metagenomes revealed microbial features associated with diets, populations, and diseases).  
Our script recovers  
・Metagenome assembled genomes (MAGs)  
・Viral genomes  
・CRISPR spacers  
from the metagenome shotgun sequencing data.  
Other scripts related to the analysis performed in the Tomofuji et al are also deposited.  

## Requirements
・bcftools (version 1.10.2)   
・beagle4.1 (27Jan18)   
・beagle5.1 (18May20)   
・bedtools (version 2.29.2)  
・bowtie2 (version 2.3.5.1)  
・fastqc (version 0.11.9)  
・GATK (version 4.1.7)   
・Picard (version 2.22.8)  
・plink (version 1.90b4.4)  
・python3 (version 3.7.6)  
・R (version 4.0.1)  
・samtools (version 2.3.5.1)   
・tidyverse (version 1.3.0)  
・Trimmomatic (version 0.39)  


# Overview (Recovery of MAGs)
<div align="center">
<img src="Figures/MAG_Pipeline.jpg" width=60%>
</div>

## Step1. Assembly
First, contigs were assembled from gut metagenome shotgun sequencing data with the script `01_assembly_and_mapping_back.sh`.  
Input file should be named as `${ID}_R1.fastq.gz` and `${ID}_R2.fastq.gz`  

Following variables are required:  
・`DIR`: Directory for analysis  
・`FASTQ_DIR`: Directory of original fastq file  
・`ID`: Sample ID   
・`THREADS`: Number of the threads  
・`LENGTH`: Minimum length of the contigs used for the subsequent analyses (2,000 bp was used in the original manuscript)  

This script outputs contigs and their coverages which were used for the subsequent binning.

## Step2. Binning
Binning with metabat2, maxbin2, and concoct was performed with the script `02_binning_and_merge.sh`.   
The outputs from three softwares were combined with the dastools.

Following variables are required:  
・`DIR`: Directory for analysis  
・`ID`: Sample ID   
・`LENGTH`: Minimum length of the contigs used for the subsequent analyses (2,000 bp was used in the original manuscript)  

This script outputs raw bins which were subjected to the subsequent quality control.

## Step3. QC and refinement  
QC (CheckM) and refinement (RefineM) of the bins made in Step2 were performed with the script `03_QC_and_refine.sh`.  

Following variables are required:  
・`DIR`: Directory for analysis  
・`ID`: Sample ID   
・`SET_NAME`: Name of the dataset (used for the file name of the QCed MAGs) 

The protein and taxonomy databases for RefineM should be given with the variables `${REFINEM_PROTEIN_DB}` and `${REFINEM_TAXONOMY_DB}`, respectively.

This script outputs QCed MAGs (`${SET_NAME}_${ID}_dastool_{1-99}.fa`) and related statistics based on the CheckM.

## Step4. Analysis on the strain-level diversity, tRNA, and rRNA   
For the QCed MAGs recovered in Step3, we calculated the strain-level diversity by inStrain (average nucleotide diversity).   
We also detected tRNA (tRNAscan-SE) and rRNA (barrnap) in the QCed MAGs.  
These analyses were performed with the script `04_strain_div_and_RNA_anno.sh`.  

Following variables are required:  
・`DIR`: Directory for analysis  
・`ID`: Sample ID   
・`FASTQ_DIR`: Directory of original fastq file  

This script outputs a single summary file which includes the summarized results of the inStrain, tRNAscan-SE, and barrnap.

## Step5. Prediction and annotation of the genes on the MAGs   
For the QCed MAGs recovered in Step3, we predicted the genes with the prodigal.   
Then, functions of these putative genes were annotated with the eggNOG-mapper.  
These analyses were performed with the script `05_Gene_annotation.sh`.  

Following variables are required:  
・`DIR`: Directory for analysis  
・`ID`: Sample ID   
・`THREADS`: Number of the threads  
・`EGG_NOG_DB_DIR`: Directory of the databases for the eggNOG-mapper  

This script outputs the result of the prodigal and eggNOG-mapper per MAGs.
Per-sample summary files for the each annoytation of the eggNOG-mapper (i.e. COG, KEGG gene, KEGG pathway, KEGG module, and CAZy) are also generated.

## Step6. Finding the CRISPR sequences in the MAGs and blast search against the viral sequences
For the QCed MAGs recovered in Step3, we predicted the CRISPR sequences with the MINCED (`06_1_CRISPR_FIND.sh`).   

Following variables are required:  
・`DIR`: Directory for analysis  
・`ID`: Sample ID   

`06_1_CRISPR_FIND.sh` outputs the spacer sequences per MAGs.

Then, CRISPR spacers were subjected to the blast search against the viral genome databases with the script `06_2_CRISPR_BLAST.sh`.

Following variables are required:  
・`DIR`: Directory for analysis  
・`ID`: Sample ID   
・`DB`: A blast database made from a viral genome database   
・`DB_NAME`: Name of the viral genome database (used for the name of the directory of the result files [`${DIR}/BIN_CRISPR/${ID}/BLAST/${DB_NAME}_${SOFT}`])  

This script outputs the result of the blast search.



# Overview (Recovery of the viral genomes)
<div align="center">
<img src="Figures/Virus_Pipeline.jpg" width=60%>
</div>

## StepV1. VirSorter and VirFinder
Viral genomes were detected by the VirSorter and VirFinder from the contigs assembled in Step1.   
This procedure is performed with the sceript `V01_VirSorter_and_VirFinder.sh`.

Following variables are required:  
・`DIR`: Directory for analysis  
・`ID`: Sample ID   
・`LENGTH`: Minimum length of the contigs used for the subsequent analyses (5,000 bp was used in the original manuscript)  
・`VIR_FINDER_DIR`: A directory which contains the VirFinder software (virfinder.R).  
・`VIR_SORTER_DB_DIR`: A directory which contains the databases for the VirFinder.  
・`SCRIPT_DIR`: A directory which contains the custom scripts used in the `V01_VirSorter_and_VirFinder.sh` (i.e. Virome_sorter_finder_merge.R).  

The results of the VirSorter and VirFinder were summarized into the `${ID}_virsorter_and_finder.txt`.

## StepV2. Extraction of the viral contigs and quality control by the CheckV
Viral genomes were extracted based on the result in the StepV1 and subjected to the quality control by the CheckV software.
This procedure is performed with the sceript `V02_contig_extraction_and_CheckV.sh`.

Following variables are required:  
・`DIR`: Directory for analysis  
・`ID`: Sample ID   
・`SET_ID`: Name of the dataset which is used for the name of the output files  
・`SCRIPT_DIR`: A directory which contains the custom scripts used in the `V02_contig_extraction_and_CheckV.sh` (i.e. Virome_ver2_metrics_summary.R and Virome_rename_ver2.py).  
・`CHECKV_DB`: A path to the database for the CheckV.  

The QCed viral genome sequences was sumamrized in the `${DIR}/VIRUS_2/${ID}/Renamed_virus.fa.gz`.

## StepV3. Prophage analysis with the recovered viral genomes
Prophage analysis was performed based on the MAGs and viral genomes recovered from the same samples.
This procedure is performed with the sceript `V03_PROPHAGE_FIND.sh`.

Following variables are required:  
・`DIR`: Directory for analysis  
・`ID`: Sample ID   

The result of the analysis was written in the `${ID}_Prophage_summary.tsv`.

## StepV4. Quantification of the viruses
Mapping-based quantification of the viruses was performed with the script `V04_virus_quantification.sh`.

Following variables are required:  
・`DIR`: Directory for analysis  
・`ID`: Sample ID   
・`FASTQ_DIR`: Directory of original fastq file  
・`NUM`: Number of the reads used for the quantification   
・`BT2_DB`: Bowtie2 index made from the viral genomes  

The result of the analysis was written in the `${ID}_${DATABASE}_coverM.tsv.gz`.





