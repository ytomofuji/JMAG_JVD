# JMAG_JVD
This is a repository of the codes used in the Tomofuji et al (prokaryotic and viral genomes recovered from 787 Japanese gut metagenomes revealed microbial features associated with diets, populations, and diseases).  
Our script recovers  
・Metagenome assembled genomes (MAGs)  
・Viral genomes  
・CRISPR spacers  
from the metagenome shotgun sequencing data.  
Other scripts related to the analysis performed in the Tomofuji et al are also deposited.  

# Overview (Recovery of MAGs)
<div align="center">
<img src="Figures/MAG_Pipeline.jpg" width=60%>
</div>

# Requirements
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

# Step1. Assembly
First, contigs were assembled from gut metagenome shotgun sequencing data with the script `01_assembly_and_mapping_back.sh`.  
Input file should be named as `${ID}_R1.fastq.gz` and `${ID}_R2.fastq.gz`  

Following variables are required:  
`DIR`: Directory for analysis  
`FASTQ_DIR`: Directory of original fastq file  
`ID`: Sample ID   
`THREADS`: Number of the threads  
`LENGTH`: Minimum length of the contigs used for the subsequent analyses (2,000 bp was used in the original manuscript)  

This script outputs contigs and their coverages which were used for the subsequent binning.

# Step2. Binning
Binning with metabat2, maxbin2, and concoct was performed with the script `02_binning_and_merge.sh`.   
The outputs from three softwares were combined with the dastools.

Following variables are required:  
`DIR`: Directory for analysis  
`ID`: Sample ID   
`LENGTH`: Minimum length of the contigs used for the subsequent analyses (2,000 bp was used in the original manuscript)  

This script outputs raw bins which were subjected to the subsequent quality control.

# Step3. QC and refinement  
QC and refinement of the bins made in Step2 were performed with the script `03_QC_and_refine.sh`.  

Following variables are required:  
`DIR`: Directory for analysis  
`ID`: Sample ID   
`SET_NAME`: Name of the dataset (used for the file name of the QCed MAGs) 

This script outputs QCed MAGs (`${SET_NAME}_${ID}_dastool_{1-99}.fa`) and related statistics based on the CheckM.

