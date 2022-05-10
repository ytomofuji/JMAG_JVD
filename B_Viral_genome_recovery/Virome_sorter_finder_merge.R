library(tidyverse)

argv=commandArgs(T)
VIRSORTER=argv[1] 
VIRFINDER=argv[2] 
SCORE=as.numeric(argv[3]) 
P_VAL=as.numeric(argv[4]) 
OUT=argv[5]

data1<-read_tsv(VIRSORTER)
data2<-read_delim(VIRFINDER,delim=" ")
data2$VirSorter_name<-gsub("\\.","_",data2$name)
data<-left_join(data2,data1,by=c("VirSorter_name"="Contig_ID"))

result<-filter(data,(Category %in% c(1,2)) | (score > SCORE & pvalue<P_VAL))
write_tsv(result,OUT)