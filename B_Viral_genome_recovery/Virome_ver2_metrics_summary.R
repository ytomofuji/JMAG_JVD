library(tidyverse)
argv=commandArgs(T)
DATA1=argv[1] 
DATA2=argv[2] 
DATA3=argv[3] 
ID=argv[4] 
SET_ID=argv[5] 

data1<-read_tsv(DATA1)%>%
dplyr::select(name,length,score,pvalue,Category,Prophage_in_virsorter)%>%
arrange(Category)%>%
dplyr::distinct(name,.keep_all=TRUE)
data2<-read_tsv(DATA2)
data3<-read_tsv(DATA3)%>%
dplyr::select(contig_id,gene_count,viral_genes,host_genes)
colnames(data3)[2:4]<-str_c("Decontaminated_",colnames(data3)[2:4])

data<-data2 %>%
left_join(data3,by="contig_id")%>%
left_join(data1,by=c("contig_id"="name"))
data$Number<-seq(1,nrow(data),1)
data<-mutate(data,New_ID=str_c(SET_ID,"_",ID,"_",Number))

write_tsv(data,"All_metrics_summary.tsv")
dplyr::select(data,contig_id,New_ID)%>%
write_tsv("Rename_table.tsv")