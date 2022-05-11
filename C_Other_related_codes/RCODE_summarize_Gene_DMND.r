library(tidyverse)

argv=commandArgs(T)
DATA1=argv[1]
DATA2=argv[2]
SEQNUM=argv[3]
OUT=argv[4]

data<-bind_rows(read_tsv(DATA1,col_names=F),read_tsv(DATA2,col_names=F))
colnames(data)<-c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")

num_data<-read_tsv(SEQNUM)

result<-tibble()

for(EVAL in c(1e-5,1e-6,1e-10,1e-15)){
for(THR in seq(90,100)){
tmp<-data %>%
filter(pident>=THR)%>%
filter(evalue<EVAL)%>%
group_by(qseqid)%>%
arrange(-bitscore)%>%
slice(1)%>%
ungroup()%>%
mutate(ALIGN_LENGTH=abs(qend-qstart)+1)%>%
group_by(sseqid)%>%
summarize(Count=n(),SeqBP=sum(ALIGN_LENGTH))%>%
mutate(Ratio=Count/sum(num_data$Number_of_reads))%>%
mutate(Depth=SeqBP/sum(num_data$Total_seq_bp))%>%
mutate(Threshold=THR)%>%
mutate(EVAL=EVAL)

result<-bind_rows(result,tmp)
}}

write_tsv(result,OUT)