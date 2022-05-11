library(tidyverse)

argv=commandArgs(T)
INPUT=argv[1]
REF=argv[2]
OUT=argv[3]


library(tidyverse)

data<-read.table(INPUT,comment="#",header=F)%>%
as_tibble()

result<-filter(data,V5<1e-5)%>%
arrange(-V6)%>%
group_by(V1,V3)%>%
slice(1)%>%
ungroup()%>%
dplyr::select(V1,V3,V5,V6,V7)
colnames(result)<-c(
    "VOG","Query","E_value","bit_score","bias"
)

ref<-read_tsv(REF)
colnames(ref)[1]<-"VOG"

result<-left_join(result,ref,by="VOG")

write_tsv(result,OUT)