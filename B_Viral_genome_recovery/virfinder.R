argv=commandArgs(T)
INPUT=argv[1]
OUTPUT=argv[2]

library(VirFinder)
predResult<-VF.pred(INPUT)
predResult<-predResult[order(predResult$pvalue),]
write.table(predResult,OUTPUT,quote=F,append=F,row.names=F)
