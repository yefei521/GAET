library("fpc")
#eps设置成单核小体+linker长度200，或半linker+半核小体=100

#H2A.Z
data<-read.delim("Z:/UTR_annotation/机器学习模型/DBSCAN-聚类/H2A.Z_group.gff3",sep="\t", header = F)
#model2 <- dbscan(data, eps=147, MinPts=5)
for (i in unique(data$V1)){  a<-subset(data,V1==i); a1<-cbind(a$V4,V2=1); colnames(a1)<-c("V1","V2"); model3 <- dbscan(a1, eps=200, MinPts=3); a1<-as.matrix(a1); a2<-cbind(model3$cluster,i,a1[,1]); write.csv(a2,file=paste("Z:/UTR_annotation/机器学习模型/DBSCAN-聚类/H2A.Z/",i,".csv"，collapse = NULL, sep=""), sep="\t") }  
data<-read.delim("Z:/UTR_annotation/机器学习模型/DBSCAN-聚类/H2A.Z_group.gff3",sep="\t", header = F)

#nucleosome
data<-read.delim("Z:/UTR_annotation/机器学习模型/DBSCAN-聚类/H2A.Z_104-160_input.srt.rmdup.bed",sep="\t", header = F) 
for (i in unique(data$V1)){  a<-subset(data,V1==i); a1<-cbind(a$V2,V2=1); colnames(a1)<-c("V1","V2"); model3 <- dbscan(a1, eps=200, MinPts=3); a1<-as.matrix(a1); a2<-cbind(model3$cluster,i,a1[,1]); write.csv(a2,file=paste("Z:/UTR_annotation/机器学习模型/DBSCAN-聚类/nucleosome/",i,".csv"，collapse = NULL, sep=""), sep="\t") }

#H3K4me3
H3K4me3<-read.delim("Z:/UTR_annotation/机器学习模型/DBSCAN-聚类/H3K4me3_group.gff3",sep="\t", header = F)
for (i in unique(data$V1)){  a<-subset(data,V1==i); a1<-cbind(a$V4,V2=1); colnames(a1)<-c("V1","V2"); model3 <- dbscan(a1, eps=200, MinPts=3); a1<-as.matrix(a1); a2<-cbind(model3$cluster,a1[,1]); write.csv(a2,file=paste("Z:/UTR_annotation/机器学习模型/DBSCAN-聚类/H3K4me3/",i,".csv"，collapse = NULL, sep=""), sep="\t") }

#6mA
data<-read.delim("Z:/UTR_annotation/机器学习模型/DBSCAN-聚类/veg_6mA.sorted.gff",sep="\t", header = F) 
for (i in unique(data$V1)){  a<-subset(data,V1==i); a1<-cbind(a$V4,V2=1); colnames(a1)<-c("V1","V2"); model3 <- dbscan(a1, eps=200, MinPts=3); a1<-as.matrix(a1); a2<-cbind(model3$cluster,i,a1[,1]); write.csv(a2,file=paste("Z:/UTR_annotation/机器学习模型/DBSCAN-聚类/6mA/",i,".csv"，collapse = NULL, sep=""), sep="\t") }





