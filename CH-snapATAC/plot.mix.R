setwd("/media/ggj/Files/scATAC/CHATAC/HumanMouseMix_20210510/")

library(tidyverse)

hb_all <- list()
for(i in c(10:40)){
  print(i)
  
  m_s1<-read.table(paste0("COL", i ,"/human.metadata.csv"),header = T,sep=",",row.names = 1)
  m_s2<-read.table(paste0("COL", i ,"/mouse.metadata.csv"),header = T,sep=",",row.names = 1)
  
  head(m_s1); head(m_s2)
  dim(m_s1)
  
  m_s1 <- m_s1[,c("UM", "TN")]
  m_s2 <- m_s2[,c("UM", "TN")]
  
  colnames(m_s1)<-c("human_UM","human_TN")
  colnames(m_s2)<-c("mouse_UM","mouse_TN")
  hb<-merge(m_s1,m_s2,by="row.names",all=TRUE)
  
  hb[is.na(hb)]<-0
  hb[,6]<-hb$human_UM+hb$mouse_UM
  hb[,7]<-hb$human_UM/hb$V6
  
  dim(hb); head(hb)
  
  hb_all[[paste0("COL", i)]] = hb
}

hb_all %>% reduce(rbind) -> hmb
hmb[is.na(hmb)] <- 0

dim(hmb); head(hmb)

hmb <- hmb[sort(hmb[,6], decreasing = T, index.return=T)$ix[1:10000], ]

hmb <- hmb[(hmb$human_UM < 30000),]
hmb <- hmb[(hmb$mouse_UM < 30000),]

hmb <- drop_na(hmb)
plot(hmb$human_UM, hmb$mouse_UM, xlab="Human UM",ylab="Mouse UM", main = "Cells",pch=20,cex=0.8)


summary(hmb$V6)
summary(hmb$V7)
summary(hmb$human_TN+hmb$mouse_TN)


result_sel <- result[sort(result[,6], decreasing = T, index.return=T)$ix[1:1000], ]
l = c(paste("Human", sum(result_sel$organism=="Human")), 
      paste("Mouse", sum(result_sel$organism=="Mouse")), 
      paste("Mixed", sum(result_sel$organism=="Mixed")))

result = hmb #hb#
colnames(result) <- c("Cell","HUMAN", "human_TN", "MOUSE", "mouse_TN", "V6", "V7")
result$organism = "Mixed"
result$organism[result$V7>0.95] = "Human"
result$organism[result$V7<0.5] = "Mouse"

dim(result); head(result)

dfHuman = result[result$organism == "Human", ]
dfMouse = result[result$organism == "Mouse", ]
dfMixed = result[result$organism == "Mixed", ]

xlim = c(0, max(result$HUMAN)); ylim = c(0, max(result$MOUSE)); point.cex = 1
colors = c("blue", "red", "grey", "purple")

pdf("mix", width = 10, height = 10)
plot(dfHuman$HUMAN, dfHuman$MOUSE, col = colors[1], pch = 16, xlim = xlim, ylim = ylim, xlab = "Human UM", ylab = "Mouse UM", cex = point.cex)
points(dfMouse$HUMAN, dfMouse$MOUSE, col = colors[2], pch = 16,cex = point.cex)
points(dfMixed$HUMAN, dfMixed$MOUSE, col = colors[3], pch = 16, cex = point.cex)
l = c(paste("Human", dim(dfHuman)[1]), paste("Mouse", dim(dfMouse)[1]), paste("Mixed", dim(dfMixed)[1]))
legend("topright", legend = l, fill = colors)
dev.off()

getwd()
#######################################################################################################
m_s1<-read.table("COL10/human.metadata.csv",header = T,sep=",",row.names = 1)
m_s2<-read.table("COL10/mouse.metadata.csv",header = T,sep=",",row.names = 1)

head(m_s1); head(m_s2)
dim(m_s1)

m_s1 <- m_s1[,c("UM", "UQ")]
m_s2 <- m_s2[,c("UM", "UQ")]

colnames(m_s1)<-c("human_UM","human_UQ")
colnames(m_s2)<-c("mouse_UM","mouse_UQ")
hb<-merge(m_s1,m_s2,by="row.names",all=TRUE)

hb[is.na(hb)]<-0
hb[,6]<-hb$human_UM+hb$mouse_UM
hb[,7]<-hb$human_UM/hb$V6

dim(hb); head(hb)
hb <- hb[sort(hb[,6], decreasing = T, index.return=T)$ix[1:1000],]

plot(hb$human_UM, hb$mouse_UM, xlab="Human UM",ylab="Mouse UM", main = "Cells",pch=20,cex=0.8)


########################################################################
m_read<-read.table("/home/ggj/Documents/YF/Mix1231/OUTPUT/Undetermined_1_filter/COL10/out_cell_readcounts.txt.gz",header = F,sep = "")
m_s<-read.table("/home/ggj/Documents/YF/Mix1231/OUTPUT/Undetermined_1_filter/COL10/_dge.summary.txt",header = T,sep="",row.names = 1)
colnames(m_read)<-"NUM_READS"
rownames(m_read)<-m_read[,2]
name<-rownames(m_s)
hb_r<-m_read[name,]
hb<-merge(hb_r,m_s,by="row.names",all=TRUE)
hbm<-hb[,-3]
rownames(hbm)<-hbm[,1]
par(mfrow=c(1,2))
plot(hbm$NUM_READS,hbm$NUM_GENES , xlab="Reads per cell",ylab="Genes per cell",
     main = "Reads to Genes of cell",pch=20,cex=0.8)
plot(hbm$NUM_READS,hbm$NUM_TRANSCRIPTS , xlab="Reads per cell",ylab="Transcripts per cell",
     main = "Reads to Transcripts of cell",pch=20,cex=0.8)
############
a=read.table("/home/ggj/Documents/YF/20190925/NN0915_NC0915_K912_C912_98-1/OUTPUT/Undetermined_1_filter/C9121/out_cell_readcounts.txt.gz",header = F,stringsAsFactors = F)
x=cumsum(a$V1)
x=x/max(x)
plot(1:length(x),x,type='l',col="blue",xlab="cell barcodes sorted by number of reads",ylab="cumulative fraction of reads",xlim=c(1,10000))
