##############################################################################################################
## plot
require(ggplot2)
background<-(theme_bw()
             +theme(plot.title = element_text(size = 16,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.title.y = element_text(size = 16,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.text.y =  element_text(size = 10,vjust = 0.5, colour = "black"))
             +theme(axis.title.x = element_text(size = 16,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.text.x = element_text(size = 12,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(strip.text  = element_text(size = 10,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.line = element_line(size = 0.5),legend.background = element_blank())
             +theme(strip.background = element_rect(fill=NA,colour = NA))
)

df <- read.table("human_gained_lost_enhancer.glm.bed",sep="\t",header = T)
df <- df[which(df$Type!="KM1"),]
df <- df[which(df$Type!="KM3"),]

KMtype <- df$Type

df[which(KMtype=="KM2"),]$Type <- 1
df[which(!(KMtype=="KM2")),]$Type <- 0

lm <- list()
for (i in 8:ncol(df)){
  lm[[i-7]] <- summary(glm(df[,i] ~ df$Type + df$peakwidth + df$conservedbppercent, df, family = binomial))
}

df_result_p <- t(data.frame(lapply(lm, function(x) x$coefficients[2,4])))
df_result_z <- t(data.frame(lapply(lm, function(x) x$coefficients[2,3])))
df_result <- cbind(as.data.frame(df_result_p),as.data.frame(df_result_z))
rownames(df_result) <- colnames(df[,8:ncol(df)])
colnames(df_result)<- c("pvalue","zvalue")

dfplot <- df_result[df_result$pvalue<0.05,]
TFinfo <- strsplit(row.names(dfplot),"_")
TFdf <- do.call(rbind.data.frame, TFinfo)
colnames(TFdf) <- c("TF","TFid")
dfplot <- cbind(dfplot,TFdf)

dfplot <- dfplot[order(dfplot$TF),]
TFfreq <- as.data.frame(table(dfplot$TF))
colnames(TFfreq) <- c("TF","Freq")

tmp1 <- dfplot[which(dfplot$TF %in% TFfreq[which(TFfreq$Freq==2),"TF"]),]
tmp1 <- tmp1[seq(1,nrow(tmp1),2),]
tmp2 <- dfplot[which(dfplot$TF %in% TFfreq[which(TFfreq$Freq==1),"TF"]),]
dfplot <- rbind(tmp1,tmp2)

dfplot <- within(dfplot,{TF <- factor(TF,levels=unique(dfplot[order(dfplot$pvalue),]$TF))})

ggplot(dfplot[dfplot$zvalue>0,],aes(x=TF,y=-log10(pvalue))) +
  geom_bar(stat="identity") +
  xlab("") + ylab("-log10(pvalue)") + labs(title = "Motifs found more often in KM2 versus KM4") +
  background +
  theme(panel.grid =element_blank(),panel.background=element_rect(fill='transparent', color="#000000"),
        axis.ticks.x = element_blank()) + coord_flip()
# c("Arx",
#   "ASCL1",
#   "Atoh1",
#   "ATOH7",
#   "BHLHA15",
#   "BHLHE22",
#   "CUX2",
#   "EOMES",
#   "EN2",
#   "GBX2",
#   "GBX1",
#   "HES6",
#   "HOXD3",
#   "LHX2",
#   "LHX6",
#   "LHX9",
#   "MEIS3",
#   "NEUROD1",
#   "NEUROG2",
#   "NFIA",
#   "OTX2",
#   "PAX3",
#   "POU6F2" ,
#   "PKNOX1",
#   "TCF12",
#   "TGIF1",
#   "TBR1",
#   "ALX3",
#   "VSX2")
ggplot(dfplot[dfplot$zvalue<0,],aes(x=TF,y=-log10(pvalue))) +
  geom_bar(stat="identity") +
  xlab("") + ylab("-log10(pvalue)") + labs(title = "Motifs found more often in KM4 versus KM2") +
  background +
  theme(panel.grid =element_blank(),panel.background=element_rect(fill='transparent', color="#000000"),
        axis.ticks.x = element_blank()) + coord_flip()

##############################################################################################################
#### TF class enrichment ####
cat /data2/yuhan/Project/3D_evolution/human_gained_lost_enhancer/TFBS/TFclass_TF.txt | cut -f1 | sort | uniq -c | awk '$1>=5{print $2}' > TFclass.5.txt
grep -f TFclass.5.txt /data2/yuhan/Project/3D_evolution/human_gained_lost_enhancer/TFBS/TFclass_TF.txt > TFclass.5.TF.txt

TFclass_5_TF <- read.table("TFclass.5.TF.txt",sep="\t",header = F)
colnames(TFclass_5_TF)<- c("TFclass","TF")

TFclass_5_TF$gained <- 0
TFclass_5_TF$lost <- 0
TFclass_5_TF[which(TFclass_5_TF$TF %in% dfplot[dfplot$zvalue>0,]$TF),]$gained <- 1
TFclass_5_TF[which(TFclass_5_TF$TF %in% dfplot[dfplot$zvalue<0,]$TF),]$lost <- 1
TFclass <- unique(TFclass_5_TF$TFclass)

#### KM4 ####
df_fisher_test <- data.frame()
for ( tfclass in TFclass) {
  df1 <- as.data.frame(table(TFclass_5_TF[which(TFclass_5_TF$TFclass==tfclass),]$lost))
  if(nrow(df1)==2)
    df_fisher_test <- rbind(df_fisher_test,data.frame(TFclass=tfclass,N=df1[df1$Var1==0,]$Freq,Y=df1[df1$Var1==1,]$Freq))
}
table(TFclass_5_TF$lost)
df_fisher_test_result <- data.frame()
for (tfclass in df_fisher_test$TFclass) {
  a <- fisher.test(matrix(c(df_fisher_test[which(df_fisher_test$TFclass==tfclass),]$Y,df_fisher_test[which(df_fisher_test$TFclass==tfclass),]$N,17,793),nrow = 2), alternative = "greater")
  df_fisher_test_result <- rbind(df_fisher_test_result,data.frame(TFclass=tfclass,pvalue=a$p.value))
}
# df_fisher_test_result <- df_fisher_test_result[which(df_fisher_test_result$pvalue<0.05),]
summary(-log10(df_fisher_test_result$pvalue))
ggplot(df_fisher_test_result,aes(x=TFclass ,y=-log10(pvalue))) +
  geom_bar(stat="identity") +
  xlab("") + ylab("-log10(pvalue)") + labs(title = "KM4 TF class enrichment") +
  background +
  theme(panel.grid =element_blank(),panel.background=element_rect(fill='transparent', color="#000000"),
        axis.ticks.x = element_blank()) + coord_flip() + scale_y_continuous(limits=c(0,4),breaks=seq(0,4,1))

#### KM2 ####
df_fisher_test <- data.frame()
for ( tfclass in TFclass) {
  df1 <- as.data.frame(table(TFclass_5_TF[which(TFclass_5_TF$TFclass==tfclass),]$gained))
  if(nrow(df1)==2)
    df_fisher_test <- rbind(df_fisher_test,data.frame(TFclass=tfclass,N=df1[df1$Var1==0,]$Freq,Y=df1[df1$Var1==1,]$Freq))
}
table(TFclass_5_TF$gained)
df_fisher_test_result <- data.frame()
for (tfclass in df_fisher_test$TFclass) {
  a <- fisher.test(matrix(c(df_fisher_test[which(df_fisher_test$TFclass==tfclass),]$Y,df_fisher_test[which(df_fisher_test$TFclass==tfclass),]$N,246,564),nrow = 2), alternative = "greater")
  df_fisher_test_result <- rbind(df_fisher_test_result,data.frame(TFclass=tfclass,pvalue=a$p.value))
}
# df_fisher_test_result <- df_fisher_test_result[which(df_fisher_test_result$pvalue<0.05),]
summary(-log10(df_fisher_test_result$pvalue))
ggplot(df_fisher_test_result,aes(x=TFclass ,y=-log10(pvalue))) +
  geom_bar(stat="identity") +
  xlab("") + ylab("-log10(pvalue)") + labs(title = "KM2 TF class enrichment") +
  background +
  theme(panel.grid =element_blank(),panel.background=element_rect(fill='transparent', color="#000000"),
        axis.ticks.x = element_blank()) + coord_flip() + scale_y_continuous(limits=c(0,15),breaks=seq(0,15,1))





