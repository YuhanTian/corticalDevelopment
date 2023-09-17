#!/bin/bash

linkINFO <- readRDS("linkINFO.rds")
write.table(linkINFO[,c("peak","KMnum")],"superficialCPN_peakRegion.txt", col.names=F, row.names=F, sep="\t", quote=F)

cat ../superficialCPN_peakRegion.txt | awk -F'[\t-]' 'BEGIN{OFS="\t";}{print $1,$2,$3,$1"-"$2"-"$3,$4,$3-$2+1;}' | sort -k1,1 -k2,2n | uniq > superficialCPN_peakRegion.bed
cat superficialCPN_peakRegion.bed > human_gained_lost_enhancer.tmp1.bed

# ##########################################################################################################################################################################
# ## hg38.phastCons100way.conservedRegion.bed
# for i in {1..22} X
# do
# zcat /data2/yuhan/Project/3D_evolution/human_gained_lost_enhancer/TFBS/hg38.100way.phastCons/chr${i}.phastCons100way.wigFix.gz | awk '{if($1=="fixedStep") print $0;else {if($0>0.4) print 0.41;else print 0;}}' > chr${i}.phastCons100way.wigFix.tmp
# wigToBigWig chr${i}.phastCons100way.wigFix.tmp /data/public/refGenome/bwa_index/hg38/hg38.fa.fai chr${i}.phastCons100way.bw
# bigWigToBedGraph chr${i}.phastCons100way.bw chr${i}.phastCons100way.bedGraph
# cat chr${i}.phastCons100way.bedGraph | awk 'BEGIN{OFS="\t";}$4!=0&&$3-$2>=20{print $1,$2+1,$3;}' > chr${i}.phastCons100way.conservedRegion.bed
# done
# cat chr*.phastCons100way.conservedRegion.bed > hg38.phastCons100way.conservedRegion.bed

# ##########################################################################################################################################################################
# ## ${PREFIX}.conservedRegion.overlapNum.txt
# cd /data3/yuhan/Project/Neuron/scMultiome/TFBS/JASPAR2022_hg38_TFBS/downloadData
# bedtools="/data/public/software/bedtools.2.25.0/bin/bedtools"
# conservedRegion=/data3/yuhan/Project/Neuron/scMultiome/TFBS/hg38.100way.phastCons/hg38.phastCons100way.conservedRegion.bed
# ls | grep ".tsv.gz" | sed 's/.tsv.gz//' > list.txt
# cat list.txt | while read PREFIX
# do
# echo "zcat ${PREFIX}.tsv.gz | cut -f1-4 | uniq | ${bedtools} intersect -a - -b ${conservedRegion} -wa -wb -f 1 | cut -f1-4 > ${PREFIX}.conservedRegion.txt" | qsub -N ${PREFIX} -l nodes=node02:ppn=1 -d ./
# done

# cd /data3/yuhan/Project/Neuron/scMultiome/TFBS/JASPAR2022_hg38_TFBS
# cat /data2/yuhan/Project/3D_evolution/human_gained_lost_enhancer/TFBS/TFclass.txt | awk '$1=="AC"||$1=="ID"{print $0;}' | sed -e 's/AC //' -e 's/ID //' > TFclass.tmp1.txt
# sed -n '1~2p' TFclass.tmp1.txt > tmp1
# sed -n '2~2p' TFclass.tmp1.txt > tmp2
# paste tmp1 tmp2 | awk '{print "mv "$1".conservedRegion.txt "$2"_"$1".conservedRegion.txt"}' > run.rename.sh
# bash run.rename.sh
# ls *.conservedRegion.txt | sed 's/.conservedRegion.txt//' > TF.name.txt

##########################################################################################################################################################################
## Differential motif enrichment analysis
## Type/peakwidth/conservedbppercent

bedtools="/data/public/software/bedtools.2.25.0/bin/bedtools"
conservedRegion=/data3/yuhan/Project/Neuron/scMultiome/TFBS/hg38.100way.phastCons/hg38.phastCons100way.conservedRegion.bed
cat human_gained_lost_enhancer.tmp1.bed | sort -k1,1 -k2,2n | ${bedtools} intersect -a - -b ${conservedRegion} -wao > human_gained_lost_enhancer.tmp2.bed
cat human_gained_lost_enhancer.tmp2.bed | awk 'BEGIN{OFS="\t";}{print $4,$10}' | awk 'BEGIN{sum=0;tmp="chr1-939846-941008";}{if($1!=tmp) {print tmp,sum;sum=0;} tmp=$1;sum=sum+$2;}END{print tmp,sum;}' > human_gained_lost_enhancer.tmp3.bed
awk 'BEGIN{OFS="\t";}NR==FNR{a[$1]=$2;next;}{print $0,a[$4]/$6;}' human_gained_lost_enhancer.tmp3.bed human_gained_lost_enhancer.tmp1.bed | sort -k1,1 -k2,2n > human_gained_lost_enhancer.tmp4.bed

## The dependent variable (TFBS) was a binary representation of whether each DA peak contained a motif of a TF or not
bedtools="/data/public/software/bedtools.2.25.0/bin/bedtools"
cp /data3/yuhan/Project/Neuron/scMultiome/5_EN_regulatory_logic/SCPNbranch_TFmotifenrichment/TF.name.txt ./
cat TF.name.txt | while read PREFIX
do
echo "${bedtools} intersect -a human_gained_lost_enhancer.tmp4.bed -b /data3/yuhan/Project/Neuron/scMultiome/TFBS/JASPAR2022_hg38_TFBS/${PREFIX}.conservedRegion.txt -wao > ${PREFIX}.conservedRegion.overlap.tmp.txt" | qsub -N ${PREFIX} -l nodes=node02:ppn=1 -d ./
done

cat TF.name.txt | while read PREFIX
do
echo "${PREFIX}" > ${PREFIX}.conservedRegion.overlapNum.txt
cat ${PREFIX}.conservedRegion.overlap.tmp.txt | awk 'BEGIN{OFS="\t";}{if($NF!=0) $NF=1; print $4,$NF;}' | uniq | cut -f2 >> ${PREFIX}.conservedRegion.overlapNum.txt
done

echo -e "chr\tstart\tend\tID\tType\tpeakwidth\tconservedbppercent" > human_gained_lost_enhancer.tmp5.bed
cat human_gained_lost_enhancer.tmp4.bed >> human_gained_lost_enhancer.tmp5.bed
paste human_gained_lost_enhancer.tmp5.bed *.conservedRegion.overlapNum.txt > human_gained_lost_enhancer.glm.bed

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
# df_result$log10pvalue <- -log10(df_result$pvalue)
# df_result$BH <- p.adjust(df_result$pvalue, method = "BH")

dfplot <- df_result[df_result$pvalue<0.05,]
# dfplot <- df_result[df_result$BH<0.05,]
TFinfo <- strsplit(row.names(dfplot),"_")
TFdf <- do.call(rbind.data.frame, TFinfo)
colnames(TFdf) <- c("TF","TFid")
dfplot <- cbind(dfplot,TFdf)
head(dfplot)
dfplot <- dfplot[dfplot$zvalue>0,]

dfplot <- dfplot[order(dfplot$pvalue),]
dfplot <- head(dfplot,22)
# dfplot <- dfplot[-which(rownames(dfplot)=="Ptf1A_MA1620.1"),]
# dfplot <- dfplot[-which(rownames(dfplot)=="Atoh1_MA0461.2"),]
dfplot <- within(dfplot,{TF <- factor(TF,levels=unique(dfplot[order(dfplot$pvalue),]$TF))})

ggplot(dfplot,aes(x=TF,y=-log10(pvalue))) +
  geom_bar(stat="identity") +
  xlab("") + ylab("-log10(pvalue)") + labs(title = "Motifs found more often in KM2 versus noKM2") +
  background +
  theme(panel.grid =element_blank(),panel.background=element_rect(fill='transparent', color="#000000"),
        axis.ticks.x = element_blank()) + coord_flip()






############################################
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

#### CTPN ####
df_fisher_test <- data.frame()
for ( tfclass in TFclass) {
  df1 <- as.data.frame(table(TFclass_5_TF[which(TFclass_5_TF$TFclass==tfclass),]$lost))
  if(nrow(df1)==2)
    df_fisher_test <- rbind(df_fisher_test,data.frame(TFclass=tfclass,N=df1[df1$Var1==0,]$Freq,Y=df1[df1$Var1==1,]$Freq))
}
table(TFclass_5_TF$lost)
df_fisher_test_result <- data.frame()
for (tfclass in df_fisher_test$TFclass) {
  a <- fisher.test(matrix(c(df_fisher_test[which(df_fisher_test$TFclass==tfclass),]$Y,df_fisher_test[which(df_fisher_test$TFclass==tfclass),]$N,3,807),nrow = 2), alternative = "greater")
  df_fisher_test_result <- rbind(df_fisher_test_result,data.frame(TFclass=tfclass,pvalue=a$p.value))
}
# df_fisher_test_result <- df_fisher_test_result[which(df_fisher_test_result$pvalue<0.05),]
summary(-log10(df_fisher_test_result$pvalue))
ggplot(df_fisher_test_result,aes(x=TFclass ,y=-log10(pvalue))) +
  geom_bar(stat="identity") +
  xlab("") + ylab("-log10(pvalue)") + labs(title = "CTPN TF class enrichment") +
  background +
  theme(panel.grid =element_blank(),panel.background=element_rect(fill='transparent', color="#000000"),
        axis.ticks.x = element_blank()) + coord_flip() + scale_y_continuous(limits=c(0,3),breaks=seq(0,3,1))

#### CSMN ####
df_fisher_test <- data.frame()
for ( tfclass in TFclass) {
  df1 <- as.data.frame(table(TFclass_5_TF[which(TFclass_5_TF$TFclass==tfclass),]$gained))
  if(nrow(df1)==2)
    df_fisher_test <- rbind(df_fisher_test,data.frame(TFclass=tfclass,N=df1[df1$Var1==0,]$Freq,Y=df1[df1$Var1==1,]$Freq))
}
table(TFclass_5_TF$gained)
df_fisher_test_result <- data.frame()
for (tfclass in df_fisher_test$TFclass) {
  a <- fisher.test(matrix(c(df_fisher_test[which(df_fisher_test$TFclass==tfclass),]$Y,df_fisher_test[which(df_fisher_test$TFclass==tfclass),]$N,24,786),nrow = 2), alternative = "greater")
  df_fisher_test_result <- rbind(df_fisher_test_result,data.frame(TFclass=tfclass,pvalue=a$p.value))
}
# df_fisher_test_result <- df_fisher_test_result[which(df_fisher_test_result$pvalue<0.05),]
ggplot(df_fisher_test_result,aes(x=TFclass ,y=-log10(pvalue))) +
  geom_bar(stat="identity") +
  xlab("") + ylab("-log10(pvalue)") + labs(title = "CSMN TF class enrichment") +
  background +
  theme(panel.grid =element_blank(),panel.background=element_rect(fill='transparent', color="#000000"),
        axis.ticks.x = element_blank()) + coord_flip() + scale_y_continuous(limits=c(0,3),breaks=seq(0,3,1))
