##########################################################################################################################################################################
# /home/yuchen/miniconda3/envs/R4.0/bin/R
library(glue)
library(Seurat)
library(future)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)
library(RColorBrewer)
library(cowplot)
library(stringr)
require("GenomicRanges")
require("Signac")
options(future.globals.maxSize = 30 * 1024^3)

LINKobj <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/3_CREgene_pair/obj.CREgene_pair.rds")
genePeak_links <- Links(LINKobj[["peaks"]])
genePeak_links <- genePeak_links[genePeak_links$score>0]
linkINFO <- data.frame(linkID=paste0("link_",c(1:length(genePeak_links))),gene=genePeak_links$gene,peak=genePeak_links$peak)

node2_branchL <- linkINFO[linkINFO$gene %in% c("GRIK3","LSAMP","SYN3","SYNE2","MIR100HG","AFDN","KIRREL3","NRXN3","STXBP6","GSG1L"),]$peak
node2_branchR <- linkINFO[linkINFO$gene %in% c("CCBE1","ADAMTS3","GALNT13","ADAMTSL1","RIPOR2","GRIK2","CDH18","LINC02082","DPP10","PCSK5"),]$peak

write.table(node2_branchL,"node2_branchL.txt", col.names=F, row.names=F, sep="\t", quote=F)
write.table(node2_branchR,"node2_branchR.txt", col.names=F, row.names=F, sep="\t", quote=F)

cat node2_branchL.txt | awk 'BEGIN{FS="-";OFS="\t";}{print $1,$2,$3,$0;}' | sort -k1,1 -k2,2n | uniq > node2_branchL_peak.bed
cat node2_branchR.txt | awk 'BEGIN{FS="-";OFS="\t";}{print $1,$2,$3,$0;}' | sort -k1,1 -k2,2n | uniq > node2_branchR_peak.bed