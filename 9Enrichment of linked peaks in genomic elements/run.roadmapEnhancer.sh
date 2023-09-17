#!/bin/bash

######################################################################################################################
bedtools="/data/public/software/bedtools.2.25.0/bin/bedtools"
inputDIR="/data1/tang/3D_evolution/DataDownload/roadmap_ChIP-seq/H3K27ac_peaks_hg38"

cat \
${inputDIR}/Brain_E067_BRN.ANG.GYR-H3K27ac.broadPeak.hg19Tohg38.bed \
${inputDIR}/Brain_E068_BRN.ANT.CAUD-H3K27ac.broadPeak.hg19Tohg38.bed \
${inputDIR}/Brain_E069_BRN.CING.GYR-H3K27ac.broadPeak.hg19Tohg38.bed \
${inputDIR}/Brain_E071_BRN.HIPP.MID-H3K27ac.broadPeak.hg19Tohg38.bed \
${inputDIR}/Brain_E072_BRN.INF.TMP-H3K27ac.broadPeak.hg19Tohg38.bed \
${inputDIR}/Brain_E073_BRN.DL.PRFRNTL.CRTX-H3K27ac.broadPeak.hg19Tohg38.bed \
${inputDIR}/Brain_E074_BRN.SUB.NIG-H3K27ac.broadPeak.hg19Tohg38.bed \
| cut -f 1-3 | sort -k1,1 -k2,2n | ${bedtools} merge -i - > tmp1.bed

cat /data3/yuhan/Project/Neuron/scMultiome/3_CREgene_pair/CRE_foldEnrichment/genePeak_links.txt | awk 'NR!=1{print $3"-"$2;}' | awk 'BEGIN{FS="-";OFS="\t";}{print $1,$2,$3,$0;}' | ${bedtools} intersect -a - -b tmp1.bed -wa -wb > tmp2.txt
cat tmp2.txt | cut -f 1-3 | sort | uniq | wc -l #10067
# 11414-10067=1347
# cat /data3/yuhan/Project/Neuron/scMultiome/3_CREgene_pair/CRE_foldEnrichment/genePeak_links.txt | cut -f3 | sort | uniq  | wc -l 11415-1=11414

require("scales")
df <- data.frame(v1_data=c(10067,1347),v2_type=c("overlapped_roadmapEnhancer","non-overlapped_roadmapEnhancer"))
df$v3_per <- percent(df$v1_data/sum(df$v1_data))
df$v4_label <- paste(df$v2_type,df$v1_data,df$v3_per,sep = ":")
pie(df$v1_data,df$v4_label)

library(VennDiagram)
overlapN=10067
nonoverlapN=1347
linkedPeak <- c(paste("A",seq(1,nonoverlapN,1),sep = ""),paste("C",seq(1,overlapN,1),sep = ""))
roadmapEnhancer <- c(paste("B",seq(1,489065-overlapN,1),sep = ""),paste("C",seq(1,overlapN,1),sep = ""))
t <- venn.diagram(list(linkedPeak=linkedPeak,roadmapEnhancer=roadmapEnhancer),filename =NULL,height=3000,width=3000,resolution=500,units="px",fill=c("#999999", "#E69F00"),
                  cat.fontfamily = "serif",cat.fontface = "bold",
                  margin = 0.05,
                  compression="lzw",cex=2, cat.cex=2)
grid.draw(t)
