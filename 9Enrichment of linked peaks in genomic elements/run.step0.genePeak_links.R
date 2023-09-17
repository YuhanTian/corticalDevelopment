#########################################################################################################################################################
# /home/yuchen/miniconda3/envs/R4.0/bin/R
library(Seurat)
library(Signac)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#########################################################################################################################################################

LINKobj <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/3_CREgene_pair/obj.CREgene_pair.rds")
genePeak_links <- Links(LINKobj[["peaks"]])
genePeak_links <- genePeak_links[genePeak_links$score>0]
write.table(as.data.frame(genePeak_links@elementMetadata),"genePeak_links.txt", col.names=T, row.names=F, sep="\t", quote=F)
# setdiff(c("TOP2A","MKI67","SOX2","PAX6","HES1","HES5","OLIG1","PTPRC","PDGFRB","GFAP","EOMES","TUBB3","NEUROD2","NEUROD6","BCL11B","LDB2","FEZF2","NR4A2","TLE4","FOXP2","CRYM","SATB2","NRP1","CUX2","RORB","LMO4","LPL","GAD1","DLX1","LHX6","NR2F2","RELN"),unique(genePeak_links$gene))

