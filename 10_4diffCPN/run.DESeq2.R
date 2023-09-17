#########################################################################################################################################################################
/home/yuchen/miniconda3/envs/R4.0/bin/R
library(Seurat)
library(Signac)
library(RColorBrewer)
library(cowplot)
library(ggplot2)
library(dplyr)
library(reshape2)

URDobj <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/4_developmental_trajectories.newVersion.cleanedcell/1_URDTree/URDobj.rds")

cell_colors <- c("RG1"="#FDAE61","RG2"="#F46D43","RG3"="#D73027","RG4"="#A50026","RG_Cyc"="#FEE090","nIPC1"="#FFFFBF","nIPC2"="#FFFFBF","nIPC3"="#FFFFBF","ImmatureN1"="#C7E9C0","ImmatureN2"="#A1D99B","SCPN"="#41AB5D","SCPN_CSMN"="#238B45","SCPN_CTPN"="#006D2C","CThPN"="#00441B","MigN"="#ABD9E9","superficialCPN"="#74ADD1","Layer4CPN"="#4575B4","deepCPN"="#313695","MGE_IN"="#54278F","CGE_IN"="#6A51A3","CR"="#225EA8","Oligo"="#253494","MG/Peric"="#081D58")

pdf(file="DimPlot.URDobj.combine_by_Combine_CellType.pdf", height=10,width=10)
DimPlot(URDobj, reduction = "umap_combine", pt.size=1.2, label=TRUE, label.size = 6, group.by='Combine_CellType') + NoLegend() + scale_color_manual(values=c(cell_colors))
dev.off()
pdf(file="DimPlot.URDobj.combine_by_sample.pdf", height=10,width=10)
DimPlot(URDobj, reduction = "umap_combine", pt.size=1.2, label=F, label.size = 6, group.by='sample') + scale_color_manual(values=c("#D52126", "#88CCEE", "#FEE52C"))
dev.off()

#########################################################################################################################################################################
superficialCPNobj <- subset(URDobj, subset=Combine_CellType %in% c("superficialCPN"))
superficialCPN_matrx <- as.data.frame(t(as.matrix(superficialCPNobj@assays$RNA@counts)))
rep1=sample(1:433,433/2)
rep2=setdiff(c(1:433),rep1)
superficialCPNCountsrep1 <- apply(superficialCPN_matrx[rep1,],2,sum)
superficialCPNCountsrep2 <- apply(superficialCPN_matrx[rep2,],2,sum)

Layer4CPNobj <- subset(URDobj, subset=Combine_CellType %in% c("Layer4CPN"))
Layer4CPN_matrx <- as.data.frame(t(as.matrix(Layer4CPNobj@assays$RNA@counts)))
rep1=sample(1:752,752/2)
rep2=setdiff(c(1:752),rep1)
Layer4CPNCountsrep1 <- apply(Layer4CPN_matrx[rep1,],2,sum)
Layer4CPNCountsrep2 <- apply(Layer4CPN_matrx[rep2,],2,sum)

deepCPNobj <- subset(URDobj, subset=Combine_CellType %in% c("deepCPN"))
deepCPN_matrx <- as.data.frame(t(as.matrix(deepCPNobj@assays$RNA@counts)))
rep1=sample(1:288,288/2)
rep2=setdiff(c(1:288),rep1)
deepCPNCountsrep1 <- apply(deepCPN_matrx[rep1,],2,sum)
deepCPNCountsrep2 <- apply(deepCPN_matrx[rep2,],2,sum)

RawCountFilt <- data.frame(GeneID=names(superficialCPNCountsrep1),
WT_REP1=superficialCPNCountsrep1,
WT_REP2=superficialCPNCountsrep2,
NKO_REP1=Layer4CPNCountsrep1,
NKO_REP2=Layer4CPNCountsrep2,
CKO_REP1=deepCPNCountsrep1,
CKO_REP2=deepCPNCountsrep2)

saveRDS(RawCountFilt, file="RawCountFilt.rds")
RawCountFilt <- readRDS("RawCountFilt.rds")


