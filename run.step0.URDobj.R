#################################################################################################################################################
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

#################################################################################################################################################
obj <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/2_cell_annotation_newVersion/obj.integratedATAC_HarmonyRNA_combine.cell_annotation_newVersion.cleanedcell.rds")
cell_type <- c("RG_Cyc","RG1","RG2","RG3","RG4","Oligo","MG/Peric","nIPC1","nIPC2","nIPC3","ImmatureN1","ImmatureN2","SCPN","SCPN_CSMN","SCPN_CTPN","CThPN","MigN","superficialCPN","Layer4CPN","deepCPN","CR","MGE_IN","CGE_IN")
cell_colors <- c("RG1"="#FDAE61","RG2"="#F46D43","RG3"="#D73027","RG4"="#A50026","RG_Cyc"="#FEE090","nIPC1"="#FFFFBF","nIPC2"="#FFFFBF","nIPC3"="#FFFFBF","ImmatureN1"="#C7E9C0","ImmatureN2"="#A1D99B","SCPN"="#41AB5D","SCPN_CSMN"="#238B45","SCPN_CTPN"="#006D2C","CThPN"="#00441B","MigN"="#ABD9E9","superficialCPN"="#74ADD1","Layer4CPN"="#4575B4","deepCPN"="#313695","MGE_IN"="#54278F","CGE_IN"="#6A51A3","CR"="#225EA8","Oligo"="#253494","MG/Peric"="#081D58")
marker <- c("TOP2A","MKI67","SOX2","PAX6","HES1","HES5","OLIG1","PTPRC","PDGFRB","GFAP","EOMES","TUBB3","NEUROD2","NEUROD6","BCL11B","LDB2","FEZF2","DIAPH3","FOXP2","TLE4","CRYM","SATB2","NRP1","CUX2","RORB","LMO4","LPL","RELN","GAD1","DLX1","LHX6","NR2F2")

pdf(file="DimPlot.combine_by_Combine_CellType.pdf", height=10,width=10)
DimPlot(obj, reduction = "umap_combine", pt.size=1.2, label=TRUE, label.size = 6, group.by='Combine_CellType') + NoLegend() + scale_color_manual(values=c(cell_colors))
dev.off()
pdf(file="DimPlot.combine_by_sample.pdf", height=10,width=10)
DimPlot(obj, reduction = "umap_combine", pt.size=1.2, label=F, label.size = 6, group.by='sample') + scale_color_manual(values=c("#D52126", "#88CCEE", "#FEE52C"))
dev.off()

## developmentalTrajectories_obj
developmentalTrajectories_obj <- subset(obj, subset=Combine_CellType %in% c("RG_Cyc","RG1","RG2","RG3","RG4","nIPC1","nIPC2","nIPC3","ImmatureN1","ImmatureN2","SCPN","SCPN_CSMN","SCPN_CTPN","CThPN","MigN","superficialCPN","Layer4CPN","deepCPN"))
tmp <- as.data.frame(developmentalTrajectories_obj@reductions$umap_combine@cell.embeddings)
developmentalTrajectories_obj <- subset(developmentalTrajectories_obj, barcodeID %in% row.names(tmp[tmp$UMAP_1>(-8)&tmp$UMAP_1<(-2)&tmp$UMAP_2>(-4)&tmp$UMAP_2<(5),]), invert=T)
pdf(file="DimPlot.combine_by_Combine_CellType.pdf", height=10,width=10)
DimPlot(developmentalTrajectories_obj, reduction = "umap_combine", pt.size=1.2, label=TRUE, label.size = 6, group.by='Combine_CellType') + NoLegend() + scale_color_manual(values=c(cell_colors))
dev.off()

## URDobj or developmentalTrajectories_obj
URDobj <- developmentalTrajectories_obj
URDobj <- RunTSNE(URDobj, reduction = "harmony", assay = "RNA", dims = 1:50, reduction.name="tsne_combine")
pdf(file="DimPlot_tsne.combine_by_Combine_CellType.pdf", height=10,width=10)
DimPlot(URDobj, reduction = "tsne_combine", pt.size=1.2, label=TRUE, label.size = 6, group.by='Combine_CellType') + NoLegend() + scale_color_manual(values=c(cell_colors))
dev.off()
pdf(file="DimPlot.URDobj.combine_by_sample.pdf", height=10,width=10)
DimPlot(URDobj, reduction = "tsne_combine", pt.size=1.2, label=F, label.size = 6, group.by='sample') + scale_color_manual(values=c("#D52126", "#88CCEE", "#FEE52C"))
dev.off()
saveRDS(URDobj, file="URDobj.rds")
