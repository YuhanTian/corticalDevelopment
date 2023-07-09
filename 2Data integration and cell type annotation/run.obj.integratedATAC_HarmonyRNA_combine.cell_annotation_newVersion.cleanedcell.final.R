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

############################################################################################################################################################################################################################################################
obj <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/2_cell_annotation_newVersion/obj.integratedATAC_HarmonyRNA_combine.cell_annotation_newVersion.cleanedcell.rds")
cell_type <- c("RG_Cyc","RG1","RG2","RG3","RG4","Oligo","MG/Peric","nIPC1","nIPC2","nIPC3","ImmatureN1","ImmatureN2","SCPN","SCPN_CSMN","SCPN_CTPN","CThPN","MigN","superficialCPN","Layer4CPN","deepCPN","CR","MGE_IN","CGE_IN")
cell_colors <- c("RG1"="#FDAE61","RG2"="#F46D43","RG3"="#D73027","RG4"="#A50026","RG_Cyc"="#FEE090","nIPC1"="#FFFFBF","nIPC2"="#FFFFBF","nIPC3"="#FFFFBF","ImmatureN1"="#C7E9C0","ImmatureN2"="#A1D99B","SCPN"="#41AB5D","SCPN_CSMN"="#238B45","SCPN_CTPN"="#006D2C","CThPN"="#00441B","MigN"="#ABD9E9","superficialCPN"="#74ADD1","Layer4CPN"="#4575B4","deepCPN"="#313695","MGE_IN"="#54278F","CGE_IN"="#6A51A3","CR"="#225EA8","Oligo"="#253494","MG/Peric"="#081D58")
marker <- c("TOP2A","MKI67","SOX2","PAX6","HES1","HES5","OLIG1","PTPRC","PDGFRB","GFAP","EOMES","TUBB3","NEUROD2","NEUROD6","BCL11B","LDB2","FEZF2","DIAPH3","FOXP2","TLE4","CRYM","SATB2","NRP1","CUX2","RORB","LMO4","LPL","RELN","GAD1","DLX1","LHX6","NR2F2")

#################################################################################################################################################
DefaultAssay(obj) <- "RNA"

## combine
obj@active.ident <- obj$combine_clusters
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
num <- length(table(Idents(obj)))+1
jBrewColors <- SpatialColors(num)
jBrewColors <- c(jBrewColors[1:ceiling(num/2)], rev(jBrewColors[(ceiling(num/2)+1):num]))
pdf(file="DimPlot_combine_by_cluster.pdf", height=10,width=10)
DimPlot(obj, reduction = "umap_combine", pt.size=1.2, label=TRUE, cols=jBrewColors, label.size = 6) + NoLegend()
dev.off()

pdf(file="integratedUMAP_combine_marker.counts.pdf", width=80, height=40)
print(FeaturePlot(obj, slot="counts",reduction = "umap_combine", ncol=8,features = marker, pt.size=1.2, order=TRUE, label=TRUE, min.cutoff = "q1", max.cutoff = "q99")+ theme(legend.position = "right"))
dev.off()

## GEX
obj@active.ident <- obj$GEX_clusters
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
num <- length(table(Idents(obj)))+1
jBrewColors <- SpatialColors(num)
jBrewColors <- c(jBrewColors[1:ceiling(num/2)], rev(jBrewColors[(ceiling(num/2)+1):num]))
pdf(file="DimPlot_GEX_by_cluster.pdf", height=10,width=10)
DimPlot(obj, reduction = "umap_harmony", pt.size=1.2, label=TRUE, cols=jBrewColors, label.size = 6) + NoLegend()
dev.off()

pdf(file="integratedUMAP_GEX_marker.counts.pdf", width=80, height=40)
print(FeaturePlot(obj, slot="counts",reduction = "umap_harmony", ncol=8,features = marker, pt.size=1.2, order=TRUE, label=TRUE, min.cutoff = "q1", max.cutoff = "q99")+ theme(legend.position = "right"))
dev.off()

## peaks
obj@active.ident <- obj$peaks_clusters
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
num <- length(table(Idents(obj)))+1
jBrewColors <- SpatialColors(num)
jBrewColors <- c(jBrewColors[1:ceiling(num/2)], rev(jBrewColors[(ceiling(num/2)+1):num]))
pdf(file="DimPlot_peaks_by_cluster.pdf", height=10,width=10)
DimPlot(obj, reduction = "umap_integrated", pt.size=1.2, label=TRUE, cols=jBrewColors, label.size = 6) + NoLegend()
dev.off()

pdf(file="integratedUMAP_peaks_marker.counts.pdf", width=80, height=40)
print(FeaturePlot(obj, slot="counts",reduction = "umap_integrated", ncol=8,features = marker, pt.size=1.2, order=TRUE, label=TRUE, min.cutoff = "q1", max.cutoff = "q99")+ theme(legend.position = "right"))
dev.off()

#################################################################################################################################################
## RNA_CellType
obj$RNA_CellType <- NA
obj$RNA_CellType[obj$GEX_clusters == "17"] <- "RG_Cyc"
obj$RNA_CellType[obj$GEX_clusters == "12"] <- "RG1"
obj$RNA_CellType[obj$GEX_clusters == "10"] <- "RG2"
obj$RNA_CellType[obj$GEX_clusters == "20"] <- "RG2"
obj$RNA_CellType[obj$GEX_clusters == "9"] <- "RG2"
obj$RNA_CellType[obj$GEX_clusters == "1"] <- "RG3"
obj$RNA_CellType[obj$GEX_clusters == "27"] <- "Oligo"
obj$RNA_CellType[obj$GEX_clusters == "28"] <- "MG"
obj$RNA_CellType[obj$GEX_clusters == "23"] <- "Peric"
obj$RNA_CellType[obj$GEX_clusters == "15"] <- "nIPC1"
obj$RNA_CellType[obj$GEX_clusters == "13"] <- "nIPC2"
obj$RNA_CellType[obj$GEX_clusters == "8"] <- "nIPC3"
obj$RNA_CellType[obj$GEX_clusters == "3"] <- "ImmatureN1"
obj$RNA_CellType[obj$GEX_clusters == "0"] <- "ImmatureN2"
obj$RNA_CellType[obj$GEX_clusters == "2"] <- "ImmatureN3"
obj$RNA_CellType[obj$GEX_clusters == "24"] <- "ImmatureN3"
obj$RNA_CellType[obj$GEX_clusters == "4"] <- "SCPN"
obj$RNA_CellType[obj$GEX_clusters == "5"] <- "SCPN"
obj$RNA_CellType[obj$GEX_clusters == "6"] <- "SCPN"
obj$RNA_CellType[obj$GEX_clusters == "7"] <- "SCPN_CSMN"
obj$RNA_CellType[obj$GEX_clusters == "25"] <- "CPN"
obj$RNA_CellType[obj$GEX_clusters == "16"] <- "CPN"
obj$RNA_CellType[obj$GEX_clusters == "19"] <- "CPN"
obj$RNA_CellType[obj$GEX_clusters == "18"] <- "CPN"
obj$RNA_CellType[obj$GEX_clusters == "26"] <- "CThPN"
obj$RNA_CellType[obj$GEX_clusters == "11"] <- "CGE_IN"
obj$RNA_CellType[obj$GEX_clusters == "14"] <- "MGE_IN"
obj$RNA_CellType[obj$GEX_clusters == "29"] <- "CR"
obj$RNA_CellType[obj$GEX_clusters == "21"] <- "UN"
obj$RNA_CellType[obj$GEX_clusters == "22"] <- "UN"
as.data.frame(table(obj$RNA_CellType))

## PEAK_CellType
obj$PEAK_CellType <- NA
obj$PEAK_CellType[obj$peaks_clusters == "12"] <- "RG_Cyc"
obj$PEAK_CellType[obj$peaks_clusters == "0"] <- "RG1"
obj$PEAK_CellType[obj$peaks_clusters == "7"] <- "RG2"
obj$PEAK_CellType[obj$peaks_clusters == "9"] <- "RG3"
obj$PEAK_CellType[obj$peaks_clusters == "11"] <- "nIPC1"
obj$PEAK_CellType[obj$peaks_clusters == "6"] <- "nIPC2"
obj$PEAK_CellType[obj$peaks_clusters == "21"] <- "Oligo"
obj$PEAK_CellType[obj$peaks_clusters == "22"] <- "MG/Peric"
obj$PEAK_CellType[obj$peaks_clusters == "16"] <- "CGE_IN"
obj$PEAK_CellType[obj$peaks_clusters == "18"] <- "CGE_IN"
obj$PEAK_CellType[obj$peaks_clusters == "15"] <- "MGE_IN"
obj$PEAK_CellType[obj$peaks_clusters == "14"] <- "MGE_IN"
obj$PEAK_CellType[obj$peaks_clusters == "1"] <- "ImmatureN1"
obj$PEAK_CellType[obj$peaks_clusters == "19"] <- "ImmatureN1"
obj$PEAK_CellType[obj$peaks_clusters == "5"] <- "ImmatureN2"
obj$PEAK_CellType[obj$peaks_clusters == "4"] <- "SCPN"
obj$PEAK_CellType[obj$peaks_clusters == "2"] <- "SCPN"
obj$PEAK_CellType[obj$peaks_clusters == "10"] <- "SCPN"
obj$PEAK_CellType[obj$peaks_clusters == "8"] <- "SCPN"
obj$PEAK_CellType[obj$peaks_clusters == "3"] <- "MigN"
# obj_integrated <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/Run_Harmony/scMultiome_data_integration/peaksCombine_CFpvalue20/obj_integratedATAC.rds")
# obj_integrated <- FindClusters(object = obj_integrated, verbose = FALSE, algorithm = 3,resolution = 1.8)
# pdf("DimPlot.obj_integrated_by_cluster.pdf", height=10,width=10)
# DimPlot(obj_integrated, reduction = "umap_integrated", pt.size=1.2, label=TRUE, label.size = 6) + NoLegend()
# dev.off()
obj$PEAK_CellType[obj$barcodeID %in% rownames(obj_integrated@meta.data)[obj_integrated$peaks_snn_res.1.8=="12"]] <- "superficialCPN"
obj$PEAK_CellType[obj$peaks_clusters == "13"] <- "Layer4CPN"
obj$PEAK_CellType[obj$peaks_clusters == "17"] <- "deepCPN"
obj$PEAK_CellType[obj$peaks_clusters == "20"] <- "CThPN"
as.data.frame(table(obj$PEAK_CellType))

## Combine_CellType


## combine
obj@active.ident <- obj$combine_clusters
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
num <- length(table(Idents(obj)))+1
jBrewColors <- SpatialColors(num)
jBrewColors <- c(jBrewColors[1:ceiling(num/2)], rev(jBrewColors[(ceiling(num/2)+1):num]))
pdf(file="DimPlot_combine_by_Combine_CellType.pdf", height=10,width=10)
DimPlot(obj, reduction = "umap_combine", pt.size=1.2, label=TRUE, cols=jBrewColors, label.size = 6, group.by='Combine_CellType') + NoLegend()
dev.off()

## GEX
obj@active.ident <- obj$GEX_clusters
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
num <- length(table(Idents(obj)))+1
jBrewColors <- SpatialColors(num)
jBrewColors <- c(jBrewColors[1:ceiling(num/2)], rev(jBrewColors[(ceiling(num/2)+1):num]))
pdf(file="DimPlot_GEX_by_RNA_CellType.pdf", height=10,width=10)
DimPlot(obj, reduction = "umap_harmony", pt.size=1.2, label=TRUE, cols=jBrewColors, label.size = 6, group.by='RNA_CellType') + NoLegend()
dev.off()

## peaks
obj@active.ident <- obj$peaks_clusters
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
num <- length(table(Idents(obj)))+1
jBrewColors <- SpatialColors(num)
jBrewColors <- c(jBrewColors[1:ceiling(num/2)], rev(jBrewColors[(ceiling(num/2)+1):num]))
pdf(file="DimPlot_peaks_by_PEAK_CellType.pdf", height=10,width=10)
DimPlot(obj, reduction = "umap_integrated", pt.size=1.2, label=TRUE, cols=jBrewColors, label.size = 6, group.by='PEAK_CellType') + NoLegend()
dev.off()


#################################################################################################################################################
pdf(file="DimPlot.combine_by_Combine_CellType.pdf", height=10,width=10)
DimPlot(obj, reduction = "umap_combine", pt.size=1.2, label=TRUE, label.size = 6, group.by='Combine_CellType') + NoLegend() + scale_color_manual(values=c("RG1"="#FDAE61","RG2"="#F46D43","RG3"="#D73027","RG4"="#A50026","RG_Cyc"="#FEE090","nIPC1"="#FFFFBF","nIPC2"="#FFFFBF","nIPC3"="#FFFFBF","ImmatureN1"="#C7E9C0","ImmatureN2"="#A1D99B","SCPN"="#41AB5D","SCPN_CSMN"="#238B45","SCPN_CTPN"="#006D2C","CThPN"="#00441B","MigN"="#ABD9E9","superficialCPN"="#74ADD1","Layer4CPN"="#4575B4","deepCPN"="#313695","MGE_IN"="#54278F","CGE_IN"="#6A51A3","CR"="#225EA8","Oligo"="#253494","MG/Peric"="#081D58"))
dev.off()
pdf(file="DimPlot.combine_by_sample.pdf", height=10,width=10)
DimPlot(obj, reduction = "umap_combine", pt.size=1.2, label=F, label.size = 6, group.by='sample') + scale_color_manual(values=c("#D52126", "#88CCEE", "#FEE52C"))
dev.off()

GEXobj <- subset(obj, RNA_CellType %in% c("UN"), invert=T)
pdf(file="DimPlot.GEX_by_RNA_CellType.pdf", height=10,width=10)
DimPlot(GEXobj, reduction = "umap_harmony", pt.size=1.2, label=TRUE, label.size = 6, group.by='RNA_CellType') + NoLegend() + scale_color_manual(values=c("RG1"="#FDAE61","RG2"="#F46D43","RG3"="#D73027","RG_Cyc"="#FEE090","nIPC1"="#FFFFBF","nIPC2"="#FFFFBF","nIPC3"="#FFFFBF","ImmatureN1"="#C7E9C0","ImmatureN2"="#A1D99B","ImmatureN3"="#41AB5D","SCPN"="#238B45","SCPN_CSMN"="#006D2C","CThPN"="#00441B","CPN"="#ABD9E9","MGE_IN"="#54278F","CGE_IN"="#6A51A3","CR"="#225EA8","Oligo"="#253494","MG"="#081D58","Peric"="#081D58"))
dev.off()
pdf(file="DimPlot.GEX_by_sample.pdf", height=10,width=10)
DimPlot(GEXobj, reduction = "umap_harmony", pt.size=1.2, label=F, label.size = 6, group.by='sample') + scale_color_manual(values=c("#D52126", "#88CCEE", "#FEE52C"))
dev.off()

pdf(file="DimPlot.peaks_by_PEAK_CellType.pdf", height=10,width=10)
DimPlot(obj, reduction = "umap_integrated", pt.size=1.2, label=TRUE, label.size = 6, group.by='PEAK_CellType') + NoLegend() + scale_color_manual(values=c("RG1"="#FDAE61","RG2"="#F46D43","RG3"="#D73027","RG_Cyc"="#FEE090","nIPC1"="#FFFFBF","nIPC2"="#FFFFBF","ImmatureN1"="#C7E9C0","ImmatureN2"="#A1D99B","SCPN"="#41AB5D","CThPN"="#00441B","MigN"="#ABD9E9","superficialCPN"="#74ADD1","Layer4CPN"="#4575B4","deepCPN"="#313695","MGE_IN"="#54278F","CGE_IN"="#6A51A3","Oligo"="#253494","MG/Peric"="#081D58"))
dev.off()
pdf(file="DimPlot.peaks_by_sample.pdf", height=10,width=10)
DimPlot(obj, reduction = "umap_integrated", pt.size=1.2, label=F, label.size = 6, group.by='sample') + scale_color_manual(values=c("#D52126", "#88CCEE", "#FEE52C"))
dev.off()

saveRDS(obj, file="obj.integratedATAC_HarmonyRNA_combine.cell_annotation_newVersion.cleanedcell.final.rds")

