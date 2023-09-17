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
cell_colors <- c("RG1"="#FDAE61","RG2"="#F46D43","RG3"="#D73027","RG4"="#A50026","RG_Cyc"="#FEE090","nIPC1"="#FFFFBF","nIPC2"="#FFFFBF","nIPC3"="#FFFFBF","ImmatureN1"="#C7E9C0","ImmatureN2"="#A1D99B","SCPN"="#41AB5D","SCPN_CSMN"="#238B45","SCPN_CTPN"="#006D2C","CThPN"="#00441B","MigN"="#ABD9E9","superficialCPN"="#74ADD1","Layer4CPN"="#4575B4","deepCPN"="#313695","MGE_IN"="#54278F","CGE_IN"="#6A51A3","CR"="#225EA8","Oligo"="#253494","MG/Peric"="#081D58")

axial_df <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/4_developmental_trajectories.newVersion.cleanedcell/1_URDTree/axial.group.ids.rds")
URDobj <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/4_developmental_trajectories.newVersion.cleanedcell/1_URDTree/URDobj.rds")

URDobj@meta.data$segment <- axial_df[rownames(URDobj@meta.data),"segment"]
URDobj@meta.data$cellPseudotime <- axial_df[rownames(URDobj@meta.data),"cellPseudotime"]

############# SCPN #############
SCPNobj <- subset(URDobj, segment %in% c(11,7))

pdf(file="DimPlot.combine_by_Combine_CellType.pdf", height=10,width=10)
DimPlot(SCPNobj, reduction = "umap_combine", pt.size=1.2, label=TRUE, label.size = 6, group.by='Combine_CellType') + NoLegend() + scale_color_manual(values=c(cell_colors))
dev.off()
saveRDS(SCPNobj, file="SCPNobj.rds")

SCPNobj@active.ident <- as.factor(SCPNobj$Combine_CellType)
cur_markers <- FindAllMarkers(SCPNobj, assays='RNA', slot='data', only.pos =TRUE, logfc.threshold =0.25)
saveRDS(cur_markers, file="SCPNmarkersGene_CellType.rds")


