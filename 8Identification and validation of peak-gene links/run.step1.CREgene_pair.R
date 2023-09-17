# /home/yuchen/miniconda3/envs/R4.0/bin/R
library(Seurat)
library(Signac)
library(presto)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)

#########################################################################################################################################################
obj <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/2_cell_annotation_newVersion/obj.integratedATAC_HarmonyRNA_combine.cell_annotation_newVersion.cleanedcell.rds")
obj <- RegionStats(obj, assay = "peaks", genome = BSgenome.Hsapiens.UCSC.hg38)
obj <- LinkPeaks(obj,
                 peak.assay = "peaks",
                 expression.assay = "RNA",
                 min.cells =1,
                 distance = 2000000,
                 genes.use = rownames(obj@assays$RNA@counts))
saveRDS(obj, file="obj.CREgene_pair.rds")
# Testing 29254 genes and 101617 peaks
# Found gene coordinates for 17151 genes