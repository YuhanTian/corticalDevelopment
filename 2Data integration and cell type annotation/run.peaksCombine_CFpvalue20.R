
# /home/yuchen/miniconda3/envs/R4.0/bin/R

library(glue)
library(Seurat)
library(future)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)
library(cowplot)
library(stringr)
require("GenomicRanges")
require("Signac")
options(future.globals.maxSize = 30 * 1024^3)

#################################################################################################################################################
NPC <- c("SOX2")
Cyc <- c("TOP2A", "MKI67")
RG <- c("SOX9", "HES1")
vRG <- c("FBXO32", "CCN2")
oRG <- c("MOXD1", "HOPX")
earlyRG <- c("NPY", "FGFR3")
lateRG <- c("CD9", "GPX3")
tRG <- c("CRYAB", "NR4A1", "FOXJ1")
mGPC <- c("SOX10", "NKX2-2", "MBP")
Oligo <- c("OLIG1","OLIG2","CLDN11")
Micoglia <- c("AIF1","CD68","PTPRC","SLC2A5","CCL3")
Astrocyte <- c("GFAP","SLC1A3","SP100","SLC1A2")
IPC <- c("EOMES", "PPP1R17", "NEUROG1")
GluN <- c("BCL11B", "SATB2", "SLC17A7","NEUROD2")
SP <- c("NR4A2", "CRYM")
INIPC <- c("DLX1", "DLX2")
IN <- c("GAD1", "GAD2")
MGE <- c("LHX6", "SST")
CGE <- c("SP8", "NR2F2")
PSB <- c("MEIS2", "ETV1")
EC <- c("CLDN5", "PECAM1")
Peric <- c("FOXC2", "PDGFRB")
VLMC <- c("COL1A1", "LUM")
RBC <- c("HEMGN")

#################################################################################################################################################
GW11 <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/9_GW11B3_FU_scMultiome_ATAC_GEX_B1/seurat.integration.rds")
GW15 <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/8_GW15B1_301_scMultiome_ATAC_GEX_B1/seurat.integration.rds")
GW20 <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/4_GW20B1_FU_scMultiome_ATAC_GEX_B1/seurat.integration.rds")

#################################################################################################################################################
## call peaks using MACS2 to produce peaksCommon/peaksUnion/peaksCombine
DefaultAssay(GW11) <- "ATAC"
DefaultAssay(GW15) <- "ATAC"
DefaultAssay(GW20) <- "ATAC"
# GW11_peaks <- CallPeaks(GW11, macs2.path = "/usr/bin/macs2")
# saveRDS(GW11_peaks, file="GW11_peaks.rds")
# GW15_peaks <- CallPeaks(GW15, macs2.path = "/usr/bin/macs2")
# saveRDS(GW15_peaks, file="GW15_peaks.rds")
# GW20_peaks <- CallPeaks(GW20, macs2.path = "/usr/bin/macs2")
# saveRDS(GW20_peaks, file="GW20_peaks.rds")
GW11_peaks <- readRDS("GW11_peaks.rds")
GW15_peaks <- readRDS("GW15_peaks.rds")
GW20_peaks <- readRDS("GW20_peaks.rds")

# peaksUnion
peaksUnion <- union(union(GW11_peaks, GW15_peaks), GW20_peaks)

# peaksCommon
ovGW11 <- findOverlaps(peaksUnion, GW11_peaks)
ovGW15 <- findOverlaps(peaksUnion, GW15_peaks)
ovGW20 <- findOverlaps(peaksUnion, GW20_peaks)
peaksCommon <- peaksUnion[intersect(intersect(queryHits(ovGW11),queryHits(ovGW15)),queryHits(ovGW20))]

# peaksCombine
peaksCombine <- read.table(file="/data3/yuhan/Project/Neuron/scMultiome/ATAC_cov_peak/callingPeak/peaksCombine_CFpvalue30.noBlackList.txt", header=F)
peaksCombine <- read.table(file="/data3/yuhan/Project/Neuron/scMultiome/ATAC_cov_peak/callingPeak/peaksCombine_CFpvalue20.noBlackList.txt", header=F)
peaksCombine <- data.frame(chr=peaksCombine$V1, start=peaksCombine$V2, end=peaksCombine$V3)
peaksCombine <- makeGRangesFromDataFrame(peaksCombine)

#################################################################################################################################################
# quantify counts in each peak
peaksUsed <- peaksCombine #peaksCommon/peaksUnion

macs2_peak_count <- function(obj,peaksUsed){
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(obj),
    features = peaksUsed,
    cells = colnames(obj))
  return(macs2_counts)
}
GW11_macs2_counts <- macs2_peak_count(GW11,peaksUsed)
GW15_macs2_counts <- macs2_peak_count(GW15,peaksUsed)
GW20_macs2_counts <- macs2_peak_count(GW20,peaksUsed)

# create a new assay using the MACS2 peak set and add it to the Seurat object
annotations <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/scMultiome_annotations/annotations.hg38.EnsDb98.rds")
GW11[["peaks"]] <- CreateChromatinAssay(counts = GW11_macs2_counts,
                                        annotation = annotations,
                                        fragments = "/data2/TangLabData/ProcessedData/scMultiome/9_GW11B3_FU_scMultiome_ATAC_GEX_B1/outs/atac_fragments.tsv.gz",
                                        genome = 'hg38')

GW15[["peaks"]] <- CreateChromatinAssay(counts = GW15_macs2_counts,
                                        annotation = annotations,
                                        fragments = "/data2/TangLabData/ProcessedData/scMultiome/8_GW15B1_301_scMultiome_ATAC_GEX_B1/outs/atac_fragments.tsv.gz",
                                        genome = 'hg38')                                             

GW20[["peaks"]] <- CreateChromatinAssay(counts = GW20_macs2_counts,
                                        annotation = annotations,
                                        fragments = "/data2/TangLabData/ProcessedData/scMultiome/4_GW20B1_FU_scMultiome_ATAC_GEX_B1/outs/atac_fragments.tsv.gz",
                                        genome = 'hg38')

scATACseq_process <- function(sampleID){
  processedObj <- get(sampleID)
  # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  library(BSgenome.Hsapiens.UCSC.hg38)
  standard_chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
  idx_standard_chroms <- which(as.character(seqnames(granges(processedObj[["peaks"]]))) %in% standard_chroms)
  processedObj[["peaks"]] <- subset(processedObj[["peaks"]],features = rownames(processedObj[["peaks"]])[idx_standard_chroms])
  seqlevels(processedObj[["peaks"]]@ranges) <- intersect(seqlevels(granges(processedObj[["peaks"]])),unique(seqnames(granges(processedObj[["peaks"]]))))
  # peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  # peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
  
  # running
  DefaultAssay(processedObj) <- "peaks"
  processedObj <- FindTopFeatures(processedObj, min.cutoff = 50)
  processedObj <- RunTFIDF(processedObj, method = 1)
  processedObj <- ScaleData(processedObj,features = row.names(processedObj))
  processedObj <- RunSVD(processedObj, n = 50)
  processedObj <- RunUMAP(processedObj,
                          reduction = "lsi",
                          dims = 2:30,
                          reduction.name = "umap_peaks",
                          reduction.key = "UMAPpeaks_")
  processedObj <- FindNeighbors(object = processedObj, reduction = 'lsi', dims = 2:30)
  processedObj <- FindClusters(object = processedObj, verbose = FALSE, algorithm = 3,resolution = 1.5)
  
  processedObj$sample <- sampleID
  processedObj <- RenameCells(object = processedObj, add.cell.id = sampleID)
  
  return(processedObj)
}
GW11_obj <- scATACseq_process("GW11")
GW15_obj <- scATACseq_process("GW15")
GW20_obj <- scATACseq_process("GW20")
saveRDS(GW11_obj, file="GW11_obj.rds")
saveRDS(GW15_obj, file="GW15_obj.rds")
saveRDS(GW20_obj, file="GW20_obj.rds")
obj_unintegrated = merge(GW11_obj, y = c(GW15_obj, GW20_obj),  merge.data = TRUE)
saveRDS(obj_unintegrated, file="obj.merge.cleanCell.rds")
# GW11_obj <- readRDS("GW11_obj.rds")
# GW15_obj <- readRDS("GW15_obj.rds")
# GW20_obj <- readRDS("GW20_obj.rds")
# obj_unintegrated <- readRDS("obj.merge.cleanCell.rds")

#################################################################################################################################################
# process the unintegrated dataset
DefaultAssay(GW11_obj) <- "peaks"
DefaultAssay(GW15_obj) <- "peaks"
DefaultAssay(GW20_obj) <- "peaks"
DefaultAssay(obj_unintegrated) <- "peaks"

obj_unintegrated <- FindTopFeatures(obj_unintegrated, min.cutoff = 50)
obj_unintegrated <- RunTFIDF(obj_unintegrated, method = 1)
obj_unintegrated <- ScaleData(obj_unintegrated,features = row.names(obj_unintegrated))
obj_unintegrated <- RunSVD(obj_unintegrated, n = 50)
obj_unintegrated <- RunUMAP(obj_unintegrated,
                            reduction = "lsi",
                            dims = 2:30,
                            reduction.name = "umap_unintegrated",
                            reduction.key = "UMAPunintegrated_")
obj_unintegrated <- FindNeighbors(object = obj_unintegrated, reduction = 'lsi', dims = 2:30)
obj_unintegrated <- FindClusters(object = obj_unintegrated, verbose = FALSE, algorithm = 3,resolution = 1.5)

pdf("DimPlot.obj_unintegrated_by_sample.pdf",height=10,width=10)
DimPlot(obj_unintegrated, reduction = "umap_unintegrated", group.by="sample", pt.size=1.2,label=F,label.box=T,repel=T, label.size = 6) + scale_color_manual(values=c("#D52126", "#88CCEE", "#FEE52C"))
dev.off()

pdf("DimPlot.obj_unintegrated_by_cluster.pdf", height=10,width=10)
DimPlot(obj_unintegrated, reduction = "umap_unintegrated", pt.size=1.2, label=TRUE, label.size = 6) + NoLegend()
dev.off()

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(GW11_obj,GW15_obj,GW20_obj),
  anchor.features = rownames(GW11_obj),
  reduction = "rlsi",
  dims = 2:30
)
# integrate LSI embeddings
obj_integratedLSI <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = obj_unintegrated[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)
saveRDS(obj_integratedLSI, file="obj_integratedLSI.rds")
# obj_integratedLSI <- readRDS("obj_integratedLSI.rds")

# create a new UMAP using the integrated embeddings
obj_integrated <- RunUMAP(obj_integratedLSI, reduction = "integrated_lsi", dims = 2:30, reduction.name = "umap_integrated",reduction.key = "UMAPintegrated_")
obj_integrated <- FindNeighbors(object = obj_integrated, reduction = 'integrated_lsi', dims = 2:30)
obj_integrated <- FindClusters(object = obj_integrated, verbose = FALSE, algorithm = 3,resolution = 1.5)
# saveRDS(obj_integrated, file="obj_integratedATAC.rds")
# obj_integrated <- readRDS("obj_integratedATAC.rds")

pdf("DimPlot.obj_integrated_by_sample.pdf",height=10,width=10)
DimPlot(obj_integrated, reduction = "umap_integrated", group.by="sample", pt.size=1.2,label=F,label.box=T,repel=T, label.size = 6) + scale_color_manual(values=c("#D52126", "#88CCEE", "#FEE52C"))
dev.off()

pdf("DimPlot.obj_integrated_by_cluster.pdf", height=10,width=10)
DimPlot(obj_integrated, reduction = "umap_integrated", pt.size=1.2, label=TRUE, label.size = 6) + NoLegend()
dev.off()

DefaultAssay(obj_integrated) <- "GeneActivity"
pdf(file="obj_integrated.FeaturePlots.GeneActivity.scale.pdf", width=23, height=10)
print(FeaturePlot(obj_integrated, slot="scale.data",ncol=4,reduction = "umap_integrated", features = c("EOMES","MKI67","HES1","GFAP","SLC1A3","OLIG1","NEUROD2","DLX2"), pt.size=1.2, order=TRUE, label=TRUE, min.cutoff = "q5", max.cutoff = "q95")+ theme(legend.position = "right"))
dev.off()
pdf(file="obj_integrated.FeaturePlots.GeneActivity.counts.pdf", width=23, height=10)
print(FeaturePlot(obj_integrated, slot="counts",ncol=4,reduction = "umap_integrated", features = c("EOMES","MKI67","HES1","GFAP","SLC1A3","OLIG1","NEUROD2","DLX2"), pt.size=1.2, order=TRUE, label=TRUE, min.cutoff = "q1", max.cutoff = "q99")+ theme(legend.position = "right"))
dev.off()

DefaultAssay(obj_integrated) <- "RNA"
pdf(file="obj_integrated.FeaturePlots.RNA.scale.pdf", width=23, height=10)
print(FeaturePlot(obj_integrated, slot="scale.data",ncol=4,reduction = "umap_integrated", features = c("EOMES","MKI67","HES1","GFAP","SLC1A3","OLIG1","NEUROD2","DLX2"), pt.size=1.2, order=TRUE, label=TRUE, min.cutoff = "q5", max.cutoff = "q95")+ theme(legend.position = "right"))
dev.off()
pdf(file="obj_integrated.FeaturePlots.RNA.counts.pdf", width=23, height=10)
print(FeaturePlot(obj_integrated, slot="counts",ncol=4,reduction = "umap_integrated", features = c("EOMES","MKI67","HES1","GFAP","SLC1A3","OLIG1","NEUROD2","DLX2"), pt.size=1.2, order=TRUE, label=TRUE, min.cutoff = "q1", max.cutoff = "q99")+ theme(legend.position = "right"))
dev.off()

#################################################################################################################################################
DefaultAssay(GW11_obj) <- "RNA"
DefaultAssay(GW15_obj) <- "RNA"
DefaultAssay(GW20_obj) <- "RNA"

pdf(file="GW11.umap_rna.noLegend.pdf", height=10,width=10)
DimPlot(GW11_obj, reduction = "umap_rna_clean2", pt.size=1.2, label=TRUE, label.size = 6, group.by='RNA_snn_res.0.8') + NoLegend()
dev.off()
pdf(file="GW15.umap_rna.noLegend.pdf", height=10,width=10)
DimPlot(GW15_obj, reduction = "umap_rna_clean", pt.size=1.2, label=TRUE, label.size = 6, group.by='RNA_snn_res.0.8') + NoLegend()
dev.off()
pdf(file="GW20.umap_rna.noLegend.pdf", height=10,width=10)
DimPlot(GW20_obj, reduction = "umap_rna", pt.size=1.2, label=TRUE, label.size = 6, group.by='RNA_snn_res.0.8') + NoLegend()
dev.off()

#################################################################################################################################################
obj <- readRDS("obj_integratedATAC.rds")
DefaultAssay(obj) <- "RNA"

# GSE162170_rna_obj <- readRDS("/data3/yuhan/Project/Neuron/GSE162170_rna/GSE162170_rna.obj.rds")
# GSE162170_rna_obj <- subset(GSE162170_rna_obj, subset=seurat_clusters!='c21' & seurat_clusters!='c22')
# objHarmony <- Run_Harmony(obj, VariableFeatures(GSE162170_rna_obj)[1:1000])

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
obj <- ScaleData(obj, features = row.names(obj))
Run_Harmony <- function(obj, used_genes) {
  obj <- RunPCA(object = obj, assay = "RNA", npcs = 50, features=used_genes)
  theta1 = 5
  pc = 50
  obj <- RunHarmony(obj,
                    group.by.vars=c("sample"),
                    assay.use="RNA",
                    reduction="pca",
                    dims.use=1:pc,
                    theta=c(theta1),
                    max.iter.harmony = 50, max.iter.cluster=200,
                    kmeans_init_nstart=2, kmeans_init_iter_max=1000,
                    return_object = TRUE,
                    plot_convergence = FALSE)
  
  obj <- RunUMAP(obj, dims = 1:30, reduction = "harmony", reduction.name = "umap_harmony", reduction.key = "UMAPharmony_")
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
  obj <- FindClusters(obj, verbose = T, resolution=1.3)
  return (obj)
}
objHarmony <- Run_Harmony(obj, VariableFeatures(obj)[1:3000])
# saveRDS(objHarmony, file="obj_integratedATAC.HarmonyRNA.rds")
# objHarmony <- readRDS("obj_integratedATAC.HarmonyRNA.rds")

library(RColorBrewer)
library(cowplot)
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
num <- length(table(Idents(objHarmony)))+1
jBrewColors <- SpatialColors(num)
jBrewColors <- c(jBrewColors[1:ceiling(num/2)], rev(jBrewColors[(ceiling(num/2)+1):num]))
pdf(file="DimPlot_by_cluster.pdf", height=10,width=10)
DimPlot(objHarmony, reduction = "umap_harmony", pt.size=1.2, label=TRUE, cols=jBrewColors, label.size = 6) + NoLegend()
dev.off()

pdf("DimPlot_by_sample.pdf",height=10,width=10)
DimPlot(objHarmony, group.by="sample", reduction = "umap_harmony", pt.size=1.2,label=F,label.box=T,repel=T, label.size = 6) + scale_color_manual(values=c("#D52126", "#88CCEE", "#FEE52C"))
dev.off()

pdf(file="harmony.FeaturePlots.counts.pdf", width=66, height=33)
print(FeaturePlot(objHarmony, slot="counts",reduction = "umap_harmony", ncol=6,features = c("PAX6","MKI67","HES1","OLIG1","PTPRC","SLC1A3","EOMES","NEUROD2","NEUROD6","BCL11B","SATB2","CRYM","NR4A2","GAD1","DLX1","LHX6","NR2F2","PDGFRB"), pt.size=1.2, order=TRUE, label=TRUE, min.cutoff = "q1", max.cutoff = "q99")+ theme(legend.position = "right"))
dev.off()

#################################################################################################################################################
# Bi-modal integrative analysis of the RNA-ATAC scMultiome data, Weighted nearest neighbor analysis
obj <- readRDS("obj_integratedATAC.HarmonyRNA.rds")
obj <- FindMultiModalNeighbors(
  object = obj,
  reduction.list = list("harmony", "integrated_lsi"), 
  dims.list = list(1:30, 2:30),
  modality.weight.name = c("RNA.weight.integration","ATAC.weight.integration"),
  k.nn=49,
  verbose = TRUE
)
obj <- RunUMAP(obj, nn.name = "weighted.nn", assay = "RNA",reduction.name="umap_combine")
obj <- FindClusters(obj, graph.name = "wsnn", resolution = 2.3)
# saveRDS(obj, file="obj_integratedATAC.HarmonyRNA.combine.rds")
# obj <- readRDS("obj_integratedATAC.HarmonyRNA.combine.rds")

library(RColorBrewer)
library(cowplot)
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
num <- length(table(Idents(obj)))+1
jBrewColors <- SpatialColors(num)
jBrewColors <- c(jBrewColors[1:ceiling(num/2)], rev(jBrewColors[(ceiling(num/2)+1):num]))

pdf("DimPlot_by_sample.pdf",height=10,width=10)
DimPlot(obj, group.by="sample", reduction = "umap_combine", pt.size=1.2,label=F,label.box=T,repel=T, label.size = 6) + scale_color_manual(values=c("#D52126", "#88CCEE", "#FEE52C"))
dev.off()

pdf(file="DimPlot_by_cluster.pdf", height=10,width=10)
DimPlot(obj, reduction = "umap_combine", pt.size=1.2, label=TRUE, cols=jBrewColors, label.size = 6) + NoLegend()
dev.off()

pdf(file="umap_combine.FeaturePlots.counts.pdf", width=66, height=33)
print(FeaturePlot(obj, slot="counts",reduction = "umap_combine", ncol=6,features = c("PAX6","MKI67","HES1","OLIG1","PTPRC","SLC1A3","EOMES","NEUROD2","NEUROD6","BCL11B","SATB2","CRYM","NR4A2","GAD1","DLX1","LHX6","NR2F2","PDGFRB"), pt.size=1.2, order=TRUE, label=TRUE, min.cutoff = "q1", max.cutoff = "q99")+ theme(legend.position = "right"))
dev.off()
