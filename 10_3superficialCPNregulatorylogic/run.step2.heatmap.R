cell_colors <- c("RG1"="#FDAE61","RG2"="#F46D43","RG3"="#D73027","RG4"="#A50026","RG_Cyc"="#FEE090","nIPC1"="#FFFFBF","nIPC2"="#FFFFBF","nIPC3"="#FFFFBF","ImmatureN1"="#C7E9C0","ImmatureN2"="#A1D99B","SCPN"="#41AB5D","SCPN_CSMN"="#238B45","SCPN_CTPN"="#006D2C","CThPN"="#00441B","MigN"="#ABD9E9","superficialCPN"="#74ADD1","Layer4CPN"="#4575B4","deepCPN"="#313695","MGE_IN"="#54278F","CGE_IN"="#6A51A3","CR"="#225EA8","Oligo"="#253494","MG/Peric"="#081D58")
marker <- c("TOP2A","MKI67","SOX2","PAX6","HES1","HES5","OLIG1","PTPRC","PDGFRB","GFAP","EOMES","TUBB3","NEUROD2","NEUROD6","BCL11B","LDB2","FEZF2","DIAPH3","FOXP2","TLE4","CRYM","SATB2","NRP1","CUX2","RORB","LMO4","LPL","RELN","GAD1","DLX1","LHX6","NR2F2")

#########################################################################################################################################################
# /home/yuchen/miniconda3/envs/R4.0/bin/R
library(Seurat)
library(Signac)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#########################################################################################################################################################
obj <- readRDS("superficialCPNobj.rds")
cur_markers <- readRDS("superficialCPNmarkersGene_CellType.rds")
markersGene <- unique(cur_markers[which(cur_markers$p_val_adj<0.01 & abs(cur_markers$avg_log2FC)>0.25),"gene"])

LINKobj <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/3_CREgene_pair/obj.CREgene_pair.rds")
genePeak_links <- Links(LINKobj[["peaks"]])
genePeak_links <- genePeak_links[genePeak_links$score>0]
genePeak_links <- genePeak_links[genePeak_links$gene %in% markersGene]
# setdiff(marker,unique(genePeak_links$gene))

#########################################################################################################################################################
metaDF <- obj@meta.data[,c("sample","Combine_CellType","segment","cellPseudotime")]
ordered_cellID <- rownames(metaDF[order(metaDF$cellPseudotime),])

pseudobulkNum <- 53
pseudobulkMem <- 50
pseudobulkCell <- matrix(ordered_cellID[1:(pseudobulkNum*pseudobulkMem)],nrow=pseudobulkNum,ncol=pseudobulkMem,byrow=T)

matrixRNACounts <- obj@assays$RNA@counts
matrixPeakCounts <- obj@assays$peaks@counts

matrixRNACounts_pseudobulk <- data.frame(geneName=rownames(obj@assays$RNA@counts))
for (cellN in c(1:pseudobulkNum)) {
  df <- as.data.frame(rowSums(matrixRNACounts[,pseudobulkCell[cellN,]]))
  colnames(df) <- paste0("pseudobulk",cellN)
  matrixRNACounts_pseudobulk <- cbind(matrixRNACounts_pseudobulk,df)
}
matrixRNACounts_pseudobulk <- matrixRNACounts_pseudobulk[,-1]
matrixRNACounts_pseudobulk <- as.matrix(matrixRNACounts_pseudobulk)
matrixRNACounts_pseudobulk <- as(matrixRNACounts_pseudobulk, 'dgCMatrix')

matrixPeakCounts_pseudobulk <- data.frame(peakCor=rownames(obj@assays$peaks@counts))
for (cellN in c(1:pseudobulkNum)) {
  df <- as.data.frame(rowSums(matrixPeakCounts[,pseudobulkCell[cellN,]]))
  colnames(df) <- paste0("pseudobulk",cellN)
  matrixPeakCounts_pseudobulk <- cbind(matrixPeakCounts_pseudobulk,df)
}
matrixPeakCounts_pseudobulk <- matrixPeakCounts_pseudobulk[,-1]
matrixPeakCounts_pseudobulk <- as.matrix(matrixPeakCounts_pseudobulk)
matrixPeakCounts_pseudobulk <- as(matrixPeakCounts_pseudobulk, 'dgCMatrix')

## Z-score
pseudobulkObj <- CreateSeuratObject(counts = matrixRNACounts_pseudobulk, assay = "RNA")
pseudobulkObj <- NormalizeData(pseudobulkObj)
pseudobulkObj <- ScaleData(pseudobulkObj, features = row.names(pseudobulkObj))

annotations <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/scMultiome_annotations/annotations.hg38.EnsDb98.rds")
pseudobulkObj[["peaks"]] <- CreateChromatinAssay(counts = matrixPeakCounts_pseudobulk,
                                                 annotation = annotations,
                                                 genome = 'hg38')
DefaultAssay(pseudobulkObj) <- "peaks"
pseudobulkObj <- RunTFIDF(pseudobulkObj, method = 1)
pseudobulkObj <- ScaleData(pseudobulkObj,features = row.names(pseudobulkObj))

## filter links
matrixRNACounts <- pseudobulkObj@assays$RNA@counts
matrixPeakCounts <- pseudobulkObj@assays$peaks@counts
matrixRNACounts <- matrixRNACounts[genePeak_links$gene,]
matrixPeakCounts <- matrixPeakCounts[genePeak_links$peak,]
rownames(matrixRNACounts) <- paste0("link_",c(1:length(genePeak_links)))
rownames(matrixPeakCounts) <- paste0("link_",c(1:length(genePeak_links)))
matrixRNACounts_sd <- apply(matrixRNACounts,1,sd)
matrixPeakCounts_sd <- apply(matrixPeakCounts,1,sd)
deleteID <- which(matrixPeakCounts_sd<1)

matrixRNAScale <- pseudobulkObj@assays$RNA@scale.data
matrixPeakScale <- pseudobulkObj@assays$peaks@scale.data
matrixRNAScale <- matrixRNAScale[genePeak_links$gene,]
matrixPeakScale <- matrixPeakScale[genePeak_links$peak,]
rownames(matrixRNAScale) <- paste0("link_",c(1:length(genePeak_links)))
rownames(matrixPeakScale) <- paste0("link_",c(1:length(genePeak_links)))
matrixRNAScale <- matrixRNAScale[-deleteID,]
matrixPeakScale <- matrixPeakScale[-deleteID,]

selectedCell_info <- metaDF[pseudobulkCell[,1],]
rownames(selectedCell_info) <- paste0("pseudobulk",c(1:pseudobulkNum))

pseudotime_colors <- c(colorRampPalette(brewer.pal(9,"Greens"))(pseudobulkNum))
names(pseudotime_colors) <- selectedCell_info$cellPseudotime
col_ha = HeatmapAnnotation(Celltype = selectedCell_info$Combine_CellType,
                           Sample = selectedCell_info$sample,
                           pseudotime = selectedCell_info$cellPseudotime,
                           col = list(Celltype = cell_colors, Sample=c(GW11="#D52126",GW15="#88CCEE",GW20="#FEE52C")),
                           height = unit(0.5, "cm"))

robust_dist = function(x, y) {
  qx = quantile(x, c(0.1, 0.9))
  qy = quantile(y, c(0.1, 0.9))
  l = x > qx[1] & x < qx[2] & y > qy[1] & y < qy[2]
  x = x[l]
  y = y[l]
  sqrt(sum((x - y)^2))
}

Peakheatmap <- ComplexHeatmap::Heatmap(
  matrixPeakScale,
  top_annotation = col_ha,
  name = "Accessibility (Z-score)",
  col = colorRamp2(c(-1, 0, 1), c(brewer.pal(11, 'PRGn' )[10], "white", brewer.pal(11, 'PRGn' )[2])),
  cluster_rows=TRUE,
  cluster_columns=F,
  column_order=rownames(selectedCell_info),
  clustering_method_rows = "ward.D2",
  clustering_distance_rows = robust_dist,
  show_row_names=F,
  show_column_names=F,
  show_row_dend=F,
  show_column_dend=F,
  row_km = 20,
  row_km_repeats = 100,
  use_raster = TRUE,show_heatmap_legend=F
)

RNAheatmap <- ComplexHeatmap::Heatmap(
  matrixRNAScale,
  top_annotation = col_ha,
  name = "Normalized gene expression (Z-score)",
  col = colorRamp2(c(-1, 0, 1), c(brewer.pal(11, 'PRGn' )[10], "white", brewer.pal(11, 'PRGn' )[2])),
  cluster_rows=F,
  cluster_columns=F,
  column_order=rownames(selectedCell_info),
  show_row_names=F,
  show_column_names=F,
  show_row_dend=F,
  show_column_dend=F,
  use_raster = TRUE,show_heatmap_legend=F
)

ht_list <- Peakheatmap + RNAheatmap
ht_list = draw(ht_list)
saveRDS(ht_list, file="superficialCPN_ht_list.rds")

pdf('superficialCPN.pseudobulk.Peak_RNA.km20.heatmap.pdf', width=13, height=10)
ht_list
dev.off()

save.image(file="superficialCPN.pseudobulk.Peak_RNA.km20.heatmap.RData")
# load(file="superficialCPN.pseudobulk.Peak_RNA.km20.heatmap.RData")

#########################################################################################################################################################
ht_list <- readRDS("superficialCPN_ht_list.rds")
RR <- row_order(ht_list)
peakOrder <- c()

for (id in as.character(c(19,15,14,16,13,17,11,10,12,20,4,3,6,7,9,8))) {
  peakOrder <- c(peakOrder,RR[[id]])
}

Peakheatmap <- ComplexHeatmap::Heatmap(
  matrixPeakScale[peakOrder,],
  top_annotation = col_ha,
  name = "Accessibility (Z-score)",
  col = colorRamp2(c(-1, 0, 1), c(brewer.pal(11, 'PRGn' )[10], "white", brewer.pal(11, 'PRGn' )[2])),
  cluster_rows=F,
  cluster_columns=F,
  column_order=rownames(selectedCell_info),
  show_row_names=F,
  show_column_names=F,
  show_row_dend=F,
  show_column_dend=F,
  use_raster = TRUE,show_heatmap_legend=F
)

RNAheatmap <- ComplexHeatmap::Heatmap(
  matrixRNAScale[peakOrder,],
  top_annotation = col_ha,
  name = "Normalized gene expression (Z-score)",
  col = colorRamp2(c(-1, 0, 1), c(brewer.pal(11, 'PRGn' )[10], "white", brewer.pal(11, 'PRGn' )[2])),
  cluster_rows=F,
  cluster_columns=F,
  column_order=rownames(selectedCell_info),
  show_row_names=F,
  show_column_names=F,
  show_row_dend=F,
  show_column_dend=F,
  use_raster = TRUE,show_heatmap_legend=F
)

ht_list_order <- Peakheatmap + RNAheatmap
ht_list_order = draw(ht_list_order)

pdf('superficialCPN_order.pseudobulk.Peak_RNA.km20.heatmap.pdf', width=13, height=10)
ht_list_order
dev.off()
