# /home/yuchen/miniconda3/envs/R4.0/bin/R

#### Section 1. Mono-modal data analysis of the scMultiome data
### Step 0. Import the required packages
library(Seurat)
library(Signac)
library(Matrix)
library(dplyr)
require("AnnotationHub")

#########################################################################################################################################################
### Step 1. Load the data and create the Seurat object
# ah <- AnnotationHub()
# ensdbs <- query(ah, c("EnsDb.Hsapiens"))
# ensdb_id <- ensdbs$ah_id[grep(paste0(" 98 EnsDb"), ensdbs$title)]
# ensdb <- ensdbs[[ensdb_id]]
# seqlevelsStyle(ensdb) <- "UCSC"
# annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
# genome(annotations) <- "hg38"
# saveRDS(annotations, file="annotations.hg38.EnsDb98.rds")

annotations <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/scMultiome_annotations/annotations.hg38.EnsDb98.rds")
counts_NC8 <- Read10X_h5("/data2/TangLabData/ProcessedData/scMultiome/4_GW20B1_FU_scMultiome_ATAC_GEX_B1/outs/filtered_feature_bc_matrix.h5")

# construct RNA assay and filter mt genes
seurat_NC8 <- CreateSeuratObject(counts = counts_NC8$`Gene Expression`,assay = "RNA")
seurat_NC8 <- PercentageFeatureSet(seurat_NC8, pattern = "^MT-", col.name = "percent.mt", assay = "RNA")
genes_use <- rownames(seurat_NC8[['RNA']]@counts)
genes_use <- genes_use[which(!grepl('^MT-', genes_use))]
genes_use <- genes_use[which(!grepl('^RPS', genes_use))]
genes_use <- genes_use[which(!grepl('^RPL', genes_use))]
seurat_NC8 <- subset(seurat_NC8, features=genes_use)

# construct ATAC assay
seurat_NC8[['ATAC']] <- CreateChromatinAssay(counts = counts_NC8$`Peaks`,
                                             annotation = annotations,
                                             fragments = "/data2/TangLabData/ProcessedData/scMultiome/4_GW20B1_FU_scMultiome_ATAC_GEX_B1/outs/atac_fragments.tsv.gz",
                                             sep = c(":", "-"),
                                             genome = 'hg38')

library(BSgenome.Hsapiens.UCSC.hg38)
standard_chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
idx_standard_chroms <- which(as.character(seqnames(granges(seurat_NC8[['ATAC']]))) %in% standard_chroms)
seurat_NC8[["ATAC"]] <- subset(seurat_NC8[["ATAC"]],features = rownames(seurat_NC8[["ATAC"]])[idx_standard_chroms])
seqlevels(seurat_NC8[['ATAC']]@ranges) <- intersect(seqlevels(granges(seurat_NC8[['ATAC']])),unique(seqnames(granges(seurat_NC8[['ATAC']]))))

### Step 2. Quality control
seurat <- seurat_NC8
seurat <- NucleosomeSignal(seurat, assay = "ATAC")
seurat <- TSSEnrichment(seurat, assay = "ATAC")
# plot
pdf(file=paste0("QC_metrics_raw.pdf"), height=10, width=50)
VlnPlot(seurat,ncol = 5,pt.size = 0,
        features = c("nFeature_RNA",
                     "percent.mt",
                     "nFeature_ATAC",
                     "TSS.enrichment",
                     "nucleosome_signal"))
dev.off()
# saveRDS(seurat, file="seurat_QualityControl.rds")
# seurat <- readRDS("seurat_QualityControl.rds")

# filter cell with high percentage of mt reads and outline
seurat <- subset(seurat,
                 subset = nFeature_RNA > 1000 &
                   nFeature_RNA < 7500 &
                   percent.mt < 30 &
                   nFeature_ATAC > 1000 &
                   nFeature_ATAC < 30000 &
                   TSS.enrichment > 1 &
                   nucleosome_signal < 2
)

### Step 3. Analysis on the RNA assay
DefaultAssay(seurat) <- "RNA"
seurat <- NormalizeData(seurat) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  CellCycleScoring(s.features = cc.genes.updated.2019$s.genes,
                   g2m.features = cc.genes.updated.2019$g2m.genes) %>%
  ScaleData(features = row.names(seurat)) %>%
  RunPCA(assay = "RNA",npcs = 50,verbose = TRUE) %>%
  RunUMAP(dims = 1:20, reduction.name = "umap_rna", reduction.key = "UMAPRNA_")
seurat <- FindNeighbors(seurat, dims = 1:50) %>%
  FindClusters(resolution = 0.8)

library(RColorBrewer)
library(cowplot)
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
num <- length(table(Idents(seurat)))+1
jBrewColors <- SpatialColors(num)
jBrewColors <- c(jBrewColors[1:ceiling(num/2)], rev(jBrewColors[(ceiling(num/2)+1):num]))

pdf(file="umap_rna.noLegend.pdf", height=10,width=10)
DimPlot(seurat, reduction = "umap_rna", pt.size=1.2, label=TRUE, cols=jBrewColors, label.size = 6) + NoLegend()
dev.off()

library(ggplot2)
markers <- c("EOMES","MKI67","GFAP","SLC1A3","OLIG1","NEUROD2","DLX2")
for (feature_input_list in c("markers")) {
  for (i in 1:length(get(feature_input_list))){
    feature_input <- get(feature_input_list)[i]
    print(feature_input)
    pdf(file=paste0("umap_rna.FeaturePlots_",feature_input_list,"_",feature_input,".scale.pdf"), width=10, height=10)
    print(FeaturePlot(seurat, slot="scale.data",reduction = "umap_rna", features = feature_input, pt.size=1.2, order=TRUE, label=TRUE, min.cutoff = "q5", max.cutoff = "q95")+ theme(legend.position = "right"))
    dev.off()
    
    pdf(file=paste0("umap_rna.FeaturePlots_",feature_input_list,"_",feature_input,".counts.pdf"), width=10, height=10)
    print(FeaturePlot(seurat, slot="counts",reduction = "umap_rna", features = feature_input, pt.size=1.2, order=TRUE, label=TRUE, min.cutoff = "q1", max.cutoff = "q99")+ theme(legend.position = "right"))
    dev.off()
  }
}

# finding top markers
# library(presto)
# DE_cl_rna <- presto::wilcoxauc(seurat, "RNA_snn_res.0.8")
# top_markers <- DE_cl_rna %>%
#   filter(logFC > log(1.2) &
#          auc > 0.7 &
#          padj < 0.01 &
#          pct_in - pct_out > 30 &
#          pct_out < 30) %>%
#   group_by(group) %>%
#   top_n(1, wt = auc)
# pdf(file="umap_rna.featurePlot.top_markers.pdf", height=30,width=30)
# FeaturePlot(seurat,features = unique(top_markers$feature),reduction="umap_rna",order = T,ncol=3) & NoAxes() & NoLegend()
# dev.off()

### Step 4. Analysis on the ATAC assay
# Step 4.1. Feature selection
DefaultAssay(seurat) <- "ATAC"
seurat <- FindTopFeatures(seurat, min.cutoff = 50)
# Step 4.2. Normalization
seurat <- RunTFIDF(seurat, method = 1)
seurat <- ScaleData(seurat,features = row.names(seurat))
# Step 4.3. Linear dimension reduction
seurat <- RunSVD(seurat, n = 50)
# Step 4.4. Non-linear dimension reduction with UMAP for visualization
seurat <- RunUMAP(seurat,
                  reduction = "lsi",
                  dims = 2:30,
                  reduction.name = "umap_atac",
                  reduction.key = "UMAPATAC_")
seurat <- FindNeighbors(object = seurat, reduction = 'lsi', dims = 2:30)
seurat <- FindClusters(object = seurat, verbose = FALSE, algorithm = 3,resolution = 1.5)

library(RColorBrewer)
library(cowplot)
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
num <- length(table(Idents(seurat)))+1
jBrewColors <- SpatialColors(num)
jBrewColors <- c(jBrewColors[1:ceiling(num/2)], rev(jBrewColors[(ceiling(num/2)+1):num]))

pdf(file="umap_atac.ATACcluster.noLegend.pdf", height=10,width=10)
DimPlot(seurat, reduction = "umap_atac", pt.size=1.2, label=TRUE, cols=jBrewColors, label.size = 6) + NoLegend()
dev.off()

# add the gene activity matrix to the Seurat object as a new assay and normalize it
gene.activities <- GeneActivity(seurat)
seurat[['GeneActivity']] <- CreateAssayObject(counts = gene.activities)
seurat <- NormalizeData(
  object = seurat,
  assay = 'GeneActivity',
  normalization.method = 'LogNormalize',
  scale.factor = median(seurat$nCount_ATAC)
)

DefaultAssay(seurat) <- 'GeneActivity'
seurat <- ScaleData(seurat,features = row.names(seurat))

markers <- c("EOMES","MKI67","GFAP","SLC1A3","OLIG1","NEUROD2","DLX2")
for (feature_input_list in c("markers")) {
  for (i in 1:length(get(feature_input_list))){
    feature_input <- get(feature_input_list)[i]
    print(feature_input)
    pdf(file=paste0("umap_atac.FeaturePlots_",feature_input_list,"_",feature_input,".scale.pdf"), width=10, height=10)
    print(FeaturePlot(seurat, slot="scale.data",reduction = "umap_atac", features = feature_input, pt.size=1.2, order=TRUE, label=TRUE, min.cutoff = "q5", max.cutoff = "q95")+ theme(legend.position = "right"))
    dev.off()
    
    # pdf(file=paste0("umap_atac.FeaturePlots_",feature_input_list,"_",feature_input,".counts.pdf"), width=10, height=10)
    # print(FeaturePlot(seurat, slot="counts",reduction = "umap_atac", features = feature_input, pt.size=1.2, order=TRUE, label=TRUE, min.cutoff = "q1", max.cutoff = "q99")+ theme(legend.position = "right"))
    # dev.off()
  }
}

# saveRDS(seurat, file="seurat.Section1.rds")
# seurat <- readRDS("seurat.Section1.rds")

#########################################################################################################################################################
#### Section 2. Bi-modal integrative analysis of the RNA-ATAC scMultiome data
# Step 1. Weighted nearest neighbor analysis
seurat <- FindMultiModalNeighbors(seurat,
                                  reduction.list = list("umap_rna", "umap_atac"),
                                  dims.list = list(1:ncol(Embeddings(seurat,"umap_rna")),
                                                   1:ncol(Embeddings(seurat,"umap_atac"))),
                                  modality.weight.name = c("RNA.weight","ATAC.weight"),
                                  verbose = TRUE)
seurat <- RunUMAP(seurat, nn.name = "weighted.nn", assay = "RNA",reduction.name="umap_combine")
seurat <- FindClusters(seurat, graph.name = "wsnn", resolution = 1.1)

library(RColorBrewer)
library(cowplot)
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
num <- length(table(Idents(seurat)))+1
jBrewColors <- SpatialColors(num)
jBrewColors <- c(jBrewColors[1:ceiling(num/2)], rev(jBrewColors[(ceiling(num/2)+1):num]))

pdf(file="umap_combine.noLegend.pdf", height=10,width=10)
DimPlot(seurat, reduction = "umap_combine", pt.size=1.2, label=TRUE, cols=jBrewColors, label.size = 6) + NoLegend()
dev.off()

DefaultAssay(seurat) <- "RNA"
markers <- c("EOMES","MKI67","GFAP","SLC1A3","OLIG1","NEUROD2","DLX2")
for (feature_input_list in c("markers")) {
  for (i in 1:length(get(feature_input_list))){
    feature_input <- get(feature_input_list)[i]
    print(feature_input)
    pdf(file=paste0("umap_combine.FeaturePlots_",feature_input_list,"_",feature_input,".scale.pdf"), width=10, height=10)
    print(FeaturePlot(seurat, slot="scale.data",reduction = "umap_combine", features = feature_input, pt.size=1.2, order=TRUE, label=TRUE, min.cutoff = "q5", max.cutoff = "q95")+ theme(legend.position = "right"))
    dev.off()
    
    # pdf(file=paste0("umap_combine.FeaturePlots_",feature_input_list,"_",feature_input,".counts.pdf"), width=10, height=10)
    # print(FeaturePlot(seurat, slot="counts",reduction = "umap_combine", features = feature_input, pt.size=1.2, order=TRUE, label=TRUE, min.cutoff = "q1", max.cutoff = "q99")+ theme(legend.position = "right"))
    # dev.off()
  }
}

# saveRDS(seurat, file="seurat.integration.rds")
# seurat <- readRDS("seurat.integration.rds")

## RNA_CellType
seurat$RNA_CellType <- NA
seurat$RNA_CellType[seurat$RNA_snn_res.0.8 == "0"] <- "EN2"
seurat$RNA_CellType[seurat$RNA_snn_res.0.8 == "1"] <- "EN3"
seurat$RNA_CellType[seurat$RNA_snn_res.0.8 == "2"] <- "EN1"
seurat$RNA_CellType[seurat$RNA_snn_res.0.8 == "3"] <- "EN4"
seurat$RNA_CellType[seurat$RNA_snn_res.0.8 == "4"] <- "IN1"
seurat$RNA_CellType[seurat$RNA_snn_res.0.8 == "5"] <- "EN5"
seurat$RNA_CellType[seurat$RNA_snn_res.0.8 == "6"] <- "EN8"
seurat$RNA_CellType[seurat$RNA_snn_res.0.8 == "7"] <- "EN10"
seurat$RNA_CellType[seurat$RNA_snn_res.0.8 == "8"] <- "Astrocyte"
seurat$RNA_CellType[seurat$RNA_snn_res.0.8 == "9"] <- "IN2"
seurat$RNA_CellType[seurat$RNA_snn_res.0.8 == "10"] <- "nIPC"
seurat$RNA_CellType[seurat$RNA_snn_res.0.8 == "11"] <- "EN9"
seurat$RNA_CellType[seurat$RNA_snn_res.0.8 == "12"] <- "EN6"
seurat$RNA_CellType[seurat$RNA_snn_res.0.8 == "13"] <- "EN7"
seurat$RNA_CellType[seurat$RNA_snn_res.0.8 == "14"] <- "RG"
seurat$RNA_CellType[seurat$RNA_snn_res.0.8 == "15"] <- "Oligo"
as.data.frame(table(seurat$RNA_CellType))

pdf(file="umap_rna.RNA_CellType.pdf", height=10,width=10)
DimPlot(seurat, group.by="RNA_CellType", reduction = "umap_rna", pt.size=1.2, label=TRUE, label.size = 6) + NoLegend() + scale_color_manual(values=c("IN1"="#FF69B4","IN2"="#FF00FF","EN1"="#2c86ef","EN2"="#0a30c6","EN3"="#006281","EN4"="#005c3b","EN5"="#256e00","EN6"="#008044","EN7"="#059876","EN8"="#589400","EN9"="#228B22","EN10"="#008000","nIPC"="#8B008B","RG"="#FFD700","Oligo"="#FFA500","Astrocyte"="#D2691E"))
dev.off()

## ATAC_CellType
seurat$ATAC_CellType <- NA
seurat$ATAC_CellType[seurat$ATAC_snn_res.1.5 == "0"] <- "EN1"
seurat$ATAC_CellType[seurat$ATAC_snn_res.1.5 == "1"] <- "EN2"
seurat$ATAC_CellType[seurat$ATAC_snn_res.1.5 == "2"] <- "EN6"
seurat$ATAC_CellType[seurat$ATAC_snn_res.1.5 == "3"] <- "EN3"
seurat$ATAC_CellType[seurat$ATAC_snn_res.1.5 == "4"] <- "IN1"
seurat$ATAC_CellType[seurat$ATAC_snn_res.1.5 == "5"] <- "nIPC"
seurat$ATAC_CellType[seurat$ATAC_snn_res.1.5 == "6"] <- "EN9"
seurat$ATAC_CellType[seurat$ATAC_snn_res.1.5 == "7"] <- "EN5"
seurat$ATAC_CellType[seurat$ATAC_snn_res.1.5 == "8"] <- "Astrocyte"
seurat$ATAC_CellType[seurat$ATAC_snn_res.1.5 == "9"] <- "EN7"
seurat$ATAC_CellType[seurat$ATAC_snn_res.1.5 == "10"] <- "IN2"
seurat$ATAC_CellType[seurat$ATAC_snn_res.1.5 == "11"] <- "EN4"
seurat$ATAC_CellType[seurat$ATAC_snn_res.1.5 == "12"] <- "EN10"
seurat$ATAC_CellType[seurat$ATAC_snn_res.1.5 == "13"] <- "EN8"
seurat$ATAC_CellType[seurat$ATAC_snn_res.1.5 == "14"] <- "Oligo"
seurat$ATAC_CellType[seurat$ATAC_snn_res.1.5 == "15"] <- "N"
as.data.frame(table(seurat$ATAC_CellType))

pdf(file="umap_atac.ATAC_CellType.pdf", height=10,width=10)
DimPlot(seurat, group.by="ATAC_CellType", reduction = "umap_atac", pt.size=1.2, label=TRUE, label.size = 6) + NoLegend() + scale_color_manual(values=c("IN1"="#FF69B4","IN2"="#FF00FF","EN1"="#2c86ef","EN2"="#0a30c6","EN3"="#006281","EN4"="#005c3b","EN5"="#256e00","EN6"="#008044","EN7"="#059876","EN8"="#589400","EN9"="#228B22","EN10"="#008000","nIPC"="#8B008B","Oligo"="#FFA500","Astrocyte"="#D2691E","N"="#A9A9A9"))
dev.off()

## Combine_CellType
seurat$Combine_CellType <- NA
seurat$Combine_CellType[seurat$wsnn_res.1.1 == "0"] <- "EN7"
seurat$Combine_CellType[seurat$wsnn_res.1.1 == "1"] <- "EN5"
seurat$Combine_CellType[seurat$wsnn_res.1.1 == "2"] <- "EN4"
seurat$Combine_CellType[seurat$wsnn_res.1.1 == "3"] <- "EN9"
seurat$Combine_CellType[seurat$wsnn_res.1.1 == "4"] <- "EN2"
seurat$Combine_CellType[seurat$wsnn_res.1.1 == "5"] <- "IN1"
seurat$Combine_CellType[seurat$wsnn_res.1.1 == "6"] <- "EN8"
seurat$Combine_CellType[seurat$wsnn_res.1.1 == "7"] <- "Astrocyte2"
seurat$Combine_CellType[seurat$wsnn_res.1.1 == "8"] <- "IN2"
seurat$Combine_CellType[seurat$wsnn_res.1.1 == "9"] <- "EN10"
seurat$Combine_CellType[seurat$wsnn_res.1.1 == "10"] <- "EN1"
seurat$Combine_CellType[seurat$wsnn_res.1.1 == "11"] <- "RG"
seurat$Combine_CellType[seurat$wsnn_res.1.1 == "12"] <- "EN6"
seurat$Combine_CellType[seurat$wsnn_res.1.1 == "13"] <- "EN3"
seurat$Combine_CellType[seurat$wsnn_res.1.1 == "14"] <- "nIPC"
seurat$Combine_CellType[seurat$wsnn_res.1.1 == "15"] <- "Astrocyte1"
seurat$Combine_CellType[seurat$wsnn_res.1.1 == "16"] <- "Oligo"
as.data.frame(table(seurat$Combine_CellType))

pdf(file="umap_combine.Combine_CellType.pdf", height=10,width=10)
DimPlot(seurat, group.by="Combine_CellType", reduction = "umap_combine", pt.size=1.2, label=TRUE, label.size = 6) + NoLegend() + scale_color_manual(values=c("IN1"="#FF69B4","IN2"="#FF00FF","EN1"="#2c86ef","EN2"="#0a30c6","EN3"="#006281","EN4"="#005c3b","EN5"="#256e00","EN6"="#008044","EN7"="#059876","EN8"="#589400","EN9"="#228B22","EN10"="#008000","nIPC"="#8B008B","RG"="#FFD700","Oligo"="#FFA500","Astrocyte1"="#F4A460","Astrocyte2"="#D2691E"))
dev.off()

# saveRDS(seurat, file="seurat.CellType.rds")
# seurat <- readRDS("seurat.CellType.rds")

#########################################################################################################################################################
# Step 2. Cell type gene/peak marker identification and visualization of the chromatin accessibility profiles
library(presto)

DefaultAssay(seurat) <- "RNA"
DE_ct <- wilcoxauc(seurat, "Combine_CellType", seurat_assay = "RNA")
top_markers_ct <- DE_ct %>%
  filter(abs(logFC) > log(1.2) &
           padj < 0.01 &
           auc > 0.65 &
           pct_in - pct_out > 30 &
           pct_out < 20) %>%
  group_by(group) %>%
  top_n(10, wt = auc)
top_markers_ct

DefaultAssay(seurat) <- "ATAC"
DA_ct <- wilcoxauc(seurat, "Combine_CellType", seurat_assay = "ATAC")
top_peaks_ct <- DA_ct %>%
  filter(abs(logFC) > log(1.1) &
           padj < 0.01 &
           auc > 0.55) %>%
  group_by(group) %>%
  top_n(100, wt = auc)
marker_peak_ct <- top_peaks_ct %>% top_n(5, wt=auc)

library(BSgenome.Hsapiens.UCSC.hg38)
seurat <- RegionStats(seurat,
                      genome = BSgenome.Hsapiens.UCSC.hg38)
seurat <- LinkPeaks(seurat,
                    peak.assay = "ATAC",
                    expression.assay = "RNA",
                    genes.use = top_markers_ct$feature)

DefaultAssay(seurat) <- "ATAC"
TopMarkerGenesID <- "TRPC4"
pdf(file=paste0("TopMarkerGenes.",TopMarkerGenesID,".pdf"), height=10,width=10)
CoveragePlot(seurat,
             region = TopMarkerGenesID,
             features = TopMarkerGenesID,
             group.by = "Combine_CellType",
             extend.upstream = 1000,
             extend.downstream = 1000)
dev.off()

TopMarkerPeaksID <- "chr1-34842652-34844740"
pdf(file=paste0("TopMarkerPeaks.",TopMarkerPeaksID,".pdf"), height=10,width=10)
CoveragePlot(seurat,
             region = TopMarkerPeaksID,
             group.by = "Combine_CellType",
             extend.upstream = 1000,
             extend.downstream = 1000)
dev.off()

# saveRDS(seurat, file="seurat.LinkPeaks.rds")
# seurat <- readRDS("seurat.LinkPeaks.rds")
# saveRDS(marker_peak_ct, file="marker_peak_ct.rds")
# marker_peak_ct <- readRDS("marker_peak_ct.rds")
# saveRDS(top_peaks_ct, file="top_peaks_ct.rds")
# top_peaks_ct <- readRDS("top_peaks_ct.rds")


# Step 3. TF binding motif enrichment analysis
conda activate scMultiome
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(JASPAR2020)
require("motifmatchr")
require('ggseqlogo')
require("ggplot2")
seurat <- readRDS("seurat.LinkPeaks.rds")
marker_peak_ct <- readRDS("marker_peak_ct.rds")
top_peaks_ct <- readRDS("top_peaks_ct.rds")

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
df_pfm <- data.frame(t(sapply(pfm, function(x)
  c(id=x@ID, name=x@name, symbol=ifelse(!is.null(x@tags$symbol),x@tags$symbol,NA)))))
seurat <- AddMotifs(seurat, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)

open_peaks <- AccessiblePeaks(seurat)
peaks_matched <- MatchRegionStats(meta.feature = seurat[['ATAC']]@meta.features[open_peaks, ],
                                  query.feature = seurat[['ATAC']]@meta.features[marker_peak_ct$feature, ],
                                  n = 50000)

motif_enrichment_RG <- FindMotifs(seurat,
                                    features = top_peaks_ct$feature[top_peaks_ct$group == "RG"],
                                    background = peaks_matched) %>%
  mutate(symbol = setNames(ifelse(is.na(df_pfm$symbol), df_pfm$name, df_pfm$symbol), df_pfm$id)[motif]) %>%
  mutate(padj = p.adjust(pvalue, method="BH"))
enriched_motif_RG <- motif_enrichment_RG %>%
  filter(padj < 0.01 & fold.enrichment > 3) %>%
  top_n(-4, wt = padj)

pdf(file="enriched_motif_RG.pdf", height=3,width=10)
MotifPlot(seurat, motifs = enriched_motif_RG$motif[1:4], ncol=4)
dev.off()

DefaultAssay(seurat) <- "RNA"
pdf(file="umap_combine.enriched_motif_RG.pdf", height=10,width=40)
print(FeaturePlot(seurat, slot="counts",reduction = "umap_combine", features = c("EMX2","GBX1","LHX1","EN2"), pt.size=1.2, order=TRUE,ncol=4, label=TRUE, 
min.cutoff = "q1", max.cutoff = "q99")+ theme(legend.position = "right"))
dev.off()

# Step 4. ChromVAR: another motif enrichment analysis
require(chromVAR)

DefaultAssay(seurat) <- "ATAC"
seurat <- RunChromVAR(seurat, genome = BSgenome.Hsapiens.UCSC.hg38)

DefaultAssay(seurat) <- "chromvar"
DA_motifs_ct <- wilcoxauc(seurat, group_by = "Combine_CellType", seurat_assay = "chromvar") %>%
  mutate(symbol = setNames(ifelse(is.na(df_pfm$symbol), df_pfm$name, df_pfm$symbol),
                           df_pfm$id)[feature])

enriched_motifs_ct <- DA_motifs_ct %>%
  filter(padj < 0.01 & auc > 0.7) %>%
  group_by(group)
top_motifs_ct <- top_n(enriched_motifs_ct, 3, wt=auc)

bluered_colscheme <- colorRampPalette(rev(c("#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1","#4575b4")))
FeaturePlot(seurat,
            features = top_motifs_ct$feature,
            cols = bluered_colscheme(30),
            reduction = "umap_combine",
            ncol = 3) & NoAxes() & NoLegend()


tfs <- read.table("extdata/Homo_sapiens_TF.txt", sep="\t", header=T)

tf_motifs_ct <- enriched_motifs_ct %>%
  filter(symbol %in% tfs$Symbol)
marker_tfs_ct <- DE_ct %>%
  filter(feature %in% tfs$Symbol &
           abs(logFC) > log(1.2) &
           padj < 0.01 &
           auc > 0.65 &
           pct_in - pct_out > 20) %>%
  inner_join(tf_motifs_ct,
             by = c("feature" = "symbol"),
             suffix = c("_tf","_motif")) %>%
  filter(group_tf == group_motif)

top_tfs_ct <- group_by(marker_tfs_ct, group_tf) %>%
  top_n(3, wt = auc_motif)

beach_colscheme <- colorRampPalette(c("#cdcdcd","#edf8b1","#7fcdbb","#41b6c4","#1d91c0","#225ea8","#0c2c84"))
DefaultAssay(seurat) <- "RNA"
p1 <- FeaturePlot(seurat,
                  top_tfs_ct$feature,
                  reduction = "umap",
                  order=T,
                  cols=beach_colscheme(30),
                  ncol=6) & NoAxes() & NoLegend()
DefaultAssay(seurat) <- "chromvar"
p2 <- FeaturePlot(seurat,
                  top_tfs_ct$feature_motif,
                  reduction = "umap",
                  order=T,
                  cols=bluered_colscheme(30),
                  ncol=6) & NoAxes() & NoLegend()


# Section 3. Gene regulatory network reconstruction

