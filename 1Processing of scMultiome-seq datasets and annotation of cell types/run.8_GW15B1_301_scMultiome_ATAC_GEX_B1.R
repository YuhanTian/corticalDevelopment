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
counts_NC8 <- Read10X_h5("/data2/TangLabData/ProcessedData/scMultiome/8_GW15B1_301_scMultiome_ATAC_GEX_B1/outs/filtered_feature_bc_matrix.h5")

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
                                             fragments = "/data2/TangLabData/ProcessedData/scMultiome/8_GW15B1_301_scMultiome_ATAC_GEX_B1/outs/atac_fragments.tsv.gz",
                                             sep = c(":", "-"),
                                             genome = 'hg38')

library(BSgenome.Hsapiens.UCSC.hg38)
standard_chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
idx_standard_chroms <- which(as.character(seqnames(granges(seurat_NC8[['ATAC']]))) %in% standard_chroms)
seurat_NC8[["ATAC"]] <- subset(seurat_NC8[["ATAC"]],features = rownames(seurat_NC8[["ATAC"]])[idx_standard_chroms])
seqlevels(seurat_NC8[['ATAC']]@ranges) <- intersect(seqlevels(granges(seurat_NC8[['ATAC']])),unique(seqnames(granges(seurat_NC8[['ATAC']]))))

### Step 2. Quality control
seurat <- seurat_NC8
DefaultAssay(seurat) <- "ATAC"
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

seurat@meta.data$nCount_RNA_log2 <- log2(seurat@meta.data$nCount_RNA)
pdf("seurat.nCount_RNA_log2.pdf", width=10, height=10)
FeaturePlot(seurat, features="nCount_RNA_log2", pt.size=1.2)+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) +
  theme(legend.position = "right")
dev.off()
# saveRDS(seurat, file="seurat_before_doubletFinder.rds")
# seurat <- readRDS("seurat_before_doubletFinder.rds")

### find doublet
library(glue)
library(Seurat)
library(harmony)
library(future)
library(dplyr)
library(ggplot2)
library(patchwork)
options(future.globals.maxSize = 30 * 1024^3)
library(DoubletFinder)

run_doubletFinder <- function(obj, cur_assum_rate){
  # cur_assum_rate = 30 / 100
  
  print('pK Identification no ground-truth ---------------------------------------------------------------------------------------')
  sweep.res <- paramSweep_v3(obj, PCs = 1:50, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  print('Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------')
  annotations = obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(cur_assum_rate*dim(obj@meta.data)[1])  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  print('Run DoubletFinder with varying classification stringencies')
  obj <- doubletFinder_v3(obj, PCs = 1:50, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  # obj <- doubletFinder_v3(obj, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
  
  meta_colnames <- colnames(obj@meta.data)
  meta_colnames[grepl('DF.classifications', meta_colnames)] <- "isDoublet"
  colnames(obj@meta.data) <- meta_colnames
  
  return(obj)
}

cur_assum_rate = 15 / 100
seurat <- run_doubletFinder(seurat, cur_assum_rate)
table(seurat@meta.data$isDoublet)
# saveRDS(seurat, file="seurat.doubletFinder.rds")
# seurat <- readRDS("seurat.doubletFinder.rds")

pdf('8_GW15B1_301_scMultiome_ATAC_GEX_B1.dimplot.pdf', h=7, w=8)
p1 <- DimPlot(seurat, reduction = "umap_rna", label = TRUE)
print(p1)
dev.off()

pdf('8_GW15B1_301_scMultiome_ATAC_GEX_B1.split.dimplot.doublet.pdf', h=7, w=15)
p1 <- DimPlot(seurat, reduction = "umap_rna", label = TRUE, split.by='isDoublet')
print(p1)
dev.off()

pdf('8_GW15B1_301_scMultiome_ATAC_GEX_B1.dimplot.doublet.pdf', h=7, w=8)
p1 <- DimPlot(seurat, reduction = "umap_rna", label = TRUE, group.by='isDoublet')
print(p1)
dev.off()

pdf("seurat.nCount_RNA.doublet.pdf", width=10, height=10)
FeaturePlot(subset(seurat, subset=isDoublet=="Doublet"), features="nCount_RNA_log2", pt.size=1.2)+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) +
  theme(legend.position = "right")
dev.off()

pdf("seurat.nCount_RNA.singlet.pdf", width=10, height=10)
FeaturePlot(subset(seurat, subset=isDoublet=="Singlet"), features="nCount_RNA_log2", pt.size=1.2)+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) +
  theme(legend.position = "right")
dev.off()

markers <- c("EOMES","MKI67","GFAP","SLC1A3","OLIG1","NEUROD2","DLX2")
for (feature_input_list in c("markers")) {
  for (i in 1:length(get(feature_input_list))){
    feature_input <- get(feature_input_list)[i]
    print(feature_input)
    pdf(file=paste0("umap_rna.FeaturePlots_",feature_input_list,"_",feature_input,".split.scale.pdf"), width=20, height=10)
    print(FeaturePlot(seurat, slot="scale.data",reduction = "umap_rna", split.by="isDoublet", features = feature_input, pt.size=1.2, order=TRUE, label=TRUE, min.cutoff = "q5", max.cutoff = "q95")+ theme(legend.position = "right"))
    dev.off()
  }
}


seurat_clean <- subset(seurat, subset=isDoublet=="Singlet")
#nCount_RNA_q0.05 <- quantile(seurat_clean@meta.data$nCount_RNA, probs=0.05)
#nCount_RNA_q0.95 <- quantile(seurat_clean@meta.data$nCount_RNA, probs=0.95)
#seurat_clean <- subset(seurat_clean, nCount_RNA <= nCount_RNA_q0.05, invert=T)
#seurat_clean <- subset(seurat_clean, nCount_RNA >= nCount_RNA_q0.95, invert=T)

seurat_clean <- FindVariableFeatures(seurat_clean, selection.method = "vst", nfeatures = 3000)
seurat_clean <- ScaleData(seurat_clean,features = row.names(seurat_clean))
seurat_clean <- RunPCA(
  object = seurat_clean,
  assay = "RNA",
  npcs = 50,
  verbose = TRUE)
seurat_clean <- RunUMAP(seurat_clean, dims = 1:50, reduction.name = "umap_rna_clean", reduction.key = "UMAPRNAclean_")
seurat_clean <- FindNeighbors(seurat_clean, dims = 1:50)
seurat_clean <- FindClusters(seurat_clean,resolution = 0.8)
table(Idents(seurat_clean))

library(RColorBrewer)
library(cowplot)
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
num <- length(table(Idents(seurat_clean)))+1
jBrewColors <- SpatialColors(num)
jBrewColors <- c(jBrewColors[1:ceiling(num/2)], rev(jBrewColors[(ceiling(num/2)+1):num]))

pdf(file="umap_rna_clean.noLegend.pdf", height=10,width=10)
DimPlot(seurat_clean, reduction = "umap_rna_clean", pt.size=1.2, label=TRUE, cols=jBrewColors, label.size = 6) + NoLegend()
dev.off()

library(ggplot2)
markers <- c("EOMES","MKI67","GFAP","SLC1A3","OLIG1","NEUROD2","DLX2")
for (feature_input_list in c("markers")) {
  for (i in 1:length(get(feature_input_list))){
    feature_input <- get(feature_input_list)[i]
    print(feature_input)
    pdf(file=paste0("umap_rna_clean.FeaturePlots_",feature_input_list,"_",feature_input,".scale.pdf"), width=10, height=10)
    print(FeaturePlot(seurat_clean, slot="scale.data",reduction = "umap_rna_clean", features = feature_input, pt.size=1.2, order=TRUE, label=TRUE, min.cutoff = "q5", max.cutoff = "q95")+ theme(legend.position = "right"))
    dev.off()
    
    # pdf(file=paste0("umap_rna_clean.FeaturePlots_",feature_input_list,"_",feature_input,".counts.pdf"), width=10, height=10)
    # print(FeaturePlot(seurat_clean, slot="counts",reduction = "umap_rna_clean", features = feature_input, pt.size=1.2, order=TRUE, label=TRUE, min.cutoff = "q1", max.cutoff = "q99")+ theme(legend.position = "right"))
    # dev.off()
  }
}

# saveRDS(seurat_clean, file="seurat.doubletFinder.final.rds")
# seurat_clean <- readRDS("seurat.doubletFinder.final.rds")

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
seurat <- seurat_clean
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
                                  reduction.list = list("umap_rna_clean", "umap_atac"),
                                  dims.list = list(1:ncol(Embeddings(seurat,"umap_rna_clean")),
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
