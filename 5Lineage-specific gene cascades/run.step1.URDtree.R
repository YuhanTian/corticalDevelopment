#################################################################################################################################################
#######################################################################################
# conda create -n R4.2
# source activate R4.2
# conda install -c conda-forge r-base
# conda install udunits2
# conda install -c conda-forge libgit2
# source("https://raw.githubusercontent.com/farrellja/URD/master/URD-Install.R")

#######################################################################################
suppressPackageStartupMessages(library(rgl))
suppressPackageStartupMessages(library(URD))
knitr::opts_chunk$set(echo = TRUE)
rgl::setupKnitr()

cell_type <- c("RG_Cyc","RG1","RG2","RG3","RG4","Oligo","MG/Peric","nIPC1","nIPC2","nIPC3","ImmatureN1","ImmatureN2","SCPN","SCPN_CSMN","SCPN_CTPN","CThPN","MigN","superficialCPN","Layer4CPN","deepCPN","CR","MGE_IN","CGE_IN")
cell_colors <- c("RG1"="#FDAE61","RG2"="#F46D43","RG3"="#D73027","RG4"="#A50026","RG_Cyc"="#FEE090","nIPC1"="#FFFFBF","nIPC2"="#FFFFBF","nIPC3"="#FFFFBF","ImmatureN1"="#C7E9C0","ImmatureN2"="#A1D99B","SCPN"="#41AB5D","SCPN_CSMN"="#238B45","SCPN_CTPN"="#006D2C","CThPN"="#00441B","MigN"="#ABD9E9","superficialCPN"="#74ADD1","Layer4CPN"="#4575B4","deepCPN"="#313695","MGE_IN"="#54278F","CGE_IN"="#6A51A3","CR"="#225EA8","Oligo"="#253494","MG/Peric"="#081D58")
marker <- c("TOP2A","MKI67","SOX2","PAX6","HES1","HES5","OLIG1","PTPRC","PDGFRB","GFAP","EOMES","TUBB3","NEUROD2","NEUROD6","BCL11B","LDB2","FEZF2","DIAPH3","FOXP2","TLE4","CRYM","SATB2","NRP1","CUX2","RORB","LMO4","LPL","RELN","GAD1","DLX1","LHX6","NR2F2")

#######################################################################################
## Create URD object
URDobj <- readRDS("URDobj.rds")
seurat.object <- URDobj
# Copy over data
ds <- new("URD")
ds@logupx.data <- as(as.matrix(seurat.object@assays$RNA@data), "dgCMatrix")
ds@count.data <- as(as.matrix(seurat.object@assays$RNA@counts[rownames(seurat.object@assays$RNA@data), colnames(seurat.object@assays$RNA@data)]), "dgCMatrix")
ds@logupx.data <- ds@logupx.data[-grep("XIST",rownames(ds@logupx.data)),]
ds@count.data <- ds@count.data[-grep("XIST",rownames(ds@count.data)),]
get.data <- as.data.frame(seurat.object@meta.data) 
get.data <- get.data[,c("nCount_RNA","nFeature_RNA","sample","combine_clusters","Combine_CellType")]
ds@meta <- get.data
ds@group.ids <- get.data[,c("sample","combine_clusters","Combine_CellType")]
# Move over tSNE projection
ds@tsne.y <- as.data.frame(seurat.object@reductions$tsne_combine@cell.embeddings)
colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
# Move over PCA results
ds@pca.load <- as.data.frame(seurat.object@reductions$harmony@feature.loadings)
ds@pca.scores <- as.data.frame(seurat.object@reductions$harmony@cell.embeddings)
ds@pca.sdev <- seurat.object@reductions$harmony@stdev
ds@pca.sig <- pcaMarchenkoPastur(M=dim(ds@pca.scores)[1], N=dim(ds@pca.load)[1], pca.sdev=ds@pca.sdev)
axial <- ds
#######################################################################################
## check tSNE (CellType and sample)
pdf(file="tSNE_CellType.pdf", width=10, height=10)
plotDim(axial, "Combine_CellType", plot.title = "tSNE: Stage") + scale_color_manual(values=c(cell_colors))
dev.off()
pdf(file="tSNE_sample.pdf", width=10, height=10)
plotDim(axial, "sample", plot.title = "tSNE: Stage")
dev.off()
#######################################################################################
## Calculate variable genes
# (Normally would do this for each stage, but there are not very many cells in this subset of the data)
# diffCV.cutoff can be varied to include more or fewer genes.
# Copy stage from @meta to @group.ids 
axial@group.ids$stage <- as.character(axial@group.ids$Combine_CellType)
# Find a list of cells from each stage.
stages <- c("RG_Cyc","RG1","RG2","RG3","RG4","nIPC1","nIPC2","nIPC3","ImmatureN1","ImmatureN2","SCPN","SCPN_CSMN","SCPN_CTPN","CThPN","MigN","superficialCPN","Layer4CPN","deepCPN")
cells.each.stage <- lapply(stages, function(stage) rownames(axial@meta)[which(axial@meta$Combine_CellType == stage)])
# Compute variable genes for each stage.
var.by.stage <- lapply(1:length(stages), function(n) findVariableGenes(axial, cells.fit = cells.each.stage[[n]], set.object.var.genes = F, diffCV.cutoff = 0.6, mean.min = 0.005, mean.max = 100, main.use = stages[[n]], do.plot = T))
# Combine the results from each group of stages into a single list of variable genes and load into the URD object
var.genes <- sort(unique(unlist(var.by.stage)))
axial@var.genes <- var.genes
#######################################################################################
## Calculate Diffusion Map
# In this case, knn=100 (larger than sqrt(n.cells)) works well because there are not many cell types.
# Sigma 16 is slightly smaller than the sigma auto-determined by using NULL parameter.
axial <- calcDM(axial, knn = 200, sigma=10)
saveRDS(axial,file="axial_calcDM.rds")
# axial <- readRDS("axial_calcDM.rds")
pdf(file="DiffusionMap.pdf", width=10, height=10)
plotDimArray(axial, reduction.use = "dm", dims.to.plot = 1:8, outer.title = "Diffusion Map (Sigma 16, 100 NNs): Stage", label="stage", plot.title="", legend=F)
dev.off()
pdf(file="tSNE_transitions.pdf", width=10, height=10)
plotDim(axial, "Combine_CellType", transitions.plot = 10000, plot.title="Developmental stage (with transitions)") + scale_color_manual(values=c(cell_colors))
dev.off()
#######################################################################################
## Calculate pseudotime
# Here we use all cells from the first stage as the root
root.cells <- cellsInCluster(axial, "stage", c("RG_Cyc"))
# Then we run 'flood' simulations
axial.floods <- floodPseudotime(axial, root.cells = root.cells, n=50, minimum.cells.flooded = 2, verbose=F)
# The we process the simulations into a pseudotime
axial <- floodPseudotimeProcess(axial, axial.floods, floods.name="pseudotime")
saveRDS(axial, file="axial_floodPseudotimeProcess.rds")
# axial <- readRDS("axial_floodPseudotimeProcess.rds")
pdf(file="Overall_Pseudotime_Stability.pdf", width=10, height=10)
pseudotimePlotStabilityOverall(axial)
dev.off()
pdf(file="tSNE_pseudotime.pdf", width=10, height=10)
plotDim(axial, "pseudotime")
dev.off()
pdf(file="Distribution_of_pseudotime.pdf", width=10, height=10)
plotDists(axial, "pseudotime", "stage", plot.title="Pseudotime by stage") + scale_color_manual(values=c(cell_colors))
dev.off()
#######################################################################################
## Find tips
# Copy cluster identities from axial.6somite object to a new clustering ("tip.clusters") in the full axial object.
axial@group.ids[which(axial@group.ids$stage=="SCPN_CSMN"),"tip.clusters"] <- "1"
axial@group.ids[which(axial@group.ids$stage=="SCPN_CTPN"),"tip.clusters"] <- "2"
axial@group.ids[which(axial@group.ids$stage=="CThPN"),"tip.clusters"] <- "3"
axial@group.ids[which(axial@group.ids$stage=="superficialCPN"),"tip.clusters"] <- "4"
axial@group.ids[which(axial@group.ids$stage=="Layer4CPN"),"tip.clusters"] <- "5"
axial@group.ids[which(axial@group.ids$stage=="deepCPN"),"tip.clusters"] <- "6"
#######################################################################################
## Biased random walks
# Determine the parameters of the logistic used to bias the transition probabilities. The procedure
# is relatively robust to this parameter, but the cell numbers may need to be modified for larger
# or smaller data sets.
axial.ptlogistic <- pseudotimeDetermineLogistic(axial, "pseudotime", optimal.cells.forward=40, max.cells.back=80, do.plot = T)
# Bias the transition matrix acording to pseudotime
axial.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(axial, "pseudotime", logistic.params=axial.ptlogistic))
# Simulate the biased random walks from each tip
axial.walks <- simulateRandomWalksFromTips(axial, tip.group.id="tip.clusters", root.cells=root.cells, transition.matrix = axial.biased.tm, n.per.tip = 100000, root.visits = 1, max.steps = 5000, verbose = F)
# Process the biased random walks into visitation frequencies
axial <- processRandomWalksFromTips(axial, axial.walks, verbose = F)
saveRDS(axial,file="axial_BiasedRandomWalks.rds")
pdf(file="tSNE_tip.pdf", width=10, height=10)
plotDim(axial, "tip.clusters", plot.title="Cells in each tip")
dev.off()
#######################################################################################
## Build tree
# Load the cells used for each tip into the URD object
axial.tree <- loadTipCells(axial, "tip.clusters")
# Build the tree
axial.tree <- buildTree(axial.tree, pseudotime = "pseudotime", save.all.breakpoint.info = T, visit.threshold = 0.7, minimum.visits = 2, bins.per.pseudotime.window = 8, cells.per.pseudotime.bin = 25, divergence.method = "preference", p.thresh = 0.01)
# Name the segments based on our previous determination of the identity of tips 1 and 2.
axial.tree <- nameSegments(axial.tree, segments=c("1","2","3","4","5","6"), segment.names = c("SCPN_CSMN","SCPN_CTPN","CThPN","superficialCPN","Layer4CPN","deepCPN"), short.names = c("SCPN_CSMN","SCPN_CTPN","CThPN","superficialCPN","Layer4CPN","deepCPN"))
saveRDS(axial.tree,file="axial_buildTree.rds")
pdf(file="Tree_CellType.pdf", width=10, height=10)
plotTree(axial.tree, "stage", title="Developmental Stage",discrete.colors=cell_colors,cell.size = 1,legend=F)
dev.off()
pdf(file="Tree_marker.BCL11B.pdf", width=10, height=10)
plotTree(axial.tree, "BCL11B", title="BCL11B")
dev.off()
pdf(file="Tree_segment.pdf", width=10, height=10)
plotTree(axial.tree, "segment", title="URD tree segment")
dev.off()
pdf(file="tSNE_segment.pdf", width=10, height=10)
plotDim(axial.tree, "segment", plot.title="URD tree segment")
dev.off()
#######################################################################################
## Force-directed layout
# Generate the force-directed layout
# axial.tree <- treeForceDirectedLayout(axial.tree, num.nn=100, cut.unconnected.segments=2, verbose=T)
# pdf(file="TreeForce_BCL11B.pdf", width=10, height=10)
# print(plotTreeForce(axial.tree, "BCL11B", title = "BCL11B", title.cex = 2, title.line=2.5))
# dev.off()
