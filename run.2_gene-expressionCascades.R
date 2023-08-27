#### Gene-expression cascades ####

#################################################################################################################################################
# conda create -n R4.2
# source activate R4.2
# conda install -c conda-forge r-base
# conda install udunits2
# conda install -c conda-forge libgit2
# source("https://raw.githubusercontent.com/farrellja/URD/master/URD-Install.R")
suppressPackageStartupMessages(library(rgl))
suppressPackageStartupMessages(library(URD))
knitr::opts_chunk$set(echo = TRUE)
rgl::setupKnitr()
# obj <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/2_cell_annotation_newVersion/obj.integratedATAC_HarmonyRNA_combine.cell_annotation_newVersion.cleanedcell.rds")
cell_type <- c("RG_Cyc","RG1","RG2","RG3","RG4","Oligo","MG/Peric","nIPC1","nIPC2","nIPC3","ImmatureN1","ImmatureN2","SCPN","SCPN_CSMN","SCPN_CTPN","CThPN","MigN","superficialCPN","Layer4CPN","deepCPN","CR","MGE_IN","CGE_IN")
cell_colors <- c("RG1"="#FDAE61","RG2"="#F46D43","RG3"="#D73027","RG4"="#A50026","RG_Cyc"="#FEE090","nIPC1"="#FFFFBF","nIPC2"="#FFFFBF","nIPC3"="#FFFFBF","ImmatureN1"="#C7E9C0","ImmatureN2"="#A1D99B","SCPN"="#41AB5D","SCPN_CSMN"="#238B45","SCPN_CTPN"="#006D2C","CThPN"="#00441B","MigN"="#ABD9E9","superficialCPN"="#74ADD1","Layer4CPN"="#4575B4","deepCPN"="#313695","MGE_IN"="#54278F","CGE_IN"="#6A51A3","CR"="#225EA8","Oligo"="#253494","MG/Peric"="#081D58")
marker <- c("TOP2A","MKI67","SOX2","PAX6","HES1","HES5","OLIG1","PTPRC","PDGFRB","GFAP","EOMES","TUBB3","NEUROD2","NEUROD6","BCL11B","LDB2","FEZF2","DIAPH3","FOXP2","TLE4","CRYM","SATB2","NRP1","CUX2","RORB","LMO4","LPL","RELN","GAD1","DLX1","LHX6","NR2F2")
axial <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/4_developmental_trajectories.newVersion.cleanedcell/1_URDTree/axial_buildTree.rds")

###################################################################
## Calculate the markers of each other population.
tips.to.run <- as.character(axial@tree$tips)
gene.markers <- list()
for (tipn in 1:length(tips.to.run)) {
  tip <- tips.to.run[tipn]
  print(paste0(Sys.time(), ": ", tip))
  markers <- aucprTestAlongTree(
    axial,
    pseudotime="pseudotime", 
    tips=tip,
    log.effect.size = 0.25,
    auc.factor = 1.25,
    max.auc.threshold = 1.5,
    frac.must.express = 0.1,
    frac.min.diff = 0.1,
    genes.use = NULL,
    root = NULL,
    segs.to.skip = NULL,
    only.return.global = F,
    must.beat.sibs = 0.6,
    report.debug = F
  )
  gene.markers[[tip]] <- markers
}
saveRDS(gene.markers, file="gene.markers.rds")
gene.markers <- readRDS("gene.markers.rds")

###################################################################
## Calculate the markers of each other population.
tips.to.run <- c("1", "2")
gene.markers <- list()
for (tipn in 1:length(tips.to.run)) {
  tip <- tips.to.run[tipn]
  print(paste0(Sys.time(), ": ", tip))
  markers <- aucprTestAlongTree(
    axial,
    pseudotime="pseudotime",
    tips=tip,
    log.effect.size = 0.25,
    auc.factor = 1, ## 1.25
    max.auc.threshold = 1, ## 1.5
    frac.must.express = 0.1,
    frac.min.diff = 0.1,
    genes.use = NULL,
    root = NULL,
    segs.to.skip = NULL,
    only.return.global = F,
    must.beat.sibs = 0.6,
    report.debug = F
  )
  gene.markers[[tip]] <- markers
}
saveRDS(gene.markers, file="gene.markers.tipSCPN.rds")
gene.markers <- readRDS("gene.markers.tipSCPN.rds")
SCPNmarkers <- unique(c(rownames(gene.markers[["2"]]),rownames(gene.markers[["1"]])))
saveRDS(SCPNmarkers, file="SCPNmarkers.rds")

###################################################################
## Generate impulse fits
gene.cascades_result <- lapply(tips.to.run, function(tip) {
  print(paste0(Sys.time(), ": Impulse Fit ", tip))
  seg.cells <- cellsAlongLineage(axial, tip, remove.root=F)
  casc <- geneCascadeProcess(axial, pseudotime='pseudotime', cells = seg.cells, genes=rownames(gene.markers[[tip]]), moving.window=5, cells.per.window=25, verbose = F)
  tip.file.name <- tip
  saveRDS(casc, file=paste0("casc_", tip.file.name, ".rds"))
  return(casc)
})
names(gene.cascades_result) <- tips.to.run

###################################################################
## Make a heatmap of every cascade in a single PDF.
gene.cascades_result <- readRDS("casc_N.rds")
pdf(file="cascades.SCPN_CSMN.pdf", width=7.5, height=10)
geneCascadeHeatmap(cascade=gene.cascades_result, color.scale = RColorBrewer::brewer.pal(9, "YlOrRd"), add.time = NULL, times.annotate = seq(0, 1, 0.1), title = "SCPN_CSMN", annotation.list = NULL, row.font.size = 1, max.font.size = 0.9)
dev.off()

## step by step
cascade=gene.cascades_result
color.scale = RColorBrewer::brewer.pal(9, "YlOrRd")
add.time = NULL
times.annotate = seq(0, 1, 0.1)
annotation.list = NULL
row.font.size = 1
max.font.size = 0.9

timing <- cascade$timing
timing[intersect(which(is.na(timing$time.on)), which(is.infinite(timing$time.off))),
       "time.on"] <- Inf
gene.order <- order(timing$time.on, timing$time.off, na.last = F)
cols <- (scales::gradient_n_pal(color.scale))(seq(0, 1, length.out = 50))
time <- as.numeric(names(cascade$pt.windows))
time.lab <- rep("", length(time))
for (annotate in times.annotate) {
  gt <- which(time >= annotate)
  if (length(gt) > 0)
    time.lab[min(gt)] <- as.character(annotate)
}

pdf(file="cascades.SCPN_CSMN.pdf", width=7.5, height=10)
gplots::heatmap.2(as.matrix(cascade$scaled.expression[gene.order,]), Rowv = F, Colv = F, dendrogram = "none", col = cols,
                  trace = "none", density.info = "none", key = F, labCol = NULL, labRow=NULL,
                  cexCol = 0.8, cexRow = 0.8, margins = c(5,8), lwid = c(0.3, 4), lhei = c(0.4, 4))
dev.off()

1"SCPN_CSMN"
2"SCPN_CTPN"
3"CThPN"
4"superficialCPN"
5"Layer4CPN"
6"deepCPN"

##################################################################
gene.markers <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/4_developmental_trajectories/2_gene-expressionCascades/gene.markers.tip1_2.rds")

df <- data.frame()
for (tip in names(gene.markers)) {
    tmp <- gene.markers[[tip]]
    tmp$geneName <- rownames(tmp)
    tmp$tip <- tip
    df <- rbind(df,tmp)
    rownames(df) <- 1:nrow(df)
}
geneBranch <- df[,c("tip","segment.maxBranch","segment.minBranch","geneName")]
geneBranch <- geneBranch[!duplicated(geneBranch),]
table(geneBranch$tip)
##################################################################
