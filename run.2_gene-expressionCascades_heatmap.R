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

axial <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/4_developmental_trajectories.newVersion.cleanedcell/1_URDTree/axial_buildTree.rds")
gene.markers <- readRDS("gene.markers.tipSCPN.rds")

###################################################################
## SCPN
gene.cascades_result <- readRDS("casc_SCPN.rds")

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

SCPNgene <- intersect(rownames(cascade$scaled.expression[gene.order,]),intersect(rownames(gene.markers[["2"]]),rownames(gene.markers[["1"]])))
SCPNcascade <- cascade

###################################################################
## CSMN
gene.cascades_result <- readRDS("casc_CSMN.rds")

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

CSMNgene <- intersect(rownames(cascade$scaled.expression[gene.order,]),setdiff(rownames(gene.markers[["1"]]),rownames(gene.markers[["2"]])))
CSMNcascade <- cascade

###################################################################
## CTPN
gene.cascades_result <- readRDS("casc_CTPN.rds")

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

CTPNgene <- intersect(rownames(cascade$scaled.expression[gene.order,]),setdiff(rownames(gene.markers[["2"]]),rownames(gene.markers[["1"]])))
CTPNcascade <- cascade

###################################################################
heatmapGene <- c(SCPNgene,CSMNgene,CTPNgene)

pdf(file="cascades.Figure2B_SCPN.pdf", width=7.5, height=10)
gplots::heatmap.2(as.matrix(SCPNcascade$scaled.expression[heatmapGene,]), Rowv = F, Colv = F, dendrogram = "none", col = cols,
                  trace = "none", density.info = "none", key = F, labCol = NULL, labRow=NULL,
                  cexCol = 0.8, cexRow = 0.8, margins = c(5,8), lwid = c(0.3, 4), lhei = c(0.4, 4))
dev.off()

pdf(file="cascades.Figure2B_CSMN.pdf", width=7.5, height=10)
gplots::heatmap.2(as.matrix(CSMNcascade$scaled.expression[heatmapGene,]), Rowv = F, Colv = F, dendrogram = "none", col = cols,
                  trace = "none", density.info = "none", key = F, labCol = NULL, labRow=NULL,
                  cexCol = 0.8, cexRow = 0.8, margins = c(5,8), lwid = c(0.3, 4), lhei = c(0.4, 4))
dev.off()

pdf(file="cascades.Figure2B_CTPN.pdf", width=7.5, height=10)
gplots::heatmap.2(as.matrix(CTPNcascade$scaled.expression[heatmapGene,]), Rowv = F, Colv = F, dendrogram = "none", col = cols,
                  trace = "none", density.info = "none", key = F, labCol = NULL, labRow=NULL,
                  cexCol = 0.8, cexRow = 0.8, margins = c(5,8), lwid = c(0.3, 4), lhei = c(0.4, 4))
dev.off()
