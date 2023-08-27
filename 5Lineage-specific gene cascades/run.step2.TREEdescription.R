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
background<-(theme_bw()
             +theme(plot.title = element_text(size = 24,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.title.y = element_text(size = 16,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.text.y =  element_text(size = 12,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.title.x = element_text(size = 16,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.text.x = element_text(size = 12,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.line = element_line(size = 0.5),legend.background = element_blank())
             +theme(panel.grid =element_blank(),panel.background=element_rect(fill='transparent', color="#000000"))
)
axial <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/4_developmental_trajectories.newVersion.cleanedcell/1_URDTree/axial_buildTree.rds")
cell_type <- c("RG_Cyc","RG1","RG2","RG3","RG4","Oligo","MG/Peric","nIPC1","nIPC2","nIPC3","ImmatureN1","ImmatureN2","SCPN","SCPN_CSMN","SCPN_CTPN","CThPN","MigN","superficialCPN","Layer4CPN","deepCPN","CR","MGE_IN","CGE_IN")
cell_colors <- c("RG1"="#FDAE61","RG2"="#F46D43","RG3"="#D73027","RG4"="#A50026","RG_Cyc"="#FEE090","nIPC1"="#FFFFBF","nIPC2"="#FFFFBF","nIPC3"="#FFFFBF","ImmatureN1"="#C7E9C0","ImmatureN2"="#A1D99B","SCPN"="#41AB5D","SCPN_CSMN"="#238B45","SCPN_CTPN"="#006D2C","CThPN"="#00441B","MigN"="#ABD9E9","superficialCPN"="#74ADD1","Layer4CPN"="#4575B4","deepCPN"="#313695","MGE_IN"="#54278F","CGE_IN"="#6A51A3","CR"="#225EA8","Oligo"="#253494","MG/Peric"="#081D58")
marker <- c("TOP2A","MKI67","SOX2","PAX6","HES1","HES5","OLIG1","PTPRC","PDGFRB","GFAP","EOMES","TUBB3","NEUROD2","NEUROD6","BCL11B","LDB2","FEZF2","DIAPH3","FOXP2","TLE4","CRYM","SATB2","NRP1","CUX2","RORB","LMO4","LPL","RELN","GAD1","DLX1","LHX6","NR2F2")

## plot segment, CellType and transitions
pdf(file="Tree_segment.pdf", width=10, height=10)
plotTree(axial, "segment", title="URD tree segment",label.segments=T,legend=F) + background
dev.off()

pdf(file="Tree_CellType.pdf", width=10, height=10)
plotTree(axial, "stage", title="Developmental trajectories",discrete.colors=cell_colors,cell.size = 1,legend=F) + background
dev.off()
for (cellType in unique(axial@group.ids$Combine_CellType)){
  axial@group.ids[which(axial@group.ids$Combine_CellType==cellType),cellType] <- cellType
  pdf(file=paste0("Tree_",cellType,".pdf"), width=10, height=10)
  print(plotTree(axial, cellType, title=cellType,discrete.colors=cell_colors,cell.size = 1,legend=F) + background)
  dev.off()
}

pdf(file="tSNE_transitions.pdf", width=10, height=10)
plotDim(axial, "Combine_CellType", transitions.plot = 10000, plot.title="tSNE_transitions",legend=F) + scale_color_manual(values=c(cell_colors)) + background
dev.off()

df <- data.frame()
for (cellType in unique(axial@group.ids$Combine_CellType)){
  tmp <- as.data.frame(table(axial@group.ids[which(axial@group.ids$Combine_CellType==cellType),"segment"]))
colnames(tmp) <- c("segment","number")
tmp$cellType <- cellType
df <- rbind(df,tmp)
}

## plot sample
pdf(file="sample.pdf", width=10, height=10)
print(plotTree(axial, "sample", title="sample") + background)
dev.off()

## plot gene expression
for (geneName in marker) {
  pdf(file=paste0("cellTypeMarker_",geneName,"_Tree.pdf"), width=10, height=10)
  print(plotTree(axial, geneName, title=geneName) + background)
  dev.off()
}

## Pseudotime
cellPseudotime <- axial@tree$pseudotime
a <- matrix(cellPseudotime[colnames(axial@count.data)],nrow=1)
rownames(a) <- "cellPseudotime"
axial@count.data <- rbind(a,axial@count.data)

a <- matrix(cellPseudotime[colnames(axial@logupx.data)],nrow=1)
rownames(a) <- "cellPseudotime"
axial@logupx.data <- rbind(a,axial@logupx.data)

pdf(file="cellPseudotime.pdf", width=10, height=10)
print(plotTree(axial, "cellPseudotime", title="cellPseudotime") + background)
dev.off()
axial@group.ids$cellPseudotime <- cellPseudotime[rownames(axial@group.ids)]

# ## Distribution of pseudotime
# df <- axial@group.ids
# df_seg9 <- df[which(df$segment==9),]
# df_seg8 <- df[which(df$segment==8),]

# pdf(file="nIPC_in_SCPN.pdf", width=11.5, height=10)
# ggplot(df_seg9[which(df_seg9$Combine_CellType %in% c("nIPC1","nIPC2","nIPC3")),], aes(x=cellPseudotime,colour=Combine_CellType)) +
#   geom_density() +
#   labs(title = "Distribution of pseudotime in SCPN") +
#   xlab("cellPseudotime") + ylab("Density") + 
#   scale_x_continuous(limits=c(0,0.8),breaks=seq(0,0.8,0.2))+
#   scale_color_manual(values=c("#D52126", "#88CCEE", "#FEE52C")) + background
# dev.off()

# pdf(file="nIPC_in_RG.pdf", width=11.5, height=10)
# ggplot(df_seg8[which(df_seg8$Combine_CellType %in% c("RG0","RG1","RG2","RG3","RG4","nIPC2")),], aes(x=cellPseudotime,colour=Combine_CellType)) +
#   geom_density() +
#   labs(title = "Distribution of pseudotime in RG") +
#   xlab("cellPseudotime") + ylab("Density") + 
#   scale_x_continuous(limits=c(0,0.4),breaks=seq(0,0.4,0.1)) +
#   scale_color_manual(values=c("#D73027", "#FDAE61",  "#FFFFBF",  "#A6D96A", "#66BD63", "#1A9850")) + 
#   background
# dev.off()

saveRDS(axial@group.ids, file="axial.group.ids.rds")


