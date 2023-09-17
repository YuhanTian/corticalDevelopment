#### branch-point-associated genes ####


#################################################################################################################################################
#### differentially expressed genes between parent and sibling branches
###################################################################
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

cell_infor <- axial@group.ids
cell_infor$pseudotime <- axial@tree$pseudotime[rownames(cell_infor)]
cell_infor <- cell_infor[order(cell_infor$segment,cell_infor$pseudotime),]

cellID_segEarlyLater <- data.frame()
for (ID in unique(cell_infor$segment)) {
  print(paste0("process seg",ID))
  df <- cell_infor[which(cell_infor$segment==ID),]
  tmp_min <- min(df$pseudotime)
  tmp_max <- max(df$pseudotime)
  print(tmp_min)
  print(tmp_max)
  df_tmp_early <- data.frame(cellID=rownames(df[which(df$pseudotime-tmp_min<0.04),]),type="early",seg=ID)
  df_tmp_later <- data.frame(cellID=rownames(df[which(tmp_max-df$pseudotime<0.04),]),type="later",seg=ID)
  cellID_segEarlyLater <- rbind(cellID_segEarlyLater,df_tmp_early,df_tmp_later)
}

background<-(theme_bw()
             +theme(plot.title = element_text(size = 24,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.title.y = element_text(size = 16,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.text.y =  element_text(size = 12,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.title.x = element_text(size = 16,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.text.x = element_text(size = 12,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.line = element_line(size = 0.5),legend.background = element_blank())
             +theme(panel.grid =element_blank(),panel.background=element_rect(fill='transparent', color="#000000"))
)


df <- as.data.frame(table(cellID_segEarlyLater$cellID))
cellID_segEarlyLater <- cellID_segEarlyLater[which(!(cellID_segEarlyLater$cellID %in% df[which(df$Freq!=1),"Var1"])),]


cellID_segEarlyLater$segEarlyLater <- paste0(cellID_segEarlyLater$seg,"_",cellID_segEarlyLater$type)
rownames(cellID_segEarlyLater) <- cellID_segEarlyLater$cellID
axial@group.ids$segEarlyLater <- cellID_segEarlyLater[rownames(axial@group.ids),"segEarlyLater"]
pdf(file="Tree_segEarlyLater.pdf", width=10, height=10)
plotTree(axial, "segEarlyLater", title="URD tree segment",label.segments=T,legend=F) + background
dev.off()
# saveRDS(cellID_segEarlyLater, file="cellID_segEarlyLater.rds")


###################################################################
/home/yuchen/miniconda3/envs/R4.0/bin/R
library(Seurat)
library(Signac)
library(RColorBrewer)
library(cowplot)
library(ggplot2)
library(dplyr)
library(reshape2)

cellID_segEarlyLater <- readRDS("cellID_segEarlyLater.rds")
URDobj <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/4_developmental_trajectories.newVersion.cleanedcell/1_URDTree/URDobj.rds")
cell_colors <- c("RG1"="#FDAE61","RG2"="#F46D43","RG3"="#D73027","RG4"="#A50026","RG_Cyc"="#FEE090","nIPC1"="#FFFFBF","nIPC2"="#FFFFBF","nIPC3"="#FFFFBF","ImmatureN1"="#C7E9C0","ImmatureN2"="#A1D99B","SCPN"="#41AB5D","SCPN_CSMN"="#238B45","SCPN_CTPN"="#006D2C","CThPN"="#00441B","MigN"="#ABD9E9","superficialCPN"="#74ADD1","Layer4CPN"="#4575B4","deepCPN"="#313695","MGE_IN"="#54278F","CGE_IN"="#6A51A3","CR"="#225EA8","Oligo"="#253494","MG/Peric"="#081D58")

pdf(file="DimPlot.URDobj.combine_by_Combine_CellType.pdf", height=10,width=10)
DimPlot(URDobj, reduction = "umap_combine", pt.size=1.2, label=TRUE, label.size = 6, group.by='Combine_CellType') + NoLegend() + scale_color_manual(values=c(cell_colors))
dev.off()
pdf(file="DimPlot.URDobj.combine_by_sample.pdf", height=10,width=10)
DimPlot(URDobj, reduction = "umap_combine", pt.size=1.2, label=F, label.size = 6, group.by='sample') + scale_color_manual(values=c("#D52126", "#88CCEE", "#FEE52C"))
dev.off()

df <- as.data.frame(table(cellID_segEarlyLater$cellID))
cellID_segEarlyLater <- cellID_segEarlyLater[which(!(cellID_segEarlyLater$cellID %in% df[which(df$Freq!=1),"Var1"])),]
cellID_segEarlyLater$compareType <- paste0(cellID_segEarlyLater$seg,cellID_segEarlyLater$type)
URDobj@meta.data$compareType <- cellID_segEarlyLater[rownames(URDobj@meta.data),]$compareType

## DEgenes
node1DEPC1 <- FindMarkers(URDobj, ident.1="7early", ident.2=c("11later"), slot='data', group.by="compareType", logfc.threshold=0.2, min.pct=0.1)
node1DEPC2 <- FindMarkers(URDobj, ident.1="10early", ident.2=c("11later"), slot='data', group.by="compareType", logfc.threshold=0.2, min.pct=0.1)
node1DEPC1gene <- rownames(node1DEPC1[which(node1DEPC1$p_val_adj<0.05 & node1DEPC1$avg_log2FC>0),])
node1DEPC2gene <- rownames(node1DEPC2[which(node1DEPC2$p_val_adj<0.05 & node1DEPC2$avg_log2FC>0),])

node2DEPC1 <- FindMarkers(URDobj, ident.1="1early", ident.2=c("7later"), slot='data', group.by="compareType", logfc.threshold=0.2, min.pct=0.1)
node2DEPC2 <- FindMarkers(URDobj, ident.1="2early", ident.2=c("7later"), slot='data', group.by="compareType", logfc.threshold=0.2, min.pct=0.1)
node2DEPC1gene <- rownames(node2DEPC1[which(node2DEPC1$p_val_adj<0.05 & node2DEPC1$avg_log2FC>0),])
node2DEPC2gene <- rownames(node2DEPC2[which(node2DEPC2$p_val_adj<0.05 & node2DEPC2$avg_log2FC>0),])

node3DEPC1 <- FindMarkers(URDobj, ident.1="9early", ident.2=c("10later"), slot='data', group.by="compareType", logfc.threshold=0.2, min.pct=0.1)
node3DEPC2 <- FindMarkers(URDobj, ident.1="4early", ident.2=c("10later"), slot='data', group.by="compareType", logfc.threshold=0.2, min.pct=0.1)
node3DEPC1gene <- rownames(node3DEPC1[which(node3DEPC1$p_val_adj<0.05 & node3DEPC1$avg_log2FC>0),])
node3DEPC2gene <- rownames(node3DEPC2[which(node3DEPC2$p_val_adj<0.05 & node3DEPC2$avg_log2FC>0),])

node4DEPC1 <- FindMarkers(URDobj, ident.1="3early", ident.2=c("9later"), slot='data', group.by="compareType", logfc.threshold=0.2, min.pct=0.1)
node4DEPC2 <- FindMarkers(URDobj, ident.1="6early", ident.2=c("9later"), slot='data', group.by="compareType", logfc.threshold=0.2, min.pct=0.1)
node4DEPC3 <- FindMarkers(URDobj, ident.1="5early", ident.2=c("9later"), slot='data', group.by="compareType", logfc.threshold=0.2, min.pct=0.1)
node4DEPC1gene <- rownames(node4DEPC1[which(node4DEPC1$p_val_adj<0.05 & node4DEPC1$avg_log2FC>0),])
node4DEPC2gene <- rownames(node4DEPC2[which(node4DEPC2$p_val_adj<0.05 & node4DEPC2$avg_log2FC>0),])
node4DEPC3gene <- rownames(node4DEPC3[which(node4DEPC3$p_val_adj<0.05 & node4DEPC3$avg_log2FC>0),])

nodeDEgene <- rbind(
  data.frame(gene=node1DEPC1gene,type="node1branchL"),
  data.frame(gene=node1DEPC2gene,type="node1branchR"),
  data.frame(gene=node2DEPC1gene,type="node2branchL"),
  data.frame(gene=node2DEPC2gene,type="node2branchR"),
  data.frame(gene=node3DEPC1gene,type="node3branchL"),
  data.frame(gene=node3DEPC2gene,type="node3branchR"),
  data.frame(gene=node4DEPC1gene,type="node4branchL"),
  data.frame(gene=node4DEPC2gene,type="node4branchI"),
  data.frame(gene=node4DEPC3gene,type="node4branchR")
)
write.table(nodeDEgene,"nodeDEgene.txt", col.names=T, row.names=F, sep="\t", quote=F)
write.csv(nodeDEgene,"nodeDEgene.csv", row.names=F)

## produce File
node1obj <- subset(URDobj, subset=compareType %in% c("7early","10early","11later"))
node1_matrx <- as.data.frame(t(as.matrix(node1obj@assays$RNA@data)))
node1_matrx <- node1_matrx[,c(node1DEPC1gene,node1DEPC2gene)]
tmp <- data.frame(class=node1obj@meta.data$compareType)
rownames(tmp) <- rownames(node1obj@meta.data)
node1Matrx_cell_gene_class <- cbind(tmp,node1_matrx)
write.csv(node1Matrx_cell_gene_class,"node1Matrx_cell_gene_class.csv", row.names=T)
write.csv(node1Matrx_cell_gene_class[,which(colnames(node1Matrx_cell_gene_class) %in% c("class",node1DEPC1gene))],"node1branchLMatrx_cell_gene_class.csv", row.names=T)
write.csv(node1Matrx_cell_gene_class[,which(colnames(node1Matrx_cell_gene_class) %in% c("class",node1DEPC2gene))],"node1branchRMatrx_cell_gene_class.csv", row.names=T)

node2obj <- subset(URDobj, subset=compareType %in% c("1early","2early","7later"))
node2_matrx <- as.data.frame(t(as.matrix(node2obj@assays$RNA@data)))
node2_matrx <- node2_matrx[,c(node2DEPC1gene,node2DEPC2gene)]
tmp <- data.frame(class=node2obj@meta.data$compareType)
rownames(tmp) <- rownames(node2obj@meta.data)
node2Matrx_cell_gene_class <- cbind(tmp,node2_matrx)
write.csv(node2Matrx_cell_gene_class,"node2Matrx_cell_gene_class.csv", row.names=T)
write.csv(node2Matrx_cell_gene_class[,which(colnames(node2Matrx_cell_gene_class) %in% c("class",node2DEPC1gene))],"node2branchLMatrx_cell_gene_class.csv", row.names=T)
write.csv(node2Matrx_cell_gene_class[,which(colnames(node2Matrx_cell_gene_class) %in% c("class",node2DEPC2gene))],"node2branchRMatrx_cell_gene_class.csv", row.names=T)

node3obj <- subset(URDobj, subset=compareType %in% c("9early","4early","10later"))
node3_matrx <- as.data.frame(t(as.matrix(node3obj@assays$RNA@data)))
node3_matrx <- node3_matrx[,c(node3DEPC1gene,node3DEPC2gene)]
tmp <- data.frame(class=node3obj@meta.data$compareType)
rownames(tmp) <- rownames(node3obj@meta.data)
node3Matrx_cell_gene_class <- cbind(tmp,node3_matrx)
write.csv(node3Matrx_cell_gene_class,"node3Matrx_cell_gene_class.csv", row.names=T)
write.csv(node3Matrx_cell_gene_class[,which(colnames(node3Matrx_cell_gene_class) %in% c("class",node3DEPC1gene))],"node3branchLMatrx_cell_gene_class.csv", row.names=T)
write.csv(node3Matrx_cell_gene_class[,which(colnames(node3Matrx_cell_gene_class) %in% c("class",node3DEPC2gene))],"node3branchRMatrx_cell_gene_class.csv", row.names=T)

node4obj <- subset(URDobj, subset=compareType %in% c("3early","6early","5early","9later"))
node4_matrx <- as.data.frame(t(as.matrix(node4obj@assays$RNA@data)))
node4_matrx <- node4_matrx[,c(node4DEPC1gene,node4DEPC2gene,node4DEPC3gene)]
tmp <- data.frame(class=node4obj@meta.data$compareType)
rownames(tmp) <- rownames(node4obj@meta.data)
node4Matrx_cell_gene_class <- cbind(tmp,node4_matrx)
write.csv(node4Matrx_cell_gene_class,"node4Matrx_cell_gene_class.csv", row.names=T)
write.csv(node4Matrx_cell_gene_class[,which(colnames(node4Matrx_cell_gene_class) %in% c("class",node4DEPC1gene))],"node4branchLMatrx_cell_gene_class.csv", row.names=T)
write.csv(node4Matrx_cell_gene_class[,which(colnames(node4Matrx_cell_gene_class) %in% c("class",node4DEPC2gene))],"node4branchIMatrx_cell_gene_class.csv", row.names=T)
write.csv(node4Matrx_cell_gene_class[,which(colnames(node4Matrx_cell_gene_class) %in% c("class",node4DEPC3gene))],"node4branchRMatrx_cell_gene_class.csv", row.names=T)

#################################################################################################################################################
####  To find the top distinguishing features/genes between cells in sibling and parent branches at a given branch-point in the development trajectory
###################################################################
# GradientBoostingClassifier





