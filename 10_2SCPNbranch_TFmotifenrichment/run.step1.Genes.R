#########################################################################################################################################################################################################
/home/yuchen/miniconda3/envs/R4.0/bin/R
library(Seurat)
library(Signac)
library(RColorBrewer)
library(cowplot)
library(ggplot2)
library(dplyr)
library(reshape2)

cellID_segEarlyLater <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/4_developmental_trajectories.newVersion.cleanedcell/3_branch-point-associatedGenes/cellID_segEarlyLater.rds")
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
node2DEPC1 <- FindMarkers(URDobj, ident.1="1early", ident.2=c("7later"), slot='data', group.by="compareType", logfc.threshold=0.2, min.pct=0.1)
node2DEPC2 <- FindMarkers(URDobj, ident.1="2early", ident.2=c("7later"), slot='data', group.by="compareType", logfc.threshold=0.2, min.pct=0.1)
node2DEPC1gene <- rownames(node2DEPC1[which(node2DEPC1$p_val_adj<0.05 & node2DEPC1$avg_log2FC>0),])
node2DEPC2gene <- rownames(node2DEPC2[which(node2DEPC2$p_val_adj<0.05 & node2DEPC2$avg_log2FC>0),])
node2DEPC_all <- unique(c(node2DEPC1gene,node2DEPC2gene))

#########################################################################################################################################################################################################
node2DE_CSMN_CTPN <- FindMarkers(URDobj, ident.1="1early", ident.2=c("2early"), slot='data', group.by="compareType", logfc.threshold=0.2, min.pct=0.1)
CSMNmarker <- rownames(node2DE_CSMN_CTPN[which(node2DE_CSMN_CTPN$p_val_adj<0.05 & node2DE_CSMN_CTPN$avg_log2FC>0),])
CTPNmarker <- rownames(node2DE_CSMN_CTPN[which(node2DE_CSMN_CTPN$p_val_adj<0.05 & node2DE_CSMN_CTPN$avg_log2FC<0),])
node2DE_PC_CSMN_CTPN <- unique(c(intersect(CSMNmarker,node2DEPC_all),intersect(CTPNmarker,node2DEPC_all)))

#########################################################################################################################################################################################################
node2obj <- subset(URDobj, subset=compareType %in% c("1early","2early","7later"))

node2_matrx <- as.data.frame(t(as.matrix(node2obj@assays$RNA@data)))
node2_matrx <- node2_matrx[,intersect(CSMNmarker,node2DEPC_all)]
tmp <- data.frame(class=node2obj@meta.data$compareType)
rownames(tmp) <- rownames(node2obj@meta.data)
node2Matrx_cell_gene_class <- cbind(tmp,node2_matrx)
write.csv(node2Matrx_cell_gene_class,"node2branchLMatrx_cell_gene_class.csv", row.names=T)

node2_matrx <- as.data.frame(t(as.matrix(node2obj@assays$RNA@data)))
node2_matrx <- node2_matrx[,intersect(CTPNmarker,node2DEPC_all)]
tmp <- data.frame(class=node2obj@meta.data$compareType)
rownames(tmp) <- rownames(node2obj@meta.data)
node2Matrx_cell_gene_class <- cbind(tmp,node2_matrx)
write.csv(node2Matrx_cell_gene_class,"node2branchRMatrx_cell_gene_class.csv", row.names=T)

#######################################################
run.step2.scikit-learn.node2_branchL.py
run.step2.scikit-learn.node2_branchR.py

#######################################################
node_branch_fun <- function(node_branch,node,branch){
  node_branch$gene <- rownames(node_branch)
  node_branch$importance <- sqrt(node_branch$importance)
  node_branch$node <- node
  node_branch$branch <- branch
  return(node_branch)
}
node2_branchL <- read.table("node2_branchL_featuresImportance.txt")
node2_branchR <- read.table("node2_branchR_featuresImportance.txt")
node2_branchL <- node_branch_fun(node2_branchL,"node2","branchL")
node2_branchR <- node_branch_fun(node2_branchR,"node2","branchR")

node_branch_df <- rbind(node2_branchL,node2_branchR)
node_branch_df$type <- paste0(node_branch_df$node,"_",node_branch_df$branch,"_",node_branch_df$gene)
rownames(node_branch_df) <- node_branch_df$type

background<-(theme_bw()
             +theme(axis.title.y = element_text(size = 20,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.text.y =  element_text(size = 10,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.title.x = element_text(size = 20,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.text.x = element_text(size = 10,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(strip.text  = element_text(size = 20,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(panel.grid =element_blank(),panel.background=element_rect(fill='transparent', color="#000000"))
             +theme(axis.line = element_line(size = 0.5),legend.background = element_blank())
             +theme(strip.background = element_rect(fill=NA,colour = NA))
)
percentage <- function(x){n=length(x);m=length(x[x>2]);res=m/n*100;return(res)}

df <- as.data.frame(table(cellID_segEarlyLater$cellID))
cellID_segEarlyLater <- cellID_segEarlyLater[which(!(cellID_segEarlyLater$cellID %in% df[which(df$Freq!=1),"Var1"])),]
cellID_segEarlyLater$segEarlyLater <- paste0(cellID_segEarlyLater$seg,cellID_segEarlyLater$type)
URDobj@meta.data$segEarlyLater <- cellID_segEarlyLater[rownames(URDobj@meta.data),]$segEarlyLater

# segEarlyLater_GeneExpression.dotplot
obj <- subset(URDobj, subset=segEarlyLater %in% c("1early","2early"))
marker <- unique(node_branch_df$gene)
GeneExpressionData <- obj@assays$RNA@counts
segEarlyLater_GeneExpression <- data.frame(segEarlyLater=obj@meta.data$segEarlyLater)
for (markerID in marker) {
  tmp <- data.frame(markerID=GeneExpressionData[markerID,])
  colnames(tmp) <- markerID
  segEarlyLater_GeneExpression <- cbind(segEarlyLater_GeneExpression,tmp)
}

segEarlyLater_GeneExpression_per <- segEarlyLater_GeneExpression %>% group_by(segEarlyLater) %>% summarise_all(percentage) %>% as.data.frame()
segEarlyLater_GeneExpression_per_df <- melt(segEarlyLater_GeneExpression_per,value.name="Percent",variable.name="Gene")
segEarlyLater_GeneExpression_per_df$type <- "node_branch"

segEarlyLater_GeneExpression_per_df[which(segEarlyLater_GeneExpression_per_df$segEarlyLater=="1early"),]$type <- "node2_branchL"
segEarlyLater_GeneExpression_per_df[which(segEarlyLater_GeneExpression_per_df$segEarlyLater=="2early"),]$type <- "node2_branchR"

segEarlyLater_GeneExpression_per_df$type <- paste0(segEarlyLater_GeneExpression_per_df$type,"_",segEarlyLater_GeneExpression_per_df$Gene)
rownames(segEarlyLater_GeneExpression_per_df) <- segEarlyLater_GeneExpression_per_df$type

node_branch_df$GEX <- segEarlyLater_GeneExpression_per_df[rownames(node_branch_df),]$Percent
# cat 1.txt | awk 'BEGIN{OFS="";ORS=",";}{print "\"",$1,"\"";}'
node_branch_df <- node_branch_df[which(node_branch_df$gene %in% c("GRIK3","LSAMP","SYN3","SYNE2","MIR100HG","AFDN","KIRREL3","NRXN3","STXBP6","GSG1L","CCBE1","ADAMTS3","GALNT13","ADAMTSL1","RIPOR2","GRIK2","CDH18","LINC02082","DPP10","PCSK5")),]

node_branch_df$branch <- factor(node_branch_df$branch,levels=c("branchL","branchR"))
node_branch_df$gene <- factor(node_branch_df$gene,levels=rev(c("GRIK3","LSAMP","SYN3","SYNE2","MIR100HG","AFDN","KIRREL3","NRXN3","STXBP6","GSG1L","CCBE1","ADAMTS3","GALNT13","ADAMTSL1","RIPOR2","GRIK2","CDH18","LINC02082","DPP10","PCSK5")))

summary(node_branch_df$importance)
pdf("CSMN.CTPN.marker.segEarlyLater_GeneExpression.dotplot.pdf",height=10,width=10)
ggplot(node_branch_df,aes(branch,gene)) +
  geom_point(aes(color=GEX,size=importance)) + background + scale_size_area(max_size=7,breaks = seq(0, 0.5, 0.1)) +
  scale_colour_viridis_c() + scale_x_discrete(labels= c("CSMN","CTPN"))
dev.off()

# #########################################################################################################################################################################################################
# ## produce File
# node2DE_PC_CSMN_CTPN <- c("RBFOX1","GRIK3","SEPTIN11","NRXN3","SYNE2","KLHL1","DSCAM","AFDN","SSBP3","LSAMP","RGS6","MIR100HG","CTNND2","CDHR3","STXBP6","GRM3","FOXP2","XKR4","GRID2","ADAMTS3","DPP10","ADAMTSL1","CARMIL1","CCBE1","CADM2","LINGO2","GALNT13","MYRIP","LINC02082","GRID1","RIPOR2","NLGN1","TENM2","PLCL1","PLXDC2","CDH9","PRKG1")
# node2obj <- subset(URDobj, subset=compareType %in% c("1early","2early"))
# node2_matrx <- as.data.frame(t(as.matrix(node2obj@assays$RNA@data)))
# node2_matrx <- node2_matrx[,node2DE_PC_CSMN_CTPN]
# tmp <- data.frame(class=node2obj@meta.data$compareType)
# rownames(tmp) <- rownames(node2obj@meta.data)
# node2Matrx_cell_gene_class <- cbind(tmp,node2_matrx)
# node2Matrx_cell_gene_class[which(node2Matrx_cell_gene_class$class=="1early"),"class"] <- "CSMN"
# node2Matrx_cell_gene_class[which(node2Matrx_cell_gene_class$class=="2early"),"class"] <- "CTPN"

# # 使用H检验找出不同成熟程度差异的features
# featureList = node2DE_PC_CSMN_CTPN
# Pval = c()
# for(ii in featureList){
#     Pval = c(Pval, kruskal.test(get(ii) ~ class, node2Matrx_cell_gene_class)$p.value)
# }
# Pval_adj = p.adjust(Pval, method="BH")
# # 取出P-values小于0.01的特征
# Pans = data.frame(Pval_adj)
# rownames(Pans) = featureList
# selfeature = rownames(Pans)[Pans$Pval_adj <= 0.01]
# length(selfeature)

# # 然后对筛选过后的特征进行层次聚类
# seldata = node2Matrx_cell_gene_class[, c(selfeature,"class")]

# # 按照TLS成熟度将所有TLS的feature进行平均，同时对所求的的结果进行横向的标准化
# avedata = c()
# for(ii in unique(seldata$class)){
#     subda = seldata[seldata$class == ii,selfeature]
#     avedata = cbind(avedata, colMeans(subda))
# }
# colnames(avedata) = unique(seldata$class)
# avedata.sorted1 <- scale(avedata)
# avedata.sorted1 <- t(avedata.sorted1)
# pdf("heatmap.v1.pdf", w=20, h=6)
# p = pheatmap(avedata.sorted1, color=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100), cluster_rows=F, cluster_cols=T, show_colnames=F, show_rownames=T, scale='none',breaks=seq(-1,1,by=0.02))
# print(p)
# dev.off()
# #########################################################################################################################################################################################################







