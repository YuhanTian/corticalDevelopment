#########################################################################################################################################################################
/home/yuchen/miniconda3/envs/R4.0/bin/R
library(Seurat)
library(Signac)
library(RColorBrewer)
library(cowplot)
library(ggplot2)
library(dplyr)
library(reshape2)

URDobj <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/4_developmental_trajectories.newVersion.cleanedcell/1_URDTree/URDobj.rds")

cell_colors <- c("RG1"="#FDAE61","RG2"="#F46D43","RG3"="#D73027","RG4"="#A50026","RG_Cyc"="#FEE090","nIPC1"="#FFFFBF","nIPC2"="#FFFFBF","nIPC3"="#FFFFBF","ImmatureN1"="#C7E9C0","ImmatureN2"="#A1D99B","SCPN"="#41AB5D","SCPN_CSMN"="#238B45","SCPN_CTPN"="#006D2C","CThPN"="#00441B","MigN"="#ABD9E9","superficialCPN"="#74ADD1","Layer4CPN"="#4575B4","deepCPN"="#313695","MGE_IN"="#54278F","CGE_IN"="#6A51A3","CR"="#225EA8","Oligo"="#253494","MG/Peric"="#081D58")

pdf(file="DimPlot.URDobj.combine_by_Combine_CellType.pdf", height=10,width=10)
DimPlot(URDobj, reduction = "umap_combine", pt.size=1.2, label=TRUE, label.size = 6, group.by='Combine_CellType') + NoLegend() + scale_color_manual(values=c(cell_colors))
dev.off()
pdf(file="DimPlot.URDobj.combine_by_sample.pdf", height=10,width=10)
DimPlot(URDobj, reduction = "umap_combine", pt.size=1.2, label=F, label.size = 6, group.by='sample') + scale_color_manual(values=c("#D52126", "#88CCEE", "#FEE52C"))
dev.off()

#########################################################################################################################################################################
## DEgenes
DECPN1 <- FindMarkers(URDobj, ident.1="superficialCPN", ident.2="Layer4CPN", slot='data', group.by="Combine_CellType", logfc.threshold=0.2, min.pct=0.1)
DECPN2 <- FindMarkers(URDobj, ident.1="superficialCPN", ident.2="deepCPN", slot='data', group.by="Combine_CellType", logfc.threshold=0.2, min.pct=0.1)
DECPN3 <- FindMarkers(URDobj, ident.1="Layer4CPN", ident.2="deepCPN", slot='data', group.by="Combine_CellType", logfc.threshold=0.2, min.pct=0.1)
# DECPN1[which(DECPN1$p_val_adj<0.05 & DECPN1$avg_log2FC>0.2),]
# DECPN2[which(DECPN2$p_val_adj<0.05 & DECPN2$avg_log2FC>0.2),]
# DECPN3[which(DECPN3$p_val_adj<0.05 & DECPN3$avg_log2FC>0.2),]
FCValue=1

# background<-(theme_bw()
#              +theme(axis.title.y = element_text(size = 16,vjust = 0.5, hjust = 0.5, colour = "black"))
#              +theme(axis.text.y =  element_text(size = 12,vjust = 0.5, hjust = 0.5, colour = "black"))
#              +theme(axis.title.x = element_text(size = 16,vjust = 0.5, hjust = 0.5, colour = "black"))
#              +theme(axis.text.x = element_text(size = 12,vjust = 0.5, hjust = 0.5, colour = "black"))
#              +theme(strip.text  = element_text(size = 10,vjust = 0.5, hjust = 0.5, colour = "black"))
#              +theme(axis.line = element_line(size = 0.5),legend.background = element_blank())
#              +theme(strip.background = element_rect(fill=NA,colour = NA))
# )

# pdf('DECPN1.pdf', width=10, height=10)
# (ggplot(DECPN1)
#  +geom_point(aes(avg_log2FC,-log10(p_val_adj),color=(p_val_adj<0.05 & abs(avg_log2FC)>FCValue)),alpha=0.5,size=0.5)
#  +geom_vline(xintercept = c(-1,1), lty=2, size=0.2)
#  +scale_color_manual(values = c("grey50","red"))
#  +xlim(-6,6)
#  +background
#  +theme(legend.position = "none")
# )
# dev.off()

upGene <- intersect(rownames(DECPN1[which(DECPN1$p_val_adj<0.05 & DECPN1$avg_log2FC>FCValue),]),rownames(DECPN2[which(DECPN2$p_val_adj<0.05 & DECPN2$avg_log2FC>FCValue),]))
downGene <- intersect(rownames(DECPN1[which(DECPN1$p_val_adj<0.05 & DECPN1$avg_log2FC<(-FCValue)),]),rownames(DECPN2[which(DECPN2$p_val_adj<0.05 & DECPN2$avg_log2FC<(-FCValue)),]))
markergeneDF <- rbind(data.frame(gene=upGene,type="up"),data.frame(gene=downGene,type="down"))
saveRDS(markergeneDF, file="markergeneDF.rds")
markergeneDF <- readRDS("markergeneDF.rds")

superficialCPNobj <- subset(URDobj, subset=Combine_CellType %in% c("superficialCPN"))
superficialCPN_matrx <- as.data.frame(t(as.matrix(superficialCPNobj@assays$RNA@data)))
superficialCPN_matrx <- superficialCPN_matrx[,markergeneDF$gene]
superficialCPNCounts <- apply(superficialCPN_matrx,2,mean)

Layer4CPNobj <- subset(URDobj, subset=Combine_CellType %in% c("Layer4CPN"))
Layer4CPN_matrx <- as.data.frame(t(as.matrix(Layer4CPNobj@assays$RNA@data)))
Layer4CPN_matrx <- Layer4CPN_matrx[,markergeneDF$gene]
Layer4CPNCounts <- apply(Layer4CPN_matrx,2,mean)

deepCPNobj <- subset(URDobj, subset=Combine_CellType %in% c("deepCPN"))
deepCPN_matrx <- as.data.frame(t(as.matrix(deepCPNobj@assays$RNA@data)))
deepCPN_matrx <- deepCPN_matrx[,markergeneDF$gene]
deepCPNCounts <- apply(deepCPN_matrx,2,mean)

df <- data.frame(superficialCPNCounts=superficialCPNCounts[markergeneDF$gene],deepCPNCounts=deepCPNCounts[markergeneDF$gene],Layer4CPNCounts=Layer4CPNCounts[markergeneDF$gene])

# library(pheatmap)
# avedata <- df
# avedata.sorted1 <- scale(avedata)
# avedata.sorted1 <- t(avedata.sorted1)
# pdf("heatmap.v1.pdf", w=20, h=6)
# p = pheatmap(avedata.sorted1, color=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100), cluster_rows=F, cluster_cols=T, show_colnames=F, show_rownames=T, scale='none',breaks=seq(-1,1,by=0.02))
# print(p)
# dev.off()

library(pheatmap)
anno_row<-data.frame(row.names = markergeneDF$gene, differential=markergeneDF$type)
pheatmapDF <- df
pdf("heatmap.v1.pdf", w=10, h=10)
pheatmap(log2(pheatmapDF+1),
         clustering_method = "ward.D",
         color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
         show_rownames = T, annotation_row = anno_row, cluster_rows = T, cluster_cols = F)
dev.off()

pdf(file="UMAP_cellType.counts.pdf", width=10, height=10)
print(FeaturePlot(URDobj, slot="counts",reduction = "umap_combine",features = c("FOXP2"), ncol=1, pt.size=1.2, order=TRUE, label=TRUE, min.cutoff = "q1", max.cutoff = "q99") + theme(legend.position = "right"))
dev.off()

ggplot(dfplot,aes(x=type ,y=log2(counts)),color=type)+
  geom_boxplot(aes(fill=factor(type)),outlier.shape = NA,notch = T) +
  ylab("Number of reads [log2]") + xlab("") +
  theme(panel.grid =element_blank(),
        panel.background=element_rect(fill='transparent', color="#000000"),
        legend.text=element_text(size=10),
        legend.position=c(0.03,0.99),
        legend.justification = c(-0.03,0.99),
        legend.key.size=unit(0.3,'cm')) +
  scale_fill_manual(values = c("#FDD692","#EC7357"))+
  scale_y_continuous(limits=c(2,8),breaks=seq(2,8,2))


