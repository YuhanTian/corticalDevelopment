#########################################################################################################################################################################################################
node_branch_fun <- function(node_branch,node,branch){
  node_branch$gene <- rownames(node_branch)
  node_branch$importance <- sqrt(node_branch$importance)
  node_branch$node <- node
  node_branch$branch <- branch
  return(node_branch)
}
node1_branchL <- read.table("/data3/yuhan/Project/Neuron/scMultiome/4_developmental_trajectories.newVersion.cleanedcell/3_branch-point-associatedGenes/node1_branchL_featuresImportance.txt")
node1_branchR <- read.table("/data3/yuhan/Project/Neuron/scMultiome/4_developmental_trajectories.newVersion.cleanedcell/3_branch-point-associatedGenes/node1_branchR_featuresImportance.txt")
node2_branchL <- read.table("/data3/yuhan/Project/Neuron/scMultiome/4_developmental_trajectories.newVersion.cleanedcell/3_branch-point-associatedGenes/node2_branchL_featuresImportance.txt")
node2_branchR <- read.table("/data3/yuhan/Project/Neuron/scMultiome/4_developmental_trajectories.newVersion.cleanedcell/3_branch-point-associatedGenes/node2_branchR_featuresImportance.txt")
node3_branchL <- read.table("/data3/yuhan/Project/Neuron/scMultiome/4_developmental_trajectories.newVersion.cleanedcell/3_branch-point-associatedGenes/node3_branchL_featuresImportance.txt")
node3_branchR <- read.table("/data3/yuhan/Project/Neuron/scMultiome/4_developmental_trajectories.newVersion.cleanedcell/3_branch-point-associatedGenes/node3_branchR_featuresImportance.txt")
node4_branchL <- read.table("/data3/yuhan/Project/Neuron/scMultiome/4_developmental_trajectories.newVersion.cleanedcell/3_branch-point-associatedGenes/node4_branchL_featuresImportance.txt")
node4_branchI <- read.table("/data3/yuhan/Project/Neuron/scMultiome/4_developmental_trajectories.newVersion.cleanedcell/3_branch-point-associatedGenes/node4_branchI_featuresImportance.txt")
node4_branchR <- read.table("/data3/yuhan/Project/Neuron/scMultiome/4_developmental_trajectories.newVersion.cleanedcell/3_branch-point-associatedGenes/node4_branchR_featuresImportance.txt")
node1_branchL <- node_branch_fun(node1_branchL,"node1","branchL")
node1_branchR <- node_branch_fun(node1_branchR,"node1","branchR")
node2_branchL <- node_branch_fun(node2_branchL,"node2","branchL")
node2_branchR <- node_branch_fun(node2_branchR,"node2","branchR")
node3_branchL <- node_branch_fun(node3_branchL,"node3","branchL")
node3_branchR <- node_branch_fun(node3_branchR,"node3","branchR")
node4_branchL <- node_branch_fun(node4_branchL,"node4","branchL")
node4_branchI <- node_branch_fun(node4_branchI,"node4","branchI")
node4_branchR <- node_branch_fun(node4_branchR,"node4","branchR")

node_branch_df <- rbind(node1_branchL,node1_branchR,node2_branchL,node2_branchR,node3_branchL,node3_branchR,node4_branchL,node4_branchI,node4_branchR)
node_branch_df$type <- paste0(node_branch_df$node,"_",node_branch_df$branch,"_",node_branch_df$gene)
rownames(node_branch_df) <- node_branch_df$type

#########################################################################################################################################################################################################
/home/yuchen/miniconda3/envs/R4.0/bin/R
library(Seurat)
library(Signac)
library(RColorBrewer)
library(cowplot)
library(ggplot2)
library(dplyr)
library(reshape2)
background<-(theme_bw()
             +theme(axis.title.y = element_text(size = 20,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.text.y =  element_text(size = 8,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.title.x = element_text(size = 20,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.text.x = element_text(size = 10,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(strip.text  = element_text(size = 20,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(panel.grid =element_blank(),panel.background=element_rect(fill='transparent', color="#000000"))
             +theme(axis.line = element_line(size = 0.5),legend.background = element_blank())
             +theme(strip.background = element_rect(fill=NA,colour = NA))
)
percentage <- function(x){n=length(x);m=length(x[x>2]);res=m/n*100;return(res)}

cellID_segEarlyLater <- readRDS("cellID_segEarlyLater.rds")
URDobj <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/4_developmental_trajectories.newVersion.cleanedcell/1_URDTree/URDobj.rds")


df <- as.data.frame(table(cellID_segEarlyLater$cellID))
cellID_segEarlyLater <- cellID_segEarlyLater[which(!(cellID_segEarlyLater$cellID %in% df[which(df$Freq!=1),"Var1"])),]
cellID_segEarlyLater$segEarlyLater <- paste0(cellID_segEarlyLater$seg,cellID_segEarlyLater$type)
URDobj@meta.data$segEarlyLater <- cellID_segEarlyLater[rownames(URDobj@meta.data),]$segEarlyLater

# segEarlyLater_GeneExpression.dotplot
obj <- subset(URDobj, subset=segEarlyLater %in% c("7early","10early","1early","2early","9early","4early","3early","6early","5early"))
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
segEarlyLater_GeneExpression_per_df[which(segEarlyLater_GeneExpression_per_df$segEarlyLater=="7early"),]$type <- "node1_branchL"
segEarlyLater_GeneExpression_per_df[which(segEarlyLater_GeneExpression_per_df$segEarlyLater=="10early"),]$type <- "node1_branchR"
segEarlyLater_GeneExpression_per_df[which(segEarlyLater_GeneExpression_per_df$segEarlyLater=="1early"),]$type <- "node2_branchL"
segEarlyLater_GeneExpression_per_df[which(segEarlyLater_GeneExpression_per_df$segEarlyLater=="2early"),]$type <- "node2_branchR"
segEarlyLater_GeneExpression_per_df[which(segEarlyLater_GeneExpression_per_df$segEarlyLater=="9early"),]$type <- "node3_branchL"
segEarlyLater_GeneExpression_per_df[which(segEarlyLater_GeneExpression_per_df$segEarlyLater=="4early"),]$type <- "node3_branchR"
segEarlyLater_GeneExpression_per_df[which(segEarlyLater_GeneExpression_per_df$segEarlyLater=="3early"),]$type <- "node4_branchL"
segEarlyLater_GeneExpression_per_df[which(segEarlyLater_GeneExpression_per_df$segEarlyLater=="6early"),]$type <- "node4_branchI"
segEarlyLater_GeneExpression_per_df[which(segEarlyLater_GeneExpression_per_df$segEarlyLater=="5early"),]$type <- "node4_branchR"
segEarlyLater_GeneExpression_per_df$type <- paste0(segEarlyLater_GeneExpression_per_df$type,"_",segEarlyLater_GeneExpression_per_df$Gene)
rownames(segEarlyLater_GeneExpression_per_df) <- segEarlyLater_GeneExpression_per_df$type

node_branch_df$GEX <- segEarlyLater_GeneExpression_per_df[rownames(node_branch_df),]$Percent
# cat 1.txt | awk 'BEGIN{OFS="";ORS=",";}{print "\"",$1,"\"";}'
node_branch_df <- node_branch_df[which(node_branch_df$gene %in% c("CACNA2D1","SRRM4","ZFPM2","EOMES","CHD7","HES6","SYNE2","MAML3","DDX5","DISP3","ADGRB3","CSMD2","SHROOM3","PTN","MOXD1","GRIK3","SYN3","NRXN3","KLHL1","CTNND2","CDH6","PLXDC2","ELMOD1","MEGF11","AFDN","LMO7","CCBE1","ADAMTS3","GALNT13","ADAMTSL1","GRIK2","MIR137HG","RBPJ","RIPOR2","CSGALNACT1","NLGN1","NAV3","GABRB1","ARPP21","AL136456.1","TENM2","DYNC1I1","UNC5C","AC008591.1","SLC26A4","CASC15","ROBO2","CASK","ATP1B3","TRIM2","LYPD6","SEZ6L","NKAIN1","PLPPR4","PPFIA2","SORCS1","CNTNAP4","AGAP1","BCL11B","ZFHX3","EPB41L4A","PARD3","SEMA3E","MEG3","LRP1B","SLC26A7","TMEM108","ZEB2","LRRC4C","CCDC85A","AL033504.1","SGCZ","NTNG2","SEMA6D","TAFA1","DAB1","FAT3","DLG2","MEF2C","CELF2","BMPR1B","LIMCH1","PCDH11X")),]

node_branch_df$node <- factor(node_branch_df$node,levels=c("node1","node2","node3","node4"))
node_branch_df$branch <- factor(node_branch_df$branch,levels=c("branchL","branchI","branchR"))
node_branch_df$gene <- factor(node_branch_df$gene,levels=rev(c("CACNA2D1","SRRM4","ZFPM2","EOMES","CHD7","HES6","SYNE2","MAML3","DDX5","DISP3","ADGRB3","CSMD2","SHROOM3","PTN","MOXD1","GRIK3","SYN3","NRXN3","KLHL1","CTNND2","CDH6","PLXDC2","ELMOD1","MEGF11","AFDN","LMO7","CCBE1","ADAMTS3","GALNT13","ADAMTSL1","GRIK2","MIR137HG","RBPJ","RIPOR2","CSGALNACT1","NLGN1","NAV3","GABRB1","ARPP21","AL136456.1","TENM2","DYNC1I1","UNC5C","AC008591.1","SLC26A4","CASC15","ROBO2","CASK","ATP1B3","TRIM2","LYPD6","SEZ6L","NKAIN1","PLPPR4","PPFIA2","SORCS1","CNTNAP4","AGAP1","BCL11B","ZFHX3","EPB41L4A","PARD3","SEMA3E","MEG3","LRP1B","SLC26A7","TMEM108","ZEB2","LRRC4C","CCDC85A","AL033504.1","SGCZ","NTNG2","SEMA6D","TAFA1","DAB1","FAT3","DLG2","MEF2C","CELF2","BMPR1B","LIMCH1","PCDH11X")))

pdf("segEarlyLater_GeneExpression.dotplot.pdf",height=10,width=10)
ggplot(node_branch_df,aes(branch,gene)) + 
  geom_point(aes(color=GEX,size=importance)) + background + scale_size_area(max_size=7,breaks = seq(0, 1, 0.25)) +
  scale_colour_viridis_c() + facet_wrap(.~node,nrow=1,scales="free_x")
dev.off()

