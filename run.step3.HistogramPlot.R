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

df <- axial@group.ids
df <- df[which(df$segment %in% c("7","11","10","9")),]
df <- df[which(df$Combine_CellType %in% c("RG1","RG2","RG3","RG4")),]

#####################################################################
dfplot <- as.data.frame(table(df[,c("Combine_CellType","segment")])/991)
pdf(file="RGsegment.pdf", width=5, height=5)
ggplot(dfplot,aes(x=factor(segment,levels = c("11","10","9","7")),y=Freq,fill=Combine_CellType)) +
  geom_bar(position="stack",stat="identity",width=0.8) +
  labs(title ="")+xlab("segment")+ylab("proportion")+
  scale_y_continuous(limits=c(0,0.55),breaks=seq(0,0.55,0.1),expand = c(0, 0)) +
  theme(panel.grid =element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background=element_rect(fill='transparent', color="#000000"),
        legend.text=element_text(size=10),
        legend.background = element_blank(),
        legend.position=c(0.8,0.9),
        legend.key.size=unit(0.3,'cm'))+
  guides(fill=guide_legend(title=NULL)) +
  scale_fill_manual(values = c("#A9D179","#BCB8D3","#F5AE6B","#84CAC0"))
dev.off()

#####################################################################
dfplot <- as.data.frame(table(df[,c("Combine_CellType","segment")]))
dfplot$perC <- 1
dfplot[which(dfplot$Combine_CellType=="RG1"),]$perC <- dfplot[which(dfplot$Combine_CellType=="RG1"),]$Freq/sum(dfplot[which(dfplot$Combine_CellType=="RG1"),]$Freq)
dfplot[which(dfplot$Combine_CellType=="RG2"),]$perC <- dfplot[which(dfplot$Combine_CellType=="RG2"),]$Freq/sum(dfplot[which(dfplot$Combine_CellType=="RG2"),]$Freq)
dfplot[which(dfplot$Combine_CellType=="RG3"),]$perC <- dfplot[which(dfplot$Combine_CellType=="RG3"),]$Freq/sum(dfplot[which(dfplot$Combine_CellType=="RG3"),]$Freq)
dfplot[which(dfplot$Combine_CellType=="RG4"),]$perC <- dfplot[which(dfplot$Combine_CellType=="RG4"),]$Freq/sum(dfplot[which(dfplot$Combine_CellType=="RG4"),]$Freq)

pdf(file="RGsegment.pdf", width=5, height=5)
ggplot(dfplot,aes(x=Combine_CellType,y=perC,fill=factor(segment,levels = c("11","10","9","7")))) +
  geom_bar(position="stack",stat="identity",width=0.8) +
  labs(title ="")+xlab("")+ylab("proportion")+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1),expand = c(0, 0)) +
  theme(panel.grid =element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background=element_rect(fill='transparent', color="#000000"),
        legend.text=element_text(size=10),
        legend.background = element_blank(),
        legend.position=c(0.8,0.9),
        legend.key.size=unit(0.3,'cm'))+
  guides(fill=guide_legend(title=NULL)) +
  scale_fill_manual(values = c("#A9D179","#BCB8D3","#F5AE6B","#84CAC0"))
dev.off()