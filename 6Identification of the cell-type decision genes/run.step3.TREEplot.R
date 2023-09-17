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
# cat 1.txt | awk 'BEGIN{OFS="";ORS=",";}{print "\"",$1,"\"";}'

## node1L
for (geneName in c("SOX11","SOX4","ST8SIA2","HELLS","KCNQ3","ELAVL2","ELAVL4","SLCO5A1")) {
  pdf(file=paste0("node1L_",geneName,"_Tree.pdf"), width=10, height=10)
  print(plotTree(axial, geneName, title=geneName) + background)
  dev.off()
}
## node1R
for (geneName in c("CCDC175","LINC01965","PTPRM","FHOD3","SHROOM3","PRDM16","NRG1","MAST4","SLC1A3","PRKG1")) {
  pdf(file=paste0("node1R_",geneName,"_Tree.pdf"), width=10, height=10)
  print(plotTree(axial, geneName, title=geneName) + background)
  dev.off()
}
## node2L
for (geneName in c("RBFOX1","GRIK3","CDH13","SORCS1","SYNE2","MEGF11","DACT1","LMO7","SEPTIN11","MIR100HG","NRXN3")) {
  pdf(file=paste0("node2L_",geneName,"_Tree.pdf"), width=10, height=10)
  print(plotTree(axial, geneName, title=geneName) + background)
  dev.off()
}
## node2R
for (geneName in c("XKR4","GRID2","ADAMTS3","DPP10","ADAMTSL1","CARMIL1","CCBE1","CADM2","MYRIP","LINGO2","GALNT13")) {
  pdf(file=paste0("node2R_",geneName,"_Tree.pdf"), width=10, height=10)
  print(plotTree(axial, geneName, title=geneName) + background)
  dev.off()
}
## node3L
for (geneName in c("MEF2C","DYNC1I1","AC245123.1","OSBPL10","CDH8","SCD5","KIAA0319","SULT4A1","DLGAP2","SLC2A13","DSCAM")) {
  pdf(file=paste0("node3L_",geneName,"_Tree.pdf"), width=10, height=10)
  print(plotTree(axial, geneName, title=geneName) + background)
  dev.off()
}
## node3R
for (geneName in c("CASC15","RBFOX3","PTPRS","HECW1","TCF4","SLCO5A1","CASK","AFDN","CPE","PBX1","FRMD4A")) {
  pdf(file=paste0("node3R_",geneName,"_Tree.pdf"), width=10, height=10)
  print(plotTree(axial, geneName, title=geneName) + background)
  dev.off()
}
## node4L
for (geneName in c("TRPM3","PDE1A","RIMS2","KCNIP1","STARD13","NEGR1","LINC01435","SORCS1","AGAP1","PCNX1","SOX5")) {
  pdf(file=paste0("node4L_",geneName,"_Tree.pdf"), width=10, height=10)
  print(plotTree(axial, geneName, title=geneName) + background)
  dev.off()
}
## node4I
for (geneName in c("LRP1B","SLC26A7","TMEM108","ADAMTSL3","SOX2-OT","ZEB2","AC011474.1","LRRC4C","SASH1","ACTN2","SAMD5")) {
  pdf(file=paste0("node4I_",geneName,"_Tree.pdf"), width=10, height=10)
  print(plotTree(axial, geneName, title=geneName) + background)
  dev.off()
}
## node4R
for (geneName in c("MEF2C","TAFA1","DLG2","DAB1","FAT3","CELF2","RFX3","PTPRD","DOK5","NKAIN2","SYBU","ABLIM1","SATB2")) {
  pdf(file=paste0("node4R_",geneName,"_Tree.pdf"), width=10, height=10)
  print(plotTree(axial, geneName, title=geneName) + background)
  dev.off()
}



