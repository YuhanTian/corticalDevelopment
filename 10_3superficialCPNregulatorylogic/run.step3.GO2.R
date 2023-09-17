load(file="superficialCPN.pseudobulk.Peak_RNA.km20.heatmap.RData")

#########################################################################################################################################################
# /home/yuchen/miniconda3/envs/R4.0/bin/R
library(Seurat)
library(Signac)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#########################################################################################################################################################
cur_markers <- readRDS("superficialCPNmarkersGene_CellType.rds")
markersGene <- unique(cur_markers[which(cur_markers$p_val_adj<0.01 & abs(cur_markers$avg_log2FC)>0.25),"gene"])

# LINKobj <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/3_CREgene_pair/obj.CREgene_pair.rds")
genePeak_links <- Links(LINKobj[["peaks"]])
genePeak_links <- genePeak_links[genePeak_links$score>0]
genePeak_links <- genePeak_links[genePeak_links$gene %in% markersGene]

# ht_list <- readRDS("superficialCPN_ht_list.rds")
RR <- row_order(ht_list)

peakOrderDF <- data.frame()
for (id in as.character(c(19,15,14,16,13,17,11,10,12,20,4,3,6,7,9,8))) {
  peakOrderDF <- rbind(peakOrderDF,data.frame(KMid=id,linkNUM=RR[[id]]))
}
rownames(peakOrderDF) <- peakOrderDF$linkNUM #hanghao

linkINFO <- data.frame(linkID=paste0("link_",c(1:length(genePeak_links))),gene=genePeak_links$gene,peak=genePeak_links$peak)
linkINFO <- linkINFO[-deleteID,]
linkINFO <- linkINFO[peakOrderDF$linkNUM,]
linkINFO$KMid <- peakOrderDF$KMid
linkINFO$heatmapRow <- c(1:nrow(linkINFO))

linkINFO$KMnum <- "KM"
linkINFO[which(linkINFO$KMid %in% c(19,15,14,16,13,17,11,10,12,20)),]$KMnum <- "KM1"
linkINFO[which(linkINFO$KMid %in% c(4,3)),]$KMnum <- "KM2"
linkINFO[which(linkINFO$KMid %in% c(6)),]$KMnum <- "KM3"
linkINFO[which(linkINFO$KMid %in% c(7,9,8)),]$KMnum <- "KM4"
intersect(c("TOP2A","MKI67","SOX2","PAX6","HES1","HES5","OLIG1","PTPRC","PDGFRB","GFAP","EOMES","TUBB3","NEUROD2","NEUROD6","BCL11B","LDB2","FEZF2","DIAPH3","FOXP2","TLE4","CRYM","SATB2","NRP1","CUX2","RORB","LMO4","LPL","RELN","GAD1","DLX1","LHX6","NR2F2"),unique(linkINFO[which(linkINFO$KMnum=="KM1"),]$gene))
intersect(c("TOP2A","MKI67","SOX2","PAX6","HES1","HES5","OLIG1","PTPRC","PDGFRB","GFAP","EOMES","TUBB3","NEUROD2","NEUROD6","BCL11B","LDB2","FEZF2","DIAPH3","FOXP2","TLE4","CRYM","SATB2","NRP1","CUX2","RORB","LMO4","LPL","RELN","GAD1","DLX1","LHX6","NR2F2"),unique(linkINFO[which(linkINFO$KMnum=="KM2"),]$gene))
intersect(c("TOP2A","MKI67","SOX2","PAX6","HES1","HES5","OLIG1","PTPRC","PDGFRB","GFAP","EOMES","TUBB3","NEUROD2","NEUROD6","BCL11B","LDB2","FEZF2","DIAPH3","FOXP2","TLE4","CRYM","SATB2","NRP1","CUX2","RORB","LMO4","LPL","RELN","GAD1","DLX1","LHX6","NR2F2"),unique(linkINFO[which(linkINFO$KMnum=="KM3"),]$gene))
intersect(c("TOP2A","MKI67","SOX2","PAX6","HES1","HES5","OLIG1","PTPRC","PDGFRB","GFAP","EOMES","TUBB3","NEUROD2","NEUROD6","BCL11B","LDB2","FEZF2","DIAPH3","FOXP2","TLE4","CRYM","SATB2","NRP1","CUX2","RORB","LMO4","LPL","RELN","GAD1","DLX1","LHX6","NR2F2"),unique(linkINFO[which(linkINFO$KMnum=="KM4"),]$gene))
[1] "TOP2A"  "MKI67"  "SOX2"   "PAX6"   "HES1"   "HES5"   "PDGFRB" "GFAP"   "EOMES"  "DIAPH3" "FOXP2"  "TLE4"  
[1] "EOMES"   "NEUROD2" "NEUROD6" "BCL11B"  "FOXP2"   "SATB2"   "NRP1"    "CUX2"   
[1] "NEUROD2" "TLE4"    "SATB2"   "NRP1"    "CUX2"    "LPL"    
[1] "LDB2"  "SATB2" "NRP1"  "CUX2"  "RORB"  "LPL"
saveRDS(linkINFO, file="linkINFO.rds")
linkINFO <- readRDS("linkINFO.rds")

#########################################################################################################################################################
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(ggplot2)

geneID_geneName <- read.table("geneID_geneName.txt",col.names=c("ID","Name"))
linkINFO <- readRDS("linkINFO.rds")

df <- data_frame()
for (KMnum in c("KM1","KM2","KM3","KM4")) {
  GOgenes <- geneID_geneName[which(geneID_geneName$Name %in% unique(linkINFO[which(linkINFO$KMnum==KMnum),]$gene)),]$ID
  gene_info <- bitr(as.character(GOgenes), fromType = "ENSEMBL", toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")
  target_gene_id <- as.character(gene_info$ENTREZID)
  length(unique(gene_info$ENSEMBL))
  length(unique(target_gene_id))
  ego_BP <- enrichGO(OrgDb="org.Hs.eg.db",
                     gene = target_gene_id,
                     qvalueCutoff = 0.05,
                     ont = "BP",
                     readable = TRUE)
  
  ego_result_BP <- as.data.frame(ego_BP)
  ego_result_BP$KM <- KMnum
  df <- rbind(df,ego_result_BP)
}
df <- df[which(df$qvalue < 0.01),]

ego_fun <- function(KMgeneName,KM){
  GOgenes <- geneID_geneName[which(geneID_geneName$Name %in% KMgeneName),]$ID
  gene_info <- bitr(as.character(GOgenes), fromType = "ENSEMBL", toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")
  target_gene_id <- as.character(gene_info$ENTREZID)
  length(unique(gene_info$ENSEMBL))
  length(unique(target_gene_id))
  ego_BP <- enrichGO(OrgDb="org.Hs.eg.db",
                     gene = target_gene_id,
                     qvalueCutoff = 0.05,
                     ont = "BP",
                     readable = TRUE)
  ego_result_BP <- as.data.frame(ego_BP)
  ego_result_BP <- ego_result_BP[which(ego_result_BP$qvalue < 0.01),]
  ego_result_BP$KM <- KM
  return(ego_result_BP)
}

KM4geneName <- setdiff(linkINFO[which(linkINFO$KMnum=="KM4"),"gene"],linkINFO[which(linkINFO$KMnum!="KM4"),"gene"])
KM4GO <- ego_fun(KM4geneName,"KM4")
KM3geneName <- setdiff(linkINFO[which(linkINFO$KMnum=="KM3"),"gene"],linkINFO[which(linkINFO$KMnum!="KM3"),"gene"])
KM3GO <- ego_fun(KM3geneName,"KM3")
KM2geneName <- setdiff(linkINFO[which(linkINFO$KMnum=="KM2"),"gene"],linkINFO[which(linkINFO$KMnum!="KM2"),"gene"])
KM2GO <- ego_fun(KM2geneName,"KM2")
KM1geneName <- setdiff(linkINFO[which(linkINFO$KMnum=="KM1"),"gene"],linkINFO[which(linkINFO$KMnum!="KM1"),"gene"])
KM1GO <- ego_fun(KM1geneName,"KM1")

df <- rbind(KM1GO,KM2GO,KM3GO,KM4GO)
write.table(df, "ego_result_BP_2.xls", quote = F, sep = "\t", col.names = T, row.names = F)

