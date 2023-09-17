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

LINKobj <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/3_CREgene_pair/obj.CREgene_pair.rds")
genePeak_links <- Links(LINKobj[["peaks"]])
genePeak_links <- genePeak_links[genePeak_links$score>0]
genePeak_links <- genePeak_links[genePeak_links$gene %in% markersGene]

ht_list <- readRDS("superficialCPN_ht_list.rds")
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
linkINFO[which(linkINFO$KMid %in% c(6,7)),]$KMnum <- "KM3"
linkINFO[which(linkINFO$KMid %in% c(9,8)),]$KMnum <- "KM4"
intersect(c("TOP2A","MKI67","SOX2","PAX6","HES1","HES5","OLIG1","PTPRC","PDGFRB","GFAP","EOMES","TUBB3","NEUROD2","NEUROD6","BCL11B","LDB2","FEZF2","DIAPH3","FOXP2","TLE4","CRYM","SATB2","NRP1","CUX2","RORB","LMO4","LPL","RELN","GAD1","DLX1","LHX6","NR2F2"),unique(linkINFO[which(linkINFO$KMnum=="KM1"),]$gene))
intersect(c("TOP2A","MKI67","SOX2","PAX6","HES1","HES5","OLIG1","PTPRC","PDGFRB","GFAP","EOMES","TUBB3","NEUROD2","NEUROD6","BCL11B","LDB2","FEZF2","DIAPH3","FOXP2","TLE4","CRYM","SATB2","NRP1","CUX2","RORB","LMO4","LPL","RELN","GAD1","DLX1","LHX6","NR2F2"),unique(linkINFO[which(linkINFO$KMnum=="KM2"),]$gene))
intersect(c("TOP2A","MKI67","SOX2","PAX6","HES1","HES5","OLIG1","PTPRC","PDGFRB","GFAP","EOMES","TUBB3","NEUROD2","NEUROD6","BCL11B","LDB2","FEZF2","DIAPH3","FOXP2","TLE4","CRYM","SATB2","NRP1","CUX2","RORB","LMO4","LPL","RELN","GAD1","DLX1","LHX6","NR2F2"),unique(linkINFO[which(linkINFO$KMnum=="KM3"),]$gene))
intersect(c("TOP2A","MKI67","SOX2","PAX6","HES1","HES5","OLIG1","PTPRC","PDGFRB","GFAP","EOMES","TUBB3","NEUROD2","NEUROD6","BCL11B","LDB2","FEZF2","DIAPH3","FOXP2","TLE4","CRYM","SATB2","NRP1","CUX2","RORB","LMO4","LPL","RELN","GAD1","DLX1","LHX6","NR2F2"),unique(linkINFO[which(linkINFO$KMnum=="KM4"),]$gene))
# [1] "TOP2A"  "MKI67"  "SOX2"   "PAX6"   "HES1"   "HES5"   "PDGFRB" "GFAP"
# [9] "EOMES"  "DIAPH3" "FOXP2"  "TLE4"
# [1] "EOMES"   "NEUROD2" "NEUROD6" "BCL11B"  "FOXP2"   "SATB2"   "NRP1"
# [8] "CUX2"
# [1] "NEUROD2" "LDB2"    "TLE4"    "SATB2"   "NRP1"    "CUX2"    "RORB"
# [8] "LPL"
# [1] "LDB2"  "SATB2" "CUX2"  "RORB"  "LPL"

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

KMlist <- list()
for (KMnum in c("KM1","KM2","KM3","KM4")) {
  GOgenes <- geneID_geneName[which(geneID_geneName$Name %in% unique(linkINFO[which(linkINFO$KMnum==KMnum),]$gene)),]$ID
  gene_info <- bitr(as.character(GOgenes), fromType = "ENSEMBL", toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")
  target_gene_id <- as.factor(gene_info$ENTREZID)
  KMlist[[KMnum]] <- target_gene_id
}

ck <- compareCluster(geneCluster = KMlist, fun = "enrichGO",OrgDb="org.Hs.eg.db",ont = "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff = 0.05)
#ck <- compareCluster(geneCluster = KMlist, fun = "enrichKEGG",organism="hsa",qvalueCutoff = 0.05)
dotplot(ck)

ck_backup <- ck

df <- ck_backup@compareClusterResult
#df[which(df$Cluster=="KM1"),"Description"]
#df[which(df$Cluster=="KM2"),"Description"]
#df[which(df$Cluster=="KM3"),"Description"]
#df[which(df$Cluster=="KM4"),"Description"]
#setdiff(df[which(df$Cluster=="KM4"),"Description"],df[which(df$Cluster=="KM2"),"Description"])

#termGO <- c("DNA replication", "nuclear division", "regulation of mitotic cell cycle phase transition", "stem cell proliferation", "cell fate specification", "pallium development", "cell fate commitment", "positive regulation of cell morphogenesis involved in differentiation", "axon extension", "neuron fate commitment", "regulation of axon guidance", "positive regulation of axonogenesis", "neural precursor cell proliferation", "axon extension involved in axon guidance", "neuron projection extension involved in neuron projection guidance", "synapse assembly", "regulation of synaptic plasticity", "regulation of cation channel activity", "ionotropic glutamate receptor signaling pathway", "positive regulation of synaptic transmission", "negative regulation of neuron projection development", "negative regulation of synapse organization", "negative regulation of neuron differentiation", "regulation of glutamate receptor signaling pathway", "regulation of neurotransmitter secretion", "regulation of synaptic transmission, glutamatergic")
termGO <- c(
  "DNA replication", "nuclear division", "regulation of mitotic cell cycle phase transition", "cell fate specification", 
  "positive regulation of cell morphogenesis involved in differentiation", "axon extension", "neuron fate commitment", "regulation of axon guidance",
  "synapse assembly", "regulation of synaptic plasticity","ionotropic glutamate receptor signaling pathway", "positive regulation of synaptic transmission")

df <- df[which(df$Description %in% termGO),]
ck@compareClusterResult<- df
dotplot(ck)

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

## PLOT
ego_result_BP <- df[which(df$ID %in% c("GO:0000280","GO:1901990","GO:0019827","GO:0050775","GO:2001222","GO:0006929","GO:0007215","GO:0048167","GO:0048846")),-8]
go_enrich_df <- data.frame(Node=ego_result_BP$KM,ID=as.character(ego_result_BP$ID),
                           Description=as.character(ego_result_BP$Description),
                           Qval = ego_result_BP$qvalue,
                           GeneNumber = ego_result_BP$Count)
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
go_enrich_df$log10Q = -log10(go_enrich_df$Qval)
go_enrich_df$padj <- ifelse(go_enrich_df$Qval < 0.05, '*', '')
background<-(theme_bw()
             +theme(axis.title.y = element_text(size = 16,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.text.y =  element_text(size = 12,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.title.x = element_text(size = 16,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.text.x = element_text(size = 12,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(strip.text  = element_text(size = 10,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.line = element_line(size = 0.5),legend.background = element_blank())
             +theme(strip.background = element_rect(fill=NA,colour = NA))
)
ggplot(data = go_enrich_df, aes(x = number, y = log10Q, fill=Node)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() +
  scale_fill_manual(values = c("#8DD3C7","#FDB462","#BEBADA","#80B1D3")) +
  scale_x_discrete(labels=rev(go_enrich_df$Description)) +
  xlab("") + ylab('-Log10(BH)') +
  background +
  theme(legend.background = element_blank(),
        legend.text=element_text(size=14),
        legend.position="top",
        legend.title=element_blank(),
        legend.key.size=unit(0.3,'cm')) +
  geom_text(aes(label = padj), vjust=0.75, hjust = -0.25, color = 1)


#########################################################################################################################################################
# /home/yuchen/miniconda3/envs/R4.0/bin/R
linkINFO <- readRDS("linkINFO.rds")
linkINFO$KMnum <- "KM"
linkINFO[which(linkINFO$KMid %in% c(12,15,18,17,19,20,14,16,13)),]$KMnum <- "KM1"
linkINFO[which(linkINFO$KMid %in% c(3,1,2,9)),]$KMnum <- "KM2"
linkINFO[which(linkINFO$KMid %in% c(5,6,4)),]$KMnum <- "KM3"
write.table(linkINFO[,c("peak","KMnum")], "linkedPeak_superficialCPN.txt", quote = F, sep = "\t", col.names = F, row.names = F)

cat linkedPeak_superficialCPN.txt | grep "KM1" | cut -f1 | sed 's/-/\t/g' > KM1_peak.bed
cat linkedPeak_superficialCPN.txt | grep "KM2" | cut -f1 | sed 's/-/\t/g' > KM2_peak.bed
cat linkedPeak_superficialCPN.txt | grep "KM3" | cut -f1 | sed 's/-/\t/g' > KM3_peak.bed
cat linkedPeak_superficialCPN.txt | cut -f1 | sed 's/-/\t/g' > universe_peak.bed

## Obtaining Necessary Data
/home/yuhan/Software/CBGR-jaspar_enrichment/bin/zenodo_fetch 5555937 -f hg38.tar.gz
/home/yuhan/Software/CBGR-jaspar_enrichment/bin/zenodo_fetch 5555932 -f sacCer3.tar.gz
tar -zxvf hg38.tar.gz

## run in Management Node
cd /home/yuhan/Software/CBGR-jaspar_enrichment
Sbed="/data3/yuhan/Project/Neuron/scMultiome/5_regulation/KM1_peak.bed"
Ubed="/data3/yuhan/Project/Neuron/scMultiome/5_regulation/universe_peak.bed"
outputdir="/data3/yuhan/Project/Neuron/scMultiome/5_regulation/JASPAREnrichmentAnalysis"
loladb_dir="/home/yuhan/Software/CBGR-jaspar_enrichment/sacCer3"
API_URL="https://jaspar.genereg.net/api/v1/matrix/"
n_cores=10
/home/yuhan/Software/CBGR-jaspar_enrichment/bin/JASPAR_enrich.sh oneSetBg ${loladb_dir} ${Sbed} ${Ubed} ${outputdir} ${API_URL} ${n_cores}




