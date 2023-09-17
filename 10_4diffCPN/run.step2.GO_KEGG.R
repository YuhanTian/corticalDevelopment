library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(ggplot2)


markergeneDF <- readRDS("markergeneDF.rds")
GENEid <- geneID_geneName[which(geneID_geneName$Name %in% markergeneDF[which(markergeneDF$type=="down"),]$gene),]$ID

## GO
gene_info <- bitr(as.character(GENEid), fromType = "ENSEMBL", toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")
target_gene_id <- as.character(gene_info$ENTREZID)
length(unique(gene_info$ENSEMBL))
length(unique(target_gene_id))
ego_BP <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_id,
                   qvalueCutoff = 0.05,
                   ont = "BP",
                   readable = TRUE)

ego_result_BP <- as.data.frame(ego_BP)


## KEGG
gene_info <- bitr(as.character(GENEid), fromType = "ENSEMBL", toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")
length(unique(gene_info$ENSEMBL))
target_gene_id <- as.character(gene_info[, 2])
length(target_gene_id)

KEGG_result <- enrichKEGG(
  gene = target_gene_id,
  keyType = "kegg",
  organism = 'hsa',
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
KEGG_result_df <- as.data.frame(KEGG_result)
KEGG_result_df[which(KEGG_result_df$qvalue<0.05),]$Description





