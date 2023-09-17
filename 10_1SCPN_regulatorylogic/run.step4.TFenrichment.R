#########################################################################################################################################################
# TF binding motif enrichment analysis
conda activate scMultiome
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(JASPAR2020)
require("motifmatchr")
require('ggseqlogo')
require("ggplot2")
library(Seurat)
library(Signac)

background<-(theme_bw()
             +theme(plot.title = element_text(size = 24,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.title.y = element_text(size = 16,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.text.y =  element_text(size = 12,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.title.x = element_text(size = 16,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.text.x = element_text(size = 12,vjust = 1, hjust = 1, colour = "black", angle = 45))
             +theme(strip.text  = element_text(size = 10,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(panel.grid =element_blank(),panel.background=element_rect(fill='transparent', color="#000000"))
             +theme(axis.line = element_line(size = 0.5),legend.background = element_blank())
             +theme(strip.background = element_rect(fill=NA,colour = NA))
)

seurat <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/5_EN_regulatory_logic.newVersion.cleanedcell/SCPNbranch_TFmotifenrichment/findMarker/genes_motif_trend/seurat_TFmotif.rds")
DefaultAssay(seurat) <- "peaks"

linkINFO <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/5_EN_regulatory_logic.newVersion.cleanedcell/SCPN_regulatorylogic/linkINFO.rds")
linkINFO$KMnum <- "KM"
linkINFO[which(linkINFO$KMid %in% c(1,2,8,7,9,10,6,5,12)),]$KMnum <- "KM1"
linkINFO[which(linkINFO$KMid %in% c(11,16)),]$KMnum <- "KM2"
linkINFO[which(linkINFO$KMid %in% c(17,13,14,15)),]$KMnum <- "KM3"

#########################################################################################################################################################
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
df_pfm <- data.frame(t(sapply(pfm, function(x)
  c(id=x@ID, name=x@name, symbol=ifelse(!is.null(x@tags$symbol),x@tags$symbol,NA)))))

# seurat <- AddMotifs(seurat, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
open_peaks <- AccessiblePeaks(seurat)
sample1000 <- sample(c(1:nrow(linkINFO)),1000)
peaks_matched <- MatchRegionStats(meta.feature = seurat[['peaks']]@meta.features[open_peaks, ],
                                  query.feature = seurat[['peaks']]@meta.features[linkINFO[sample1000,"peak"], ],
                                  n = 50000)

## KM1
motif_enrichment_KM1 <- FindMotifs(seurat,
                                    features = linkINFO$peak[linkINFO$KMnum == "KM1"],
                                    background = peaks_matched) %>%
  mutate(symbol = setNames(ifelse(is.na(df_pfm$symbol), df_pfm$name, df_pfm$symbol), df_pfm$id)[motif]) %>%
  mutate(padj = p.adjust(pvalue, method="BH"))
motif_enrichment_KM1[which(motif_enrichment_KM1$padj < 0.01 & motif_enrichment_KM1$fold.enrichment > 2),"motif.name"]
## KM2
motif_enrichment_KM2 <- FindMotifs(seurat,
                                    features = linkINFO$peak[linkINFO$KMnum == "KM2"],
                                    background = peaks_matched) %>%
  mutate(symbol = setNames(ifelse(is.na(df_pfm$symbol), df_pfm$name, df_pfm$symbol), df_pfm$id)[motif]) %>%
  mutate(padj = p.adjust(pvalue, method="BH"))
motif_enrichment_KM2[which(motif_enrichment_KM2$padj < 0.01 & motif_enrichment_KM2$fold.enrichment > 2),"motif.name"]
## KM3
motif_enrichment_KM3 <- FindMotifs(seurat,
                                    features = linkINFO$peak[linkINFO$KMnum == "KM3"],
                                    background = peaks_matched) %>%
  mutate(symbol = setNames(ifelse(is.na(df_pfm$symbol), df_pfm$name, df_pfm$symbol), df_pfm$id)[motif]) %>%
  mutate(padj = p.adjust(pvalue, method="BH"))
motif_enrichment_KM3[which(motif_enrichment_KM3$padj < 0.01 & motif_enrichment_KM3$fold.enrichment > 2),"motif.name"]

#########################################################################################################################################################
dfplot <- rbind(data.frame(head(motif_enrichment_KM1[which(motif_enrichment_KM1$padj < 0.01 & motif_enrichment_KM1$fold.enrichment > 2),c("motif.name","pvalue")],10),tpye="KM1"),
data.frame(head(motif_enrichment_KM2[which(motif_enrichment_KM2$padj < 0.01 & motif_enrichment_KM2$fold.enrichment > 2),c("motif.name","pvalue")],15),tpye="KM2"),
data.frame(head(motif_enrichment_KM3[which(motif_enrichment_KM3$padj < 0.01 & motif_enrichment_KM3$fold.enrichment > 2),c("motif.name","pvalue")],11),tpye="KM3"))

dfplot$TFenrichment <- -log10(dfplot$pvalue)
dfplot <- dfplot[which(!(rownames(dfplot)=="MA0826.1")),]
dfplot <- dfplot[which(!(rownames(dfplot)=="MA0607.1")),]
dfplot <- dfplot[which(!(rownames(dfplot)=="MA0461.2")),]
dfplot <- dfplot[which(!(rownames(dfplot)=="MA0678.1")),]
dfplot <- dfplot[which(!(rownames(dfplot)=="MA0827.1")),]
dfplot <- dfplot[which(!(rownames(dfplot)=="MA0623.21")),]

pdf("dotplot.pdf",height=10,width=10)
(ggplot(dfplot,aes(factor(tpye,levels = c("KM1","KM2","KM3")),factor(motif.name,levels = rev(dfplot$motif.name))))
  +geom_point(aes(color=tpye,size=TFenrichment),na.rm=T)
  +scale_size_area(max_size=7,breaks = seq(0, 100, 25),limits=c(0,100),na.value=7)
  +scale_color_manual(values=c("KM1"="#F7FCB9","KM2"="#A1D99B","KM3"="#74ADD1"))
  +theme_bw(base_size=16)+labs(x="Enrichment -log10(pvalue)",y="",title = "TF motif enrichment (linked peaks)") + background
)
dev.off()

rm(list=setdiff(ls(),c("background","dfplot","motif_enrichment_KM1","motif_enrichment_KM2","motif_enrichment_KM3","sample1000")))
save.image(file="figure5B.TFmotif.RData")
load(file="figure5B.TFmotif.RData")
