#!/bin/bash

cat /data3/yuhan/Project/Neuron/scMultiome/3_CREgene_pair/CRE_foldEnrichment/genePeak_links.txt | awk 'NR!=1{print $3"-"$2;}' | awk 'BEGIN{FS="-";OFS="\t";}{print $1,$2,$3,$0;}' | grep -v "chrY" > genePeak_links.txt

#########################################################################################################################################
featureFile=genePeak_links.txt
stateFile=E081_25_imputed12marks_hg38lift_dense.bed

[(number of bases in state AND overlap feature)/(number of bases in genome)]/[(number of bases overlap feature)/(number of bases in genome) * (number of bases in state)/(number of bases in genome)]
(NUMbases_stateANDfeature/3031042417)/((16719920/3031042417)*(state/3031042417))

############################################################
zcat /data1/yuhan/BasicBrowser/RoadMap_ChromatinState/${stateFile}.gz | awk 'BEGIN{OFS="\t";}$1~"chr"{print $1,$2,$3,$4;}' | sed "s/'//" | grep -v "chrY" > ${stateFile}
cat /data2/yuhan/Project/3D_evolution/human_gained_lost_enhancer/forebrain/forebrain_enhancer.hg38.bed | grep -v "chrY" | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,"26_forebrainEnh"}' >> ${stateFile}
cat /data2/yuhan/Project/3D_evolution/human_gained_lost_enhancer/forebrain/nonforebrain_enhancer.hg38.bed | grep -v "chrY" | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,"27_nonforebrainEnh"}' >> ${stateFile}

## number of bases in state
for state in {1..27}
do
cat ${stateFile} | awk -v STATEnum=${state} 'BEGIN{SUM=0;OFS="\t";}{split($4,a,"_");if(a[1]==STATEnum) {SUM=SUM+($3-$2);STATEid=$4;}} END{print STATEid,SUM;}' >> NUMbases_state.txt
done

## number of bases overlap feature: 16719920
cat ${featureFile} | awk 'BEGIN{SUM=0;}{SUM=SUM+($3-$2);}END{print SUM;}'

## number of bases in genome: 3031042417
cat /data/public/refGenome/bwa_index/hg38/ChromInfo.txt | grep -v "chrY" | awk 'BEGIN{SUM=0;}{SUM=SUM+$2;}END{print SUM;}'

bedtools="/data/public/software/bedtools.2.25.0/bin/bedtools"
for state in {1..27}
do
cat ${stateFile} | awk -v STATEnum=${state} '{split($4,a,"_");if(a[1]==STATEnum) print $0;}' | ${bedtools} intersect -a ${featureFile} -b - -wo | awk 'BEGIN{SUM=0;OFS="\t";}{SUM=SUM+$9;STATEid=$8;}END{print STATEid,SUM;}' >> NUMbases_stateANDfeature.txt
done

paste NUMbases_stateANDfeature.txt NUMbases_state.txt | awk 'BEGIN{OFS="\t";}{print $1,($2/3031042417)/((16719920/3031042417)*($4/3031042417))}' > fold_enrichment_of_overlap_ChromatinState.txt

############################################################
require(ggplot2)
background<-(theme_bw()
             +theme(plot.title = element_text(size = 16,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.title.y = element_text(size = 16,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.text.y =  element_text(size = 12,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.title.x = element_text(size = 16,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.text.x = element_text(size = 12,vjust = 0.5, hjust = 0.5, colour = "black",angle = 90))
             +theme(strip.text  = element_text(size = 10,vjust = 0.5, hjust = 0.5, colour = "black"))
             +theme(axis.line = element_line(size = 0.5),legend.background = element_blank())
             +theme(strip.background = element_rect(fill=NA,colour = NA))
)

dfplot <- read.table("fold_enrichment_of_overlap_ChromatinState.txt",sep="\t")
names(dfplot) <- c("state","foldEnrichment")
dfplot <- within(dfplot,{state <- factor(state,levels=dfplot$state)})

summary(log10(dfplot$foldEnrichment))
ggplot(dfplot,aes(x=state,y=log10(foldEnrichment))) +
  geom_bar(stat="identity") +
  xlab("") + ylab("log10(fold enrichment)") + labs(title = "linked peak") +
  scale_y_continuous(limits=c(-1.5,2),breaks=seq(-1.5,2,0.5)) +
  background +
  theme(panel.grid =element_blank(),panel.background=element_rect(fill='transparent', color="#000000"),
        axis.ticks.x = element_blank())

