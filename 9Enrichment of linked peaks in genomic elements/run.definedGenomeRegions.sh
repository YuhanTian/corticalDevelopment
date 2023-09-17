#!/bin/bash

#########################################################################################################################################
bedtools="/data/public/software/bedtools.2.25.0/bin/bedtools"
inputGTF=/data1/tang/3D_evolution/primates_annotation/hg38/Homo_sapiens.GRCh38.97_withchr.gtf

cat ${inputGTF} | grep -v "^#" | awk -F'[\t ]' 'BEGIN{OFS="\t";}$1~/chr[0-9X]+/&&$3=="gene"{print $1,$4,$5;}' | sort -k1,1 -k2,2n | ${bedtools} merge -i - | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,"gene";}' > gene.txt
cat ${inputGTF} | grep -v "^#" | awk -F'[\t ]' 'BEGIN{OFS="\t";}$1~/chr[0-9X]+/&&$3=="gene"&&$18=="\"protein_coding\";"{print $1,$4,$5;}' | sort -k1,1 -k2,2n | ${bedtools} merge -i - | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,"protein_coding";}' > protein_coding.txt
cat ${inputGTF} | grep -v "^#" | awk -F'[\t ]' 'BEGIN{OFS="\t";}$1~/chr[0-9X]+/&&$3=="gene"&&$18=="\"lncRNA\";"{print $1,$4,$5;}' | sort -k1,1 -k2,2n | ${bedtools} merge -i - | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,"lncRNA";}' > lncRNA.txt
cat ${inputGTF} | grep -v "^#" | awk -F'[\t ]' 'BEGIN{OFS="\t";}$1~/chr[0-9X]+/&&$3=="gene"&&$18=="\"miRNA\";"{print $1,$4,$5;}' | sort -k1,1 -k2,2n | ${bedtools} merge -i - | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,"miRNA";}' > miRNA.txt
cat ${inputGTF} | grep -v "^#" | awk -F'[\t ]' 'BEGIN{OFS="\t";}$1~/chr[0-9X]+/&&$3=="exon"{print $1,$4,$5;}' | sort -k1,1 -k2,2n | ${bedtools} merge -i - | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,"exon";}' > exon.txt

cat /data/public/refGenome/bwa_index/hg38/ChromInfo.txt | awk 'BEGIN{OFS="\t";}$1~/chr[0-9X]+/{print $0;}' | sort -k1,1 -k2,2n > chromSize.txt
${bedtools} complement -i gene.txt -g chromSize.txt > intergenic.tmp.txt
cat intergenic.tmp.txt | awk '$2!=0{print $0;}' | awk 'NR==FNR{a[$1]=$2;next;}{if(a[$1]!=$3) print $0}' chromSize.txt - | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,"intergenic";}' > intergenic.txt

${bedtools} complement -i exon.txt -g chromSize.txt > intron.tmp.txt
${bedtools} subtract -a intron.tmp.txt -b intergenic.tmp.txt | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,"intron";}' > intron.txt

cat ${inputGTF} | grep -v "^#" | awk -F'[\t ]' 'BEGIN{OFS="\t";}$1~/chr[0-9X]+/&&$3=="gene"{if($7=="+") print $1,$4-1000,$4+1000;else print $1,$5-1000,$5+1000;}' | sort -k1,1 -k2,2n | ${bedtools} merge -i - | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,"TSS";}' > TSS.txt
cat ${inputGTF} | grep -v "^#" | awk -F'[\t ]' 'BEGIN{OFS="\t";}$1~/chr[0-9X]+/&&$3=="gene"{if($7=="+") print $1,$4-2000,$4+1000;else print $1,$5-2000,$5+1000;}' | sort -k1,1 -k2,2n | ${bedtools} merge -i - | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,"promoter";}' > promoter.txt
cat ${inputGTF} | grep -v "^#" | awk -F'[\t ]' 'BEGIN{OFS="\t";}$1~/chr[0-9X]+/&&$3=="gene"{if($7=="-") print $1,$4-1000,$4+1000;else print $1,$5-1000,$5+1000;}' | sort -k1,1 -k2,2n | ${bedtools} merge -i - | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,"TES";}' > TES.txt

cat ${inputGTF} | grep -v "^#" | awk -F'[\t ]' 'BEGIN{OFS="\t";}$1~/chr[0-9X]+/&&$3=="five_prime_utr"{print $1,$4,$5;}' | sort -k1,1 -k2,2n | ${bedtools} merge -i - | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,"fiveUTR";}' > fiveUTR.txt
cat ${inputGTF} | grep -v "^#" | awk -F'[\t ]' 'BEGIN{OFS="\t";}$1~/chr[0-9X]+/&&$3=="three_prime_utr"{print $1,$4,$5;}' | sort -k1,1 -k2,2n | ${bedtools} merge -i - | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,"threeUTR";}' > threeUTR.txt

#########################################################################################################################################
cat /data3/yuhan/Project/Neuron/scMultiome/3_CREgene_pair/CRE_foldEnrichment/genePeak_links.txt | awk 'NR!=1{print $3"-"$2;}' | awk 'BEGIN{FS="-";OFS="\t";}{print $1,$2,$3,$0;}' | grep -v "chrY" > genePeak_links.txt

featureFile=genePeak_links.txt

[(number of bases in state AND overlap feature)/(number of bases in genome)]/[(number of bases overlap feature)/(number of bases in genome) * (number of bases in state)/(number of bases in genome)]
(NUMbases_stateANDfeature/3031042417)/((16719920/3031042417)*(state/3031042417))

############################################################
## number of bases in state
for state in gene protein_coding lncRNA miRNA exon intergenic intron TSS promoter TES fiveUTR threeUTR
do
cat ${state}.txt | awk -v STATEid=${state} 'BEGIN{SUM=0;OFS="\t";} {SUM=SUM+($3-$2);} END{print STATEid,SUM;}' >> NUMbases_state.txt
done

## number of bases overlap feature: 16719920
cat ${featureFile} | awk 'BEGIN{SUM=0;}{SUM=SUM+($3-$2);}END{print SUM;}'

## number of bases in genome: 3031042417
cat /data/public/refGenome/bwa_index/hg38/ChromInfo.txt | grep -v "chrY" | awk 'BEGIN{SUM=0;}{SUM=SUM+$2;}END{print SUM;}'

bedtools="/data/public/software/bedtools.2.25.0/bin/bedtools"
for state in gene protein_coding lncRNA miRNA exon intergenic intron TSS promoter TES fiveUTR threeUTR
do
cat ${state}.txt | ${bedtools} intersect -a ${featureFile} -b - -wo | awk 'BEGIN{SUM=0;OFS="\t";}{SUM=SUM+$9;STATEid=$8;}END{print STATEid,SUM;}' >> NUMbases_stateANDfeature.txt
done

paste NUMbases_stateANDfeature.txt NUMbases_state.txt | awk 'BEGIN{OFS="\t";}{print $1,($2/3031042417)/((16719920/3031042417)*($4/3031042417))}' > fold_enrichment_of_overlap_definedGenomeRegions.txt

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

setwd("C:/Users/dell/Desktop/Neuron/plot")
dfplot <- read.table("fold_enrichment_of_overlap_definedGenomeRegions.txt",sep="\t")
names(dfplot) <- c("state","foldEnrichment")
dfplot <- within(dfplot,{state <- factor(state,levels=c("gene", "protein_coding", "lncRNA", "miRNA", "exon", "intergenic", "intron", "TSS", "promoter", "TES", "fiveUTR", "threeUTR"))})

summary(log10(dfplot$foldEnrichment))
ggplot(dfplot,aes(x=state,y=log10(foldEnrichment))) +
  geom_bar(stat="identity") +
  xlab("") + ylab("log10(fold enrichment)") + labs(title = "linked peak") +
  scale_y_continuous(limits=c(-0.2,1.3),breaks=seq(-0.2,1.3,0.2)) +
  background +
  theme(panel.grid =element_blank(),panel.background=element_rect(fill='transparent', color="#000000"),
        axis.ticks.x = element_blank())

