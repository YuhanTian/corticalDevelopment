# /home/yuchen/miniconda3/envs/R4.0/bin/R
library(Seurat)
library(Signac)
library(presto)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)

obj <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/3_CREgene_pair/obj.CREgene_pair.rds")
CREgene_pair <- Links(obj[["peaks"]])
CREgene_pair <- as.data.frame(CREgene_pair)
# > dim(CREgene_pair)
# [1] 17000    10
# > dim(CREgene_pair[which(CREgene_pair$score>0),])
# [1] 15269    10
# > dim(CREgene_pair[which(CREgene_pair$score<0),])
# [1] 1731   10
write.table(CREgene_pair,"CREgene_pair.txt", col.names=F, row.names=F, sep="\t", quote=F)

#########################################################################################################################################################
zcat /home/yuhan/Software/cellranger-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz | awk -F'[ \t]' 'BEGIN{OFS="\t";}$3=="gene"{print $1"-"$4"-"$5,$7,$16;}' | sed -e 's/"//g' -e 's/;//g' | awk 'BEGIN{OFS="\t";}{if($2=="+") $2="P";else $2="N";print $0;}' > tmp1.txt
cat CREgene_pair.txt | cut -f6,7,8 > tmp2.txt

awk 'BEGIN{OFS="\t";}NR==FNR{a[$3]=$0;next}{if($2 in a) print a[$2],$1,$3;}' tmp1.txt tmp2.txt | awk '{$3=0;print $0;}' | sed 's/-/\t/g' | awk 'BEGIN{OFS="\t";}{if($4=="P") print $1,$2,$2+1,$7,$8,$9,$6;else print $1,$3-1,$3,$7,$8,$9,$6;}' | awk -F'[ \t]' 'BEGIN{OFS="\t";}{print $1,$2,$3,$4":"$5"-"$6","$7;print $4,$5,$6,$1":"$2"-"$3","$7;}' | sort -k1,1 -k2,2n > CREgene_pair.bedpe
bgzip CREgene_pair.bedpe
tabix -p bed CREgene_pair.bedpe.gz

awk 'BEGIN{OFS="\t";}NR==FNR{a[$3]=$1;next}{if($2 in a) print a[$2],$2;}' tmp1.txt tmp2.txt | sed 's/-/\t/g' | sort -k1,1 -k2,2n | uniq > CREgene_pair.gene.bed
bgzip CREgene_pair.gene.bed
tabix -p bed CREgene_pair.gene.bed.gz

cat CREgene_pair.txt | cut -f8 | sed 's/-/\t/g' | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,$1"-"$2"-"$3;}' | sort -k1,1 -k2,2n | uniq > CREgene_pair.CRE.bed
bgzip CREgene_pair.CRE.bed
tabix -p bed CREgene_pair.CRE.bed.gz


zcat /home/yuhan/Software/cellranger-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz | awk -F'[ \t]' 'BEGIN{OFS="\t";}$3=="gene"&&($14~"protein_coding"||$14~"lncRNA"){print $16,$14;}' | sed -e 's/"//g' -e 's/;//g' | grep "" > tmp3.txt
