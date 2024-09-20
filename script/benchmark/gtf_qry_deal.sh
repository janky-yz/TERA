#!/bin/bash
sample=$1
ref=$2
ref_TE=$3

basepath=$(cd `dirname $0`; pwd)

source ~/miniconda3/bin/activate gff
gffread --sort-alpha -T -o ${sample}.sorted.gtf ${sample}.gtf

grep exon ${sample}.sorted.gtf | awk '{if($9=="transcript_id"){print $1,$4-1,$5,$10,$12,$7}else{print $1,$4-1,$5,$12,$10,$7}}' OFS='\t' | tr -d '"' | tr -d ';' | sort -k1V,1 -k2n,2 -k3n,3 >${sample}_exon.bed
Rscript ${basepath}/exon_deal.R ${sample}

bedtools intersect -a ${sample}_exon.bed -b ${ref_TE} -s -wo | awk '{if($NF>20){print $0}}' >${sample}_TE_exon.txt
cut -f1-7 ${sample}_TE_exon.txt | uniq >${sample}_TE_exon.bed
Rscript ${basepath}/TE_exon_deal.R ${sample}

bedtools intersect -a ${sample}_TE_exon.bed -b ${ref}_TE_exon.bed -s -wo | awk '{if($NF>20 && $7==$14){print $0}}' >${sample}_ref_TE_exon_overlap.txt
bedtools intersect -a ${sample}_TE_exon.bed -b ${sample}_TE_exon.bed -s -wo | awk '{if($NF>20 && $7==$14){print $0}}' >${sample}_TE_exon_self_overlap.txt
bedtools intersect -a ${sample}_TE_all_exon.bed -b ${ref}_TE_all_exon.bed -s -wo | awk '{if($8==$16 && $7==$15){print $0}}' >${sample}_ref_TE_transcript_overlap.txt
bedtools intersect -a ${sample}_TE_all_exon.bed -b ${sample}_TE_all_exon.bed -s -wo | awk '{if($8==$16 && $7==$15){print $0}}' >${sample}_TE_transcript_self_overlap.txt
cut -f4 ${sample}_TE_exon.bed | sort | uniq >${sample}_TET_ID.list

Rscript ${basepath}/exon_class.R ${sample} 5
Rscript ${basepath}/tran_class.R ${sample} 5

Rscript ${basepath}/Extract_TE_gtf.R ${sample}
grep exon ${sample}.TE.gtf | cut -d '"' -f2 | uniq -u >${sample}.TE.SE.list

rm ${sample}_*self_overlap.txt
