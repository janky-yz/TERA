#!/bin/bash
sample=$1
ref=$2
ref_TE=$3

basepath=$(cd `dirname $0`; pwd)

source ~/miniconda3/bin/activate gff
gffread --sort-alpha -T -o ${sample}.sorted.gtf ${sample}.gtf

grep exon ${sample}.sorted.gtf | awk '{if($9=="transcript_id"){print $1,$4-1,$5,$10,$12,$7}else{print $1,$4-1,$5,$12,$10,$7}}' OFS='\t' | tr -d '"' | tr -d ';' | sort -k1V,1 -k2n,2 -k3n,3 >${sample}_exon.bed
cp ${sample}_exon.bed ${sample}.TE.exon.sorted.bed

Rscript ${basepath}/exon_deal.R ${sample}

cp ${sample}_TE_exon.bed ${sample}.TE.exon.bed
cut -f 1-7 ${sample}.TE.exon.bed >${sample}_TE_exon.bed
cut -f 1-6 ${sample}_TE_exon.bed >${sample}.TE.exon.bed

bedtools intersect -a ${sample}_TE_exon.bed -b ${ref_TE} -s -wo | awk '{if($NF>20){print $0}}' >${sample}_TE_exon.txt
#Rscript ${basepath}/TE_exon_deal.R ${sample}

bedtools intersect -a ${sample}_TE_exon.bed -b ${ref}_TE_exon.bed -s -wo | awk '{if($NF>20 && $7==$14){print $0}}' >${sample}_ref_TE_exon_overlap.txt
bedtools intersect -a ${sample}_TE_exon.bed -b ${sample}_TE_exon.bed -s -wo | awk '{if($NF>20 && $7==$14){print $0}}' >${sample}_TE_exon_self_overlap.txt
#bedtools intersect -a ${sample}_TE_all_exon.bed -b ${ref}_TE_all_exon.bed -s -wo | awk '{if($8==$16 && $7==$15){print $0}}' >${sample}_ref_TE_transcript_overlap.txt
#bedtools intersect -a ${sample}_TE_all_exon.bed -b ${sample}_TE_all_exon.bed -s -wo | awk '{if($8==$16 && $7==$15){print $0}}' >${sample}_TE_transcript_self_overlap.txt
cut -f4 ${sample}_TE_exon.bed | sort | uniq >${sample}_TET_ID.list

Rscript ${basepath}/Extract_TE_gtf.R ${sample}
Rscript ${basepath}/exon_class.R ${sample} 5
#Rscript ${basepath}/tran_class.R ${sample} 5

grep exon ${sample}.TE.gtf | cut -d '"' -f2 | uniq -u >${sample}.TE.SE.list

rm ${sample}_*self_overlap.txt

grep -w transcript ${sample}.TE.gtf | awk '{print $1,$4-1,$5,$10,$12,$7}' OFS='\t' | tr -d '"' | tr -d ';' | sort -k1V,1 -k2n,2 -k3n,3 >${sample}.TE.tran.sorted.bed
bedtools intersect -a ${sample}.TE.tran.sorted.bed -b /media/data6/sjq/HERV/HEK293T/T2T/chm13v2.0.gene.sort.bed -v >${sample}.TE.intergenic.bed
bedtools intersect -a ${sample}.TE.tran.sorted.bed -b /media/data6/sjq/HERV/HEK293T/T2T/chm13v2.0.gene.sort.bed -s -v >${sample}.TE.intergenic.antisense.bed
cat ${sample}.TE.intergenic.bed ${sample}.TE.intergenic.antisense.bed | sort | uniq -u >${sample}.TE.antisense.bed
bedtools intersect -a ${sample}.TE.tran.sorted.bed -b /media/data6/sjq/HERV/HEK293T/T2T/chm13v2.0.allexon.collapse.bed -s -v >${sample}.TE.intergenic.antisense.intronic.bed
cat ${sample}.TE.intergenic.antisense.bed ${sample}.TE.intergenic.antisense.intronic.bed | sort | uniq -u >${sample}.TE.intronic.bed

bedtools intersect -a /media/data6/sjq/HERV/HEK293T/T2T/chm13v2.nrph.TE.anno.bed -b ${sample}.TE.exon.sorted.bed -s -wo | awk '{if($NF>20){print $0}}' >${sample}.TE.overlap.txt
Rscript ${basepath}/overlap.R ${sample}

bedtools intersect -a ${sample}.TE.overlap.pos.bed -b /media/data6/sjq/HERV/HEK293T/T2T/chm13v2.0.allexon.collapse.bed -s -f 1 -v >${sample}.TE.chimeric.bed
Rscript ${basepath}/anno.R ${sample}
