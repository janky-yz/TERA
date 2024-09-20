#!/bin/bash
basepath=$(cd `dirname $0`; pwd)
sample=$1
grep -w transcript ${sample}.TE.gtf | awk '{print $1,$4-1,$5,$10,$12,$7}' OFS='\t' | tr -d '"' | tr -d ';' | sort -k1V,1 -k2n,2 -k3n,3 >${sample}.TE.tran.sorted.bed
grep -w exon ${sample}.TE.gtf | awk '{print $1,$4-1,$5,$10,$12,$7}' OFS='\t' | tr -d '"' | tr -d ';' | sort -k1V,1 -k2n,2 -k3n,3 >${sample}.TE.exon.sorted.bed
bedtools intersect -a ${sample}.TE.tran.sorted.bed -b /media/data6/sjq/HERV/HEK293T/T2T/chm13v2.0.gene.sort.bed -v >${sample}.TE.intergenic.bed
bedtools intersect -a ${sample}.TE.tran.sorted.bed -b /media/data6/sjq/HERV/HEK293T/T2T/chm13v2.0.gene.sort.bed -s -v >${sample}.TE.intergenic.antisense.bed
cat ${sample}.TE.intergenic.bed ${sample}.TE.intergenic.antisense.bed | sort | uniq -u >${sample}.TE.antisense.bed
bedtools intersect -a ${sample}.TE.tran.sorted.bed -b /media/data6/sjq/HERV/HEK293T/T2T/chm13v2.0.allexon.collapse.bed -s -v >${sample}.TE.intergenic.antisense.intronic.bed
cat ${sample}.TE.intergenic.antisense.bed ${sample}.TE.intergenic.antisense.intronic.bed | sort | uniq -u >${sample}.TE.intronic.bed
#rm U251.TE.intergenic.antisense.intronic.bed  U251.TE.intergenic.antisense.bed

bedtools intersect -a /media/data6/sjq/HERV/HEK293T/T2T/chm13v2.nrph.TE.anno.bed -b ${sample}.TE.exon.sorted.bed -s -wo | awk '{if($NF>20){print $0}}' >${sample}.TE.overlap.txt
Rscript ${basepath}/overlap.R ${sample}

bedtools intersect -a ${sample}.TE.overlap.pos.bed -b /media/data6/sjq/HERV/HEK293T/T2T/chm13v2.0.allexon.collapse.bed -s -f 1 -v >${sample}.TE.chimeric.bed
Rscript ${basepath}/anno.R ${sample}
