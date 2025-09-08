args<-commandArgs(T)
qry_id=args[1]

exon_bed <- read.table(paste0(qry_id, "_gene_exon.bed"), stringsAsFactors=F)
TE_gene <- read.table(paste0(qry_id, "_TE_gene_ID.list"), stringsAsFactors=F)

TE_exon_bed <- exon_bed[exon_bed[,4] %in% TE_gene[,1],]
exon_num <- as.data.frame(table(TE_exon_bed[,4]))
out <- merge(TE_exon_bed, exon_num, by.x="V4", by.y="Var1", sort=F)
out <- out[,c(2:4,1,5:8)]
out <- out[order(out[,1], out[,2], out[,3]),]
write.table(out, file=paste0(qry_id, "_TE_gene_exon.bed"), sep='\t', col.names=F, row.names=F, quote=F)
