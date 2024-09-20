args<-commandArgs(T)
sample=args[1]

overlap_file = paste0(sample, ".TE.overlap.txt")
overlap = read.table(overlap_file, stringsAsFactors=F)
overlap$TE_ID <- paste(overlap[,1], overlap[,2], overlap[,3], overlap[,4], overlap[,6], sep='_')

intergenic_file = paste0(sample, ".TE.intergenic.bed")
intronic_file = paste0(sample, ".TE.intronic.bed")
antisense_file = paste0(sample, ".TE.antisense.bed")

intergenic = read.table(intergenic_file, stringsAsFactors=F)
intronic = read.table(intronic_file, stringsAsFactors=F)
antisense = read.table(antisense_file, stringsAsFactors=F)

chimeric_file = paste0(sample, ".TE.chimeric.bed")
chimeric = read.table(chimeric_file, stringsAsFactors=F)
chimeric$TE_ID <- paste(chimeric[,7], chimeric[,8], chimeric[,9], chimeric[,10], chimeric[,12], sep='_')
chimeric_TE_ID <- unique(chimeric$TE_ID)

intergenic_TE_ID <- unique(overlap[overlap[,12] %in% intergenic[,4],]$TE_ID)
intronic_TE_ID <- unique(overlap[overlap[,12] %in% intronic[,4],]$TE_ID)
antisense_TE_ID <- unique(overlap[overlap[,12] %in% antisense[,4],]$TE_ID)

TE_element <- unique(overlap[,c(1:8,16)])
names(TE_element) <- c("chrom", "start", "end", "element", "family", "strand", "superfamily", "class", "TE_ID")

TE_element$type = "exonic"
TE_element[TE_element$TE_ID %in% chimeric_TE_ID,]$type = "chimeric"
TE_element[TE_element$TE_ID %in% intronic_TE_ID,]$type = "intronic"
TE_element[TE_element$TE_ID %in% intergenic_TE_ID,]$type = "intergenic"
TE_element[TE_element$TE_ID %in% antisense_TE_ID,]$type = "antisense"

element_file = paste0(sample, ".TE.element.anno.txt")
write.table(TE_element, file=element_file, sep='\t', col.names=T, row.names=F, quote=F)

overlap$exon_ID <- paste0(overlap[,9], ":", overlap[,10], "-", overlap[,11], ":", overlap[,12], "(", overlap[,14], ")")

chimeric$exon_ID <- paste0(chimeric[,15], ":", chimeric[,16], "-", chimeric[,17], ":", chimeric[,18], "(", chimeric[,20], ")")
chimeric_exon_ID <- unique(chimeric$exon_ID)

intergenic_exon_ID <- unique(overlap[overlap[,12] %in% intergenic[,4],]$exon_ID)
intronic_exon_ID <- unique(overlap[overlap[,12] %in% intronic[,4],]$exon_ID)
antisense_exon_ID <- unique(overlap[overlap[,12] %in% antisense[,4],]$exon_ID)


TE_exon <- unique(overlap[,c(9:14,17)])
names(TE_exon) <- c("chrom", "start", "end", "tid", "gid", "strand", "exon_ID")

TE_exon$type = "exonic"
TE_exon[TE_exon$exon_ID %in% chimeric_exon_ID,]$type = "chimeric"
TE_exon[TE_exon$exon_ID %in% intronic_exon_ID,]$type = "intronic"
TE_exon[TE_exon$exon_ID %in% intergenic_exon_ID,]$type = "intergenic"
TE_exon[TE_exon$exon_ID %in% antisense_exon_ID,]$type = "antisense"

exon_cluster <- read.table(paste0(sample, "_exon_cluster.txt"), header=T, stringsAsFactors=F)
TE_exon <- merge(TE_exon, exon_cluster, by="exon_ID")

TE_exon$SM <- "ME"
if(file.info(paste0(sample, ".TE.SE.list"))$size != 0){
	SE <- read.table(paste0(sample, ".TE.SE.list"), stringsAsFactors=F)
	TE_exon[TE_exon$tid %in% SE[,1],]$SM <- "SE"
}

exon_file = paste0(sample, ".TE.exon.anno.txt")
write.table(TE_exon, file=exon_file, sep='\t', col.names=T, row.names=F, quote=F)

TE_tran = read.table(paste0(sample, ".TE.tran.sorted.bed"), stringsAsFactors=F)
names(TE_tran) <- c("chrom", "start", "end", "tid", "gid", "strand")

TE_tran$type = "exonic"
TE_tran[TE_tran[,4] %in% unique(chimeric[,18]),]$type = "chimeric"
TE_tran[TE_tran[,4] %in% intronic[,4],]$type = "intronic"
TE_tran[TE_tran[,4] %in% intergenic[,4],]$type = "intergenic"
TE_tran[TE_tran[,4] %in% antisense[,4],]$type = "antisense"

tran_cluster <- read.table(paste0(sample, "_tran_cluster.txt"), header=T, stringsAsFactors=F)
TE_tran <- merge(TE_tran, tran_cluster, by.x="tid", by.y="tran_ID")

TE_tran$SM <- "ME"
if(file.info(paste0(sample, ".TE.SE.list"))$size != 0){
	TE_tran[TE_tran$tid %in% SE[,1],]$SM <- "SE"
}

tran_file = paste0(sample, ".TE.tran.anno.txt")
write.table(TE_tran, file=tran_file, sep='\t', col.names=T, row.names=F, quote=F)
