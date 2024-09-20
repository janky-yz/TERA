args<-commandArgs(T)
sample=args[1]

bed_file <- paste0(sample, "_exon.bed")
bed <- read.table(bed_file, sep='\t', stringsAsFactors=F)
bed$transcript_id <- bed[,4]

exon <- aggregate(bed[,2], by=list(bed[,4]), min)
names(exon) <- c("transcript_id", "min")
exon$max <- aggregate(bed[,3], by=list(bed[,4]), max)[,2]
bed_anno <- merge(bed, exon, by="transcript_id", sort=F)
bed_anno <- bed_anno[,-1]
bed_anno$type <- "middle"
bed_anno[bed_anno[,2]==bed_anno$min & bed_anno[,6]=="+",]$type <- "start"
bed_anno[bed_anno[,3]==bed_anno$max & bed_anno[,6]=="-",]$type <- "start"
bed_anno[bed_anno[,3]==bed_anno$max & bed_anno[,6]=="+",]$type <- "end"
bed_anno[bed_anno[,2]==bed_anno$min & bed_anno[,6]=="-",]$type <- "end"
bed_anno[bed_anno[,2]==bed_anno$min & bed_anno[,3]==bed_anno$max,]$type <- "SE"
bed_anno <- bed_anno[order(bed_anno[,1], bed_anno[,2], bed_anno[,3]),]

write.table(bed_anno[,c(1:6,9)], file=bed_file, sep='\t', col.names=F, row.names=F, quote=F)
