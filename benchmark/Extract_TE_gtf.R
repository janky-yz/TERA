args<-commandArgs(T)
sample=args[1]

gtf_file <- paste0(sample, ".sorted.gtf")
gtf <- read.table(gtf_file, stringsAsFactors=F, fill=T)

ID_file <- paste0(sample, "_TET_ID.list")
ID <- read.table(ID_file, stringsAsFactors=F)

TE_gtf <- gtf[gtf[,10] %in% ID[,1],]
TE_gtf[,9] <- paste0('transcript_id "', TE_gtf[,10], '"; gene_id "', TE_gtf[,13], '";')

out_file <- paste0(sample, ".TE.gtf")
write.table(TE_gtf[,1:9], file=out_file, sep='\t', col.names=F, row.names=F, quote=F)
