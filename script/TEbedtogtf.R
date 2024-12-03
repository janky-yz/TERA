args<-commandArgs(T)
prefix=args[1]

bed_file = paste0(prefix, ".bed")
a=read.table(bed_file, stringsAsFactors=F)
a$TE_ID <- paste0(a[,4], "_", a[,1], "_", a[,2]+1, "_", a[,3])
a$TE_name <- paste0(a[,5], "_", a[,1], "_", a[,2]+1, "_", a[,3])
a$des <- paste0('transcript_id "', a$TE_ID, '"; gene_id "', a$TE_name, '"; family_id "', a[,7], '"; class_id "', a[,8], '";')
a$attr <- "TEtranscripts"
a$entry <- "exon"
a$score <- "."
out <- data.frame(a[,1], a$attr, a$entry, a[,2]+1, a[,3], a$score, a[,6], a$score, a$des)

gtf_file = paste0(prefix, ".gtf")
write.table(out, file=gtf_file, sep='\t', col.names=F, row.names=F, quote=F)
