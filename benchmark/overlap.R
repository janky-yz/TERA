args<-commandArgs(T)
sample=args[1]

overlap_file = paste0(sample, ".TE.overlap.txt")
intergenic_file = paste0(sample, ".TE.intergenic.bed")
intronic_file = paste0(sample, ".TE.intronic.bed")
antisense_file = paste0(sample, ".TE.antisense.bed")

overlap = read.table(overlap_file, stringsAsFactors=F)
intergenic = read.table(intergenic_file, stringsAsFactors=F)
intronic = read.table(intronic_file, stringsAsFactors=F)
antisense = read.table(antisense_file, stringsAsFactors=F)

overlap$start <- as.numeric(apply(overlap,1,function(x){max(c(x[2],x[10]))}))
overlap$end <- as.numeric(apply(overlap,1,function(x){min(c(x[3],x[11]))}))
out <- overlap[!(overlap[,12] %in% c(intergenic[,4], intronic[,4], antisense[,4])),c(1,16,17,4,12,6,1:14)]

out_file = paste0(sample, ".TE.overlap.pos.bed")
write.table(out, file=out_file, sep='\t', col.names=F, row.names=F, quote=F)
