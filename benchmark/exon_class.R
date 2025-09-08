args<-commandArgs(T)
qry_id=args[1]
diff=as.numeric(args[2])

qry_self_file <- paste0(qry_id, "_TE_exon_self_overlap.txt")
qry_self <- read.table(qry_self_file, stringsAsFactors=F, sep='\t')

qry_self$qry_id <- paste0(qry_self[,1], ":", qry_self[,2], "-", qry_self[,3], ":", qry_self[,4], "(", qry_self[,6], ")")
qry_self$qry_id2 <- paste0(qry_self[,8], ":", qry_self[,9], "-", qry_self[,10], ":", qry_self[,11], "(", qry_self[,13], ")")


exon_keep <- function(exon_overlap){
        keep_middle <- exon_overlap[,7]=="middle" & abs(exon_overlap[,2]-exon_overlap[,9])<=diff & abs(exon_overlap[,3]-exon_overlap[,10])<=diff
        keep_start_rev <- exon_overlap[,7]=="start" & abs(exon_overlap[,2]-exon_overlap[,9])<=diff & exon_overlap[,6]=="-"
        keep_start_for <- exon_overlap[,7]=="start" & abs(exon_overlap[,3]-exon_overlap[,10])<=diff & exon_overlap[,6]=="+"
        keep_end_rev <- exon_overlap[,7]=="end" & abs(exon_overlap[,3]-exon_overlap[,10])<=diff & exon_overlap[,6]=="-"
        keep_end_for <- exon_overlap[,7]=="end" & abs(exon_overlap[,2]-exon_overlap[,9])<=diff & exon_overlap[,6]=="+"
        keep_SE <- exon_overlap[,7]=="SE"
        keep <- keep_middle | keep_start_rev | keep_start_for | keep_end_rev | keep_end_for | keep_SE
        return(keep)
}

filt <- qry_self[qry_self$qry_id != qry_self$qry_id2,]
filt <- filt[exon_keep(filt),]

filt$min <- apply(filt[,16:17],1,function(x){sort(x)[1]})
filt$max <- apply(filt[,16:17],1,function(x){sort(x)[2]})
test <- unique(filt[,18:19])

exon_ID <- c()
exon_class <- c()
n <- 0
for(i in 1:nrow(test)){
	if(!(test[i,1] %in% exon_ID) & !(test[i,2] %in% exon_ID)){exon_ID <- c(exon_ID,test[i,1], test[i,2]); n <- n+1; exon_class <- c(exon_class, n, n)}
	if((test[i,1] %in% exon_ID) & !(test[i,2] %in% exon_ID)){exon_ID <- c(exon_ID, test[i,2]); m <- exon_class[which(exon_ID==test[i,1])]; exon_class <- c(exon_class, m)}
	if(!(test[i,1] %in% exon_ID) & (test[i,2] %in% exon_ID)){exon_ID <- c(exon_ID, test[i,1]); m <- exon_class[which(exon_ID==test[i,2])]; exon_class <- c(exon_class, m)}
}

data1 <- data.frame(exon_ID, exon_class)

all_ID <- unique(qry_self$qry_id)
remain_ID <- all_ID[!(all_ID %in% data1$exon_ID)]
data2 <- data.frame(exon_ID=remain_ID, exon_class=(n+1):(n+length(remain_ID)))

data <- rbind(data1, data2)

write.table(data, file=paste0(qry_id, "_exon_cluster.txt"), sep='\t', col.names=T, row.names=F, quote=F)
