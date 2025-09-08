args<-commandArgs(T)
qry_id=args[1]
diff=as.numeric(args[2])

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

self <- read.table(paste0(qry_id, "_TE_transcript_self_overlap.txt"), stringsAsFactors=F)
filt <- self[self[,4]!=self[,12],]
filt <- filt[,-8]
filt <- filt[exon_keep(filt),]
filt$ID <- paste0(filt[,4], ":", filt[,11], ":", filt[,15])

count1 <- as.data.frame(table(filt$ID))
count2 <- unique(filt[,c(4,11,15,17)])
names(count2) <- c("ref_id", "qry_id", "exon_num", "ID")
names(count1) <- c("ID", "Freq")
count <- merge(count2, count1, by="ID")
count <- count[count$exon_num==count$Freq,]



count$min <- apply(count[,2:3],1,function(x){sort(x)[1]})
count$max <- apply(count[,2:3],1,function(x){sort(x)[2]})
test <- unique(count[,c("min", "max")])

tran_ID <- c()
tran_class <- c()
n <- 0
for(i in 1:nrow(test)){
	if(!(test[i,1] %in% tran_ID) & !(test[i,2] %in% tran_ID)){tran_ID <- c(tran_ID, test[i,1], test[i,2]); n <- n+1; tran_class <- c(tran_class, n, n)}
	if((test[i,1] %in% tran_ID) & !(test[i,2] %in% tran_ID)){tran_ID <- c(tran_ID, test[i,2]); m <- tran_class[which(tran_ID==test[i,1])]; tran_class <- c(tran_class, m)}
	if(!(test[i,1] %in% tran_ID) & (test[i,2] %in% tran_ID)){tran_ID <- c(tran_ID, test[i,1]); m <- tran_class[which(tran_ID==test[i,2])]; tran_class <- c(tran_class, m)}
}

data1 <- data.frame(tran_ID, tran_class)
tid <- unique(self[,4])
tid <- tid[!(tid %in% data1$tran_ID)]
data2 <- data.frame(tran_ID=tid, tran_class=(n+1):(n+length(tid)))
data <- rbind(data1, data2)

write.table(data, file=paste0(qry_id, "_tran_cluster.txt"), sep='\t', col.names=T, row.names=F, quote=F)
