args <- commandArgs(T);
if (length(args) != 5 ){
    cat("\033[0;32m
    USAGE: 
        category.txt go.background.txt kegg.backgroud.txt gene.list outprefix
    \033[0m\n"
       );
    q();
}


category <- read.delim(args[1], header=F, sep="\t", row.names=1)
go.bg <- read.delim(args[2], header=F, sep="\t")
kegg.bg <- read.delim(args[3], header=F, sep="\t")

map <- function(diff.bg) {
	id <- strsplit(as.character(diff.bg[2]), ",")[[1]]
	term <- strsplit(as.character(diff.bg[3]), "|", fix=T)[[1]]
	return(as.data.frame(cbind(diff.bg[1], id, term)))
}

enrichment <- function (diff, bg) {	
	bg.total <- unlist(lapply(as.vector(bg[[2]]), function(x) strsplit(x, ",") ))
	bg.total <- as.data.frame(table(bg.total))
	diff.bg <- bg[!is.na(match(bg[,1], diff[,1])),]
	if(nrow(diff.bg)==0 ){ return(NULL) }

	d <- apply(diff.bg, 1, map)
	d <- Reduce(rbind, d)
	enrich <- as.data.frame( table(d[,2]) )
	colnames(enrich)[c(1,2)] <- c("id", "ListHits")

	for(i in 1:nrow(enrich)) {
		enrich[i, "term"] <- d[which(d[,2]==enrich[i,1]), 3][1]
		enrich[i, "Gene"] <- paste(d[which(d[,2]==enrich[i,1]), 1], collapse="; ")
		enrich[i, "PopHits"] <- bg.total[which(bg.total[,1]==as.vector(enrich[i,1])), 2]
	}
	enrich["ListTotal"] <- nrow(diff.bg)
	enrich["PopTotal"] <- nrow(bg)
	enrich["pval"] <- phyper(enrich[,"ListHits"], enrich[,"PopHits"], 
		enrich[,"PopTotal"]-enrich[,"PopHits"], enrich[,"ListTotal"], lower.tail=F)

	enrich["padj"] <- p.adjust(enrich[,"pval"], method="fdr")
	enrich["Enrichment_score"] <- 
		(enrich["ListHits"]*enrich["PopTotal"])/(enrich["ListTotal"]*enrich["PopHits"])
	enrich <- enrich[order(enrich["pval"]), ]
	return(enrich[c(1,3,2,6,5,7,8,9,10,4)])
}

#for(i in 4:length(args)) {
	f <- args[4]
	diff <- read.delim(f, header=F, sep="\t")
	d <- enrichment(diff, go.bg)
	d["category"] <- category[as.vector(d[,1]), ]
	d <- d[, c(1, 2, ncol(d), 3:(ncol(d)-1))]

	if(!is.null(d)) { 
		write.table(d, paste0(args[5],"/enrichment-go.txt"), 
		sep="\t", row.names=F, quote=F) }

	d <- enrichment(diff, kegg.bg)
	if(!is.null(d)) { 
		write.table(d, paste0(args[5],"/enrichment-kegg.txt"), 
		sep="\t", row.names=F, quote=F) }
#}
