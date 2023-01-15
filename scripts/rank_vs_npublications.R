setwd("C:/Users/Tania_2/Internship/R/Results/ranked.lists")
## NOTE!: - We could have done a simplier code (rank_vs_nºpublications_memmoryproblem.R), but due to a memmory problem we had to use the present script instead
##        - Here Matthias run on his computer a script to get all pubmed ids with term 'stem cell' either only in Title/Abstract or on all article
##          Then those pubmed ids can be crossed with pubmed ids for publications with a certain gene

# the following part was run on Matthias computer because of memory error in mine
# so i just opened the file rankpublicationhuman.RData
# StemCellTA is the vector with pubmed ids for 'stem cell' in Title/Abstract:
# library(XML)
# library(rentrez)
# stemhist <- entrez_search(db="pubmed", term="\"stem cell\"[Title/Abstract]", use_history=TRUE)
# maxRecs <- stemhist$count
# tmp <- ""
# seqRecs <- c(0, seq(10000,maxRecs,10000),maxRecs)
# for (i in 1:(length(seqRecs)-1)){
#         
#         recs <- entrez_fetch(db="pubmed", web_history=stemhist$web_history,
#                              rettype="xml", retmax=seqRecs[i+1]-seqRecs[i], retstart=seqRecs[i],parsed=TRUE)
#         
#         cat(seqRecs[i+1],"references downloaded\r")
#         recsList <- xmlToList(recs)
#         
#         for (i in 1:length(recsList)) tmp<- c(tmp,recsList[[i]]$Medline$PMID$text)
# }
# tmp2 <- tmp[-1]
# tmp2 <- unique(tmp2)
# StemCellTA <- tmp2

# just opened the file StemCellPubmed2.RData
# StemCellWR is the vector with pubmed ids for 'stem cell' in all article:
# library(XML)
# library(rentrez)
# stemhist <- entrez_search(db="pubmed", term="\"stem cell\"", use_history=TRUE)
# maxRecs <- stemhist$count
# tmp <- ""
# seqRecs <- c(0, seq(10000,maxRecs,10000),maxRecs)
# for (i in 1:(length(seqRecs)-1)){
#         
#         recs <- entrez_fetch(db="pubmed", web_history=stemhist$web_history,
#                              rettype="xml", retmax=seqRecs[i+1]-seqRecs[i], retstart=seqRecs[i],parsed=TRUE)
#         
#         cat(seqRecs[i+1],"references downloaded\r")
#         recsList <- xmlToList(recs)
#         
#         for (i in 1:length(recsList)) tmp<- c(tmp,recsList[[i]]$Medline$PMID$text)
# }
# tmp2 <- tmp[-1]
# tmp2 <- unique(tmp2)
# StemCellWR <- tmp2

# ## HUMAN
# # retrieve unyfying signature and overlap with articles previously selected
# library(rentrez)
# library(httr)
# set_config(use_proxy(url="10.3.100.207",port=8080))
# rank_genes <- read.csv("C:/Users/Tania_2/Internship/R/Results/ranked.lists/ranked.list.human.stem.signature.csv", header=T, stringsAsFactors = F)
# genes <- rank_genes[1:164,1] # use gene IDs to count for alias
# HitsTA <- list()
# Hitsall <- list()
# NumberHitsTA <- character(length(genes))
# NumberHitsall <- character(length(genes))
# for (i in 1:length(genes)){
#         pubmedlistid <- entrez_link(dbfrom='gene', id=genes[i], db='pubmed')
#         a <- pubmedlistid[['links']]$gene_pubmed
#         b <- pubmedlistid[['links']]$gene_pubmed_citedinomim
#         c <- pubmedlistid[['links']]$gene_pubmed_pmc_nucleotide
#         d <- pubmedlistid[['links']]$gene_pubmed_rif
#         pubmedids <- unique(c(a, b, c, d))
#         AbstractsTA <- pubmedids[pubmedids%in%StemCellTA] # pubmedIDs for a gene with 'stem cell' in title/abstract
#         NumberTA <- sum(pubmedids%in%StemCellTA) # number of articles for a gene with 'stem cell' in title/abstract
#         Abstractsall <- pubmedids[pubmedids%in%StemCellWR] # pubmedIDs for a gene with 'stem cell' in all article
#         NumberWR <- sum(pubmedids%in%StemCellWR) # number of articles for a gene with 'stem cell' in all article
#         HitsTA[[i]] <- AbstractsTA
#         NumberHitsTA[i] <- NumberTA
#         Hitsall[[i]] <- Abstractsall
#         NumberHitsall[i] <- NumberWR
# }
load('C:/Users/Tania_2/Internship/R/article/scripts_used_in_article/rankpublicationhuman.RData')
# add gene names
gene_name <- rank_genes[1:164,2]
names(HitsTA) <- gene_name
names(Hitsall) <- gene_name
names(NumberHitsTA) <- gene_name
names(NumberHitsall) <- gene_name
# add rank of unifying signature:
rank_uni <- rank_genes[1:164,3]

# table with info all together:
togetherTA <- data.frame(genes=gene_name, rank=rank_uni, publication_number=as.numeric(NumberHitsTA), stringsAsFactors = F)
togetherall <- data.frame(genes=gene_name, rank=rank_uni, publication_number=as.numeric(NumberHitsall), stringsAsFactors = F)

# graphs Gene Stemness Score vs Gene number of publications
# 1st graph with all genes of integrative stemness signature and log scale in number of publications
 # NOTE! Genes without publications (publication_number=0) are not in this plot as it is in log scale!

svg(file="C:/Users/Tania_2/Internship/R/Results/rank.vs.publication/human/human_integrativesig_SteminTitleorAbstract.svg")#export file to png
set.seed(1)
plot(jitter(togetherTA[,"publication_number"]),togetherTA[,"rank"],main="Rank vs Number of Publications", xlab= "Number of Publications", ylab= "Gene Score", type="p", col= "blue", bg= "blue", pch = 21, lwd = 0.5, ylim=c(3.5,12.5),yaxt = "n", log="x", xaxt = "n", xlim = c(1,1000), cex.main=2, cex.lab=2)
ticks = seq(0,3,by=1)
labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
axis(1, at=c(1, 10, 100, 1000), labels=labels, cex.axis = 1.4)
axis(side = 2, at = seq(from = 4, to = 12, by = 1), cex.axis=1.4)
text(jitter(togetherTA[,"publication_number"]),togetherTA[,"rank"], labels=gene_name, cex= 1.3, pos=3)
dev.off()

svg(file="C:/Users/Tania_2/Internship/R/Results/rank.vs.publication/human/human_integrativesig_SteminAllArticle.svg")#export file to svg
set.seed(1)
plot(jitter(togetherall[,"publication_number"]),togetherall[,"rank"],main="Rank vs Number of Publications", xlab= "Number of Publications", ylab= "Gene Score", type="p", col= "blue", bg= "blue", pch = 21, lwd = 0.5, ylim=c(3.5,12.5),yaxt = "n", log="x", xaxt = "n", xlim = c(1,1000), cex.main=2, cex.lab=2)
ticks = seq(0,3,by=1)
labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
axis(1, at=c(1, 10, 100, 1000), labels=labels, cex.axis = 1.4)
axis(side = 2, at = seq(from = 4, to = 12, by = 1), cex.axis=1.4)
text(jitter(togetherall[,"publication_number"]),togetherall[,"rank"], labels=gene_name, cex= 1.3, pos=3)
dev.off()

# selects those with publications < 4
togetherTA_sub <- subset(togetherTA, togetherTA[,"publication_number"] < 4)
togetherall_sub <- subset(togetherall, togetherall[,"publication_number"] < 4)

svg(file="C:/Users/Tania_2/Internship/R/Results/rank.vs.publication/human/human_integrativesig_SteminTitleorAbstract_until3publications_withoutlog.svg")#export file to png
plot(togetherTA_sub[,"publication_number"],togetherTA_sub[,"rank"],main="Rank vs Number of Publications", xlab= "Number of Publications", ylab= "Gene Score", type="p", col= "blue",bg= "blue", pch = 21, lwd = 0.5, xlim=c(-0.5,3.5), ylim=c(3.5,7.5), xaxt = "n", yaxt = "n", cex.main=2, cex.lab=2)
axis(side = 1, at = seq(from = 0, to = 5, by = 1),cex.axis = 1.4)
axis(side = 2, at = seq(from = 4, to = 7, by = 1), cex.axis = 1.4)
text(togetherTA_sub[,"publication_number"], togetherTA_sub[,"rank"], labels=togetherTA_sub[,"genes"], cex=1.1, pos=3)
dev.off()

svg(file="C:/Users/Tania_2/Internship/R/Results/rank.vs.publication/human/human_integrativesig_SteminAllArticle_until3publications_withoutlog.svg")#export file to png
plot(togetherall_sub[,"publication_number"],togetherall_sub[,"rank"],main="Rank vs Number of Publications", xlab= "Number of Publications", ylab= "Gene Score", type="p", col= "blue",bg= "blue", pch = 21, lwd = 0.5, xlim=c(-0.5,3.5), ylim=c(3.5,7.5), xaxt = "n", yaxt = "n", cex.main=2, cex.lab=2)
axis(side = 1, at = seq(from = 0, to = 5, by = 1),cex.axis = 1.4)
axis(side = 2, at = seq(from = 4, to = 7, by = 1),cex.axis = 1.4)
text(togetherall_sub[,"publication_number"], togetherall_sub[,"rank"], labels=togetherall_sub[,"genes"], cex=1.1, pos=3)
dev.off()

# make table with those never published :
pub0TA <- subset(togetherTA_sub, togetherTA_sub[,"publication_number"] == 0)
write.table(pub0TA, "C:/Users/Tania_2/Internship/R/Results/rank.vs.publication/human/human_integrativesig_SteminTitleorAbstract_neverpublished.txt", sep = "\t", row.names = F, quote = F)

pub0all <- subset(togetherall_sub, togetherall_sub[,"publication_number"] == 0)
write.table(pub0all, "C:/Users/Tania_2/Internship/R/Results/rank.vs.publication/human/human_integrativesig_SteminAllArticle_neverpublished.txt", sep = "\t", row.names = F, quote = F)

### MOUSE
## UNIFYING SIGNATURE
# retrieve unyfying signature and overlap with articles previously selected
# I just run the environment file rankpublicationmouse.RData
# library(rentrez)
# rank_genes <- read.csv("C:/Users/Tania_2/Internship/R/Results/ranked.lists/ranked.list.mouse.stem.signature.csv", header=T, stringsAsFactors = F)
# genes <- rank_genes[1:115,1]
# HitsTA <- list()
# Hitsall <- list()
# NumberHitsTA <- character(length(genes))
# NumberHitsall <- character(length(genes))
# for (i in 1:length(genes)){
#         pubmedlistid <- entrez_link(dbfrom='gene', id=genes[i], db='pubmed')
#         a <- pubmedlistid[['links']]$gene_pubmed
#         b <- pubmedlistid[['links']]$gene_pubmed_citedinomim
#         c <- pubmedlistid[['links']]$gene_pubmed_pmc_nucleotide
#         d <- pubmedlistid[['links']]$gene_pubmed_rif
#         pubmedids <- unique(c(a, b, c, d))
#         AbstractsTA <- pubmedids[pubmedids%in%StemCellTA] # pubmedIDs for a gene with 'stem cell' in title/abstract
#         NumberTA <- sum(pubmedids%in%StemCellTA) # number of articles for a gene with 'stem cell' in title/abstract
#         Abstractsall <- pubmedids[pubmedids%in%StemCellWR] # pubmedIDs for a gene with 'stem cell' in all article
#         NumberWR <- sum(pubmedids%in%StemCellWR) # number of articles for a gene with 'stem cell' in all article
#         HitsTA[[i]] <- AbstractsTA
#         NumberHitsTA[i] <- NumberTA
#         Hitsall[[i]] <- Abstractsall
#         NumberHitsall[i] <- NumberWR
# }
load('C:/Users/Tania_2/Internship/R/article/scripts_used_in_article/rankpublicationmouse.RData')
gene_name <- rank_genes[1:115,2]
names(HitsTA) <- gene_name
names(Hitsall) <- gene_name
names(NumberHitsTA) <- gene_name
names(NumberHitsall) <- gene_name
# add rank of unifying signature:
rank_uni <- rank_genes[1:115,3]

# table with info all together:
togetherTA <- data.frame(genes=gene_name, rank=rank_uni, publication_number=as.numeric(NumberHitsTA), stringsAsFactors = F)
togetherall <- data.frame(genes=gene_name, rank=rank_uni, publication_number=as.numeric(NumberHitsall), stringsAsFactors = F)

# graphs Gene Stemness Score vs Gene number of publications
# 1st graph with all genes of integrative stemness signature and log scale in number of publications
# NOTE! Genes without publications (publication_number=0) are not in this plot as it is in log scale!

svg(file="C:/Users/Tania_2/Internship/R/Results/rank.vs.publication/mouse/mouse_integrativesig_SteminTitleorAbstract.svg")
set.seed(1)
plot(jitter(togetherTA[,"publication_number"]),togetherTA[,"rank"],main="Rank vs Number of Publications", xlab= "Number of Publications", ylab= "Gene Score", type="p", col= "blue", bg= "blue", pch = 21, lwd = 0.5, ylim=c(6.5,13.5),yaxt = "n", log="x", xaxt = "n", xlim = c(1,1000), cex.main=2, cex.lab=2)
ticks = seq(0,3,by=1)
labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
axis(1, at=c(1, 10, 100, 1000), labels=labels, cex.axis = 1.4)
axis(side = 2, at = seq(from = 7, to = 13, by = 1), cex.axis = 1.4)
text(jitter(togetherTA[,"publication_number"]),togetherTA[,"rank"], labels=gene_name, cex= 1.3, pos=3)
dev.off()

svg(file="C:/Users/Tania_2/Internship/R/Results/rank.vs.publication/mouse/mouse_integrativesig_SteminAllArticle.svg")#export file to png
set.seed(1)
plot(jitter(togetherall[,"publication_number"]),togetherall[,"rank"],main="Rank vs Number of Publications", xlab= "Number of Publications", ylab= "Gene Score", type="p", col= "blue", bg= "blue", pch = 21, lwd = 0.5, ylim=c(6.5,13.5),yaxt = "n", log="x", xaxt = "n", xlim = c(1,1000), cex.main=2, cex.lab=2)
ticks = seq(0,3,by=1)
labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
axis(1, at=c(1, 10, 100, 1000), labels=labels, cex.axis = 1.4)
axis(side = 2, at = seq(from = 7, to = 13, by = 1), cex.axis = 1.4)
text(jitter(togetherall[,"publication_number"]),togetherall[,"rank"], labels=gene_name, cex= 1.3, pos=3)
dev.off()

# selects those with publications < 4
togetherTA_sub <- subset(togetherTA, togetherTA[,"publication_number"] < 4)
togetherall_sub <- subset(togetherall, togetherall[,"publication_number"] < 4)

svg(file="C:/Users/Tania_2/Internship/R/Results/rank.vs.publication/mouse/mouse_integrativesig_SteminTitleorAbstract_until3publications_withoutlog.svg")#export file to png
plot(togetherTA_sub[,"publication_number"],togetherTA_sub[,"rank"],main="Rank vs Number of Publications", xlab= "Number of Publications", ylab= "Gene Score", type="p", col= "blue",bg= "blue", pch = 21, lwd = 0.5, xlim=c(-0.5,3.5), ylim=c(6.5,10.5), xaxt = "n", yaxt = "n", cex.main=2, cex.lab=2)
axis(side = 1, at = seq(from = 0, to = 5, by = 1), cex.axis = 1.4)
axis(side = 2, at = seq(from = 7, to = 10, by = 1), cex.axis = 1.4)
text(togetherTA_sub[,"publication_number"], togetherTA_sub[,"rank"], labels=togetherTA_sub[,"genes"], cex=1.3, pos=3)
dev.off()

svg(file="C:/Users/Tania_2/Internship/R/Results/rank.vs.publication/mouse/mouse_integrativesig_SteminAllArticle_until3publications_withoutlog.svg")#export file to png
plot(togetherall_sub[,"publication_number"],togetherall_sub[,"rank"],main="Rank vs Number of Publications", xlab= "Number of Publications", ylab= "Gene Score", type="p", col= "blue",bg= "blue", pch = 21, cex = 0.7, lwd = 0.5, xlim=c(-0.5,3.5), ylim=c(6.5,9.5), xaxt = "n", yaxt = "n",cex.main=2, cex.lab=2)
axis(side = 1, at = seq(from = 0, to = 5, by = 1), cex.axis = 1.4)
axis(side = 2, at = seq(from = 7, to = 9, by = 1), cex.axis = 1.4)
text(togetherall_sub[,"publication_number"], togetherall_sub[,"rank"], labels=togetherall_sub[,"genes"], cex=1.1, pos=3)
dev.off()

# make table with those never published :
pub0TA <- subset(togetherTA_sub, togetherTA_sub[,"publication_number"] == 0)
write.table(pub0TA, "C:/Users/Tania_2/Internship/R/Results/rank.vs.publication/mouse/mouse_integrativesig_SteminTitleorAbstract_neverpublished.txt", sep = "\t", row.names = F, quote = F)

pub0all <- subset(togetherall_sub, togetherall_sub[,"publication_number"] == 0)
write.table(pub0all, "C:/Users/Tania_2/Internship/R/Results/rank.vs.publication/mouse/mouse_integrativesig_SteminAllArticle_neverpublished.txt", sep = "\t", row.names = F, quote = F)
