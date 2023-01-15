###### TABLES
### Human
## REACTOME IDs and corresponding genes
library(reactome.db)
library(ReactomePA)
library(org.Hs.eg.db)

# get reactome ids of reactome only:
xx.reactome <- as.list(reactomeEXTID2PATHID) # a list in which each element name is a reactomeID and each element is a vector of all gene entrez IDs corresponding to that reactomeID
reactomegene <- unique(names(xx.reactome)) # all reactome genes
reactomeid <- unique(unlist(xx.reactome))# all reactome terms
library(biomaRt)
ensembl = useMart(biomart ="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="may2017.archive.ensembl.org")
d <- getBM(attributes = "entrezgene",values = "*",mart = ensembl)
h <- as.character(unique(d[,1])) # gets all human genes
allgenes <- h[h%in%reactomegene]# gets all reactome genes that are human!

### top genes (the ones with score equal or higher than rank.limit) from stemness gene sets
setwd("C:/Users/Tania/Internship/R/Results/ranked.lists")
rank.genes <- unique(read.csv("ranked.list.human.stem.signature.csv", header=T, stringsAsFactors=F)) #create a list with the data from all files
rank.limit <- 4
top.rank <- as.character(subset(rank.genes$gene_id, rank.genes$score >= rank.limit))
# NOTE: if we want to work with the top 100 genes instead(not good idea cause this splits genes of the same ranking (rank 4) into the group of the selected ones and not selected ones):
#top.rank <- rank.genes[1:100]

treact <- enrichPathway(gene=top.rank, organism="human", pvalueCutoff=2, qvalueCutoff = 2, pAdjustMethod="fdr", universe=allgenes, readable=T, minGSSize=1)
treact2 <- treact@result
signifh <- subset(treact2, treact2$p.adjust < 0.05)
write.csv(signifh, "C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/Reactome/reactome.human.top.stemness.csv", row.names=F)

### Mouse
## REACTOME IDs and corresponding genes
library(biomaRt)
ensembl = useMart(biomart ="ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host="may2017.archive.ensembl.org")
d <- getBM(attributes = "entrezgene",values = "*",mart = ensembl)
m <- as.character(unique(d[,1])) # gets all mouse genes
allgenes <- m[m%in%reactomegene]# gets all reactome genes that are mouse!

### top genes (the ones with rank equal or higher than rank.limit) from stemness gene sets
setwd("C:/Users/Tania/Internship/R/Results/ranked.lists")
rank.genes <- unique(read.csv("ranked.list.mouse.stem.signature.csv", header=T, stringsAsFactors=F)) #create a list with the data from all files
### top genes (the ones with score equal or higher than rank.limit) from stemness gene sets
rank.limit <- 7
top.rank <- as.character(subset(rank.genes$gene_id, rank.genes$score >= rank.limit))
# NOTE: if we want to work with the top 100 genes instead(not good idea cause this splits genes of the same ranking (rank 4) into the group of the selected ones and not selected ones):
#top.rank <- rank.genes[1:100]

treact <- enrichPathway(gene=top.rank, organism="mouse", pvalueCutoff=2, qvalueCutoff = 2, pAdjustMethod="fdr", universe=allgenes, readable=T, minGSSize=1)
treact2 <- treact@result
signifm <- subset(treact2, treact2$p.adjust < 0.05)
write.csv(signifm, "C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/Reactome/reactome.mouse.top.stemness.csv", row.names=F)

### Bar plots (human and mouse together in same barplot)
setwd("C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/Reactome")
rank.human <- read.csv("Reactome.human.top.stemness.csv", header=T, stringsAsFactors=F)
rank.human <- rank.human[,c(2,6)]
rank.mouse <- read.csv("Reactome.mouse.top.stemness.csv", header=T, stringsAsFactors=F)
rank.mouse <- rank.mouse[,c(2,6)]
# put column FDR with the indication if it's from mouse or human
names(rank.mouse)[2] <- "Mouse"
names(rank.human)[2] <- "Human"
# Enriched in both mouse and human - put human and mouse in the same matrix:
# NOTE!!!: as there are some Reactome terms for mouse that don't exist in human and vice-versa
#          it doesn't make sense to input all terms of mouse or human. we could in the case a term was signifficant for mouse for example to put the p-value column of human as 1, but in that case the reader would interpret that that term wasn't enriched in human, while in fact that term didn't exist in human
#          SO WE LOOK FOR THE ENRICHED TERMS IN BOTH!
both <- merge(rank.mouse, rank.human, by.x=1, by.y=1, all=F)
# give 1st column of matrix as row.names:
row.names(both) <- both[,1]
# remove 1st column:
both <- both[,-1]
# transpose:
botha <- t(both)
# transform adjusted p-value into -log(p-value):
botha <- -log10(botha)
# do the graph:
# las - changes angle of names.arg rotation/ x -label location
# horiz = T - changes the orientation of the plot horizontaly
# names.arg - gives names to the bars
# cex.axis - x label size
# par(mar=c(3.5, 3.5, 2, 1)) - bottom, left, top, right margins of plot
# args.legend - allows to add aditional arguments for legend
# cex - gives the sixe of axis and corresponding numbers
# xlim - limits for the x axis
# png(file = "C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/Reactome/reactome.enrich.bothmousehuman.top.stemness.png", width = 80*300, height = 65*300, res = 300)#export file to png
# par(mar=c(20, max(nchar(colnames(botha)))+135, 10,5))
# barplot(botha, horiz = T, beside=T, las=1, space = c(6/ncol(botha),12/ncol(botha)), col=colorRampPalette(c("forestgreen","red"))(n=2), cex.names=7, legend.text=T, args.legend=list(cex = 3), axes = F, xlim = c(0,10)) # both axis have to be set false so that then we can overlaigh a axis with parameters we want
# # to use main title with its corresponding parameters (for ex. size) we have to set the title appart from barplot:
# # title(main = list("Reactome Pathways - both mouse and human", cex = 5))
# axis(side=1, at=seq(from=0,to = 10, by=2), lwd=4, tck = -0.01, labels=NA, line = -5) # this is for the ticks only in x axis. that's why lables=NA.side - an integer indicating the side of the graph to draw the axis (1=bottom, 2=left, 3=top, 4=right); at - a numeric vector indicating where tick marks should be drawn; lwd - line width;
# axis(side=1, at=seq(from=0,to = 10, by=2), lwd=0, cex=5, cex.axis=7, line = 3) # this is for the numbering only in x axis. that's why lwd=0.
# title(xlab=expression(-log[10](`adj.p-value`)), cex.lab = 9, line = 17) # adds x axis label
# dev.off()

png(file = "C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/Reactome/reactome.enrich.bothmousehuman.top.stemness.png", width = 35*600, height = 30*600, res = 600)#export file to png
par(mar=c(10, max(nchar(colnames(botha)))+70, 15,10))
#barplot(botha, horiz = T, beside=T, las=1, space = c(15/ncol(botha),45/ncol(botha)), col=colorRampPalette(c("forestgreen","red"))(n=2), cex.names=4.9, legend.text=T, args.legend=list(x=0, y=0, cex = 3), axes = F, xlim = c(0,20)) # both axis have to be set false so that then we can overlaigh a axis with parameters we want
barplot(botha, horiz = T, beside=T, las=1, space = c(15/ncol(botha),45/ncol(botha)), col=colorRampPalette(c("forestgreen","red"))(n=2), cex.names=4.9, axes = F, xlim = c(0,10)) # both axis have to be set false so that then we can overlaigh a axis with parameters we want
# to use main title with its corresponding parameters (for ex. size) we have to set the title appart from barplot:
# title(main = list("Reactome Pathways - both mouse and human", cex = 5))
#par(cex.axis=1, cex.lab=3, cex.main=500)
axis(side=1, at=seq(from=0,to = 10, by=5), lwd=5, tck = -0.02, labels=NA, line = -2) # this is for the ticks only in x axis. that's why lables=NA.side - an integer indicating the side of the graph to draw the axis (1=bottom, 2=left, 3=top, 4=right); at - a numeric vector indicating where tick marks should be drawn; lwd - line width;
axis(side=1, at=seq(from=0,to = 10, by=5), lwd=0, cex.axis=5, line = 0.5) # this is for the numbering only in x axis. that's why lwd=0.
title(xlab=expression(-log[10](`adj.p-value`)), cex.lab =5.5, line = 9) # adds x axis label
dev.off()

tiff(file = "C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/Reactome/reactome.enrich.bothmousehuman.top.stemness.tif", width = 35*600, height = 30*600, res = 600)#export file to png
par(mar=c(10, max(nchar(colnames(botha)))+70, 15,10))
#barplot(botha, horiz = T, beside=T, las=1, space = c(15/ncol(botha),45/ncol(botha)), col=colorRampPalette(c("forestgreen","red"))(n=2), cex.names=4.9, legend.text=T, args.legend=list(x=0, y=0, cex = 3), axes = F, xlim = c(0,20)) # both axis have to be set false so that then we can overlaigh a axis with parameters we want
barplot(botha, horiz = T, beside=T, las=1, space = c(15/ncol(botha),45/ncol(botha)), col=colorRampPalette(c("forestgreen","red"))(n=2), cex.names=4.9, axes = F, xlim = c(0,10)) # both axis have to be set false so that then we can overlaigh a axis with parameters we want
# to use main title with its corresponding parameters (for ex. size) we have to set the title appart from barplot:
# title(main = list("Reactome Pathways - both mouse and human", cex = 5))
#par(cex.axis=1, cex.lab=3, cex.main=500)
axis(side=1, at=seq(from=0,to = 10, by=5), lwd=5, tck = -0.02, labels=NA, line = -2) # this is for the ticks only in x axis. that's why lables=NA.side - an integer indicating the side of the graph to draw the axis (1=bottom, 2=left, 3=top, 4=right); at - a numeric vector indicating where tick marks should be drawn; lwd - line width;
axis(side=1, at=seq(from=0,to = 10, by=5), lwd=0, cex.axis=5, line = 0.5) # this is for the numbering only in x axis. that's why lwd=0.
title(xlab=expression(-log[10](`adj.p-value`)), cex.lab =5.5, line = 9) # adds x axis label
dev.off()


# Enriched only in human and not in mouse:
humanonly <- subset(rank.human, !(rank.human$Description%in%rank.mouse$Description))
 # remove 1st column:
humanonlya <- humanonly[,-1]
# transform matrix into vector:
humanonlya <- as.vector(humanonlya)
# give terms as names of the vector:
names(humanonlya) <- humanonly[,1]
# reverse order of vector (so that later on the terms most enriched are shown first in barplot):
humanonlya <- rev(humanonlya)
# change rows by columns:
humanonlya <- t(humanonlya)
 # transform adjusted p-value into -log(p-value):
humanonlya <- -log10(humanonlya)
 # do graph:

png(file = "C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/Reactome/reactome.enrich.humanonly.top.stemness.png", width = 45*450, height = 58*450, res = 450)#export file to png
par(mar=c(10, max(nchar(colnames(humanonlya)))+65, 5,15))
barplot(humanonlya, horiz = T, beside=T, las=1, space = c(45/ncol(humanonlya),45/ncol(humanonlya)), col=colorRampPalette(c("red"))(n=1), cex.names=4.1, axes=F, xlim=c(0,12)) # both axis have to be set false so that then we can overlaigh a axis with parameters we want
# to use main title with its corresponding parameters (for ex. size) we have to set the title appart from barplot:
# title(main = list("Reactome Pathways - human only", cex = 5))
axis(side=1, at=seq(from=0,to = 12, by=6), lwd=7, tck = -0.04, labels=NA, line = -6) # this is for the ticks only in x axis. that's why lables=NA.side - an integer indicating the side of the graph to draw the axis (1=bottom, 2=left, 3=top, 4=right); at - a numeric vector indicating where tick marks should be drawn; lwd - line width;
axis(side=1, at=seq(from=0,to = 12, by=6), lwd=0, cex.axis=5, line = -2) # this is for the numbering only in x axis. that's why lwd=0.
title(xlab=expression(-log[10](`adj.p-value`)), cex.lab = 5.5, line = 6) # adds x axis label
dev.off()

tiff(file = "C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/Reactome/reactome.enrich.humanonly.top.stemness.tif", width = 45*450, height = 58*450, res = 450)#export file to png
par(mar=c(10, max(nchar(colnames(humanonlya)))+65, 5,15))
barplot(humanonlya, horiz = T, beside=T, las=1, space = c(45/ncol(humanonlya),45/ncol(humanonlya)), col=colorRampPalette(c("red"))(n=1), cex.names=4.1, axes=F, xlim=c(0,12)) # both axis have to be set false so that then we can overlaigh a axis with parameters we want
# to use main title with its corresponding parameters (for ex. size) we have to set the title appart from barplot:
# title(main = list("Reactome Pathways - human only", cex = 5))
axis(side=1, at=seq(from=0,to = 12, by=6), lwd=7, tck = -0.04, labels=NA, line = -6) # this is for the ticks only in x axis. that's why lables=NA.side - an integer indicating the side of the graph to draw the axis (1=bottom, 2=left, 3=top, 4=right); at - a numeric vector indicating where tick marks should be drawn; lwd - line width;
axis(side=1, at=seq(from=0,to = 12, by=6), lwd=0, cex.axis=5, line = -2) # this is for the numbering only in x axis. that's why lwd=0.
title(xlab=expression(-log[10](`adj.p-value`)), cex.lab = 5.5, line = 6) # adds x axis label
dev.off()

# tiff(file = "C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/Reactome/reactome.enrich.humanonly.top.stemness.tiff", width = 70*300, height = 85*300, res = 300)#export file to png
# par(mar=c(10, max(nchar(colnames(humanonlya)))+180, 5, 5))
# barplot(humanonlya, horiz = T, beside=T, las=1, space = c(24/ncol(humanonlya),24/ncol(humanonlya)), col=colorRampPalette(c("red"))(n=1), cex.names=6.5, axes=F, xlim=c(0,12)) # both axis have to be set false so that then we can overlaigh a axis with parameters we want
# # to use main title with its corresponding parameters (for ex. size) we have to set the title appart from barplot:
# # title(main = list("Reactome Pathways - human only", cex = 5))
# axis(side=1, at=seq(from=0,to = 12, by=4), lwd=7, tck = -0.02, labels=NA, line = -10) # this is for the ticks only in x axis. that's why lables=NA.side - an integer indicating the side of the graph to draw the axis (1=bottom, 2=left, 3=top, 4=right); at - a numeric vector indicating where tick marks should be drawn; lwd - line width;
# axis(side=1, at=seq(from=0,to = 12, by=4), lwd=0, cex=5, cex.axis=5, line = -4) # this is for the numbering only in x axis. that's why lwd=0.
# title(xlab=expression(-log[10](`adj.p-value`)), cex.lab = 7, line = 5.5) # adds x axis label
# dev.off()

# Enriched only in mouse and not in human:
mouseonly <- subset(rank.mouse, !(rank.mouse$Description%in%rank.human$Description))
# remove 1st column:
mouseonlya <- mouseonly[,-1]
# transform matrix into vector:
mouseonlya <- as.vector(mouseonlya)
# give terms as names of the vector:
names(mouseonlya) <- mouseonly[,1]
# reverse order of vector (so that later on the terms most enriched are shown first in barplot):
mouseonlya <- rev(mouseonlya)
# change rows by columns:
mouseonlya <- t(mouseonlya)
# transform adjusted p-value into -log(p-value):
mouseonlya <- -log10(mouseonlya)
# do graph:

png(file = "C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/Reactome/reactome.enrich.mouseonly.top.stemness.png", width = 47*450, height = 40*450, res = 450)#export file to png
par(mar=c(15, max(nchar(colnames(mouseonlya)))+80, 5,15))
barplot(mouseonlya, horiz = T, beside=T, las=1, space = c(15/ncol(mouseonlya),15/ncol(mouseonlya)), col=colorRampPalette(c("forestgreen"))(n=1), cex.names=5, axes=F, xlim=c(0,2)) # both axis have to be set false so that then we can overlaigh a axis with parameters we want
# to use main title with its corresponding parameters (for ex. size) we have to set the title appart from barplot:
# title(main = list("Reactome Pathways - mouse only", cex = 5))
axis(side=1, at=seq(from=0,to = 2, by=1), lwd=8, tck = -0.05, labels=NA, line = -2) # this is for the ticks only in x axis. that's why lables=NA.side - an integer indicating the side of the graph to draw the axis (1=bottom, 2=left, 3=top, 4=right); at - a numeric vector indicating where tick marks should be drawn; lwd - line width;
axis(side=1, at=seq(from=0,to = 2, by=1), lwd=0, cex.axis=5.5, line = 2) # this is for the numbering only in x axis. that's why lwd=0.
title(xlab=expression(-log[10](`adj.p-value`)), cex.lab = 6, line = 10) # adds x axis label
dev.off()

tiff(file = "C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/Reactome/reactome.enrich.mouseonly.top.stemness.tif", width = 47*450, height = 40*450, res = 450)#export file to png
par(mar=c(15, max(nchar(colnames(mouseonlya)))+80, 5,15))
barplot(mouseonlya, horiz = T, beside=T, las=1, space = c(15/ncol(mouseonlya),15/ncol(mouseonlya)), col=colorRampPalette(c("forestgreen"))(n=1), cex.names=5, axes=F, xlim=c(0,2)) # both axis have to be set false so that then we can overlaigh a axis with parameters we want
# to use main title with its corresponding parameters (for ex. size) we have to set the title appart from barplot:
# title(main = list("Reactome Pathways - mouse only", cex = 5))
axis(side=1, at=seq(from=0,to = 2, by=1), lwd=8, tck = -0.05, labels=NA, line = -2) # this is for the ticks only in x axis. that's why lables=NA.side - an integer indicating the side of the graph to draw the axis (1=bottom, 2=left, 3=top, 4=right); at - a numeric vector indicating where tick marks should be drawn; lwd - line width;
axis(side=1, at=seq(from=0,to = 2, by=1), lwd=0, cex.axis=5.5, line = 2) # this is for the numbering only in x axis. that's why lwd=0.
title(xlab=expression(-log[10](`adj.p-value`)), cex.lab = 6, line = 10) # adds x axis label
dev.off()

# tiff(file = "C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/Reactome/reactome.enrich.mouseonly.top.stemness.tiff", width = 30*300, height = 30*300, res = 300)#export file to png
# par(mar=c(15, max(nchar(colnames(humanonlya)))+20, 15, 10))
# barplot(mouseonlya, horiz = T, beside=T, las=1, space = c(7/ncol(mouseonlya),7/ncol(mouseonlya)), col=colorRampPalette(c("forestgreen"))(n=1), cex.names=4, axes=F, xlim = c(0,2)) # both axis have to be set false so that then we can overlaigh a axis with parameters we want
# # to use main title with its corresponding parameters (for ex. size) we have to set the title appart from barplot:
# #title(main = list("Reactome Pathways - mouse only", cex = 5))
# axis(side=1, at=seq(from=0,to = 2, by=1), lwd=7, tck = -0.02, labels=NA, line = -2) # this is for the ticks only in x axis. that's why lables=NA.side - an integer indicating the side of the graph to draw the axis (1=bottom, 2=left, 3=top, 4=right); at - a numeric vector indicating where tick marks should be drawn; lwd - line width;
# axis(side=1, at=seq(from=0,to = 2, by=1), lwd=0, cex=5, cex.axis=3.5, line = -1) # this is for the numbering only in x axis. that's why lwd=0.
# title(xlab=expression(-log[10](`adj.p-value`)), cex.lab = 4, line = 5.5) # adds x axis label
# dev.off()


####### OTHER STUFF: in case we want to put the terms that are enriched in human and not in mouse and vice-versa together with those present in both
# both <- merge(rank.human,rank.mouse, by.x=1, by.y=1, all.x= T,all.y=T)
# # replace NAs (terms not signifficant in human but signifficant in mouse and vice-versa) by the correct adj.p-value, so that these terms are also shown in the graph
# # find terms of NAs
# NAs.human <- both$Term[is.na(both$Human)]
# NAs.mouse <- both$Term[is.na(both$Mouse)]
# # find p-values of NAs
# extra.human <- subset(onto.overReactome.Reactome.h1[,c(8,3)],onto.overReactome.Reactome.h1$Term%in%NAs.human) 
# extra.mouse <- subset(onto.overReactome.Reactome.m1[,c(8,3)],onto.overReactome.Reactome.m1$Term%in%NAs.mouse) 
# # put column FDR with the indication if it's from mouse or human
# names(extra.mouse)[2] <- "Mouse"
# names(extra.human)[2] <- "Human"
# # p-values to data.frame
# for(i in 1:nrow(both)){
#         if (is.na(both[i,2])){
#                 both[i,2] <- subset(extra.human[,2], extra.human[,1]%in%both[i,1])
#                 
#         }
# }
# 
# 
# 
# 
# apply(both, 1, function(x) if(is.na(x[2])){x[2] = subset(extra.human[,2], extra.human[,1]%in%x[1])})
