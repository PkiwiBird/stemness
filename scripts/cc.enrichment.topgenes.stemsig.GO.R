min.count <- 0
min.size <- 0
###### TABLES
### Human
## GO IDs and corresponding genes
library(org.Hs.eg.db)
# get go ids of CC only:
go <- as.data.frame(org.Hs.egGO2ALLEGS)
ccid <- go$go_id[go$Ontology=="CC"]
ccid <- unique(ccid)
xx.GO <- as.list(org.Hs.egGO2ALLEGS) # a list in which each element name is a GOID and each element is a vector of all gene entrez IDs corresponding to that GOID
                                     # we could have used: org.Hs.egGOEG instead of org.Hs.egGO2ALLEGS ; but org.Mm.egGO2ALLEGS besides providing mappings between a given GO identifier and
                                     # all of the Entrez Gene identifiers annotated at that GO term also provides annotation between a given GOID and GENIDS OF IT'S CHILD NODES
                                     # in the GO ontology. Thus, this mapping is much larger and more inclusive than org.Mm.egGO2EG.

# list with human cc ids and corresponding genes:
xx.cc <- xx.GO[ccid]
# to get gene universe for Hs CC
allgenes <- unique(unlist(xx.cc))

### top genes (the ones with score equal or higher than rank.limit) from stemness gene sets
setwd("C:/Users/Tania/Internship/R/Results/ranked.lists")
rank.genes <- unique(read.csv("ranked.list.human.stem.signature.csv", header=T, stringsAsFactors=F)) #create a list with the data from all files
rank.limit <- 4
top.rank <- as.character(subset(rank.genes$gene_id, rank.genes$score >= rank.limit))
# NOTE: if we want to work with the top 100 genes instead(not good idea cause this splits genes of the same ranking (rank 4) into the group of the selected ones and not selected ones):
#top.rank <- rank.genes[1:100]

## Hypergeometric test with conditioning
# NOTE: Often an analysis for GO term associations results in the identification of directly related
# GO terms with considerable overlap of genes. This is because each GO term inherits all
# annotations from its more specific descendants. To alleviate this problem, we have implemented
# a method which conditions on all child terms that are themselves significant at a specified pvalue
# cutoff. Given a subgraph of one of the three GO ontologies, we test the leaves of the
# graph, that is, those terms with no child terms. Before testing the terms whose children have
# already been tested, we remove all genes annotated at significant children from the parent's
# gene list. This continues until all terms have been tested.
# The advantage of using GOstats library functions to calculate p-value  over calculating it by ourselves (as shown in the script....) is precisely the possibility to do conditioning
# To activate the conditioning, when creating a GOHyperGParams class use the option conditional=T and pvalueCutoff=0.05
# For a conditional analysis, the cutoff is used during the computation to perform the conditioning: child terms with a p-value less than pvalueCutoff are conditioned out of the test for their parent term.
#  When the test being performed is non-conditional, pvalueCutoff is only used as a default value for printing and suHsarizing the results.

# "GOHyperGParams" is a parameter class (group of parameters) for analyzing for the GO.
# There're others: the KEGGHyperGParams and PFAMHyperGParams
# in geneIds use "rank.genes[rank.genes %in% allgenes]" to use from the genes we want to analyze those that have a corresponding GO CC Hs slim term (are annotated)
# annotation - annotation package used(for the function to know what organism we are using)
# ontology - specifies the GO ontology: is it biological process (CC), molecular function (CC) or celullar compartment (CC) we're analyzing
# test direction:  detects over or under represented GO terms

library(GOstats)
x.2GO <- org.Hs.egSYMBOL 
mapped_genesGO <- mappedkeys(x.2GO) 
entrez2symGO <- as.vector(as.list(x.2GO[mapped_genesGO])) #it's still a list
paramsGO <- new("GOHyperGParams", geneIds =top.rank[top.rank %in% allgenes], universeGeneIds = allgenes, annotation = "org.Hs.eg.db", ontology = "CC", testDirection = "over", conditional=T, pvalueCutoff=0.05, categorySubsetIds = ccid) # all parameters are explained just above
hypGO <- hyperGTest(paramsGO)
tgo <- summary(hypGO, pvalue = 1.1) # only terms with p-value less than the cutoff specified will be shown in the table.
                                    # we use 1.1 to include p-value 1. we include p-value 1 cause 1st we need to calculate adjusted p-value before applying the cutoff and for using these tables for the heatmap
                                    # note: table shows GO terms with size 0: I think it's because some child terms are not considered with conditional=T (when we do the analysis with conditional F and without pvalueCutoff parameter set we don't get GO terms with size 0)                                  
                                    # selection of signifficant terms for the outputed tables is done downstream
go1 <- subset(tgo,tgo$GOCCID%in%ccid) # to have only go.slim terms in the table
go1 <- data.frame(go1[, 1:2], FDR = p.adjust(go1[, 2], method = "fdr", n=length(ccid)),go1[, 3:7]) # !IMPORTANT: there are some no significant GO IDs that are not shown even after selecting p-value threshold as 1.1. So the remaining GO terms although not in the table will have to be included in the p-value adjustment calculation: so we added n=length(ccid)
# we may have a small number of FDR (significant) and still have a small number of counts(small number of genes in the gene set with that GO term)       
# because p-values (which are afterwards transformed by FDR into q-values)compare the difference between two values with the variance (i.e. the range over which the values fall), even with a small number of counts we may have a significant difference when the variance is small 
# besides having a significant value (small FDR) we may also want only the GO terms with a high number of counts, to further narrow the number of GO terms found
go1 <- go1[go1$Count > min.count, ]
# another way of narrowing the number of GO terms found: consider an enriched GO term that only has in total 10 genes associated with it and 9 are in our gene set list, if the same GO term had 1000 genes in total would we still say that those 9 genes show enrichment?
# lists of low sizes do not give much confidence in the enrichment analysis 
go1 <- go1[go1$Size > min.size, ]
        # if number of rows is below 0 (there're no GO terms) then we do not do the next operation
        # for each GO term of the data frame
        # tmp: find the genes in our gene set that are associated with each GO term. [xx.go[[go[i, 1]]] %in% top.rank] is a logical vector with all genes of data set that are associated with a GO term
        # entrez.genes[i]: paste with collapse just to separate the genes by ','
        # sym.genes[i]: convert to gene symbol
if (nrow(go1) > 0) {
        # create empty vectors:
        entrez.genesGO <- character()
        sym.genesGO <- character()
         for (i in 1:nrow(go1)) {
                        tmpGO <- unique(xx.GO[[go1[i, 1]]][xx.GO[[go1[i, 1]]] %in% top.rank])
                        entrez.genesGO[i] <- paste(tmpGO, collapse = ",")
                        sym.genesGO[i] <- paste(entrez2symGO[tmpGO], collapse = ",")
        }
}
# merge onto.over with gene ids and gene symbols for each GO term
onto.overGO.cc.h <- data.frame(go1,entrez.id = entrez.genesGO,gene.symbol = sym.genesGO)
# find correponding term explanation for GOIDs
GOIDs<-go1[,1]
GOIDs.Description<-select(GO.db, keys=GOIDs, columns="DEFINITION", keytype="GOID")
# merge onto.over.h with GO term explanation
onto.overGO.cc.h<-data.frame(onto.overGO.cc.h[,1:8],GO.description=GOIDs.Description[,2],onto.overGO.cc.h[,9:ncol(onto.overGO.cc.h)])
# order data frame from most signifficant term (lower FDR) to less signifficant term:
ord <- order(onto.overGO.cc.h$FDR)      
onto.overGO.cc.h1 <- onto.overGO.cc.h[ord,]
# to select the terms with signifficant FDR:
signifh <- onto.overGO.cc.h1[onto.overGO.cc.h1$FDR < 0.05,]
write.csv(signifh, "C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/GOCC/cc.human.top.stemness.csv", row.names=F)

### Mouse
## GO IDs and corresponding genes
library(org.Mm.eg.db)
# get go ids of CC only:
go <- as.data.frame(org.Mm.egGO2ALLEGS)
ccid <- go$go_id[go$Ontology=="CC"]
ccid <- unique(ccid)
xx.GO <- as.list(org.Mm.egGO2ALLEGS) # a list in which each element name is a GOID and each element is a vector of all gene entrez IDs corresponding to that GOID
# we could have used: org.Mm.egGOEG instead of org.Mm.egGO2ALLEGS ; but org.Mm.egGO2ALLEGS besides providing mappings between a given GO identifier and
# all of the Entrez Gene identifiers annotated at that GO term also provides annotation between a given GOID and GENIDS OF IT'S CHILD NODES
# in the GO ontology. Thus, this mapping is much larger and more inclusive than org.Mm.egGO2EG.

# list with human cc ids and corresponding genes:
xx.cc <- xx.GO[ccid]
# to get gene universe for Hs CC
allgenes <- unique(unlist(xx.cc))

### top genes (the ones with score equal or higher than rank.limit) from stemness gene sets
setwd("C:/Users/Tania/Internship/R/Results/ranked.lists")
rank.genes <- unique(read.csv("ranked.list.mouse.stem.signature.csv", header=T, stringsAsFactors=F)) #create a list with the data from all files
rank.limit <- 7
top.rank <- as.character(subset(rank.genes$gene_id, rank.genes$score >= rank.limit))
# NOTE: if we want to work with the top 100 genes instead(not good idea cause this splits genes of the same ranking (rank 4) into the group of the selected ones and not selected ones):
#top.rank <- rank.genes[1:100]

## Hypergeometric test with conditioning
# NOTE: Often an analysis for GO term associations results in the identification of directly related
# GO terms with considerable overlap of genes. This is because each GO term inherits all
# annotations from its more specific descendants. To alleviate this problem, we have implemented
# a method which conditions on all child terms that are themselves significant at a specified pvalue
# cutoff. Given a subgraph of one of the three GO ontologies, we test the leaves of the
# graph, that is, those terms with no child terms. Before testing the terms whose children have
# already been tested, we remove all genes annotated at significant children from the parent's
# gene list. This continues until all terms have been tested.
# The advantage of using GOstats library functions to calculate p-value  over calculating it by ourselves (as shown in the script....) is precisely the possibility to do conditioning
# To activate the conditioning, when creating a GOHyperGParams class use the option conditional=T and pvalueCutoff=0.05
# For a conditional analysis, the cutoff is used during the computation to perform the conditioning: child terms with a p-value less than pvalueCutoff are conditioned out of the test for their parent term.
#  When the test being performed is non-conditional, pvalueCutoff is only used as a default value for printing and suMmarizing the results.

# "GOHyperGParams" is a parameter class (group of parameters) for analyzing for the GO.
# There're others: the KEGGHyperGParams and PFAMHyperGParams
# in geneIds use "rank.genes[rank.genes %in% allgenes]" to use from the genes we want to analyze those that have a corresponding GO CC Mm slim term (are annotated)
# annotation - annotation package used(for the function to know what organism we are using)
# ontology - specifies the GO ontology: is it biological process (CC), molecular function (CC) or celullar compartment (CC) we're analyzing
# test direction:  detects over or under represented GO terms

library(GOstats)
x.2GO <- org.Mm.egSYMBOL 
mapped_genesGO <- mappedkeys(x.2GO) 
entrez2symGO <- as.vector(as.list(x.2GO[mapped_genesGO])) #it's still a list
paramsGO <- new("GOHyperGParams", geneIds =top.rank[top.rank %in% allgenes], universeGeneIds = allgenes, annotation = "org.Mm.eg.db", ontology = "CC", testDirection = "over", conditional=T, pvalueCutoff=0.05, categorySubsetIds = ccid) # all parameters are explained just above
hypGO <- hyperGTest(paramsGO)
tgo <- summary(hypGO, pvalue = 1.1) # only terms with p-value less than the cutoff specified will be shown in the table.
# we use 1.1 to include p-value 1. we include p-value 1 cause 1st we need to calculate adjusted p-value before applying the cutoff and for using these tables for the heatmap
# note: table shows GO terms with size 0: I think it's because some child terms are not considered with conditional=T (when we do the analysis with conditional F and without pvalueCutoff parameter set we don't get GO terms with size 0)                                  
# selection of signifficant terms for the outputed tables is done downstream
go1 <- subset(tgo,tgo$GOCCID%in%ccid) # to have only go.slim terms in the table
go1 <- data.frame(go1[, 1:2], FDR = p.adjust(go1[, 2], method = "fdr", n=length(ccid)),go1[, 3:7]) # !IMPORTANT: there are some no significant GO IDs that are not shown even after selecting p-value threshold as 1.1. So the remaining GO terms although not in the table will have to be included in the p-value adjustment calculation: so we added n=length(ccid)
# we may have a small number of FDR (significant) and still have a small number of counts(small number of genes in the gene set with that GO term)       
# because p-values (which are afterwards transformed by FDR into q-values)compare the difference between two values with the variance (i.e. the range over which the values fall), even with a small number of counts we may have a significant difference when the variance is small 
# besides having a significant value (small FDR) we may also want only the GO terms with a high number of counts, to further narrow the number of GO terms found
go1 <- go1[go1$Count > min.count, ]
# another way of narrowing the number of GO terms found: consider an enriched GO term that only has in total 10 genes associated with it and 9 are in our gene set list, if the same GO term had 1000 genes in total would we still say that those 9 genes show enrichment?
# lists of low sizes do not give much confidence in the enrichment analysis 
go1 <- go1[go1$Size > min.size, ]
# if number of rows is below 0 (there're no GO terms) then we do not do the next operation
# for each GO term of the data frame
# tmp: find the genes in our gene set that are associated with each GO term. [xx.go[[go[i, 1]]] %in% top.rank] is a logical vector with all genes of data set that are associated with a GO term
# entrez.genes[i]: paste with collapse just to separate the genes by ','
# sym.genes[i]: convert to gene symbol
if (nrow(go1) > 0) {
        # create empty vectors:
        entrez.genesGO <- character()
        sym.genesGO <- character()
        for (i in 1:nrow(go1)) {
                tmpGO <- unique(xx.GO[[go1[i, 1]]][xx.GO[[go1[i, 1]]] %in% top.rank])
                entrez.genesGO[i] <- paste(tmpGO, collapse = ",")
                sym.genesGO[i] <- paste(entrez2symGO[tmpGO], collapse = ",")
        }
}
# merge onto.over with gene ids and gene symbols for each GO term
onto.overGO.cc.m <- data.frame(go1,entrez.id = entrez.genesGO,gene.symbol = sym.genesGO)
# find correponding term explanation for GOIDs
GOIDs<-go1[,1]
GOIDs.Description<-select(GO.db, keys=GOIDs, columns="DEFINITION", keytype="GOID")
# merge onto.over.m with GO term explanation
onto.overGO.cc.m<-data.frame(onto.overGO.cc.m[,1:8],GO.description=GOIDs.Description[,2],onto.overGO.cc.m[,9:ncol(onto.overGO.cc.m)])
# order data frame from most signifficant term (lower FDR) to less signifficant term:
ord <- order(onto.overGO.cc.m$FDR)      
onto.overGO.cc.m1 <- onto.overGO.cc.m[ord,]
# to select the terms with signifficant FDR:
signifm <- onto.overGO.cc.m1[onto.overGO.cc.m1$FDR < 0.05,]
write.csv(signifm, "C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/GOCC/cc.mouse.top.stemness.csv", row.names=F)

### Bar plots (human and mouse together in same barplot)
setwd("C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/GOCC")
rank.human <- read.csv("cc.human.top.stemness.csv", header=T, stringsAsFactors=F)
rank.human <- rank.human[,c(8,3)]
rank.mouse <- read.csv("cc.mouse.top.stemness.csv", header=T, stringsAsFactors=F)
rank.mouse <- rank.mouse[,c(8,3)]
# put column FDR with the indication if it's from mouse or human
names(rank.mouse)[2] <- "Mouse"
names(rank.human)[2] <- "Human"
# Enriched in both mouse and human - put human and mouse in the same matrix:
# NOTE!!!: as there are some CC terms for mouse that don't exist in human and vice-versa
#          it doesn't make sense to input all terms of mouse or human. we could in the case a term was signifficant for mouse for example to put the p-value column of human as 1, but in that case the reader would interpret that that term wasn't enriched in human, while in fact that term didn't exist in human
#          SO WE LOOK FOR THE ENRICHED TERMS IN BOTH!
both <- merge(rank.mouse, rank.human, by.x=1, by.y=1, all=F)
# give 1st column of matrix as row.names:
row.names(both) <- both[,1]
# remove 1st column:
both <- both[,-1]
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
# png(file = "C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/GOCC/cc.enrich.bothmousehuman.top.stemness.png", width = 35*300, height = 30*300, res = 300)#export file to png
# par(mar=c(20, max(nchar(colnames(botha)))+20, 20,5))
# barplot(botha, horiz = T, beside=T, las=1, space = c(3/ncol(botha),6/ncol(botha)), col=colorRampPalette(c("forestgreen","red"))(n=2), cex.names=4, legend.text=T, args.legend=list(cex = 3), axes = F, xlim = c(0,18)) # both axis have to be set false so that then we can overlaigh a axis with parameters we want
# # to use main title with its corresponding parameters (for ex. size) we have to set the title appart from barplot:
# # title(main = list("Celular Component - both mouse and human", cex = 5))
# axis(side=1, at=seq(from=0,to = 18, by=2), lwd=4, tck = -0.01, labels=NA, line = -0.5) # this is for the ticks only in x axis. that's why lables=NA.side - an integer indicating the side of the graph to draw the axis (1=bottom, 2=left, 3=top, 4=right); at - a numeric vector indicating where tick marks should be drawn; lwd - line width;
# axis(side=1, at=seq(from=0,to = 18, by=2), lwd=0, cex=5, cex.axis=3.5, line = 2) # this is for the numbering only in x axis. that's why lwd=0.
# title(xlab=expression(-log[10](`adj.p-value`)), cex.lab = 4.5, line = 9) # adds x axis label
# dev.off()

png(file = "C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/GOCC/cc.enrich.bothmousehuman.top.stemness.png", width = 25*600, height = 15*600, res = 600)#export file to png
par(mar=c(14, max(nchar(colnames(botha)))+30, 15,10))
#barplot(botha, horiz = T, beside=T, las=1, space = c(15/ncol(botha),45/ncol(botha)), col=colorRampPalette(c("forestgreen","red"))(n=2), cex.names=4.9, legend.text=T, args.legend=list(x=0, y=0, cex = 3), axes = F, xlim = c(0,20)) # both axis have to be set false so that then we can overlaigh a axis with parameters we want
barplot(botha, horiz = T, beside=T, las=1, space = c(1/ncol(botha),7.5/ncol(botha)), col=colorRampPalette(c("forestgreen","red"))(n=2), cex.names=4.9, axes = F, xlim = c(0,18)) # both axis have to be set false so that then we can overlaigh a axis with parameters we want
# to use main title with its corresponding parameters (for ex. size) we have to set the title appart from barplot:
# title(main = list("Celular Component - both mouse and human", cex = 5))
#par(cex.axis=1, cex.lab=3, cex.main=500)
axis(side=1, at=seq(from=0,to = 18, by=4), lwd=5, tck = -0.02, labels=NA, line = 2) # this is for the ticks only in x axis. that's why lables=NA.side - an integer indicating the side of the graph to draw the axis (1=bottom, 2=left, 3=top, 4=right); at - a numeric vector indicating where tick marks should be drawn; lwd - line width;
axis(side=1, at=seq(from=0,to = 18, by=4), lwd=0, cex.axis=5, line = 5) # this is for the numbering only in x axis. that's why lwd=0.
title(xlab=expression(-log[10](`adj.p-value`)), cex.lab =5.5, line = 13) # adds x axis label
dev.off()

tiff(file = "C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/GOCC/cc.enrich.bothmousehuman.top.stemness.tif", width = 25*600, height = 15*600, res = 600)#export file to png
par(mar=c(14, max(nchar(colnames(botha)))+30, 15,10))
#barplot(botha, horiz = T, beside=T, las=1, space = c(15/ncol(botha),45/ncol(botha)), col=colorRampPalette(c("forestgreen","red"))(n=2), cex.names=4.9, legend.text=T, args.legend=list(x=0, y=0, cex = 3), axes = F, xlim = c(0,20)) # both axis have to be set false so that then we can overlaigh a axis with parameters we want
barplot(botha, horiz = T, beside=T, las=1, space = c(1/ncol(botha),7.5/ncol(botha)), col=colorRampPalette(c("forestgreen","red"))(n=2), cex.names=4.9, axes = F, xlim = c(0,18)) # both axis have to be set false so that then we can overlaigh a axis with parameters we want
# to use main title with its corresponding parameters (for ex. size) we have to set the title appart from barplot:
# title(main = list("Celular Component - both mouse and human", cex = 5))
#par(cex.axis=1, cex.lab=3, cex.main=500)
axis(side=1, at=seq(from=0,to = 18, by=4), lwd=5, tck = -0.02, labels=NA, line = 2) # this is for the ticks only in x axis. that's why lables=NA.side - an integer indicating the side of the graph to draw the axis (1=bottom, 2=left, 3=top, 4=right); at - a numeric vector indicating where tick marks should be drawn; lwd - line width;
axis(side=1, at=seq(from=0,to = 18, by=4), lwd=0, cex.axis=5, line = 5) # this is for the numbering only in x axis. that's why lwd=0.
title(xlab=expression(-log[10](`adj.p-value`)), cex.lab =5.5, line = 13) # adds x axis label
dev.off()


# Note: it gives error in this case cause there are no common enriched terms

# Enriched only in human and not in mouse:
humanonly <- subset(rank.human, !(rank.human$Term%in%rank.mouse$Term))
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
png(file = "C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/GOCC/cc.enrich.humanonly.top.stemness.png", width = 15*600, height = 7*600, res = 600)#export file to png
par(mar=c(10, max(nchar(colnames(humanonlya)))+5, 5,15))
barplot(humanonlya, horiz = T, beside=T, las=1, space = c(7/ncol(humanonlya),7/ncol(humanonlya)), col=colorRampPalette(c("red"))(n=1), cex.names=2, axes=F, xlim=c(0,14)) # both axis have to be set false so that then we can overlaigh a axis with parameters we want
# to use main title with its corresponding parameters (for ex. size) we have to set the title appart from barplot:
# title(main = list("Celular Component - human only", cex = 5))
axis(side=1, at=seq(from=0,to = 14, by=7), lwd=6, tck = -0.02, labels=NA, line = 0) # this is for the ticks only in x axis. that's why lables=NA.side - an integer indicating the side of the graph to draw the axis (1=bottom, 2=left, 3=top, 4=right); at - a numeric vector indicating where tick marks should be drawn; lwd - line width;
axis(side=1, at=seq(from=0,to = 14, by=7), lwd=0, cex.axis=2, line = 1) # this is for the numbering only in x axis. that's why lwd=0.
title(xlab=expression(-log[10](`adj.p-value`)), cex.lab = 2, line = 6) # adds x axis label
dev.off()

tiff(file = "C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/GOCC/cc.enrich.humanonly.top.stemness.tif", width = 15*600, height = 7*600, res = 600)#export file to png
par(mar=c(10, max(nchar(colnames(humanonlya)))+5, 5,15))
barplot(humanonlya, horiz = T, beside=T, las=1, space = c(7/ncol(humanonlya),7/ncol(humanonlya)), col=colorRampPalette(c("red"))(n=1), cex.names=2, axes=F, xlim=c(0,14)) # both axis have to be set false so that then we can overlaigh a axis with parameters we want
# to use main title with its corresponding parameters (for ex. size) we have to set the title appart from barplot:
# title(main = list("Celular Component - human only", cex = 5))
axis(side=1, at=seq(from=0,to = 14, by=7), lwd=6, tck = -0.02, labels=NA, line = 0) # this is for the ticks only in x axis. that's why lables=NA.side - an integer indicating the side of the graph to draw the axis (1=bottom, 2=left, 3=top, 4=right); at - a numeric vector indicating where tick marks should be drawn; lwd - line width;
axis(side=1, at=seq(from=0,to = 14, by=7), lwd=0, cex.axis=2, line = 1) # this is for the numbering only in x axis. that's why lwd=0.
title(xlab=expression(-log[10](`adj.p-value`)), cex.lab = 2, line = 6) # adds x axis label
dev.off()

# tiff(file = "C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/GOCC/cc.enrich.humanonly.top.stemness.tiff", width = 35*300, height = 15*300, res = 300)#export file to png
# par(mar=c(23, max(nchar(colnames(humanonlya)))+20, 10, 5))
# barplot(humanonlya, horiz = T, beside=T, las=1, space = c(5/ncol(humanonlya),5/ncol(humanonlya)), col=colorRampPalette(c("red"))(n=1), cex.names=3.5, axes=F, xlim=c(0,14)) # both axis have to be set false so that then we can overlaigh a axis with parameters we want
# # to use main title with its corresponding parameters (for ex. size) we have to set the title appart from barplot:
# # title(main = list("Celular Component - human only", cex = 5))
# axis(side=1, at=seq(from=0,to = 14, by=2), lwd=7, tck = -0.02, labels=NA, line = 1) # this is for the ticks only in x axis. that's why lables=NA.side - an integer indicating the side of the graph to draw the axis (1=bottom, 2=left, 3=top, 4=right); at - a numeric vector indicating where tick marks should be drawn; lwd - line width;
# axis(side=1, at=seq(from=0,to = 14, by=2), lwd=0, cex=2.5, cex.axis=4, line = 4) # this is for the numbering only in x axis. that's why lwd=0.
# title(xlab=expression(-log[10](`adj.p-value`)), cex.lab = 4, line = 12) # adds x axis label
# dev.off()

# Enriched only in mouse and not in human:
mouseonly <- subset(rank.mouse, !(rank.mouse$Term%in%rank.human$Term))
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

png(file = "C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/GOCC/cc.enrich.mouseonly.top.stemness.png", width = 15*600, height = 6*600, res = 600)#export file to png
par(mar=c(10, max(nchar(colnames(mouseonlya)))+5, 5,10))
barplot(mouseonlya, horiz = T, beside=T, las=1, space = c(7/ncol(mouseonlya),7/ncol(mouseonlya)), col=colorRampPalette(c("forestgreen"))(n=1), cex.names=1.5, axes=F, xlim=c(0,4)) # both axis have to be set false so that then we can overlaigh a axis with parameters we want
# to use main title with its corresponding parameters (for ex. size) we have to set the title appart from barplot:
# title(main = list("Celular Component - human only", cex = 5))
axis(side=1, at=seq(from=0,to = 4, by=2), lwd=4, tck = -0.02, labels=NA, line = 0) # this is for the ticks only in x axis. that's why lables=NA.side - an integer indicating the side of the graph to draw the axis (1=bottom, 2=left, 3=top, 4=right); at - a numeric vector indicating where tick marks should be drawn; lwd - line width;
axis(side=1, at=seq(from=0,to = 4, by=2), lwd=0, cex.axis=1.5, line = 0) # this is for the numbering only in x axis. that's why lwd=0.
title(xlab=expression(-log[10](`adj.p-value`)), cex.lab = 1.5, line = 3) # adds x axis label
dev.off()

tiff(file = "C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/GOCC/cc.enrich.mouseonly.top.stemness.tif", width = 15*600, height = 6*600, res = 600)#export file to png
par(mar=c(10, max(nchar(colnames(mouseonlya)))+5, 5,10))
barplot(mouseonlya, horiz = T, beside=T, las=1, space = c(7/ncol(mouseonlya),7/ncol(mouseonlya)), col=colorRampPalette(c("forestgreen"))(n=1), cex.names=1.5, axes=F, xlim=c(0,4)) # both axis have to be set false so that then we can overlaigh a axis with parameters we want
# to use main title with its corresponding parameters (for ex. size) we have to set the title appart from barplot:
# title(main = list("Celular Component - human only", cex = 5))
axis(side=1, at=seq(from=0,to = 4, by=2), lwd=4, tck = -0.02, labels=NA, line = 0) # this is for the ticks only in x axis. that's why lables=NA.side - an integer indicating the side of the graph to draw the axis (1=bottom, 2=left, 3=top, 4=right); at - a numeric vector indicating where tick marks should be drawn; lwd - line width;
axis(side=1, at=seq(from=0,to = 4, by=2), lwd=0, cex.axis=1.5, line = 0) # this is for the numbering only in x axis. that's why lwd=0.
title(xlab=expression(-log[10](`adj.p-value`)), cex.lab = 1.5, line = 3) # adds x axis label
dev.off()

# tiff(file = "C:/Users/Tania/Internship/R/Results/enrichment.top/stemness_signatures/GOCC/cc.enrich.mouseonly.top.stemness.tiff", width = 50*300, height = 15*300, res = 300)#export file to png
# par(mar=c(17, max(nchar(colnames(mouseonlya)))+53, 5, 5))
# barplot(mouseonlya, horiz = T, beside=T, las=1, space = c(3/ncol(mouseonlya),3/ncol(mouseonlya)), col=colorRampPalette(c("forestgreen"))(n=1), cex.names=5, axes=F, xlim = c(0,4)) # both axis have to be set false so that then we can overlaigh a axis with parameters we want
# # to use main title with its corresponding parameters (for ex. size) we have to set the title appart from barplot:
# #title(main = list("Celular Component - mouse only", cex = 5))
# axis(side=1, at=seq(from=0,to = 4, by=1), lwd=7, tck = -0.02, labels=NA, line = 2) # this is for the ticks only in x axis. that's why lables=NA.side - an integer indicating the side of the graph to draw the axis (1=bottom, 2=left, 3=top, 4=right); at - a numeric vector indicating where tick marks should be drawn; lwd - line width;
# axis(side=1, at=seq(from=0,to = 4, by=1), lwd=0, cex=5, cex.axis=3.5, line = 4) # this is for the numbering only in x axis. that's why lwd=0.
# title(xlab=expression(-log[10](`adj.p-value`)), cex.lab =5, line = 12) # adds x axis label
# dev.off()

####### OTHER STUFF: in case we want to put the terms that are enriched in human and not in mouse and vice-versa
# both <- merge(rank.human,rank.mouse, by.x=1, by.y=1, all.x= T,all.y=T)
# # replace NAs (terms not signifficant in human but signifficant in mouse and vice-versa) by the correct adj.p-value, so that these terms are also shown in the graph
# # find terms of NAs
# NAs.human <- both$Term[is.na(both$Human)]
# NAs.mouse <- both$Term[is.na(both$Mouse)]
# # find p-values of NAs
# extra.human <- subset(onto.overGO.cc.h1[,c(8,3)],onto.overGO.cc.h1$Term%in%NAs.human) 
# extra.mouse <- subset(onto.overGO.cc.m1[,c(8,3)],onto.overGO.cc.m1$Term%in%NAs.mouse) 
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
