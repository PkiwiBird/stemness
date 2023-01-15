# !!! THIS SCRIPT WAS RUN IN MY NEW ASUS CAUSE THE r VERSION INSTALLED IN THE TOSHIBA DOESN'T ALLOW FOR THE COMMAND ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "uswest.ensembl.org", ensemblRedirect = F) TO WORK
setwd('/home/tania/Desktop/stemness/')
## REACTOME IDs and corresponding genes
library(reactome.db)
library(ReactomePA)
library(org.Mm.eg.db)
library(xlsx)

# get reactome ids of reactome human only:
xx.reactome <- as.list(reactomeEXTID2PATHID) # a list in which each element name is a entrez gene id and each element is a vector of all reactome ids corresponding to that gene
reactomegene <- unique(names(xx.reactome)) # all reactome genes
reactomeid <- unique(unlist(xx.reactome))# all reactome ids
library(biomaRt)
# ensembl = useMart(biomart ="ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host="ensembl.org")
# when ensembl main site is down use instead the mirrors uswest.ensembl.org or useast.ensembl.org: 
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "uswest.ensembl.org", ensemblRedirect = F)
d <- getBM(attributes = "entrezgene",values = "*",mart = ensembl)
h <- as.character(unique(d[,1])) # gets all mouse genes
allgenes <- h[h%in%reactomegene]# gets all reactome genes that are mouse!
# get gene sets:
files_list1 <- list.files("data/clusters_forenrichment_data/mouse", full.names=TRUE)   #list files in the folder
g <- lapply(files_list1, function(x) as.character(read.xlsx(x,1, as.data.frame=T, header=F, colClasses = c("character"))[,1]))
datasetnames<-vector()
for (i in 1:length(files_list1)){
        datasetnames[i]<-unlist(strsplit(files_list1[i],'[/.]'))[4]   #parse the names for the list of data
}
names(g) <- datasetnames
# convert gene symbol into entrez id:
# to see list of atributes: l <- listAttributes(mart=ensembl)
# to see list of filters: l <- listAttributes(mart=ensembl)

genesets <- list()
for (i in 1:length(g)){
        conv <- getBM(attributes=c('mgi_symbol','entrezgene'),filters = 'mgi_symbol',values = g[[i]] ,mart = ensembl)
        conv2 <- conv[,2] 
        genesets[[i]] <- conv2
}
names(genesets) <- datasetnames
# small number of duplicates: continue
# do enrichment for reactome:
signifh <- list()
for (j in 1:length(genesets)){
        if (j != 3){ # because the 3rd geneset didn't contain any genes in reactome and so was giving error and preventing the code to run
                treact <- enrichPathway(gene=as.character(genesets[[j]]), organism="mouse", pvalueCutoff=2, qvalueCutoff = 2, pAdjustMethod="fdr", universe=allgenes, readable=T, minGSSize=1)
                treact2 <- treact@result
                signifh[[j]] <- treact2
        }
}
names(signifh) <- datasetnames
signifh2 <- lapply(signifh, function(x) subset(x, x$p.adjust < 0.01))
lapply(names(signifh2), function (x) write.csv(signifh2[[x]], paste(paste("/home/tania/Desktop/stemness/results/clusters_forenrichment_results/mouse",x,sep="/"),"csv",sep="."), row.names=F))


