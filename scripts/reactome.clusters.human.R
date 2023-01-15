# !!! THIS SCRIPT WAS RUN IN MY NEW ASUS CAUSE THE r VERSION INSTALLED IN THE TOSHIBA DOESN'T ALLOW FOR THE COMMAND ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "uswest.ensembl.org", ensemblRedirect = F) TO WORK
setwd('/home/tania/Desktop/stemness/')
# NOTE! when installing packages with R version 3.2 use biocLite('packagename') instead of install.packages()
## REACTOME IDs and corresponding genes
library(reactome.db)
library(ReactomePA)
library(org.Hs.eg.db)
library(xlsx)

# get reactome ids of reactome human only:
xx.reactome <- as.list(reactomeEXTID2PATHID) # a list in which each element name is a entrez gene id and each element is a vector of all reactome ids corresponding to that gene
reactomegene <- unique(names(xx.reactome)) # all reactome genes
reactomeid <- unique(unlist(xx.reactome))# all reactome ids
library(biomaRt)
# ensembl = useMart(biomart ="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="ensembl.org")
# when ensembl main site is down use instead the mirrors uswest.ensembl.org or useast.ensembl.org: 
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "uswest.ensembl.org", ensemblRedirect = F)
d = getBM(attributes = 'entrezgene', values = '*', mart = ensembl)
h <- as.character(unique(d[,1])) # gets all human genes
allgenes <- h[h%in%reactomegene]# gets all reactome genes that are human!
# get gene sets:
files_list1 <- list.files("data/clusters_forenrichment_data/human", full.names=TRUE)   #list files in the folder
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
conv <- getBM(attributes=c('hgnc_symbol','entrezgene'),filters = 'hgnc_symbol',values = g[[i]] ,mart = ensembl)
conv2 <- conv[,2] 
genesets[[i]] <- conv2
}
names(genesets) <- datasetnames
 # small number of duplicates: continue
# do enrichment for reactome:
signifh <- list()
for (j in 1:length(genesets)){
        treact <- enrichPathway(gene=genesets[[j]], organism="human", pvalueCutoff=2, qvalueCutoff = 2, pAdjustMethod="fdr", universe=allgenes, readable=T, minGSSize=1)
        treact2 <- treact@result
        signifh[[j]] <- treact2
}
names(signifh) <- datasetnames
signifh2 <- lapply(signifh, function(x) subset(x, x$p.adjust < 0.01))
lapply(names(signifh2), function (x) write.csv(signifh2[[x]], paste(paste("/home/tania/Desktop/stemness/results/clusters_forenrichment_results/human",x,sep="/"),"csv",sep="."), row.names=F))


#### DON'T USE CODE BELOW!
#### CODE BELOW IS FOR OTHER THINGS LIKE HEATMAPS AND NEED TO BE ADJUSTED TO THIS ANALYSIS

## Heatmap gene sets vs gene sets for Reactome # note: we don't want to have to run everything from the beginning if something changes so we're importing the tables outputed above
setwd("C:/Users/Tania_2/Internship/R/Results/enrichment.each.list/stemness_signatures/Reactome")
files_list <- list.files("human", full.names=TRUE)   #list files in the folder
term <- lapply(files_list, function(x) as.character(unique(read.csv(x, header=T, stringsAsFactors=F)[,2]))) # list of gene sets with enriched kegg  human terms inside
datasetnames <- vector()
for (i in 1:length(files_list)){
        datasetnames[i]<-unlist(strsplit(files_list[i],'[/.]'))[2]   #parse the names for the list of data
}
# correspond gene set names to the genes
names(term) <- datasetnames

q <- matrix(0,nrow=length(term),ncol=length(term))
dimnames(q) <- list(names(term),names(term))

for (i in 1:nrow(q)){
        for (j in 1:ncol(q)){
                q[i,j] <- sum(unique(term[[i]]) %in% unique(term[[j]])) # counting the number of common genes
        }
}

m <- matrix(0,nrow=length(term),ncol=length(term))
dimnames(m) <- list(names(term),names(term))
# in this case we consider the universe all reactome human terms:
react <- as.data.frame(reactomePATHNAME2ID) # dataframe with reactome ids in 1st column and reactome terms/names in the 2nd
hreactometerm <- unique(react$path_name[grepl("Homo sapiens:",react$path_name)]) # all reactome human terms
# remove leading "Homo sapiens: ":
hreactometerm <- unique(gsub("Homo sapiens: ", "", hreactometerm))
for (i in 1:nrow(m)){
        for (j in 1:ncol(m)){
                m[i,j] <- sum(unique(term[[i]]) %in% hreactometerm)
        }
}

n <- matrix(0,nrow=length(term),ncol=length(term))
dimnames(n) <- list(names(term),names(term))
for (i in 1:nrow(n)){
        for (j in 1:ncol(n)){
                n[i,j] <- length(hreactometerm) - m[[i]] # counting the number of common genes
        }
}

k <- matrix(0,nrow=length(term),ncol=length(term))
dimnames(k) <- list(names(term),names(term))
for (i in 1:nrow(k)){
        for (j in 1:ncol(k)){
                k[i,j] <- sum(unique(term[[j]]) %in% hreactometerm) # counting the number of common genes
        }
}

hyp.test <- matrix(0,nrow=length(term),ncol=length(term))
dimnames(hyp.test) <- list(names(term),names(term))
for (i in 1:nrow(hyp.test)){
        for (j in 1:ncol(hyp.test)){
                hyp.test[i,j] <- phyper(q[i,j]-1,m[i,j],n[i,j],k[i,j], lower.tail = F)
        }
}

hyp.test.adj <- matrix(0,nrow=length(term),ncol=length(term))
dimnames(hyp.test.adj) <- list(names(term),names(term))
for (i in 1: nrow(hyp.test.adj)){
        for (j in 1:ncol(hyp.test.adj)){
                hyp.test.adj[i,j] <- p.adjust(hyp.test[i,j], method = "bonferroni")
        }
}
# Some adjusted p-values get above one, but probabilities can't be above one. So some people force the values to be one...
hyp.test.adj[hyp.test.adj > 1] <- 1
# if we wanted to use logarythm:
# before doing the log we have to replace the 0's by 10E-16 because otherwise we would get -Inf values after logaritmizing and the heatmap function wouldn't be able to plot -inf values (there would be an error).
hyp.test.adj[hyp.test.adj ==0] <- 10E-16 # we could had have also replaced by the smallest value in the matrix (smallest p-value)
log.hyp.test.adj <- log10(hyp.test.adj)
log.hyp.test.adj.thres <- log.hyp.test.adj
# p-values lower than 10E-16  are set to 10E-16:The reasons are the following:
# from the statistical point of view, there is not much difference between a p-value of  10E-16 and 10E-72 (both are highly significant) , but much more between 10E-1 and 10E-5.
# However, since the heatmap.2 uses in its default version a eucleadian clustering distance, the 10E-16 and 10E-72 is much more different, than 10E-1 and 10E-5.
# this is used when in the p-value matrix there're very small values, so the smallest will determine the clustering.
# To prevent this, we can set a min. threshold for significance.
# In this case we don't have very small values in the matrix, but since this data is going to be compared with gene set clustering heatmaps with genes which had this step of "seting a min. threshold for significance" we do the same here. so it doesn't affect the comparison
log.hyp.test.adj.thres[log.hyp.test.adj.thres < -16] <- -16
library(gplots)
my_palette <- colorRampPalette(c("red","beige","cornflowerblue"))(n = 1499)
breaks = c(seq(log10(10E-16),log10(0.049),length=1000),seq(log10(0.05),log10(0.25),length=200),seq(log10(0.251),log10(1),length=300)) # number of breaks should be number of colors +1 (length(breaks) = length(my_palette) + 1)
lmat = rbind(c(0,4,0),c(0,3,0),c(2,1,0))
lwid = c(0.5,4,0.5)
lhei = c(0.9,0.5,4)
png(file="C:/Users/Tania_2/Internship/R/Results/clustering.genesets.byterms/Reactome/Hs_stemness_up.png", width = 8*300, height = 9*300, res = 300)#export file to png
heatmap.2(log.hyp.test.adj.thres,scale="none",col=my_palette,breaks=breaks,density.info="none",trace="none", margins =c(13,10),  lmat = lmat, lwid = lwid, lhei = lhei, key.xlab="", key = T, key.title=F, symkey=F, key.par=list(mar=c(4,1,4,10)), sepcolor="black", colsep=1:ncol(log.hyp.test.adj.thres), rowsep=1:nrow(log.hyp.test.adj.thres), sepwidth=c(0.01,0.01),revC=T)
dev.off()

## Heatmap
# get entrez gene ids of each gene set:
setwd("C:/Users/Tania_2/Internship/R/Data")
# get gene sets:
files_list1 <- list.files("Hs_up_without_TF", full.names=TRUE)   #list files in the folder
genesets <- lapply(files_list1, function(x) as.character(unique(read.csv(x, header=F, stringsAsFactors=F)[,1]))) #create a list with the data from all files
datasetnames<-vector()
for (i in 1:length(files_list1)){
        datasetnames[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(genesets) <- datasetnames

# do enrichment for reactome:
datah <- list()
for (i in 1:length(genesets)){
        treact <- enrichPathway(gene=genesets[[i]], organism="human", pvalueCutoff=2, qvalueCutoff = 2, pAdjustMethod="fdr", universe=allgenes, readable=T, minGSSize=1)
        treact2 <- treact@result
        # select only the term and adj. p-value:
        reac1 <- treact2[,c(2,6)]
        # to select the reactome human terms that are not in each table (remember that terms with p-value =1 and Count=0 are not included and to build the heatmap we need all terms):
        extra <- unique(hreactometerm[!(hreactometerm%in%reac1$Description)])
        # build dataframes with reactome terms that are not in each table and the correponding adj p-value (which is 1)
        extra3 <- data.frame(Description=extra, p.adjust=1)
        # bind with reactome terms that exist in each table:
        frame <- rbind(reac1,extra3) 
        datah[[i]] <- frame
}
names(datah) <- names(genesets)

# put name of gene set as name of the p.adjust column:
for (i in 1:length(datah)){
        colnames(datah[[i]]) <- c("TERM", names(datah[i]))
}
# merge all dataframes in the list by the term
f <- function(x,y) merge(x,y, by.x="TERM",by.y="TERM", all.x=T, all.y=T)
m <- unique(Reduce (f,datah)) # reduce applies the function 1st to 1st and 2nd elements, than applies the result of merging 1st and 2nd elemnts to the 3rd elemnt and so on. does an iteration
#put 1st column (names of pathways/terms) as names of matrix
rownames(m) <- m[,1]
# remove 1st column:
m <- m[,-1]
# remove rows that have no enrichment for all columns (pathways with no enrichment for any of the sets)
m2 <- m[apply(m, 1, function(y) !all(y >= 0.01)),]
# Some adjusted p-values get above one, but probabilities can't be above one. So some people force the values to be one...
m2[m2 > 1] <- 1
# if we wanted to use logarythm:
# before doing the log we have to replace the 0's by 10E-16 because otherwise we would get -Inf values after logaritmizing and the heatmap function wouldn't be able to plot -inf values (there would be an error).
m2[m2==0] <- 10E-16 # we could had have also replaced by the smallest value in the matrix (smallest p-value)
m3 <- log10(m2)
# remove Nas:
m3 <- m3[complete.cases(m3),]
# convert into a matrix:
m3 <- as.matrix(m3)

# p-values lower than 10E-16  are set to 10E-16:The reasons are the following:
# from the statistical point of view, there is not much difference between a p-value of  10E-16 and 10E-72 (both are highly significant) , but much more between 10E-1 and 10E-5.
# However, since the heatmap.2 uses in its default version a eucleadian clustering distance, the 10E-16 and 10E-72 is much more different, than 10E-1 and 10E-5.
# this is used when in the p-value matrix there're very small values, so the smallest will determine the clustering.
# To prevent this, we can set a min. threshold for significance.
# In this case we don't have very small values in the matrix, but since this data is going to be compared with gene set clustering heatmaps with genes which had this step of "seting a min. threshold for significance" we do the same here. so it doesn't affect the comparison
#m4 <- m3
#m4[m4 < -16] <- -16
library(gplots)
my_palette <- colorRampPalette(c("red","white","forestgreen"))(n = 1499)
breaks = c(seq(log10(10E-16),log10(0.049),length=1000),seq(log10(0.05),log10(0.25),length=200),seq(log10(0.251),log10(1),length=300)) # number of breaks should be number of colors +1 (length(breaks) = length(my_palette) + 1)
lmat = rbind(c(0,4,0),c(0,3,0),c(2,1,0))               
lwid = c(0.5,4,0.5)
lhei = c(0.3,0.5,4)
png(file="C:/Users/Tania_2/Internship/R/Results/heatmap.4.enrichment/Reactome/Hs_stemness_up.png", width = 16*300, height = 23*300, res = 300)#export file to png
heatmap.2(m3,scale="none",col=my_palette,breaks=breaks,density.info="none",trace="none", margins =c(25,30), cexRow = 1.5, cexCol = 1.5,  lmat = lmat, lwid = lwid, lhei = lhei, key.xlab="", key = T, key.title=F, symkey=F, key.par=list(mar=c(4,1,4,10)), sepcolor="black", colsep=1:ncol(m3), rowsep=1:nrow(m3), sepwidth=c(0.01,0.01))
dev.off()
