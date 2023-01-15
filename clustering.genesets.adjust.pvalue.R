# !!!NOTE: Here we consider as universe the human annotated genes or all human genes that exist (for m1, n1, k1)
#          a more appropriate aproach for future analysis is to use as universe the union of genes tested in of all gene sets experiments 
setwd("D:/Tania_2/Internship/R/Data")
#####################################################
#### Hs stemness gene sets
files_list1 <- list.files("Hs_up_without_TF", full.names=TRUE)   #list files in the folder
data_files_list <- lapply(files_list1, function(x) as.character(unique(read.csv(x, header=F, stringsAsFactors=F)[,1]))) #create a list with the data from all files
datasetnames<-vector()
for (i in 1:length(files_list1)){
        datasetnames[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list) <- datasetnames   #give names to each list element

# consider one gene set(the one in the row) to be sucesses and the other gene set (the one in the column) to be the sample:
# q is the number of successes in the sample: number of overlapping genes between the two gene sets
q <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
dimnames(q) <- list(names(data_files_list),names(data_files_list))

for (i in 1:nrow(q)){
        for (j in 1:ncol(q)){
                q[i,j] <- sum(unique(data_files_list[[i]]) %in% unique(data_files_list[[j]])) # counting the number of common genes  
        }
}

# m is the number of sucesses in the population: the number of genes of the first gene set that is in the genome
library(org.Hs.eg.db)
m <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
dimnames(m) <- list(names(data_files_list),names(data_files_list))
allgenes <- unique(mappedkeys(org.Hs.egGO))
for (i in 1:nrow(m)){
        for (j in 1:ncol(m)){
                m[i,j] <- sum(unique(data_files_list[[i]]) %in% allgenes)
        }
}
##NOTE! parameters with 1 in front (m1, n1 and k1) are for the calculations considering the universe as all genes of the organism, not only the annotated genes
# m1 is the number of sucesses in the population: the number of genes of the first gene set that is in the genome
#m1 <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
#dimnames(m1) <- list(names(data_files_list),names(data_files_list))
#library(biomaRt)
#ensembl = useMart(biomart ="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="may2017.archive.ensembl.org")
# d <- getBM(attributes = "entrezgene",values = "*",mart = ensembl)
# allgenes1 <- unique(d[,1]) 
# for (i in 1:nrow(m1)){
        #for (j in 1:ncol(m1)){
                #m1[i,j] <- sum(unique(data_files_list[[i]]) %in% allgenes1)
        #}
#}
# n is the number of insucesses in the population = size of population - number of sucesses in the population: genome size - m
n <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
dimnames(n) <- list(names(data_files_list),names(data_files_list))
for (i in 1:nrow(n)){
        for (j in 1:ncol(n)){
                n[i,j] <- length(allgenes) - m[[i]] # counting the number of common genes  
        }
}
# n1 is the number of insucesses in the population = size of population - number of sucesses in the population: genome size - m
#n1 <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
#dimnames(n1) <- list(names(data_files_list),names(data_files_list))
#for (i in 1:nrow(n1)){
        #for (j in 1:ncol(n1)){
                #n1[i,j] <- length(allgenes1) - m1[[i]] # counting the number of common genes  
        #}
#}

# k is the size of the sample: size of the 2nd gene set (the one in the column)
k <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
dimnames(k) <- list(names(data_files_list),names(data_files_list))
for (i in 1:nrow(k)){
        for (j in 1:ncol(k)){
                k[i,j] <- sum(unique(data_files_list[[j]]) %in% allgenes) # counting the number of common genes  
        }
}
# k1 is the size of the sample: size of the 2nd gene set (the one in the column)
#k1 <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
#dimnames(k1) <- list(names(data_files_list),names(data_files_list))
#for (i in 1:nrow(k1)){
        #for (j in 1:ncol(k1)){
                #k1[i,j] <- sum(unique(data_files_list[[j]]) %in% allgenes1) # counting the number of common genes  
        #}
#}
# hypergeometric test uses the above parameters: q,m,n,k
# we do 1-hyper(...) because we want the probability of getting that number or higher of sucesses and not that number or lower of sucesses. Instead, we could have used: phyper(...,lower.tail=FALSE)
# hyper computes P[X>value] so if we want to include the value, P[X>=value], we need to substract one from number of sucesses in the sample: phyper(q-1) 
# log.p=T gives logaritmic value as output
hyp.test <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
dimnames(hyp.test) <- list(names(data_files_list),names(data_files_list))
for (i in 1:nrow(hyp.test)){
        for (j in 1:ncol(hyp.test)){
                hyp.test[i,j] <- phyper(q[i,j]-1,m[i,j],n[i,j],k[i,j], lower.tail = F)
        }
}
#hyp.test1 <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
#dimnames(hyp.test1) <- list(names(data_files_list),names(data_files_list))
#for (i in 1:nrow(hyp.test1)){
        #for (j in 1:ncol(hyp.test1)){
                #hyp.test1[i,j] <- phyper(q[i,j]-1,m1[i,j],n1[i,j],k1[i,j], lower.tail = F)
        #}
#}
# Bonferroni correction to ajust the p-value for multiple testing
hyp.test.adj <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
dimnames(hyp.test.adj) <- list(names(data_files_list),names(data_files_list))
for (i in 1: nrow(hyp.test.adj)){
        for (j in 1:ncol(hyp.test.adj)){
                hyp.test.adj[i,j] <- p.adjust(hyp.test[i,j], method = "bonferroni")
        }
}

#hyp.test.adj1 <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
#dimnames(hyp.test.adj1) <- list(names(data_files_list),names(data_files_list))
#for (i in 1: nrow(hyp.test.adj1)){
        #for (j in 1:ncol(hyp.test.adj1)){
                #hyp.test.adj1[i,j] <- p.adjust(hyp.test1[i,j], method = "bonferroni")
        #}
#}
# Some adjusted p-values get above one, but probabilities can't be above one. So some people force the values to be one...
hyp.test.adj[hyp.test.adj > 1] <- 1
#hyp.test.adj1[hyp.test.adj1 > 1] <- 1
# if we wanted to use logarythm:
# before doing the log we have to replace the 0's by 10E-16 because otherwise we would get -Inf values after logaritmizing and the heatmap function wouldn't be able to plot -inf values (there would be an error).
hyp.test.adj[hyp.test.adj ==0] <- 10E-16 # we could had have also replaced by the smallest value in the matrix (smallest p-value)
#hyp.test.adj1[hyp.test.adj1 ==0] <- 10E-16
log.hyp.test.adj <- log10(hyp.test.adj)
#log.hyp.test.adj1 <- log10(hyp.test.adj1)
# p-values lower than 10E-16  are set to 10E-16:The reasons are the following:
    # from the statistical point of view, there is not much difference between a p-value of  10E-16 and 10E-72 (both are highly significant) , but much more between 10E-1 and 10E-5.
    # However, since we heatmap.2 uses in its default version a eucleadian clustering distance, the 10E-16 and 10E-72 is much more different, than 10E-1 and 10E-5.
    # I noticed in the p-value matrix are very small values, so the smallest will determine the clustering.
    # To prevent this, we can set a min. threshold for significance. 
log.hyp.test.adj.thres <- log.hyp.test.adj
log.hyp.test.adj.thres[log.hyp.test.adj.thres < -16] <- -16
library(gplots)
# !note:
# to assign colors by their name we could use library(RColorBrewer)
# and the argument (n = 299) lets us define how many individuals colors we want to have in our palette. the higher the number of individual colors, the smoother the transition will be
# overview of the different color names in R can be found at http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
# the function would be:
my_palette <- colorRampPalette(c("red","beige","cornflowerblue"))(n = 1499)
# but instead we chose to use colors from:
# load("D:/Tania_2/Internship/R/Packages_extras/rwb.RData")
# my_palette <- rwb[20:50]
# sometimes we want to have a little skewed color range depending on the data we are analyzing
# Let's assume that our example data set ranges from -1 to 1, and we are particularly interested in samples with relatively high values (0.8 to 1.0)
# We want to highlight these samples in our heat map by only showing values from 0.8 to 1 in green
# In this case:
# col_breaks = c(seq(-1,0,length=100), # for green
#                seq(0,0.8,length=100), # for yellow
#                seq(0.8,1,length=100)) # for red
# NOTE: number of breaks should be number of colors + 1
#       see: https://stackoverflow.com/questions/10749367/heat-map-adjusting-color-range
# here:
breaks = c(seq(log10(10E-16),log10(0.049),length=1000),seq(log10(0.05),log10(0.25),length=200),seq(log10(0.251),log10(1),length=300)) # number of breaks should be number of colors +1 (length(breaks) = length(my_palette) + 1)
# position of each element in the heatmap.2 plot can be controlled using the lmat, lhei and lwid parameters:
# to understand this go to http://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
lmat = rbind(c(0,4,0),c(0,3,0),c(2,1,0))               
lwid = c(0.5,4,0.5)
lhei = c(0.9,0.5,4)
##########
# density.info="none",  # turns off density plot inside color legend
# trace="none",         # turns off trace lines inside the heat map
# margins =c(12,9),     # numeric vector of length 2 containing the margins for column and row names, respectively.
# breaks=...,           # enable color transition at specified limits
# cexRow or cezCol      # changes font size of the row or column labels: cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc) 
# key.xlab=""           # removes the key label
# key.title = F         # removes key title
# key.par=list(mar=c(bottom, left, top, right)) # controls the position of color key within the row defined by lmat
# revC=T # does the simetry for the columns (puts the 1st column as the last and vice-versa and for all other columns) maintining the order of the clustering. Important for the order of the gene set names from top to bottom to be the same as from left to right
# cellnote # given a matrix writes matrix values inside heatmap squares
# notecol # specifies color of cellnote
# notecex # numeric scaling factor for cellnote items

tiff(file="D:/Tania_2/Internship/R/Results/clustering.genesets/by_log.p-value/heatmapgenecluster_Hs_stemness_up.tif", width = 8*600, height = 9*600, res = 600)
heatmap.2(log.hyp.test.adj.thres, cellnote = q, notecol = 'black', notecex = 0.85 , scale="none",col=my_palette,breaks=breaks,density.info="none",trace="none", margins =c(17,11),  lmat = lmat, lwid = lwid, lhei = lhei, key.xlab="", key = T, key.title=F, symkey=F, keysize=10, key.par=list(mar=c(4,1,4,10)), sepcolor="black", colsep=1:ncol(log.hyp.test.adj.thres), rowsep=1:nrow(log.hyp.test.adj.thres), sepwidth=c(0.01,0.01),revC=T, cexRow=1.7, cexCol = 1.7)
dev.off()

png(file="D:/Tania_2/Internship/R/Results/clustering.genesets/by_log.p-value/heatmapgenecluster_Hs_stemness_up.png", width = 8*600, height = 9*600, res = 600)
heatmap.2(log.hyp.test.adj.thres,cellnote = q, notecol = 'black', notecex = 0.85, scale="none",col=my_palette,breaks=breaks,density.info="none",trace="none", margins =c(17,11),  lmat = lmat, lwid = lwid, lhei = lhei, key.xlab="", key = T, key.title=F, symkey=F, keysize=10, key.par=list(mar=c(4,1,4,10)), sepcolor="black", colsep=1:ncol(log.hyp.test.adj.thres), rowsep=1:nrow(log.hyp.test.adj.thres), sepwidth=c(0.01,0.01),revC=T, cexRow=1.7, cexCol = 1.7)
dev.off()

# get the number of signifficant comparisons:
signif <- log.hyp.test.adj.thres < log10(0.05)
signif2 <- sum(signif, na.rm=TRUE)
signif2 # 254 comparisons are signifficnat
#####################################################
#### Mm stemness gene sets
files_list1 <- list.files("Mm_up_without_TF", full.names=TRUE)   #list files in the folder
data_files_list<-lapply(files_list1, function(x) read.csv(x, header=F, stringsAsFactors=F)[,1]) #create a list with the data from all files
# data_files_list<-lapply(files_list1, read.table)
datasetnames<-vector()
for (i in 1:length(files_list1)){
        datasetnames[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list) <- datasetnames   #give names to each list element

q <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
dimnames(q) <- list(names(data_files_list),names(data_files_list))

for (i in 1:nrow(q)){
        for (j in 1:ncol(q)){
                q[i,j] <- sum(unique(data_files_list[[i]]) %in% unique(data_files_list[[j]])) # counting the number of common genes  
        }
}

library(org.Mm.eg.db)
m <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
dimnames(m) <- list(names(data_files_list),names(data_files_list))
allgenes <- unique(mappedkeys(org.Mm.egGO))
for (i in 1:nrow(m)){
        for (j in 1:ncol(m)){
                m[i,j] <- length(unique(data_files_list[[i]])%in%allgenes)
        }
}

n <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
dimnames(n) <- list(names(data_files_list),names(data_files_list))
for (i in 1:nrow(n)){
        for (j in 1:ncol(n)){
                n[i,j] <- length(allgenes) - m[[i]] # counting the number of common genes  
        }
}

k <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
dimnames(k) <- list(names(data_files_list),names(data_files_list))
for (i in 1:nrow(k)){
        for (j in 1:ncol(k)){
                k[i,j] <- sum(unique(data_files_list[[j]]) %in% allgenes) # counting the number of common genes  
        }
}

hyp.test <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
dimnames(hyp.test) <- list(names(data_files_list),names(data_files_list))
for (i in 1:nrow(hyp.test)){
        for (j in 1:ncol(hyp.test)){
                hyp.test[i,j] <- phyper(q[i,j]-1,m[i,j],n[i,j],k[i,j], lower.tail = F)
        }
}

hyp.test.adj <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
dimnames(hyp.test.adj) <- list(names(data_files_list),names(data_files_list))
for (i in 1: nrow(hyp.test.adj)){
        for (j in 1:ncol(hyp.test.adj)){
                hyp.test.adj[i,j] <- p.adjust(hyp.test[i,j], method = "bonferroni")
        }
}

hyp.test.adj[hyp.test.adj > 1] <- 1
hyp.test.adj[hyp.test.adj ==0] <- 10E-16
log.hyp.test.adj <- log10(hyp.test.adj)
log.hyp.test.adj.thres <- log.hyp.test.adj
log.hyp.test.adj.thres[log.hyp.test.adj.thres < -16] <- -16

library(gplots)
my_palette <- colorRampPalette(c("red","beige","cornflowerblue"))(n = 1499)
breaks = c(seq(log10(10E-16),log10(0.049),length=1000),seq(log10(0.05),log10(0.25),length=200),seq(log10(0.251),log10(1),length=300)) # number of breaks should be number of colors +1 (length(breaks) = length(my_palette) + 1)
lmat = rbind(c(0,4,0),c(0,3,0),c(2,1,0))               
lwid = c(0.5,4,0.5)
lhei = c(0.9,0.5,4)

tiff(file="D:/Tania_2/Internship/R/Results/clustering.genesets/by_log.p-value/heatmapgenecluster_Mm_stemness_up.tif", width = 8*600, height = 9*600, res = 600)
heatmap.2(log.hyp.test.adj.thres, cellnote = q, notecol = 'black', notecex = 0.85, scale="none", col=my_palette, breaks=breaks, density.info="none", trace="none", margins =c(17,11),  lmat = lmat, lwid = lwid, lhei = lhei, key.xlab="", key = T, key.title=F, symkey=F, keysize=10, key.par=list(mar=c(4,1,4,10)), sepcolor="black", colsep=1:ncol(log.hyp.test.adj.thres), rowsep=1:nrow(log.hyp.test.adj.thres), sepwidth=c(0.01,0.01),revC=T, cexRow=1.7, cexCol = 1.7)
dev.off()

png(file="D:/Tania_2/Internship/R/Results/clustering.genesets/by_log.p-value/heatmapgenecluster_Mm_stemness_up.png", width = 8*600, height = 9*600, res = 600)
heatmap.2(log.hyp.test.adj.thres, cellnote = q, notecol = 'black', notecex = 0.85, scale="none", col=my_palette, breaks=breaks, density.info="none", trace="none", margins =c(17,11),  lmat = lmat, lwid = lwid, lhei = lhei, key.xlab="", key = T, key.title=F, symkey=F, keysize=10, key.par=list(mar=c(4,1,4,10)), sepcolor="black", colsep=1:ncol(log.hyp.test.adj.thres), rowsep=1:nrow(log.hyp.test.adj.thres), sepwidth=c(0.01,0.01),revC=T, cexRow=1.7, cexCol = 1.7)
dev.off()

#svg(file="C:/Users/Tania_2/Internship/R/Results/clustering.genesets/by_log.p-value/heatmapgenecluster_Mm_stemness_up.svg")#export file to png
#heatmap.2(log.hyp.test.adj,scale="none",notecol="black",cellnote=round(log.hyp.test.adj, digits = 2),notecex=0.7,col=my_palette,breaks=breaks,density.info="none",trace="none", margins =c(13,6),  lmat = lmat, lwid = lwid, lhei = lhei, key.xlab="", key = T, key.title=F)
#dev.off()

#pdf(file="C:/Users/Tania_2/Internship/R/Results/clustering.genesets/log.-value/heatmapgenecluster_Mm_stemness_up.pdf")#export file to png
#heatmap.2(log.hyp.test.adj,scale="none",notecol="black",cellnote=round(log.hyp.test.adj, digits = 2),notecex=0.7,col=my_palette,breaks=breaks,density.info="none",trace="none", margins =c(13,6),  lmat = lmat, lwid = lwid, lhei = lhei, key.xlab="", key = T, key.title=F)
#dev.off()

# get the number of signifficant comparisons:
signif <- log.hyp.test.adj.thres < log10(0.05)
signif2 <- sum(signif, na.rm=TRUE)
signif2 # 317 comparisons are signifficnat


# #####################################################
# #### Hs TF gene sets
# files_list1 <- list.files("Hs_up_TF", full.names=TRUE)   #list files in the folder
# data_files_list <- lapply(files_list1, function(x) as.character(unique(read.csv(x, header=F, stringsAsFactors=F)[,1]))) #create a list with the data from all files
# datasetnames<-vector()
# for (i in 1:length(files_list1)){
#         datasetnames[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
# }
# names(data_files_list) <- datasetnames   #give names to each list element
# 
# q <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
# dimnames(q) <- list(names(data_files_list),names(data_files_list))
# 
# for (i in 1:nrow(q)){
#         for (j in 1:ncol(q)){
#                 q[i,j] <- sum(unique(data_files_list[[i]]) %in% unique(data_files_list[[j]])) # counting the number of common genes  
#         }
# }
# 
# library(org.Hs.eg.db)
# m <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
# dimnames(m) <- list(names(data_files_list),names(data_files_list))
# allgenes <- unique(mappedkeys(org.Hs.egGO))
# for (i in 1:nrow(m)){
#         for (j in 1:ncol(m)){
#                 m[i,j] <- sum(unique(data_files_list[[i]]) %in% allgenes)
#         }
# }
# 
# n <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
# dimnames(n) <- list(names(data_files_list),names(data_files_list))
# for (i in 1:nrow(n)){
#         for (j in 1:ncol(n)){
#                 n[i,j] <- length(allgenes) - m[[i]] # counting the number of common genes  
#         }
# }
# 
# k <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
# dimnames(k) <- list(names(data_files_list),names(data_files_list))
# for (i in 1:nrow(k)){
#         for (j in 1:ncol(k)){
#                 k[i,j] <- sum(unique(data_files_list[[j]]) %in% allgenes) # counting the number of common genes  
#         }
# }
# 
# hyp.test <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
# dimnames(hyp.test) <- list(names(data_files_list),names(data_files_list))
# for (i in 1:nrow(hyp.test)){
#         for (j in 1:ncol(hyp.test)){
#                 hyp.test[i,j] <- phyper(q[i,j]-1,m[i,j],n[i,j],k[i,j], lower.tail = F)
#         }
# }
# 
# hyp.test.adj <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
# dimnames(hyp.test.adj) <- list(names(data_files_list),names(data_files_list))
# for (i in 1: nrow(hyp.test.adj)){
#         for (j in 1:ncol(hyp.test.adj)){
#                 hyp.test.adj[i,j] <- p.adjust(hyp.test[i,j], method = "bonferroni")
#         }
# }
# 
# hyp.test.adj[hyp.test.adj > 1] <- 1
# hyp.test.adj[hyp.test.adj ==0] <- 10E-16
# log.hyp.test.adj <- log10(hyp.test.adj)
# log.hyp.test.adj.thres <- log.hyp.test.adj
# log.hyp.test.adj.thres[log.hyp.test.adj.thres < -16] <- -16
# 
# library(gplots)
# my_palette <- colorRampPalette(c("red","beige","cornflowerblue"))(n = 1499)
# breaks = c(seq(log10(10E-16),log10(0.049),length=1000),seq(log10(0.05),log10(0.25),length=200),seq(log10(0.251),log10(1),length=300)) # number of breaks should be number of colors +1 (length(breaks) = length(my_palette) + 1)
# lmat = rbind(c(0,4,0),c(0,3,0),c(2,1,0))               
# lwid = c(0.5,4,0.5)
# lhei = c(0.9,0.5,4)
# 
# png(file="C:/Users/Tania_2/Internship/R/Results/clustering.genesets/by_log.p-value/heatmapgenecluster_Hs_TF_up.png", width = 8*300, height = 9*300, res = 300)#export file to png
# heatmap.2(log.hyp.test.adj.thres,scale="none",col=my_palette,breaks=breaks,density.info="none",trace="none", margins =c(13,10),  lmat = lmat, lwid = lwid, lhei = lhei, key.xlab="", key = T, key.title=F, symkey=F, key.par=list(mar=c(4,1,4,10)), sepcolor="black", colsep=1:ncol(log.hyp.test.adj.thres), rowsep=1:nrow(log.hyp.test.adj.thres), sepwidth=c(0.01,0.01),revC=T)
# dev.off()
# 
# #svg(file="C:/Users/Tania_2/Internship/R/Results/clustering.genesets/by_log.p-value/heatmapgenecluster_Hs_TF_up.svg")#export file to png
# #heatmap.2(log.hyp.test.adj,scale="none",notecol="black",cellnote=round(log.hyp.test.adj, digits = 2),notecex=0.7,col=my_palette,breaks=breaks1,density.info="none",trace="none", margins =c(13,6),  lmat = lmat, lwid = lwid, lhei = lhei, key.xlab="", key = T, key.title=F)
# #dev.off()
# 
# #pdf(file="C:/Users/Tania_2/Internship/R/Results/clustering.genesets/by_log.p-value/heatmapgenecluster_Hs_TF_up.pdf")#export file to png
# #heatmap.2(log.hyp.test.adj,scale="none",notecol="black",cellnote=round(log.hyp.test.adj, digits = 2),notecex=0.7,col=my_palette,breaks=breaks1,density.info="none",trace="none", margins =c(13,6),  lmat = lmat, lwid = lwid, lhei = lhei, key.xlab="", key = T, key.title=F)
# #dev.off()

# #####################################################
# #### Mm TF gene sets
# files_list1 <- list.files("Mm_up_TF", full.names=TRUE)   #list files in the folder
# data_files_list<-lapply(files_list1, function(x) read.csv(x, header=F, stringsAsFactors=F)[,1]) #create a list with the data from all files
# # data_files_list<-lapply(files_list1, read.table)
# datasetnames<-vector()
# for (i in 1:length(files_list1)){
#         datasetnames[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
# }
# names(data_files_list) <- datasetnames   #give names to each list element
# 
# q <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
# dimnames(q) <- list(names(data_files_list),names(data_files_list))
# 
# for (i in 1:nrow(q)){
#         for (j in 1:ncol(q)){
#                 q[i,j] <- sum(unique(data_files_list[[i]]) %in% unique(data_files_list[[j]])) # counting the number of common genes  
#         }
# }
# 
# library(org.Mm.eg.db)
# m <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
# dimnames(m) <- list(names(data_files_list),names(data_files_list))
# allgenes <- unique(mappedkeys(org.Mm.egGO))
# for (i in 1:nrow(m)){
#         for (j in 1:ncol(m)){
#                 m[i,j] <- length(unique(data_files_list[[i]])%in%allgenes)
#         }
# }
# 
# n <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
# dimnames(n) <- list(names(data_files_list),names(data_files_list))
# for (i in 1:nrow(n)){
#         for (j in 1:ncol(n)){
#                 n[i,j] <- length(allgenes) - m[[i]] # counting the number of common genes  
#         }
# }
# 
# k <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
# dimnames(k) <- list(names(data_files_list),names(data_files_list))
# for (i in 1:nrow(k)){
#         for (j in 1:ncol(k)){
#                 k[i,j] <- sum(unique(data_files_list[[j]]) %in% allgenes) # counting the number of common genes  
#         }
# }
# 
# hyp.test <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
# dimnames(hyp.test) <- list(names(data_files_list),names(data_files_list))
# for (i in 1:nrow(hyp.test)){
#         for (j in 1:ncol(hyp.test)){
#                 hyp.test[i,j] <- phyper(q[i,j]-1,m[i,j],n[i,j],k[i,j], lower.tail = F)
#         }
# }
# 
# hyp.test.adj <- matrix(0,nrow=length(data_files_list),ncol=length(data_files_list))
# dimnames(hyp.test.adj) <- list(names(data_files_list),names(data_files_list))
# for (i in 1:nrow(hyp.test.adj)){
#         for (j in 1:ncol(hyp.test.adj)){
#                 hyp.test.adj[i,j] <- p.adjust(hyp.test[i,j], method = "bonferroni")
#         }
# }
# 
# hyp.test.adj[hyp.test.adj > 1] <- 1
# hyp.test.adj[hyp.test.adj ==0] <- 10E-16
# log.hyp.test.adj <- log10(hyp.test.adj)
# log.hyp.test.adj.thres <- log.hyp.test.adj
# log.hyp.test.adj.thres[log.hyp.test.adj.thres < -16] <- -16
# 
# library(gplots)
# my_palette <- colorRampPalette(c("red","beige","cornflowerblue"))(n = 1499)
# breaks = c(seq(log10(10E-16),log10(0.049),length=1000),seq(log10(0.05),log10(0.25),length=200),seq(log10(0.251),log10(1),length=300)) # number of breaks should be number of colors +1 (length(breaks) = length(my_palette) + 1)
# lmat = rbind(c(0,4,0),c(0,3,0),c(2,1,0))               
# lwid = c(0.5,4,0.5)
# lhei = c(0.9,0.5,4)
# 
# png(file="C:/Users/Tania_2/Internship/R/Results/clustering.genesets/by_log.p-value/heatmapgenecluster_Mm_TF_up.png", width = 8*300, height = 9*300, res = 300)#export file to png
# heatmap.2(log.hyp.test.adj.thres,scale="none",col=my_palette,breaks=breaks,density.info="none",trace="none", margins =c(13,10),  lmat = lmat, lwid = lwid, lhei = lhei, key.xlab="", key = T, key.title=F, symkey=F, key.par=list(mar=c(4,1,4,10)), sepcolor="black", colsep=1:ncol(log.hyp.test.adj.thres), rowsep=1:nrow(log.hyp.test.adj.thres), sepwidth=c(0.01,0.01), cexRow=0.6, cexCol=0.6,revC=T)
# dev.off()
# 
# #svg(file="C:/Users/Tania_2/Internship/R/Results/clustering.genesets/by_log.p-value/heatmapgenecluster_Mm_TF_up.svg")#export file to png
# #heatmap.2(log.hyp.test.adj,scale="none",notecol="black",cellnote=round(log.hyp.test.adj, digits = 2),notecex=0.7,col=my_palette,breaks=breaks,density.info="none",trace="none", margins =c(13,6),  lmat = lmat, lwid = lwid, lhei = lhei, key.xlab="", key = T, key.title=F)
# #dev.off()
# 
# #pdf(file="C:/Users/Tania_2/Internship/R/Results/clustering.genesets/by_log.p-value/heatmapgenecluster_Mm_TF_up.pdf")#export file to png
# #heatmap.2(log.hyp.test.adj,scale="none",notecol="black",cellnote=round(log.hyp.test.adj, digits = 2),notecex=0.7,col=my_palette,breaks=breaks,density.info="none",trace="none", margins =c(13,6),  lmat = lmat, lwid = lwid, lhei = lhei, key.xlab="", key = T, key.title=F)
# #dev.off()

