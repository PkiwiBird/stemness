setwd("C:/Users/Tania_2/Internship/R/Data")
homo <- read.table("/Users/Tania_2/Internship/R/Data/Homo_sapiens.gene_info.txt", header = T,  stringsAsFactors = F, sep="\t", fill=T, quote="")
mouse <- read.table("/Users/Tania_2/Internship/R/Data/Mus_musculus.gene_info.txt", header = T,  stringsAsFactors = F, sep="\t", fill=T, quote="")

###################################################################################################################
###### Human stemness signatures
files_list <- list.files("Hs_up_without_TF", full.names=TRUE) 
datasetnames<-vector()
for (i in 1:length(files_list)){
        datasetnames[i]<-unlist(strsplit(files_list[i],'[/.]'))[2]
}
genes<-character()
list.datasets <- list()
for (i in 1:length(files_list)){
        tmp <-  as.character(read.csv(files_list[i], fileEncoding="latin1", header=F)[,1])
        for (j in 1:length(tmp)) {tmp[j] <- gsub("[^0-9]","",tmp[j]) # when writting the entrez id names in the files of each gene set, same weird spaces were introduced that were affecting the ranking of genes
                                                                     # so this step is to get rid of those spaces. [^0-9] regular expression means find any character that is not a digit
        }
        list.datasets[[i]] <- unique(tmp) # to guarantee that each set doesn't have repeated genes that could be counted twice in the ranking
        genes<-unique(c(genes,tmp)) # to get a vector with all genes of all datasets
}
B<-matrix(0,length(genes),length(datasetnames))
for (i in 1:length(genes)){
        for (j in 1:length(datasetnames)){
                B[i,j] <- sum(list.datasets[[j]]%in%genes[i])   
        }
}
dimnames(B)<-list(genes,datasetnames)
score<-vector()
for (i in 1:length(genes)){
        score[i]<-sum(B[i,])
}
M<-data.frame(genes,score)
ord<-order(M[,2],decreasing=TRUE) # vector with indices starting by the indice of the highest values and then the second highest value and so on
m<-M[ord,]
# arrange m to get 3 columns (1st=geneid,2nd=genesymbol,3rd=score)
m1<-data.frame(m[,1],character(length=nrow(m)),m[,2]) # character(length=nrow(m)) opens space for a new column
colnames(m1)<-c("gene_id","gene_symbol","score")
# change 1st and 2nd columns of m1 from factor class to charcater class:
m1[,1] <- as.character(m1[,1])
m1[,2] <- as.character(m1[,2])
# if there's gene symbols in the entrez2sym.h for genes in our ranked list we get those
# we don't get the gene symbol for some genes of our ranked list, so we have to get those latter by hand
# library to retrieve object org.Hs.eg.db and associated functions:
library(org.Hs.eg.db)
# get list in which each element is a symbol for human gene and is named with corresponding entrez gene id:
x <- org.Hs.egSYMBOL
mapped_genes<- mappedkeys(x)
entrez2sym.h <- as.list(x[mapped_genes])

for (i in 1:nrow(m1)){
        if(m1[i,1]%in%names(entrez2sym.h)){
                m1[i,2]<-entrez2sym.h[[m1[i,1]]]
        } else if (m1[i,1]%in%homo$GeneID) {
                m1[i,2]<-subset(homo$Symbol,homo$GeneID==m1[i,1])
        }
}

write.csv(m1,"C:/Users/Tania_2/Internship/R/Results/ranked.lists/ranked.list.human.stem.signature.csv",row.names=F)
C<-data.frame(B,score)
ord1<-order(C[,ncol(C)],decreasing=TRUE)
c<-C[ord1,]
c1 <- data.frame(gene_symbol=m1[,2],gene_id=m1[,1],c[,], stringsAsFactors = F)
write.csv(c1,"C:/Users/Tania_2/Internship/R/Results/ranked.lists/ranked.list.human.stem.signature.genes.csv",row.names=F)

#########################################################################
######### Mouse stemness signature
files_list <- list.files("Mm_up_without_TF", full.names=TRUE) 
datasetnames<-vector()
for (i in 1:length(files_list)){
        datasetnames[i]<-unlist(strsplit(files_list[i],'[/.]'))[2]
}
genes<-character()
list.datasets <- list()
for (i in 1:length(files_list)){
        tmp <-  as.character(read.csv(files_list[i], fileEncoding="latin1", header=F)[,1])
        for (j in 1:length(tmp)) {tmp[j] <- gsub("[^0-9]","",tmp[j]) 
        }
        list.datasets[[i]] <- unique(tmp)
        genes<-unique(c(genes,tmp))
}
B<-matrix(0,length(genes),length(datasetnames))
for (i in 1:length(genes)){
        for (j in 1:length(datasetnames)){
                B[i,j] <- sum(list.datasets[[j]]%in%genes[i])   
        }
}
dimnames(B)<-list(genes,datasetnames)
score<-vector()
for (i in 1:length(genes)){
        score[i]<-sum(B[i,])
}
M<-data.frame(genes,score)
ord<-order(M[,2],decreasing=TRUE)
m<-M[ord,]
# library to retrieve object org.Mm.eg.db and associated functions:
library(org.Mm.eg.db)
# arrange m to get 3 rows (1st=geneid,2nd=genesymbol,3rd=score)
m1<-data.frame(m[,1],character(length=nrow(m)),m[,2])
colnames(m1)<-c("gene_id","gene_symbol","score")
# change 1st and 2nd columns of m1 from factor class to charcater class:
m1[,1] <- as.character(m1[,1])
m1[,2] <- as.character(m1[,2])
# get list in which each element is a symbol for mouse gene and is named with corresponding entrez gene id:
x <- org.Mm.egSYMBOL
mapped_genes<- mappedkeys(x)
entrez2sym.m <- as.list(x[mapped_genes])
# if there's gene symbols in the entrez2sym.m for genes in our ranked list we get those
# we don't get the gene symbol for some genes of our ranked list, so we have to get those latter by hand
for (i in 1:nrow(m1)){
        if(m1[i,1]%in%names(entrez2sym.m)){
                m1[i,2]<-entrez2sym.m[[m1[i,1]]]
        } else if (m1[i,1]%in%mouse$GeneID) {
                m1[i,2]<-subset(mouse$Symbol,mouse$GeneID==m1[i,1])
        }
}
write.csv(m1,"C:/Users/Tania_2/Internship/R/Results/ranked.lists/ranked.list.mouse.stem.signature.csv",row.names=F)
C<-data.frame(B,score)
ord1<-order(C[,ncol(C)],decreasing=TRUE)
c<-C[ord1,]
c1 <- data.frame(gene_symbol=m1[,2],gene_id=m1[,1],c[,], stringsAsFactors = F)
write.csv(c1,"C:/Users/Tania_2/Internship/R/Results/ranked.lists/ranked.list.mouse.stem.signature.genes.csv",row.names=F)

###################################################################################################################
###### Human ESC signatures
files_list <- list.files("Hs_ESC", full.names=TRUE) 
datasetnames<-vector()
for (i in 1:length(files_list)){
        datasetnames[i]<-unlist(strsplit(files_list[i],'[/.]'))[2]
}
genes<-character()
list.datasets <- list()
for (i in 1:length(files_list)){
        tmp <-  as.character(read.csv(files_list[i], fileEncoding="latin1", header=F)[,1])
        for (j in 1:length(tmp)) {tmp[j] <- gsub("[^0-9]","",tmp[j]) 
        }
        list.datasets[[i]] <- unique(tmp)
        genes<-unique(c(genes,tmp))
}
B<-matrix(0,length(genes),length(datasetnames))
for (i in 1:length(genes)){
        for (j in 1:length(datasetnames)){
                B[i,j] <- sum(list.datasets[[j]]%in%genes[i])   
        }
}
dimnames(B)<-list(genes,datasetnames)
score<-vector()
for (i in 1:length(genes)){
        score[i]<-sum(B[i,])
}
M<-data.frame(genes,score)
ord<-order(M[,2],decreasing=TRUE)
m<-M[ord,]
# arrange m to get 3 rows (1st=geneid,2nd=genesymbol,3rd=score)
m1<-data.frame(m[,1],character(length=nrow(m)),m[,2])
colnames(m1)<-c("gene_id","gene_symbol","score")
# change 1st and 2nd columns of m1 from factor class to charcater class:
m1[,1] <- as.character(m1[,1])
m1[,2] <- as.character(m1[,2])
# if there's gene symbols in the entrez2sym.h for genes in our ranked list we get those
# we don't get the gene symbol for some genes of our ranked list, so we have to get those latter by hand
# library to retrieve object org.Hs.eg.db and associated functions:
library(org.Hs.eg.db)
# get list in which each element is a symbol for human gene and is named with corresponding entrez gene id:
x <- org.Hs.egSYMBOL
mapped_genes<- mappedkeys(x)
entrez2sym.h <- as.list(x[mapped_genes])

for (i in 1:nrow(m1)){
        if(m1[i,1]%in%names(entrez2sym.h)){
                m1[i,2]<-entrez2sym.h[[m1[i,1]]]
        } else if (m1[i,1]%in%homo$GeneID) {
                m1[i,2]<-subset(homo$Symbol,homo$GeneID==m1[i,1])
        }
}

write.csv(m1,"C:/Users/Tania_2/Internship/R/Results/ranked.lists/ranked.list.human.ESC.signature.csv",row.names=F)
C<-data.frame(B,score)
ord1<-order(C[,ncol(C)],decreasing=TRUE)
c<-C[ord1,]
c1 <- data.frame(gene_symbol=m1[,2],gene_id=m1[,1],c[,], stringsAsFactors = F)
write.csv(c1,"C:/Users/Tania_2/Internship/R/Results/ranked.lists/ranked.list.human.ESC.signature.genes.csv",row.names=F)

###################################################################################################################
###### Human HSC signatures
files_list <- list.files("Hs_HSC", full.names=TRUE) 
datasetnames<-vector()
for (i in 1:length(files_list)){
        datasetnames[i]<-unlist(strsplit(files_list[i],'[/.]'))[2]
}
genes<-character()
list.datasets <- list()
for (i in 1:length(files_list)){
        tmp <-  as.character(read.csv(files_list[i], fileEncoding="latin1", header=F)[,1])
        for (j in 1:length(tmp)) {tmp[j] <- gsub("[^0-9]","",tmp[j]) 
        }
        list.datasets[[i]] <- unique(tmp)
        genes<-unique(c(genes,tmp))
}
B<-matrix(0,length(genes),length(datasetnames))
for (i in 1:length(genes)){
        for (j in 1:length(datasetnames)){
                B[i,j] <- sum(list.datasets[[j]]%in%genes[i])   
        }
}
dimnames(B)<-list(genes,datasetnames)
score<-vector()
for (i in 1:length(genes)){
        score[i]<-sum(B[i,])
}
M<-data.frame(genes,score)
ord<-order(M[,2],decreasing=TRUE)
m<-M[ord,]
# arrange m to get 3 rows (1st=geneid,2nd=genesymbol,3rd=score)
m1<-data.frame(m[,1],character(length=nrow(m)),m[,2])
colnames(m1)<-c("gene_id","gene_symbol","score")
# change 1st and 2nd columns of m1 from factor class to charcater class:
m1[,1] <- as.character(m1[,1])
m1[,2] <- as.character(m1[,2])
# if there's gene symbols in the entrez2sym.h for genes in our ranked list we get those
# we don't get the gene symbol for some genes of our ranked list, so we have to get those latter by hand
# library to retrieve object org.Hs.eg.db and associated functions:
library(org.Hs.eg.db)
# get list in which each element is a symbol for human gene and is named with corresponding entrez gene id:
x <- org.Hs.egSYMBOL
mapped_genes<- mappedkeys(x)
entrez2sym.h <- as.list(x[mapped_genes])

for (i in 1:nrow(m1)){
        if(m1[i,1]%in%names(entrez2sym.h)){
                m1[i,2]<-entrez2sym.h[[m1[i,1]]]
        } else if (m1[i,1]%in%homo$GeneID) {
                m1[i,2]<-subset(homo$Symbol,homo$GeneID==m1[i,1])
        }
}

write.csv(m1,"C:/Users/Tania_2/Internship/R/Results/ranked.lists/ranked.list.human.HSC.signature.csv",row.names=F)
C<-data.frame(B,score)
ord1<-order(C[,ncol(C)],decreasing=TRUE)
c<-C[ord1,]
c1 <- data.frame(gene_symbol=m1[,2],gene_id=m1[,1],c[,], stringsAsFactors = F)
write.csv(c1,"C:/Users/Tania_2/Internship/R/Results/ranked.lists/ranked.list.human.HSC.signature.genes.csv",row.names=F)

###################################################################################################################
###### Human NSC signatures
files_list <- list.files("Hs_NSC", full.names=TRUE) 
datasetnames<-vector()
for (i in 1:length(files_list)){
        datasetnames[i]<-unlist(strsplit(files_list[i],'[/.]'))[2]
}
genes<-character()
list.datasets <- list()
for (i in 1:length(files_list)){
        tmp <-  as.character(read.csv(files_list[i], fileEncoding="latin1", header=F)[,1])
        for (j in 1:length(tmp)) {tmp[j] <- gsub("[^0-9]","",tmp[j]) 
        }
        list.datasets[[i]] <- unique(tmp)
        genes<-unique(c(genes,tmp))
}
B<-matrix(0,length(genes),length(datasetnames))
for (i in 1:length(genes)){
        for (j in 1:length(datasetnames)){
                B[i,j] <- sum(list.datasets[[j]]%in%genes[i])   
        }
}
dimnames(B)<-list(genes,datasetnames)
score<-vector()
for (i in 1:length(genes)){
        score[i]<-sum(B[i,])
}
M<-data.frame(genes,score)
ord<-order(M[,2],decreasing=TRUE)
m<-M[ord,]
# arrange m to get 3 rows (1st=geneid,2nd=genesymbol,3rd=score)
m1<-data.frame(m[,1],character(length=nrow(m)),m[,2])
colnames(m1)<-c("gene_id","gene_symbol","score")
# change 1st and 2nd columns of m1 from factor class to charcater class:
m1[,1] <- as.character(m1[,1])
m1[,2] <- as.character(m1[,2])
# if there's gene symbols in the entrez2sym.h for genes in our ranked list we get those
# we don't get the gene symbol for some genes of our ranked list, so we have to get those latter by hand
# library to retrieve object org.Hs.eg.db and associated functions:
library(org.Hs.eg.db)
# get list in which each element is a symbol for human gene and is named with corresponding entrez gene id:
x <- org.Hs.egSYMBOL
mapped_genes<- mappedkeys(x)
entrez2sym.h <- as.list(x[mapped_genes])

for (i in 1:nrow(m1)){
        if(m1[i,1]%in%names(entrez2sym.h)){
                m1[i,2]<-entrez2sym.h[[m1[i,1]]]
        } else if (m1[i,1]%in%homo$GeneID) {
                m1[i,2]<-subset(homo$Symbol,homo$GeneID==m1[i,1])
        }
}

write.csv(m1,"C:/Users/Tania_2/Internship/R/Results/ranked.lists/ranked.list.human.NSC.signature.csv",row.names=F)
C<-data.frame(B,score)
ord1<-order(C[,ncol(C)],decreasing=TRUE)
c<-C[ord1,]
c1 <- data.frame(gene_symbol=m1[,2],gene_id=m1[,1],c[,], stringsAsFactors = F)
write.csv(c1,"C:/Users/Tania_2/Internship/R/Results/ranked.lists/ranked.list.human.NSC.signature.genes.csv",row.names=F)

###################################################################################################################
###### MOuse ESC signatures
files_list <- list.files("Mm_ESC", full.names=TRUE) 
datasetnames<-vector()
for (i in 1:length(files_list)){
        datasetnames[i]<-unlist(strsplit(files_list[i],'[/.]'))[2]
}
genes<-character()
list.datasets <- list()
for (i in 1:length(files_list)){
        tmp <-  as.character(read.csv(files_list[i], fileEncoding="latin1", header=F)[,1])
        for (j in 1:length(tmp)) {tmp[j] <- gsub("[^0-9]","",tmp[j]) 
        }
        list.datasets[[i]] <- unique(tmp)
        genes<-unique(c(genes,tmp))
}
B<-matrix(0,length(genes),length(datasetnames))
for (i in 1:length(genes)){
        for (j in 1:length(datasetnames)){
                B[i,j] <- sum(list.datasets[[j]]%in%genes[i])   
        }
}
dimnames(B)<-list(genes,datasetnames)
score<-vector()
for (i in 1:length(genes)){
        score[i]<-sum(B[i,])
}
M<-data.frame(genes,score)
ord<-order(M[,2],decreasing=TRUE)
m<-M[ord,]
# arrange m to get 3 rows (1st=geneid,2nd=genesymbol,3rd=score)
m1<-data.frame(m[,1],character(length=nrow(m)),m[,2])
colnames(m1)<-c("gene_id","gene_symbol","score")
# change 1st and 2nd columns of m1 from factor class to charcater class:
m1[,1] <- as.character(m1[,1])
m1[,2] <- as.character(m1[,2])
# if there's gene symbols in the entrez2sym.h for genes in our ranked list we get those
# we don't get the gene symbol for some genes of our ranked list, so we have to get those latter by hand
# library to retrieve object org.Hs.eg.db and associated functions:
library(org.Mm.eg.db)
# get list in which each element is a symbol for human gene and is named with corresponding entrez gene id:
x <- org.Mm.egSYMBOL
mapped_genes<- mappedkeys(x)
entrez2sym.m <- as.list(x[mapped_genes])

for (i in 1:nrow(m1)){
        if(m1[i,1]%in%names(entrez2sym.m)){
                m1[i,2]<-entrez2sym.m[[m1[i,1]]]
        } else if (m1[i,1]%in%mouse$GeneID) {
                m1[i,2]<-subset(mouse$Symbol,mouse$GeneID==m1[i,1])
        }
}

write.csv(m1,"C:/Users/Tania_2/Internship/R/Results/ranked.lists/ranked.list.mouse.ESC.signature.csv",row.names=F)
C<-data.frame(B,score)
ord1<-order(C[,ncol(C)],decreasing=TRUE)
c<-C[ord1,]
c1 <- data.frame(gene_symbol=m1[,2],gene_id=m1[,1],c[,], stringsAsFactors = F)
write.csv(c1,"C:/Users/Tania_2/Internship/R/Results/ranked.lists/ranked.list.mouse.ESC.signature.genes.csv",row.names=F)

###################################################################################################################
###### MOuse HSC signatures
files_list <- list.files("Mm_HSC", full.names=TRUE) 
datasetnames<-vector()
for (i in 1:length(files_list)){
        datasetnames[i]<-unlist(strsplit(files_list[i],'[/.]'))[2]
}
genes<-character()
list.datasets <- list()
for (i in 1:length(files_list)){
        tmp <-  as.character(read.csv(files_list[i], fileEncoding="latin1", header=F)[,1])
        for (j in 1:length(tmp)) {tmp[j] <- gsub("[^0-9]","",tmp[j]) 
        }
        list.datasets[[i]] <- unique(tmp)
        genes<-unique(c(genes,tmp))
}
B<-matrix(0,length(genes),length(datasetnames))
for (i in 1:length(genes)){
        for (j in 1:length(datasetnames)){
                B[i,j] <- sum(list.datasets[[j]]%in%genes[i])   
        }
}
dimnames(B)<-list(genes,datasetnames)
score<-vector()
for (i in 1:length(genes)){
        score[i]<-sum(B[i,])
}
M<-data.frame(genes,score)
ord<-order(M[,2],decreasing=TRUE)
m<-M[ord,]
# arrange m to get 3 rows (1st=geneid,2nd=genesymbol,3rd=score)
m1<-data.frame(m[,1],character(length=nrow(m)),m[,2])
colnames(m1)<-c("gene_id","gene_symbol","score")
# change 1st and 2nd columns of m1 from factor class to charcater class:
m1[,1] <- as.character(m1[,1])
m1[,2] <- as.character(m1[,2])
# if there's gene symbols in the entrez2sym.h for genes in our ranked list we get those
# we don't get the gene symbol for some genes of our ranked list, so we have to get those latter by hand
# library to retrieve object org.Hs.eg.db and associated functions:
library(org.Mm.eg.db)
# get list in which each element is a symbol for human gene and is named with corresponding entrez gene id:
x <- org.Mm.egSYMBOL
mapped_genes<- mappedkeys(x)
entrez2sym.m <- as.list(x[mapped_genes])

for (i in 1:nrow(m1)){
        if(m1[i,1]%in%names(entrez2sym.m)){
                m1[i,2]<-entrez2sym.m[[m1[i,1]]]
        } else if (m1[i,1]%in%mouse$GeneID) {
                m1[i,2]<-subset(mouse$Symbol,mouse$GeneID==m1[i,1])
        }
}

write.csv(m1,"C:/Users/Tania_2/Internship/R/Results/ranked.lists/ranked.list.mouse.HSC.signature.csv",row.names=F)
C<-data.frame(B,score)
ord1<-order(C[,ncol(C)],decreasing=TRUE)
c<-C[ord1,]
c1 <- data.frame(gene_symbol=m1[,2],gene_id=m1[,1],c[,], stringsAsFactors = F)
write.csv(c1,"C:/Users/Tania_2/Internship/R/Results/ranked.lists/ranked.list.mouse.HSC.signature.genes.csv",row.names=F)

###################################################################################################################
###### MOuse NSC signatures
files_list <- list.files("Mm_NSC", full.names=TRUE) 
datasetnames<-vector()
for (i in 1:length(files_list)){
        datasetnames[i]<-unlist(strsplit(files_list[i],'[/.]'))[2]
}
genes<-character()
list.datasets <- list()
for (i in 1:length(files_list)){
        tmp <-  as.character(read.csv(files_list[i], fileEncoding="latin1", header=F)[,1])
        for (j in 1:length(tmp)) {tmp[j] <- gsub("[^0-9]","",tmp[j]) 
        }
        list.datasets[[i]] <- unique(tmp)
        genes<-unique(c(genes,tmp))
}
B<-matrix(0,length(genes),length(datasetnames))
for (i in 1:length(genes)){
        for (j in 1:length(datasetnames)){
                B[i,j] <- sum(list.datasets[[j]]%in%genes[i])   
        }
}
dimnames(B)<-list(genes,datasetnames)
score<-vector()
for (i in 1:length(genes)){
        score[i]<-sum(B[i,])
}
M<-data.frame(genes,score)
ord<-order(M[,2],decreasing=TRUE)
m<-M[ord,]
# arrange m to get 3 rows (1st=geneid,2nd=genesymbol,3rd=score)
m1<-data.frame(m[,1],character(length=nrow(m)),m[,2])
colnames(m1)<-c("gene_id","gene_symbol","score")
# change 1st and 2nd columns of m1 from factor class to charcater class:
m1[,1] <- as.character(m1[,1])
m1[,2] <- as.character(m1[,2])
# if there's gene symbols in the entrez2sym.h for genes in our ranked list we get those
# we don't get the gene symbol for some genes of our ranked list, so we have to get those latter by hand
# library to retrieve object org.Hs.eg.db and associated functions:
library(org.Hs.eg.db)
# get list in which each element is a symbol for human gene and is named with corresponding entrez gene id:
x <- org.Mm.egSYMBOL
mapped_genes<- mappedkeys(x)
entrez2sym.m <- as.list(x[mapped_genes])

for (i in 1:nrow(m1)){
        if(m1[i,1]%in%names(entrez2sym.m)){
                m1[i,2]<-entrez2sym.m[[m1[i,1]]]
        } else if (m1[i,1]%in%mouse$GeneID) {
                m1[i,2]<-subset(mouse$Symbol,mouse$GeneID==m1[i,1])
        }
}

write.csv(m1,"C:/Users/Tania_2/Internship/R/Results/ranked.lists/ranked.list.mouse.NSC.signature.csv",row.names=F)
C<-data.frame(B,score)
ord1<-order(C[,ncol(C)],decreasing=TRUE)
c<-C[ord1,]
c1 <- data.frame(gene_symbol=m1[,2],gene_id=m1[,1],c[,], stringsAsFactors = F)
write.csv(c1,"C:/Users/Tania_2/Internship/R/Results/ranked.lists/ranked.list.mouse.NSC.signature.genes.csv",row.names=F)
