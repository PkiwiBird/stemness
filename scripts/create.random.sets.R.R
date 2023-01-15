setwd("D:/Tania_2/Internship/R/Data")
R <- 100000 # number of random draws

### Human
## get histogram of number of genes that appear a specific number of times in the different human stemness gene sets:
humanISS <- read.table("D:/Tania_2/Internship/R/Results/ranked.lists/ranked.list.human.stem.signature.csv", sep=",", header= TRUE)
CountHs <- hist(humanISS$score,breaks=seq(0,21,1),plot=FALSE)$counts # number of genes that show up a specific number of times in human stemness signatures. each position the vector is a number of times (it starts with genes that show up 1 time, then genes that show up 2 times,etc...)
                                                                     # 21 is number of human stemness gene sets
tiff(file="D:/Tania_2/Internship/R/Results/randomization_procedure/human_stemsets_score_distribution.tif", width = 8*600, height = 8*600, res = 600)
barplot(CountHs[1:max(humanISS$score)],col="black",xlab="Score",ylab="Genes",main="Stemness signatures",names.arg=seq(1,max(humanISS$score),1), ylim = c(0,3000)) # we just show distribution until 12 times for human, cause that's max score we get in human stemness list
dev.off()

png(file="D:/Tania_2/Internship/R/Results/randomization_procedure/human_stemsets_score_distribution.png", width = 8*600, height = 8*600, res = 600)
barplot(CountHs[1:max(humanISS$score)],col="black",xlab="Score",ylab="Genes",main="Stemness signatures",names.arg=seq(1,max(humanISS$score),1), ylim = c(0,3000)) # we just show distribution until 12 times for human, cause that's max score we get in human stemness list
dev.off()
                                                                                                                              # max(humanISS$score) is 12
## get histogram of number of genes that appear a specific number of times in the different randomly generated stemness gene sets:
# retrieve human stemness gene sets names:
files_list <- list.files("Hs_up_without_TF", full.names=TRUE) 
datasetnames <- vector()
for (i in 1:length(files_list)){
        datasetnames[i]<-unlist(strsplit(files_list[i],'[/.]'))[1]
}

# retrieve lengths of human stemness gene sets:
list.lengths <- vector()
for (i in 1:length(files_list)){
        tmp <-  as.character(read.csv(files_list[i], fileEncoding="latin1", header=F)[,1])
        for (j in 1:length(tmp)) {tmp[j] <- gsub("[^0-9]","",tmp[j]) # when writting the entrez id names in the files of each gene set, same weird spaces were introduced that were affecting the ranking of genes
        # so this step is to get rid of those spaces. [^0-9] regular expression means find any character that is not a digit
        }
        list.lengths[i] <- length(unique(tmp)) # to get the sizes of the sets. to guarantee that each set doesn't have repeated genes that could be counted twice in the ranking we use unique
}
# get list with all possible human genes that can be included in random generated human gene sets:
library(org.Hs.eg.db)
x <- org.Hs.egSYMBOL
mapped_genes<- mappedkeys(x)
entrez2sym.h <- as.list(x[mapped_genes])
allgenes <- unique(unlist(entrez2sym.h))
# for each random draw (each time that we randomly split genes between gene sets with equal size to human stemness gene sets) we count how many genes appear how many times:
RandomCount <- rep(0,length(list.lengths)) # starts as a vector of 0's, but will be a matrix in which each row is a random draw
for (i in 1:R){
        list.datasetsR <- sapply(list.lengths, function(x) sample(allgenes,x)) # for each random draw: randomly select human genes for each dataset. the size (# of genes) of each dataset is same as in gene signatures.
                                                                               # each set of genes is saved into a list
        names(list.datasetsR) <- datasetnames # names each gene set in the list with name of stemness gene signature
        N <- hist(table(unlist(list.datasetsR)),breaks=seq(0,length(list.datasetsR),1),plot=FALSE) # unlist(list.datasetsR) creates a vector with all genes of all genesets (non-unique)
                                                                                                   # table gives a vector in which each element name is a gene present at least one time in the randomly generated gene sets
                                                                                                   # and each vector element is number of times that gene shows up in different gene sets
                                                                                                   # hist works as table of table, gives vector with each element is number of genes that show a certain number of times in gene sets
                                                                                                   # N$counts gives numeric result of hist. plot=FALSE avoids plotting
        RandomCount <- rbind(RandomCount,N$counts)
}

RandomCountHs <- RandomCount[-1,] # removes the first row which is an empty starting row.
                                  # gives matrix, each row corresponds to a random draw (each time that we randomly split genes between gene sets with equal size to human stemness gene sets)
                                  # each row is a vector and the index of each element of that vector corresponds to a number of times genes appear in random sets generated
                                  # each elemnt of the vector is number of genes that appear that number of times

CountHsR <- apply(RandomCountHs,2,sum)/R # sums elements of each column of matrix and divides by number of random draws
                                         # get average number of genes that show up a certain number of times considering all random draws                    
tiff(file="D:/Tania_2/Internship/R/Results/randomization_procedure/human_randomsets_score_distribution.tif", width = 8*600, height = 8*600, res = 600)
barplot(CountHsR[1:max(humanISS$score)],col="gray",xlab="Score",ylab="Genes",main="Random gene lists",names.arg=seq(1,max(humanISS$score),1), ylim = c(0,5000))
dev.off()

png(file="D:/Tania_2/Internship/R/Results/randomization_procedure/human_randomsets_score_distribution.png", width = 8*600, height = 8*600, res = 600)
barplot(CountHsR[1:max(humanISS$score)],col="gray",xlab="Score",ylab="Genes",main="Random gene lists",names.arg=seq(1,max(humanISS$score),1), ylim = c(0,5000))
dev.off()

## empirical FDR
obs <- rep(0,length(list.datasetsR))
expected <- rep(0,length(list.datasetsR))
FDR <- rep(NA,length(list.datasetsR))
for (i in length(list.datasetsR):1){ # goes from 21 to 1, 1 by 1
  obs[i] <- sum(CountHs[i:length(list.datasetsR)]) # 1st CountHs[21:21], then CountHs[20:21], then CountHs[19:21], etc
                                                   # for stemness signatures: number of genes with score 21, number of genes with score 20 or more, number of genes with score 19 or more
                                                   # score is number of times a gene shows up in the different human stemness sigrantures
  expected[i] <- sum(CountHsR[i:length(list.datasetsR)]) # 1st CountHsR[21:21], then CountHsR[20:21], then CountHsR[19:21], etc
                                                         # for randomly generated gene sets: average number of genes with score 21, average number of genes with score 20 or more, average genes with score 19 or more, etc
  if (obs[i] != 0) FDR[i] <- expected[i]/obs[i] # when score is 12 or lower, FDR is calculated for each score as:
                                                # number of genes with that score or higher in the random generated gene sets dividing by number of gene with that score or higher in stemness signatures
}

# FDR[FDR < 10^(-5)] <- 10^(-5) # as we plot FDR as log10(FDR) and log10(0) gives -INf wich cannot be plotted
                              # we have to turn 0 values into a very small value, and this is the solution
                              # we don't do FDR[FDR == 0] <- 10^(-5), cause there may be FDR values that are non-zero and smaller than 10 ^(-5)
FDR[FDR == 0] <- sort(unique(FDR))[2] # as we plot FDR as log10(FDR) and log10(0) gives -INf wich cannot be plotted
                                      # we have to turn 0 values into a very small value
                                      # we can replace 0 by the next smaller FDR value
FDR[FDR > 1] <- 1 # as FDR represents a p-value and p-values greater than 1 are contra-intuitive, we set all FDR > 1 to be 1.
tiff(file="D:/Tania_2/Internship/R/Results/randomization_procedure/human_empiricalFDR.tif", width = 8*600, height = 8*600, res = 600)
plot(log10(FDR)[1:max(humanISS$score)],type="b",col="red",pch=20,xlab="Score",xaxt="n", main = "Empirical FDR",font.lab=2,ylab="log10(FDR)")
axis(1, at=seq(1, max(humanISS$score), 1), las=1)
dev.off()

png(file="D:/Tania_2/Internship/R/Results/randomization_procedure/human_empiricalFDR.png", width = 8*600, height = 8*600, res = 600)
plot(log10(FDR)[1:max(humanISS$score)],type="b",col="red",pch=20,xlab="Score",xaxt="n", main = "Empirical FDR",font.lab=2,ylab="log10(FDR)")
axis(1, at=seq(1, max(humanISS$score), 1), las=1)
dev.off()

##### mouse 
## get histogram of number of genes that appear a specific number of times in the different human stemness gene sets:
mouseISS <- read.table("D:/Tania_2/Internship/R/Results/ranked.lists/ranked.list.mouse.stem.signature.csv", sep=",", header= TRUE)
CountMm <- hist(mouseISS$score,breaks=seq(0,21,1),plot=FALSE)$counts # number of genes that show up a specific number of times in mouse stemness signatures. each position the vector is a number of times (it starts with genes that show up 1 time, then genes that show up 2 times,etc...)
                                                                     # 21 is number of mouse stemness gene sets
tiff(file="D:/Tania_2/Internship/R/Results/randomization_procedure/mouse_stemsets_score_distribution.tif", width = 8*600, height = 8*600, res = 600)
barplot(CountMm[1:max(mouseISS$score)],col="black",xlab="Score",ylab="Genes",main="Stemness signatures",names.arg=seq(1,max(mouseISS$score),1), ylim = c(0,3500)) # we just show distribution until 13 times for mouse, cause that's max score we get in mouse stemness list
dev.off()

png(file="D:/Tania_2/Internship/R/Results/randomization_procedure/mouse_stemsets_score_distribution.png", width = 8*600, height = 8*600, res = 600)
barplot(CountMm[1:max(mouseISS$score)],col="black",xlab="Score",ylab="Genes",main="Stemness signatures",names.arg=seq(1,max(mouseISS$score),1), ylim = c(0,3500)) # we just show distribution until 13 times for mouse, cause that's max score we get in mouse stemness list
dev.off()

## get histogram of number of genes that appear a specific number of times in the different randomly generated stemness gene sets:
# retrieve mouse stemness gene sets names:
files_list <- list.files("Mm_up_without_TF", full.names=TRUE) 
datasetnames <- vector()
for (i in 1:length(files_list)){
  datasetnames[i]<-unlist(strsplit(files_list[i],'[/.]'))[1]
}

# retrieve lengths of mouse stemness gene sets:
list.lengths <- vector()
for (i in 1:length(files_list)){
  tmp <-  as.character(read.csv(files_list[i], fileEncoding="latin1", header=F)[,1])
  for (j in 1:length(tmp)) {tmp[j] <- gsub("[^0-9]","",tmp[j]) # when writting the entrez id names in the files of each gene set, same weird spaces were introduced that were affecting the ranking of genes
  # so this step is to get rid of those spaces. [^0-9] regular expression means find any character that is not a digit
  }
  list.lengths[i] <- length(unique(tmp)) # to get the sizes of the sets. to guarantee that each set doesn't have repeated genes that could be counted twice in the ranking we use unique
}

# get list with all possible mouse genes that can be included in random generated mouse gene sets:
library(org.Mm.eg.db) # library to retrieve object org.Mm.eg.db and associated functions:
x <- org.Mm.egSYMBOL
mapped_genes<- mappedkeys(x)
entrez2sym.m <- as.list(x[mapped_genes])
allgenes <- unique(unlist(entrez2sym.m))
# for each random draw (each time that we randomly split genes between gene sets with equal size to mouse stemness gene sets) we count how many genes appear how many times:
RandomCount <- rep(0,length(list.lengths)) # starts as a vector of 0's, but will be a matrix in which each row is a random draw
for (i in 1:R){
  list.datasetsR <- sapply(list.lengths, function(x) sample(allgenes,x)) # for each random draw: randomly select mouse genes for each dataset. the size (# of genes) of each dataset is same as in gene signatures.
  # each set of genes is saved into a list
  names(list.datasetsR) <- datasetnames # names each gene set in the list with name of stemness gene signature
  N <- hist(table(unlist(list.datasetsR)),breaks=seq(0,length(list.datasetsR),1),plot=FALSE) # unlist(list.datasetsR) creates a vector with all genes of all genesets (non-unique)
  # table gives a vector in which each element name is a gene present at least one time in the randomly generated gene sets
  # and each vector element is number of times that gene shows up in different gene sets
  # hist works as table of table, gives vector with each element is number of genes that show a certain number of times in gene sets
  # N$counts gives numeric result of hist. plot=FALSE avoids plotting
  RandomCount <- rbind(RandomCount,N$counts)
}

RandomCountMm <- RandomCount[-1,] # removes the first row which is an empty starting row.
# gives matrix, each row corresponds to a random draw (each time that we randomly split genes between gene sets with equal size to mouse stemness gene sets)
# each row is a vector and the index of each element of that vector corresponds to a number of times genes appear in random sets generated
# each elemnt of the vector is number of genes that appear that number of times

CountMmR <- apply(RandomCountMm,2,sum)/R # sums elements of each column of matrix and divides by number of random draws
# get average number of genes that show up a certain number of times considering all random draws                    
tiff(file="D:/Tania_2/Internship/R/Results/randomization_procedure/mouse_randomsets_score_distribution.tif", width = 8*600, height = 8*600, res = 600)
barplot(CountMmR[1:max(mouseISS$score)],col="gray",xlab="Score",ylab="Genes",main="Random gene lists",names.arg=seq(1,max(mouseISS$score),1), ylim = c(0,5000))
dev.off()

png(file="D:/Tania_2/Internship/R/Results/randomization_procedure/mouse_randomsets_score_distribution.png", width = 8*600, height = 8*600, res = 600)
barplot(CountMmR[1:max(mouseISS$score)],col="gray",xlab="Score",ylab="Genes",main="Random gene lists",names.arg=seq(1,max(mouseISS$score),1), ylim = c(0,5000))
dev.off()

## empirical FDR
obs <- rep(0,length(list.datasetsR))
expected <- rep(0,length(list.datasetsR))
FDR <- rep(NA,length(list.datasetsR))
for (i in length(list.datasetsR):1){ # goes from 21 to 1, 1 by 1
  obs[i] <- sum(CountMm[i:length(list.datasetsR)]) # 1st CountMm[21:21], then CountMm[20:21], then CountMm[19:21], etc
  # for stemness signatures: number of genes with score 21, number of genes with score 20 or more, number of genes with score 19 or more
  # score is number of times a gene shows up in the different mouse stemness sigrantures
  expected[i] <- sum(CountMmR[i:length(list.datasetsR)]) # 1st CountMmR[21:21], then CountMmR[20:21], then CountMmR[19:21], etc
  # for randomly generated gene sets: average number of genes with score 21, average number of genes with score 20 or more, average genes with score 19 or more, etc
  if (obs[i] != 0) FDR[i] <- expected[i]/obs[i] # when score is 13 or lower, FDR is calculated for each score as:
  # number of genes with that score or higher in the random generated gene sets dividing by number of gene with that score or higher in stemness signatures
}

# FDR[FDR < 10^(-5)] <- 10^(-5) # as we plot FDR as log10(FDR) and log10(0) gives -INf wich cannot be plotted
# we have to turn 0 values into a very small value, and this is the solution
# we don't do FDR[FDR == 0] <- 10^(-5), cause there may be FDR values that are non-zero and smaller than 10 ^(-5)
FDR[FDR == 0] <- sort(unique(FDR))[2] # as we plot FDR as log10(FDR) and log10(0) gives -INf wich cannot be plotted
# we have to turn 0 values into a very small value
# we can replace 0 by the next smaller FDR value
FDR[FDR > 1] <- 1 # as FDR represents a p-value and p-values greater than 1 are contra-intuitive, we set all FDR > 1 to be 1.
tiff(file="D:/Tania_2/Internship/R/Results/randomization_procedure/mouse_empiricalFDR.tif", width = 8*600, height = 8*600, res = 600)
plot(log10(FDR)[1:max(mouseISS$score)],type="b",col="red",pch=20,xlab="Score",xaxt="n", main = "Empirical FDR",font.lab=2,ylab="log10(FDR)")
axis(1, at=seq(1, max(mouseISS$score), 1), las=1)
dev.off()

png(file="D:/Tania_2/Internship/R/Results/randomization_procedure/mouse_empiricalFDR.png", width = 8*600, height = 8*600, res = 600)
plot(log10(FDR)[1:max(mouseISS$score)],type="b",col="red",pch=20,xlab="Score",xaxt="n", main = "Empirical FDR",font.lab=2,ylab="log10(FDR)")
axis(1, at=seq(1, max(mouseISS$score), 1), las=1)
dev.off()





