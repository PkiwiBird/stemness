setwd("C:/Users/Tania_2/Internship/R/Data")

a <- read.table("C:/Users/Tania_2/Internship/R/Data/gene_history4.txt", header = T,  stringsAsFactors = F, comment.char="")
withdraws <- as.character(subset(a$Discontinued_GeneID,a$GeneID == "-"))
replaced <- subset(a,a$GeneID != "-")
# read.table with pseudos in human is adapted because of the different row length given by the fact that some gene names have tabs inside them:
b <- read.table("/Users/Tania_2/Internship/R/Data/Homo_sapiens.gene_info.txt", header = T, stringsAsFactors = F, sep="\t", fill=T, quote="")
# put columns as characters:
for(i in 1:ncol(b)){
        b[,i] <- as.character(b[,i])
}
# retrieving pseudogenes in homosapiens
pseudo <- unique(subset(b$GeneID,b$type_of_gene=="pseudo"))

# read.table with pseudos in mouse is adapted because of the different row length given by the fact that some gene names have tabs inside them:
b2 <- read.table("/Users/Tania_2/Internship/R/Data/Mus_musculus.gene_info.txt", header = T, stringsAsFactors = F, sep="\t", fill=T, quote="")
# put columns as characters:
for(i in 1:ncol(b2)){
        b2[,i] <- as.character(b2[,i])
}
# retrieving pseudogenes in mouse
pseudo2 <- unique(subset(b2$GeneID,b2$type_of_gene=="pseudo"))


# put withdraws and pseudos together in the same vector:
withdraws.pseudos <- unique(c(withdraws, pseudo, pseudo2))

############
### Hs_down
# lists all files in that folder:
files_list <- list.files("Hs_down", full.names=TRUE)
data_files_list <- lapply(files_list, function(x) as.character(unique(read.csv(x, header=F, stringsAsFactors=F)[,1])))
# creates dataset names vector:
datasetnames<-vector()
for (i in 1:length(files_list)){
        datasetnames[i]<-unlist(strsplit(files_list[i],'[/.]'))[2]
}
names(data_files_list) <- datasetnames
data_files_list1 <- lapply(data_files_list, function(x) x[!(x%in%withdraws.pseudos)])
lapply(names(data_files_list1), function (x) write.table(data_files_list1[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/Hs_down",x,sep="/"),"csv",sep="."), row.names=F, col.names=F))

############
### Hs_up_TF
files_list <- list.files("Hs_up_TF", full.names=TRUE)
data_files_list <- lapply(files_list, function(x) as.character(unique(read.csv(x, header=F, stringsAsFactors=F)[,1])))
datasetnames<-vector()
for (i in 1:length(files_list)){
        datasetnames[i]<-unlist(strsplit(files_list[i],'[/.]'))[2]
}
names(data_files_list) <- datasetnames
data_files_list1 <- lapply(data_files_list, function(x) x[!(x%in%withdraws.pseudos)])
lapply(names(data_files_list1), function (x) write.table(data_files_list1[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/Hs_up_TF",x,sep="/"),"csv",sep="."), row.names=F, col.names=F))

############
### Hs_up_without_TF
files_list <- list.files("Hs_up_without_TF", full.names=TRUE)
data_files_list <- lapply(files_list, function(x) as.character(unique(read.csv(x, header=F, stringsAsFactors=F)[,1])))
datasetnames<-vector()
for (i in 1:length(files_list)){
        datasetnames[i]<-unlist(strsplit(files_list[i],'[/.]'))[2]
}
names(data_files_list) <- datasetnames
data_files_list1 <- lapply(data_files_list, function(x) x[!(x%in%withdraws.pseudos)])
lapply(names(data_files_list1), function (x) write.table(data_files_list1[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/Hs_up_without_TF",x,sep="/"),"csv",sep="."), row.names=F, col.names=F))

############
### Mm_down
files_list <- list.files("Mm_down", full.names=TRUE)
data_files_list <- lapply(files_list, function(x) as.character(unique(read.csv(x, header=F, stringsAsFactors=F)[,1])))
datasetnames<-vector()
for (i in 1:length(files_list)){
        datasetnames[i]<-unlist(strsplit(files_list[i],'[/.]'))[2]
}
names(data_files_list) <- datasetnames
data_files_list1 <- lapply(data_files_list, function(x) x[!(x%in%withdraws.pseudos)])
lapply(names(data_files_list1), function (x) write.table(data_files_list1[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/Mm_down",x,sep="/"),"csv",sep="."), row.names=F, col.names=F))

############
### Mm_up_TF
# lists all files in that folder:
files_list <- list.files("Mm_up_TF", full.names=TRUE)
data_files_list <- lapply(files_list, function(x) as.character(unique(read.csv(x, header=F, stringsAsFactors=F)[,1])))
datasetnames<-vector()
for (i in 1:length(files_list)){
        datasetnames[i]<-unlist(strsplit(files_list[i],'[/.]'))[2]
}
names(data_files_list) <- datasetnames
data_files_list1 <- lapply(data_files_list, function(x) x[!(x%in%withdraws.pseudos)])
lapply(names(data_files_list1), function (x) write.table(data_files_list1[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/Mm_up_TF",x,sep="/"),"csv",sep="."), row.names=F, col.names=F))

############
### Mm_up_without_TF
files_list <- list.files("Mm_up_without_TF", full.names=TRUE)
data_files_list <- lapply(files_list, function(x) as.character(unique(read.csv(x, header=F, stringsAsFactors=F)[,1])))
datasetnames<-vector()
for (i in 1:length(files_list)){
        datasetnames[i]<-unlist(strsplit(files_list[i],'[/.]'))[2]
}
names(data_files_list) <- datasetnames
data_files_list1 <- lapply(data_files_list, function(x) x[!(x%in%withdraws.pseudos)])
lapply(names(data_files_list1), function (x) write.table(data_files_list1[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/Mm_up_without_TF",x,sep="/"),"csv",sep="."), row.names=F, col.names=F))

