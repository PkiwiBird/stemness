setwd("C:/Users/Tania_2/Internship/R/Data")
mus.musculus <- read.table("/Users/Tania_2/Internship/R/Data/Mus_musculus.gene_info.txt", header = T, stringsAsFactors = F, sep="\t", fill=T, quote="")
# put columns as characters:
for(i in 1:ncol(mus.musculus)){
        mus.musculus[,i] <- as.character(mus.musculus[,i])
}
homo.sap <- read.table("/Users/Tania_2/Internship/R/Data/Homo_sapiens.gene_info.txt", header = T, stringsAsFactors = F, sep="\t", fill=T, quote="")
# put columns as characters:
for(i in 1:ncol(homo.sap)){
        homo.sap[,i] <- as.character(homo.sap[,i])
}

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
data_files_list1 <- lapply(data_files_list, function(x) x[x%in%homo.sap$GeneID])
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
data_files_list1 <- lapply(data_files_list, function(x) x[x%in%homo.sap$GeneID])
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
data_files_list1 <- lapply(data_files_list, function(x) x[x%in%homo.sap$GeneID])
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
data_files_list1 <- lapply(data_files_list, function(x) x[x%in%mus.musculus$GeneID])
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
data_files_list1 <- lapply(data_files_list, function(x) x[x%in%mus.musculus$GeneID])
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
data_files_list1 <- lapply(data_files_list, function(x) x[x%in%mus.musculus$GeneID])
lapply(names(data_files_list1), function (x) write.table(data_files_list1[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/Mm_up_without_TF",x,sep="/"),"csv",sep="."), row.names=F, col.names=F))



