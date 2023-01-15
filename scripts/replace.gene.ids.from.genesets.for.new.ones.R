setwd("C:/Users/Tania_2/Internship/R/Data")
a <- read.table("C:/Users/Tania_2/Internship/R/Data/gene_history4.txt", header = T,  stringsAsFactors = F)
replaced <- subset(a,a$GeneID != "-")
replace.ids <- data.frame(replaced$Discontinued_GeneID, replaced$GeneID, stringsAsFactors=F)
# change 1st column from int to chr:
for (i in 1:nrow(replace.ids)){
        replace.ids[i,1] <- as.character(replace.ids[i,1])
}
# change 2nd column from int to chr:
for (i in 1:nrow(replace.ids)){
        replace.ids[i,2] <- as.character(replace.ids[i,2])
}
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
data_files_list1 <- lapply(data_files_list, function(x) c(x[!(x%in%replace.ids[,1])],replace.ids[,2][replace.ids[,1]%in%x]))
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
data_files_list1 <- lapply(data_files_list, function(x) c(x[!(x%in%replace.ids[,1])],replace.ids[,2][replace.ids[,1]%in%x]))
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
data_files_list1 <- lapply(data_files_list, function(x) c(x[!(x%in%replace.ids[,1])],replace.ids[,2][replace.ids[,1]%in%x]))
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
data_files_list1 <- lapply(data_files_list, function(x) c(x[!(x%in%replace.ids[,1])],replace.ids[,2][replace.ids[,1]%in%x]))
lapply(names(data_files_list1), function (x) write.table(data_files_list1[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/Mm_down",x,sep="/"),"csv",sep="."), row.names=F, col.names=F))


############
### Mm_up_TF
files_list <- list.files("Mm_up_TF", full.names=TRUE)
data_files_list <- lapply(files_list, function(x) as.character(unique(read.csv(x, header=F, stringsAsFactors=F)[,1])))
datasetnames<-vector()
for (i in 1:length(files_list)){
        datasetnames[i]<-unlist(strsplit(files_list[i],'[/.]'))[2]
}
names(data_files_list) <- datasetnames
data_files_list1 <- lapply(data_files_list, function(x) c(x[!(x%in%replace.ids[,1])],replace.ids[,2][replace.ids[,1]%in%x]))
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
data_files_list1 <- lapply(data_files_list, function(x) c(x[!(x%in%replace.ids[,1])],replace.ids[,2][replace.ids[,1]%in%x]))
lapply(names(data_files_list1), function (x) write.table(data_files_list1[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/Mm_up_without_TF",x,sep="/"),"csv",sep="."), row.names=F, col.names=F))

