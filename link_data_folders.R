setwd("C:/Users/Tania_2/Internship/R/Data")
#########
# link Hs_up_TF and Hs_up_without_TF to Hs_up_with&without_TF
files_list1 <- list.files("Hs_up_TF", full.names=TRUE)   #list files in the folder
data_files_list1<-lapply(files_list1, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames1<-vector()
for (i in 1:length(files_list1)){
        datasetnames1[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list1) <- datasetnames1
files_list2 <- list.files("Hs_up_without_TF", full.names=TRUE)   #list files in the folder
data_files_list2<-lapply(files_list2, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames2<-vector()
for (i in 1:length(files_list2)){
        datasetnames2[i]<-unlist(strsplit(files_list2[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list2) <- datasetnames2
lapply(names(data_files_list1), function (x) write.table(data_files_list1[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/Hs_up_with&without_TF",x, sep="/"),"csv",sep="."), col.names=F, row.names=F))
lapply(names(data_files_list2), function (x) write.table(data_files_list2[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/Hs_up_with&without_TF",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))
#####
###link Mm_up_TF and Mm_up_without_TF to Mm_up_with&without_TF
files_list1 <- list.files("Mm_up_TF", full.names=TRUE)   #list files in the folder
data_files_list1<-lapply(files_list1, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames1<-vector()
for (i in 1:length(files_list1)){
        datasetnames1[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list1) <- datasetnames1
files_list2 <- list.files("Mm_up_without_TF", full.names=TRUE)   #list files in the folder
data_files_list2<-lapply(files_list2, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames2<-vector()
for (i in 1:length(files_list2)){
        datasetnames2[i]<-unlist(strsplit(files_list2[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list2) <- datasetnames2
lapply(names(data_files_list1), function (x) write.table(data_files_list1[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/Mm_up_with&without_TF",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))
lapply(names(data_files_list2), function (x) write.table(data_files_list2[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/Mm_up_with&without_TF",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))

# link Hs_up_with&without_TF, and Mm_up_with&without_TF to allup
files_list1 <- list.files("Hs_up_with&without_TF", full.names=TRUE)   #list files in the folder
data_files_list1<-lapply(files_list1, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames1<-vector()
for (i in 1:length(files_list1)){
        datasetnames1[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list1) <- datasetnames1
files_list2 <- list.files("Mm_up_with&without_TF", full.names=TRUE)   #list files in the folder
data_files_list2<-lapply(files_list2, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames2<-vector()
for (i in 1:length(files_list2)){
        datasetnames2[i]<-unlist(strsplit(files_list2[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list2) <- datasetnames2
lapply(names(data_files_list1), function (x) write.table(data_files_list1[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/allup",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))
lapply(names(data_files_list2), function (x) write.table(data_files_list2[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/allup",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))

# link Hs_down, and Mm_down to alldown
files_list1 <- list.files("Hs_down", full.names=TRUE)   #list files in the folder
data_files_list1<-lapply(files_list1, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames1<-vector()
for (i in 1:length(files_list1)){
        datasetnames1[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list1) <- datasetnames1
files_list2 <- list.files("Mm_down", full.names=TRUE)   #list files in the folder
data_files_list2<-lapply(files_list2, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames2<-vector()
for (i in 1:length(files_list2)){
        datasetnames2[i]<-unlist(strsplit(files_list2[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list2) <- datasetnames2
lapply(names(data_files_list1), function (x) write.table(data_files_list1[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/alldown",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))
lapply(names(data_files_list2), function (x) write.table(data_files_list2[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/alldown",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))


# link Mm_down, and Mm_up_with&without_TF to allMm
files_list1 <- list.files("Mm_down", full.names=TRUE)   #list files in the folder
data_files_list1<-lapply(files_list1, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames1<-vector()
for (i in 1:length(files_list1)){
        datasetnames1[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list1) <- datasetnames1
files_list2 <- list.files("Mm_up_with&without_TF", full.names=TRUE)   #list files in the folder
data_files_list2<-lapply(files_list2, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames2<-vector()
for (i in 1:length(files_list2)){
        datasetnames2[i]<-unlist(strsplit(files_list2[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list2) <- datasetnames2
lapply(names(data_files_list1), function (x) write.table(data_files_list1[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/allMm",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))
lapply(names(data_files_list2), function (x) write.table(data_files_list2[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/allMm",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))

# link Hs_down, and Hs_up_with&without_TF to allHs
files_list1 <- list.files("Hs_down", full.names=TRUE)   #list files in the folder
data_files_list1<-lapply(files_list1, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames1<-vector()
for (i in 1:length(files_list1)){
        datasetnames1[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list1) <- datasetnames1
files_list2 <- list.files("Hs_up_with&without_TF", full.names=TRUE)   #list files in the folder
data_files_list2<-lapply(files_list2, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames2<-vector()
for (i in 1:length(files_list2)){
        datasetnames2[i]<-unlist(strsplit(files_list2[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list2) <- datasetnames2
lapply(names(data_files_list1), function (x) write.table(data_files_list1[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/allHs",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))
lapply(names(data_files_list2), function (x) write.table(data_files_list2[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/allHs",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))

# link Hs_up_without_TF, Hs_down, Mm_up_without_TF and Mm_down to allgenesets_without_TF
files_list1 <- list.files("Hs_up_without_TF", full.names=TRUE)   #list files in the folder
data_files_list1<-lapply(files_list1, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames1<-vector()
for (i in 1:length(files_list1)){
        datasetnames1[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list1) <- datasetnames1
files_list2 <- list.files("Mm_up_without_TF", full.names=TRUE)   #list files in the folder
data_files_list2<-lapply(files_list2, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames2<-vector()
for (i in 1:length(files_list2)){
        datasetnames2[i]<-unlist(strsplit(files_list2[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list2) <- datasetnames2
files_list3 <- list.files("Hs_down", full.names=TRUE)   #list files in the folder
data_files_list3<-lapply(files_list3, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames3<-vector()
for (i in 1:length(files_list3)){
        datasetnames3[i]<-unlist(strsplit(files_list3[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list3) <- datasetnames3
files_list4 <- list.files("Mm_down", full.names=TRUE)   #list files in the folder
data_files_list4<-lapply(files_list4, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames4<-vector()
for (i in 1:length(files_list4)){
        datasetnames4[i]<-unlist(strsplit(files_list4[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list4) <- datasetnames4
lapply(names(data_files_list1), function (x) write.table(data_files_list1[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/allgenesets_without_TF",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))
lapply(names(data_files_list2), function (x) write.table(data_files_list2[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/allgenesets_without_TF",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))
lapply(names(data_files_list3), function (x) write.table(data_files_list3[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/allgenesets_without_TF",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))
lapply(names(data_files_list4), function (x) write.table(data_files_list4[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/allgenesets_without_TF",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))

# link Hs_up_TF, and Mm_up_TF to allTF
files_list1 <- list.files("Hs_up_TF", full.names=TRUE)   #list files in the folder
data_files_list1<-lapply(files_list1, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames1<-vector()
for (i in 1:length(files_list1)){
        datasetnames1[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list1) <- datasetnames1
files_list2 <- list.files("Mm_up_TF", full.names=TRUE)   #list files in the folder
data_files_list2<-lapply(files_list2, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames2<-vector()
for (i in 1:length(files_list2)){
        datasetnames2[i]<-unlist(strsplit(files_list2[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list2) <- datasetnames2
lapply(names(data_files_list1), function (x) write.table(data_files_list1[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/allTF",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))
lapply(names(data_files_list2), function (x) write.table(data_files_list2[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/allTF",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))


# link allHs, and allMm to allgenesets
files_list1 <- list.files("allHs", full.names=TRUE)   #list files in the folder
data_files_list1<-lapply(files_list1, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames1<-vector()
for (i in 1:length(files_list1)){
        datasetnames1[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list1) <- datasetnames1
files_list2 <- list.files("allMm", full.names=TRUE)   #list files in the folder
data_files_list2<-lapply(files_list2, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames2<-vector()
for (i in 1:length(files_list2)){
        datasetnames2[i]<-unlist(strsplit(files_list2[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list2) <- datasetnames2
lapply(names(data_files_list1), function (x) write.table(data_files_list1[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/allgenesets",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))
lapply(names(data_files_list2), function (x) write.table(data_files_list2[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/allgenesets",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))

# link Hs_up_without_TF and Hs_down to allHs_without_TF. Note: this isn't the same as allHs because allHs includes the TFs
files_list1 <- list.files("Hs_up_without_TF", full.names=TRUE)   #list files in the folder
data_files_list1<-lapply(files_list1, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames1<-vector()
for (i in 1:length(files_list1)){
        datasetnames1[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list1) <- datasetnames1
files_list2 <- list.files("Hs_down", full.names=TRUE)   #list files in the folder
data_files_list2<-lapply(files_list2, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames2<-vector()
for (i in 1:length(files_list2)){
        datasetnames2[i]<-unlist(strsplit(files_list2[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list2) <- datasetnames2
lapply(names(data_files_list1), function (x) write.table(data_files_list1[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/allHs_without_TF",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))
lapply(names(data_files_list2), function (x) write.table(data_files_list2[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/allHs_without_TF",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))

# link Mm_up_without_TF and Mm_down to allMm_without_TF. Note: this isn't the same as allMm because allMm includes the TFs
files_list1 <- list.files("Mm_up_without_TF", full.names=TRUE)   #list files in the folder
data_files_list1<-lapply(files_list1, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames1<-vector()
for (i in 1:length(files_list1)){
        datasetnames1[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list1) <- datasetnames1
files_list2 <- list.files("Mm_down", full.names=TRUE)   #list files in the folder
data_files_list2<-lapply(files_list2, function(x) read.csv(x, header=F, stringsAsFactors=F))   #create a list with the data from all files
datasetnames2<-vector()
for (i in 1:length(files_list2)){
        datasetnames2[i]<-unlist(strsplit(files_list2[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list2) <- datasetnames2
lapply(names(data_files_list1), function (x) write.table(data_files_list1[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/allMm_without_TF",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))
lapply(names(data_files_list2), function (x) write.table(data_files_list2[[x]], paste(paste("C:/Users/Tania_2/Internship/R/Data/allMm_without_TF",x,sep="/"),"csv",sep="."), col.names=F, row.names=F))

