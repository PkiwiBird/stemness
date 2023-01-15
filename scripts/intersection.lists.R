########################################################################
######All gene sets
setwd("C:/Users/Tania_2/Internship/R/Data")
files_list1 <- list.files("allgenesets", full.names=TRUE)   #list files in the folder
data_files_list<-lapply(files_list1, function(x) as.character(read.table(x)[,1]))   #create a list with the data from all files
datasetnames<-vector()
for (i in 1:length(files_list1)){
        datasetnames[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list) <- datasetnames

int<-Reduce(intersect,data_files_list) # Reduce applies a function to elements of a vector or list
int
write.table(int,"C:/Users/Tania_2/Internship/R/Results/intersected.lists/intersected.all.lists.csv",row.names=F, col.names=F)
########################################################################
######All gene up-regulated
setwd("C:/Users/Tania_2/Internship/R/Data")
files_list1 <- list.files("allup", full.names=TRUE)   #list files in the folder
data_files_list<-lapply(files_list1, function(x) as.character(read.table(x)[,1]))   #create a list with the data from all files
datasetnames<-vector()
for (i in 1:length(files_list1)){
        datasetnames[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list) <- datasetnames

int<-Reduce(intersect,data_files_list) # Reduce applies a function to elements of a vector or list
int
write.table(int,"C:/Users/Tania_2/Internship/R/Results/intersected.lists/intersected.all.up.csv",row.names=F, col.names=F)
########################################################################
######All gene down-regulated
setwd("C:/Users/Tania_2/Internship/R/Data")
files_list1 <- list.files("alldown", full.names=TRUE)   #list files in the folder
data_files_list<-lapply(files_list1, function(x) as.character(read.table(x)[,1]))   #create a list with the data from all files
datasetnames<-vector()
for (i in 1:length(files_list1)){
        datasetnames[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list) <- datasetnames

int<-Reduce(intersect,data_files_list) # Reduce applies a function to elements of a vector or list
int
write.table(int,"C:/Users/Tania_2/Internship/R/Results/intersected.lists/intersected.all.down.csv",row.names=F, col.names=F)
########################################################################
######All human
setwd("C:/Users/Tania_2/Internship/R/Data")
files_list1 <- list.files("allHs", full.names=TRUE)   #list files in the folder
data_files_list<-lapply(files_list1, function(x) as.character(read.table(x)[,1]))   #create a list with the data from all files
datasetnames<-vector()
for (i in 1:length(files_list1)){
        datasetnames[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list) <- datasetnames

int<-Reduce(intersect,data_files_list) # Reduce applies a function to elements of a vector or list
int
write.table(int,"C:/Users/Tania_2/Internship/R/Results/intersected.lists/intersected.all.human.csv",row.names=F, col.names=F)
########################################################################
######All mouse
setwd("C:/Users/Tania_2/Internship/R/Data")
files_list1 <- list.files("allMm", full.names=TRUE)   #list files in the folder
data_files_list<-lapply(files_list1, function(x) as.character(read.table(x)[,1]))   #create a list with the data from all files
datasetnames<-vector()
for (i in 1:length(files_list1)){
        datasetnames[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list) <- datasetnames

int<-Reduce(intersect,data_files_list) # Reduce applies a function to elements of a vector or list
int
write.table(int,"C:/Users/Tania_2/Internship/R/Results/intersected.lists/intersected.all.mouse.csv",row.names=F, col.names=F)
########################################################################
######Human down
setwd("C:/Users/Tania_2/Internship/R/Data")
files_list1 <- list.files("Hs_down", full.names=TRUE)   #list files in the folder
data_files_list<-lapply(files_list1, function(x) as.character(read.table(x)[,1]))   #create a list with the data from all files
datasetnames<-vector()
for (i in 1:length(files_list1)){
        datasetnames[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list) <- datasetnames

int<-Reduce(intersect,data_files_list) # Reduce applies a function to elements of a vector or list
int
write.table(int,"C:/Users/Tania_2/Internship/R/Results/intersected.lists/intersected.human.down.csv",row.names=F, col.names=F)
########################################################################
######Human up
setwd("C:/Users/Tania_2/Internship/R/Data")
files_list1 <- list.files("Hs_up_with&without_TF", full.names=TRUE)   #list files in the folder
data_files_list<-lapply(files_list1, function(x) as.character(read.table(x)[,1]))   #create a list with the data from all files
datasetnames<-vector()
for (i in 1:length(files_list1)){
        datasetnames[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list) <- datasetnames

int<-Reduce(intersect,data_files_list) # Reduce applies a function to elements of a vector or list
int
write.table(int,"C:/Users/Tania_2/Internship/R/Results/intersected.lists/intersected.human.up.csv",row.names=F, col.names=F)
########################################################################
######Mouse down
setwd("C:/Users/Tania_2/Internship/R/Data")
files_list1 <- list.files("Mm_down", full.names=TRUE)   #list files in the folder
data_files_list<-lapply(files_list1, function(x) as.character(read.table(x)[,1]))   #create a list with the data from all files
datasetnames<-vector()
for (i in 1:length(files_list1)){
        datasetnames[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list) <- datasetnames

int<-Reduce(intersect,data_files_list) # Reduce applies a function to elements of a vector or list
int
write.table(int,"C:/Users/Tania_2/Internship/R/Results/intersected.lists/intersected.mouse.down.csv",row.names=F,col.names=F)
########################################################################
######Mouse up
setwd("C:/Users/Tania_2/Internship/R/Data")
files_list1 <- list.files("Mm_up_with&without_TF", full.names=TRUE)   #list files in the folder
data_files_list<-lapply(files_list1, function(x) as.character(read.table(x)[,1]))   #create a list with the data from all files
datasetnames<-vector()
for (i in 1:length(files_list1)){
        datasetnames[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list) <- datasetnames

int<-Reduce(intersect,data_files_list) # Reduce applies a function to elements of a vector or list
int
write.table(int,"C:/Users/Tania_2/Internship/R/Results/intersected.lists/intersected.mouse.up.csv",row.names=F, col.names=F)

########################################################################
######Human up without TF
setwd("C:/Users/Tania_2/Internship/R/Data")
files_list1 <- list.files("Hs_up_without_TF", full.names=TRUE)   #list files in the folder
data_files_list<-lapply(files_list1, function(x) as.character(read.table(x)[,1]))   #create a list with the data from all files
datasetnames<-vector()
for (i in 1:length(files_list1)){
        datasetnames[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list) <- datasetnames
int <- Reduce(intersect,data_files_list) # Reduce applies a function to elements of a vector or list
int
write.table(int,"C:/Users/Tania_2/Internship/R/Results/intersected.lists/intersected.human.up.without.TF.csv",row.names=F, col.names=F)
########################################################################
######Mouse up without TF
setwd("C:/Users/Tania_2/Internship/R/Data")
files_list1 <- list.files("Mm_up_without_TF", full.names=TRUE)   #list files in the folder
data_files_list<-lapply(files_list1, function(x) as.character(read.table(x)[,1]))   #create a list with the data from all files
datasetnames<-vector()
for (i in 1:length(files_list1)){
        datasetnames[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list) <- datasetnames

int<-Reduce(intersect,data_files_list) # Reduce applies a function to elements of a vector or list
int
write.table(int,"C:/Users/Tania_2/Internship/R/Results/intersected.lists/intersected.mouse.up.without.TF.csv",row.names=F, col.names=F)
########################################################################
######Human up only TF
setwd("C:/Users/Tania_2/Internship/R/Data")
files_list1 <- list.files("Hs_up_TF", full.names=TRUE)   #list files in the folder
data_files_list<-lapply(files_list1, function(x) as.character(read.table(x)[,1]))   #create a list with the data from all files
datasetnames<-vector()
for (i in 1:length(files_list1)){
        datasetnames[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list) <- datasetnames

int<-Reduce(intersect,data_files_list) # Reduce applies a function to elements of a vector or list
int
write.table(int,"C:/Users/Tania_2/Internship/R/Results/intersected.lists/intersected.human.up.only.TF.csv",row.names=F, col.names=F)
########################################################################
######Mouse up only TF
setwd("C:/Users/Tania_2/Internship/R/Data")
files_list1 <- list.files("Mm_up_TF", full.names=TRUE)   #list files in the folder
data_files_list<-lapply(files_list1, function(x) as.character(read.table(x)[,1]))   #create a list with the data from all files
datasetnames<-vector()
for (i in 1:length(files_list1)){
        datasetnames[i]<-unlist(strsplit(files_list1[i],'[/.]'))[2]   #parse the names for the list of data
}
names(data_files_list) <- datasetnames

int<-Reduce(intersect,data_files_list) # Reduce applies a function to elements of a vector or list
int
write.table(int,"C:/Users/Tania_2/Internship/R/Results/intersected.lists/intersected.mouse.up.only.TF.csv",row.names=F, col.names=F)

