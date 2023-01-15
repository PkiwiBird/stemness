## Human
setwd ("C:/Users/Tania_2/Internship/R/cytoscape_graph_analysis/human")
library(xlsx)
cluster_genes <- read.xlsx("human_cluster.xlsx", sheetName= "humancluster", as.data.frame=T, header=F, colClasses = c("character", "character"))
cluster_genes[,1] <- as.character(cluster_genes[,1]) # turn factor variable as character
cluster <- sapply(cluster_genes[,1], function(x) unlist(strsplit(x," "))[2]) # remove the cluster word
cluster_genes[,1] <- cluster
names(cluster_genes) <- c("Cluster","Gene")
# remove non-significant clusters (we saw with clusterone of cytoscape that until cluster 11 is p value below 0.05)
cluster_genes[,1] <- as.integer(cluster_genes[,1])# 1st convert cluster into number
cluster_genes <- subset(cluster_genes, cluster_genes[,1] < 12)
# get the ranking of each gene:
rank <- read.xlsx("C:/Users/Tania_2/Internship/R/Results/ranked.lists/ranked.list.human.stem.signature.xlsx", sheetName= "ranked.list.human.stem.signatur", as.data.frame = T, header=T, colClasses = c("character","character","integer"))
rank <- rank[,2:3]
# merge the info:
m <- merge(cluster_genes,rank, by.x = 2, by.y = 1)
# calculate mean score:
as <- list()
for (i in 1:max(m$Cluster)){
        
        sub <- subset(m, m$Cluster==i)
        mean1 <- mean(sub$Score)
        total <- sum(sub$Score)
        mean_c <- data.frame(paste("Cluster", i, ""),mean1,total)
        as <- rbind(as,mean_c)
}
colnames(as) <- c("Cluster", "Mean_score","Total_score")
# order by total score
order(as)
write.table(as, "C:/Users/Tania_2/Internship/R/cytoscape_graph_analysis/human/human_cluster_mean.txt", col.names = T, sep='\t', row.names = F, quote = F)

## Mouse
setwd ("C:/Users/Tania_2/Internship/R/cytoscape_graph_analysis/mouse")
library(xlsx)
cluster_genes <- read.xlsx("C:/Users/Tania_2/Internship/R/cytoscape_graph_analysis/mouse/mouse_cluster.xlsx", sheetName= "mousecluster", as.data.frame=T, header=F, colClasses = c("character", "character"))
cluster_genes[,1] <- as.character(cluster_genes[,1]) # turn factor variable as character
cluster <- sapply(cluster_genes[,1], function(x) unlist(strsplit(x," "))[2]) # remove the cluster word
cluster_genes[,1] <- cluster
names(cluster_genes) <- c("Cluster","Gene")
# remove non-significant clusters (we saw with clusterone of cytoscape that until cluster 17 is p value below 0.05)
cluster_genes[,1] <- as.integer(cluster_genes[,1])# 1st convert cluster into number
cluster_genes <- subset(cluster_genes, cluster_genes[,1] < 17)
# get the ranking of each gene:
rank <- read.xlsx("C:/Users/Tania_2/Internship/R/Results/ranked.lists/ranked.list.mouse.stem.signature.xlsx", sheetName= "ranked.list.mouse.stem.signatur", as.data.frame = T, header=T, colClasses = c("character","character","integer"))
rank <- rank[,2:3]
# merge the info:
m <- merge(cluster_genes,rank, by.x = 2, by.y = 1)
# calculate mean score:
as <- list()
for (i in 1:max(m$Cluster)){
        
        sub <- subset(m, m$Cluster==i)
        mean1 <- mean(sub$Score)
        total <- sum(sub$Score)
        mean_c <- data.frame(paste("Cluster", i, ""),mean1,total)
        as <- rbind(as,mean_c)
}
colnames(as) <- c("Cluster", "Mean_score","Total_score")
# order by total score
order(as)
write.table(as, "C:/Users/Tania_2/Internship/R/cytoscape_graph_analysis/mouse/mouse_cluster_mean.txt", col.names = T, sep='\t', row.names = F, quote = F)






