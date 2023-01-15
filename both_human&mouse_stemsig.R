setwd("C:/Users/Tania_2/Internship/R/Results/ranked.lists")
library(DBI)
library("hom.Hs.inp.db")
# list available mappings:
ls("package:hom.Hs.inp.db")
# have a look at the mapping between human and mouse:
as.list(hom.Hs.inpMUSMU[1:4])
# as it gives mappings from human to mouse orthologs with ensembl protein id, we need to convert from entrez to ensembl protein ID.
 # 1st, get integrative stemness signature for human:
humanlist <- read.table("ranked.list.human.stem.signature.csv",header=T,sep=",", stringsAsFactors = F)
human <- as.character(subset(humanlist[,1], humanlist[,"score"] >= 4))
 # convert to ensembl protein ID:
library(org.Hs.eg.db)
entrez2ensemblprot <- select(org.Hs.eg.db, keys=human, columns=c("ENTREZID","ENSEMBLPROT"), keytype="ENTREZID")
 # note: select' resulted in 1:many mapping between keys and return rows.
 # get orthologs:
mouseortho <- select(hom.Hs.inp.db, keys=entrez2ensemblprot[,2], columns="MUS_MUSCULUS", keytype="HOMO_SAPIENS")
 # convert mouse ensembl protein back to mouse entrez id:
library(org.Mm.eg.db)
ensemblprot2entrez <- select(org.Mm.eg.db, keys=mouseortho[,2], columns=c("ENTREZID","SYMBOL"), keytype="ENSEMBLPROT")
orthologs <- unique(ensemblprot2entrez[,2])
# get integrative signature for mouse:
mouselist <- read.table("ranked.list.mouse.stem.signature.csv",header=T,sep=",", stringsAsFactors = F)
mouse <- as.character(subset(mouselist[,1], mouselist[,"score"] >= 7))
# intersect orthologs with mouse integrative signature:
int <- mouse[mouse%in%orthologs]
int
# retrieve the scores and gene symbols for mouse:
intersectionmouse <- subset(mouselist,mouselist[,1]%in%int)
# retrieve the scores and gene symbols for human:
 # convert lower case symbol from mouse into uppercase symbol from human:
symb <- toupper(intersectionmouse[,"gene_symbol"])
intersectionhuman <- subset(humanlist,(humanlist[,"gene_symbol"])%in%symb)
# merge two dataframes by gene symbol:
intersectionmouse[,2] <- symb
all <- merge.data.frame(intersectionhuman, intersectionmouse, by=2, all = T) 
all2 <- data.frame(all$gene_symbol,tolower(all$gene_symbol),all[,2:ncol(all)], stringsAsFactors = F)
colnames(all2) <- c("human_symbol", "mouse_symbol", "human_id", "human_score", "mouse_id", "mouse_score")
# add column with sum of the two scores:
all3 <- data.frame(all2,"Overalscore"= all2$human_score + all2$mouse_score, stringsAsFactors = F)
# order dataframe:
final <- all3[order(all3$Overalscore, decreasing = T),]
write.table(final,"both_humanmouse_stemsig.csv",sep=",", col.names = T, row.names = F)
