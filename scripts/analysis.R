# Analysis for revision of manuscript
library(stringr)
library(readxl)
library(cowplot)
library(clusterProfiler)

setwd("C:/Users/Matthias/OneDrive - Universidade do Algarve/Ambiente de Trabalho/Manus/stemness/Genes/analysis")

# READING DATA
human_data <- read_xlsx("Human.xlsx")
mouse_data <- read_xlsx("Mouse.xlsx")
sets <- read_xlsx("Sets.xlsx")
sets$Set <- str_trim(sets$Set)

human <- as.data.frame(human_data[,4:dim(human_data)[[2]]]) 
dimnames(human)[[1]] <- as.data.frame(human_data[,2])[,1]
mouse <- as.data.frame(mouse_data[,4:dim(mouse_data)[[2]]]) 
dimnames(mouse)[[1]] <- as.data.frame(mouse_data[,2])[,1]

human_scores <- data.frame(Overall = human_data$Score)
dimnames(human_scores)[[1]] <- as.data.frame(human_data[,2])[,1]
mouse_scores <- data.frame(Overall = mouse_data$Score)
dimnames(mouse_scores)[[1]] <- as.data.frame(mouse_data[,2])[,1]

human_scores_n <- data.frame(Overall = human_data$Score/dim(human)[[2]])
dimnames(human_scores_n)[[1]] <- as.data.frame(human_data[,2])[,1]
mouse_scores_n <- data.frame(Overall = mouse_data$Score/dim(mouse)[[2]])
dimnames(mouse_scores_n)[[1]] <- as.data.frame(mouse_data[,2])[,1]

############################################################################
# CLASSIFICATIONS 
# pluripotent vs multipotent

pluri_human <- sets$Set[sets$Class =="Pluripotent" & sets$Organism=="Human"] 
pluri_mouse <- sets$Set[sets$Class =="Pluripotent" & sets$Organism=="Mouse"] 
multi_human <- sets$Set[sets$Class =="Multipotent" & sets$Organism=="Human"]
npluri_human <- sets$Set[sets$Class !="Pluripotent" & sets$Organism=="Human"]
multi_mouse <- sets$Set[sets$Class =="Multipotent" & sets$Organism=="Mouse"]

nsc_mouse <- sets$Set[sets$Type2 =="NSC" & sets$Organism=="Mouse"]
hsc_mouse <- sets$Set[sets$Type2 =="HSC" & sets$Organism=="Mouse"]


# expression vs non expression 
ex_human <- sets$Set[sets$Evidence =="Expression" & sets$Organism=="Human"]
ex_mouse <- sets$Set[sets$Evidence =="Expression" & sets$Organism=="Mouse"]
nex_human <- sets$Set[sets$Evidence !="Expression" & sets$Organism=="Human"]
nex_mouse <- sets$Set[sets$Evidence !="Expression" & sets$Organism=="Mouse"]

# pluripotent expression 
pluri_ex_human <- sets$Set[sets$Class =="Pluripotent" & sets$Evidence =="Expression" & sets$Organism=="Human"]
pluri_ex_mouse <- sets$Set[sets$Class =="Pluripotent" & sets$Evidence =="Expression" & sets$Organism=="Mouse"]
pluri_nex_human <- sets$Set[sets$Class =="Pluripotent" & sets$Evidence !="Expression" & sets$Organism=="Human"]
pluri_nex_mouse <- sets$Set[sets$Class =="Pluripotent" & sets$Evidence !="Expression" & sets$Organism=="Mouse"]

#####################################################
## SCORES
human_scores <- data.frame(human_scores, 
                           Pluri =  apply(human[, pluri_human], 1, sum), 
                           Multi =  apply(human[, multi_human], 1, sum), 
                           NPluri =  apply(human[, npluri_human], 1, sum), 
                           Ex =   apply(human[, ex_human], 1, sum),
                           Nex =   apply(human[, nex_human], 1, sum),
                           Pluri_ex = apply(human[, pluri_ex_human], 1, sum), 
                           Pluri_nex = apply(human[, pluri_nex_human], 1, sum))


human_scores_n <- data.frame(human_scores_n, 
                           Pluri =  apply(human[, pluri_human], 1, sum)/length(pluri_human), 
                           Multi =  apply(human[, multi_human], 1, sum)/length(multi_human), 
                           NPluri =  apply(human[, npluri_human], 1, sum)/length(npluri_human), 
                           Ex =   apply(human[, ex_human], 1, sum)/length(ex_human),
                           Nex =   apply(human[, nex_human], 1, sum)/length(nex_human),
                           Pluri_ex = apply(human[, pluri_ex_human], 1, sum)/length(pluri_ex_human), 
                           Pluri_nex = apply(human[, pluri_nex_human], 1, sum)/length(pluri_nex_human))

######

mouse_scores <- data.frame(mouse_scores, 
                           Pluri =  apply(mouse[, pluri_mouse], 1, sum), 
                           Multi =  apply(mouse[, multi_mouse], 1, sum), 
                           NSC  =  apply(mouse[, nsc_mouse], 1, sum), 
                           HSC  =  apply(mouse[, hsc_mouse], 1, sum), 
                           Ex =   apply(mouse[, ex_mouse], 1, sum),
                           Nex =   apply(mouse[, nex_mouse], 1, sum),
                           Pluri_ex = apply(mouse[, pluri_ex_mouse], 1, sum), 
                           Pluri_nex = apply(mouse[, pluri_nex_mouse], 1, sum))

mouse_scores_n <- data.frame(mouse_scores_n, 
                             Pluri =  apply(mouse[, pluri_mouse], 1, sum)/length(pluri_mouse), 
                             Multi =  apply(mouse[, multi_mouse], 1, sum)/length(multi_mouse), 
                             NSC  =  apply(mouse[, nsc_mouse], 1, sum)/length(nsc_mouse), 
                             HSC  =  apply(mouse[, hsc_mouse], 1, sum)/length(hsc_mouse), 
                             Ex =   apply(mouse[, ex_mouse], 1, sum)/length(ex_mouse),
                             Nex =   apply(mouse[, nex_mouse], 1, sum)/length(nex_mouse),
                             Pluri_ex = apply(mouse[, pluri_ex_mouse], 1, sum)/length(pluri_ex_mouse), 
                             Pluri_nex = apply(mouse[, pluri_nex_mouse], 1, sum)/length(pluri_nex_mouse))


##############################################################################################################
### CLUSTERING OF  STEMNESS SIGNATURES  

library(pheatmap)
#pheatmap(mouse) # no clear structures
#pheatmap(human) # no clear structures 

## human

tmp <- merge(data.frame(Set=dimnames(human)[[2]]),sets, by = "Set", sort= FALSE, all.x =TRUE)
sample_human = data.frame(Class = tmp$Class, Evidence = tmp$Evidence)
row.names(sample_human) <- dimnames(human)[[2]]
ann_colors_h = list(Class = c(Mixed="mediumpurple1", Multipotent="red1", Pluripotent="royalblue2"), 
                  Evidence = c(Computational="gold1",Literature="goldenrod1", Expression="goldenrod3", RNAi="goldenrod4"))

indexS_human <- order(sample_human$Class)
indexS_human <- c(indexS_human[6:11],indexS_human[1:5], indexS_human[12:21])

h1 <- pheatmap(head(human[,indexS_human],20), annotation_col = sample_human, cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_colors = ann_colors_h, legend = FALSE)

h1pluri <- pheatmap(head(human[order(-human_scores$Pluri),indexS_human],30), cluster_rows = FALSE,  cluster_cols = FALSE,
         annotation_col = sample_human, annotation_colors = ann_colors_h, legend =  FALSE)

pheatmap(head(human[order(-human_scores$NPluri),indexS_human],20), cluster_rows = FALSE,  cluster_cols = FALSE,
         annotation_col = sample_human, annotation_colors = ann_colors_h, legend =  FALSE)


h1multi <- pheatmap(head(human[order(-human_scores$Multi),indexS_human],30), cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_col = sample_human, annotation_colors = ann_colors_h, legend =  FALSE)

index <- human_scores_n$Pluri > 0.5 & human_scores_n$Multi < 0.5

indexG <- human_scores_n$Pluri > 0.3 & human_scores$Multi < 1
sum(indexG) # 61 genes 
entrezPs <- human_data$Gene_id[indexG]

tmp <- human_scores[indexG,]
tmp <- tmp[order(-tmp$Pluri),]

pheatmap(head(human[dimnames(tmp)[[1]],indexS_human],30),annotation_col = sample_human, cluster_rows = FALSE, 
         annotation_colors = ann_colors_h,
         cluster_cols = FALSE, legend = FALSE)

dimnames(human_scores)[[1]][indexG]

indexG <- human_scores$Pluri < 1 & human_scores_n$Multi > 0.3
sum(indexG) # 100 genes 
entrezMs <- human_data$Gene_id[indexG]

tmp <- human_scores[indexG,]
tmp <- tmp[order(-tmp$Multi),]

pheatmap(head(human[dimnames(tmp)[[1]],indexS_human],30),annotation_col = sample_human, cluster_rows = FALSE, 
         annotation_colors = ann_colors_h,
         cluster_cols = FALSE, legend = FALSE)

dimnames(human_scores)[[1]][indexG]


pheatmap(head(human[indexG,indexS_human],30),annotation_col = sample_human, cluster_rows = FALSE, 
         annotation_colors = ann_colors_h,
         cluster_cols = FALSE, legend = FALSE)

# m <- dimnames(human[order(-human_scores$Multi),])[[1]][1:100]
# p <- dimnames(human[order(-human_scores$Pluri),])[[1]][1:100]
# 
# v <- venn.diagram(list(Multipotent = m, Pluripotent=p),
#                   fill = c("red1", "royalblue2"),
#                   alpha = c(0.5, 0.5), cat.cex = 0, cex=2,
#                   filename=NULL)#, area.vector = c(1,1,0.5), direct.area = TRUE )
# 
# 
# grid.newpage()
# grid.draw(v)
# 
# 
# v[[5]]$label <- v[[6]]$label <- "91 genes"
# 
# v[[8]]$label <- v[[9]]$label <- ""
# 
# v[[7]]$label <- paste(intersect(m, p), collapse="\n")  
# 
# grid.newpage()
# grid.draw(v)

####### mouse 
tmp <- merge(data.frame(Set=dimnames(mouse)[[2]]),sets, by = "Set", sort= FALSE, all.x =TRUE)
samples_mouse = data.frame(Class = tmp$Class, Evidence = tmp$Evidence)
row.names(samples_mouse) <- dimnames(mouse)[[2]]
ann_colors_mouse = list(Class = c(Multipotent="red1", Pluripotent="royalblue2", Unipotent = "royalblue4"), 
                  Evidence = c(Literature="goldenrod1", Expression="goldenrod3", RNAi="goldenrod4"))

indexS_mouse <- order(samples_mouse$Class)
indexS_mouse <- c(indexS_mouse[21],indexS_mouse[1:20])


m1 <- pheatmap(head(mouse[, indexS_mouse],20), annotation_col = samples_mouse,  
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_colors = ann_colors_mouse, legend = FALSE)

m3 <- pheatmap(head(mouse[order(-mouse_scores$Multi), indexS_mouse],30), annotation_col = samples_mouse,  
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_colors = ann_colors_mouse, legend = FALSE)


m2 <- pheatmap(head(mouse[order(-mouse_scores$Pluri), indexS_mouse],30), annotation_col = samples_mouse,
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_colors = ann_colors_mouse, legend = FALSE)


pheatmap(head(mouse[order(-mouse_scores$NSC), indexS_mouse],20), annotation_col = samples_mouse,
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_colors = ann_colors_mouse, legend = FALSE)


pheatmap(head(mouse[order(-mouse_scores$HSC), indexS_mouse],20), annotation_col = samples_mouse,
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_colors = ann_colors_mouse, legend = FALSE)


#pluripotency specific mouse genes 
indexG <- mouse_scores_n$Pluri > 0.3 & mouse_scores$Multi < 1
entrezPms <- mouse_data$Gene_id[indexG]
sum(indexG) #34
tmp <- mouse_scores[indexG,]
tmp <- tmp[order(-tmp$Pluri),]


pheatmap(head(mouse[dimnames(tmp)[[1]], indexS_mouse],30), annotation_col = samples_mouse,
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_colors = ann_colors_mouse, legend = FALSE)

#multipotency specific mouse genes 
indexG <- mouse_scores$Pluri <1 & mouse_scores_n$Multi > 0.3
entrezMms <- mouse_data$Gene_id[indexG]
tmp <- mouse_scores[indexG,]
tmp <- tmp[order(-tmp$Multi),]
sum(indexG) # 175

pheatmap(head(mouse[dimnames(tmp)[[1]], indexS_mouse],30), annotation_col = samples_mouse,
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_colors = ann_colors_mouse, legend = FALSE)




cowplot::plot_grid(NULL, h1$gtable,NULL,NULL,  m1$gtable, rel_widths = c(0.05,1,0.05,0.05,1),
                   labels = c('A', '','','B',''), nrow = 1, label_size = 20)

# 
# 
# m <- dimnames(mouse[order(-mouse_scores$Multi),])[[1]][1:100]
# p <- dimnames(mouse[order(-mouse_scores$Pluri),])[[1]][1:100]
# 
# v <- venn.diagram(list(Multipotent = m, Pluripotent=p),
#                   fill = c("red1", "royalblue2"),
#                   alpha = c(0.5, 0.5), cat.cex = 0, cex=2,
#                   filename=NULL) #, area.vector = c(2,2,0.5), direct.area = TRUE )
# 
# 
# grid.newpage()
# grid.draw(v)
# 
# v[[5]]$label <- v[[6]]$label <- "81 genes"
# 
# v[[8]]$label <- v[[9]]$label <- ""
# 
# v[[7]]$label <- paste(intersect(m, p), collapse="\n")  
# 
# grid.newpage()
# grid.draw(v)
# 
# grid.newpage()
# grid.draw(v)


## COMPARATIVE FUNCTIONAL ENRICHMENT ANALYSIS 

# HUMAN 
## highest ranking genes
N <- 200
indexG <- order(-human_scores_n$Pluri)
entrezP <- human_data$Gene_id[indexG[1:N]]
indexG <- order(-human_scores_n$Multi)
entrezM <- human_data$Gene_id[indexG[1:N]]


genell <- list(Pluripotent = entrezP, Multipotent = entrezM)
xxKh <- compareCluster(genell,fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)
dotplot(xxKh, title = "KEGG enrichment analysis - Human genes")

xxGh <- compareCluster(genell,fun="enrichGO",
                     OrgDb= "org.Hs.eg.db", ont="BP", pvalueCutoff=0.05)
dotplot(xxGh, title= "GO enrichment analysis - Human genes")

### SPECIFIC 
genell <- list(Pluripotent = entrezPs, Multipotent = entrezMs)
xxKhs <- compareCluster(genell,fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)
dotplot(xxKhs,  title = "KEGG enrichment analysis - Human genes")

xxGhs <- compareCluster(genell,fun="enrichGO",
                     OrgDb= "org.Hs.eg.db", ont="BP", pvalueCutoff=0.05)
dotplot(xxGhs,  title = "GO enrichment analysis -  Human genes")

# Mouse
N <- 200
indexG <- order(-mouse_scores_n$Pluri)
entrezP <- mouse_data$Gene_id[indexG[1:N]]
indexG <- order(-mouse_scores_n$Multi)
entrezM <- mouse_data$Gene_id[indexG[1:N]]


genell <- list(Pluripotent = entrezP, Multipotent = entrezM)
xxKm <- compareCluster(genell,fun="enrichKEGG",
                     organism="mmu", pvalueCutoff=0.05)
dotplot(xxKm, title = "KEGG enrichment analysis - Murine genes" )

xxGm <- compareCluster(genell,fun="enrichGO",
                     OrgDb= "org.Mm.eg.db", ont="BP", pvalueCutoff=0.05)
dotplot(xxGm,  title = "GO enrichment analysis - Murine genes")

### SPECIFIC 
genell <- list(Pluripotent = entrezPms, Multipotent = entrezMms)
xxKms <- compareCluster(genell,fun="enrichKEGG",
                     organism="mmu", pvalueCutoff=0.05)
dotplot(xxKms, title = "KEGG enrichment analysis - Murine genes")

xxGms <- compareCluster(genell,fun="enrichGO",
                     OrgDb= "org.Mm.eg.db", ont="BP", pvalueCutoff=0.05)
dotplot(xxGms, title="GO functional enrichment - Murine genes")


#### POLYCOMB COMPLEX 

pg_human <- c("PHC1", "PHC2", "PHC3", "CBX2", "CBX4","CBX6","CBX7", "CBX8", "SCMH1", "SCML2", 
"RING1","RING1A", "RING1B",  
"PCGF1", "PCGF2", "PCGF3", "PCGF4","PCGF5", "PCGF6",  
"EZH1","EZH2", 
"EED", "SUZ12")

pg_mouse <- str_to_title(pg_human)

pg_human <- pg_human[pg_human %in% dimnames(human_scores)[[1]]]
pg_mouse <- pg_mouse[pg_mouse %in% dimnames(mouse_scores)[[1]]]

pheatmap(head(human[pg_human,indexS_human],20), annotation_col = sample_human, cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_colors = ann_colors_h, legend = FALSE)

pheatmap(head(mouse[pg_mouse , indexS_mouse],20), annotation_col = samples_mouse,  
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_colors = ann_colors_mouse, legend = FALSE)
