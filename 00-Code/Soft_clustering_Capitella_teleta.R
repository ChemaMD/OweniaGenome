library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(corrplot)
library(RColorBrewer)
library(Mfuzz)
library(ComplexHeatmap)

#######################################################################
########################  1. Fuzzy clustering #########################
#######################################################################

df_norm <- read.table("02-Capitella_teleta_DESeq2_replicates.txt", header = TRUE)
colnames(df_norm) <- c("Gene_ID","64_cells_1","64_cells_2","gastrula_1","gastrula_2",
                      "st4tt_1","st4tt_2","st5_1","st5_2","st7_1","st7_2",
                      "pre-competent_1","pre-competent_2","pre-competent_3",
                      "pre-competent_4","competent_1","competent_2",
                      "competent_3","competent4")

df_average_raw <- df_norm[-c(1)]
df_average_raw$blastula <- rowMeans(df_average_raw[,c(1,2)], na.rm = TRUE)
df_average_raw$gastrula <- rowMeans(df_average_raw[,c(3,4)], na.rm = TRUE)
df_average_raw$st4tt <- rowMeans(df_average_raw[,c(5,6)], na.rm = TRUE)
df_average_raw$st5 <- rowMeans(df_average_raw[,c(7,8)], na.rm = TRUE)
df_average_raw$st7 <- rowMeans(df_average_raw[,c(9,10)], na.rm = TRUE)
df_average_raw$precompetent <- rowMeans(df_average_raw[,c(11:14)], na.rm = TRUE)
df_average_raw$competent <- rowMeans(df_average_raw[,c(15:18)], na.rm = TRUE)
df_average <- df_average_raw[,c(19:ncol(df_average_raw))]

df_average_to_export <- df_average
rownames(df_average_to_export) <- t(df_norm[1])
write.table(df_average_to_export, "02-Capitella_teleta_DESeq2_average.txt", sep = '\t', quote = F)

mat_average <- as.matrix(df_average)
eset_rna <- new('ExpressionSet', exprs = mat_average)
# 41,221 genes
clean_eset_rna <- filter.std(eset_rna,min.std=0) 
# 1,407 genes with null expression at all stages: remaining 39,814
norm_eset_rna <- standardise(clean_eset_rna)
fuzzifier <- round(mestimate(norm_eset_rna),2)
# fuzzifier: 1.54
centroids <- Dmin(norm_eset_rna, fuzzifier, crange = seq(4,20,1), repeats=3, visu=TRUE)
# optimal centroid with 12 clusters
cluster_number <- 12
# Run the two following lines together to ensure reproducibility
set.seed(123)
clusters <- mfuzz(norm_eset_rna, c = 12, m = fuzzifier)
clusters$size 
# 2976 2889 5377 2770 6137 2315 2846 2206 2243 4830 2477 2748

df_names <- df_average
rownames(df_names) <- t(df_norm[1])
df_names_notnull <- df_names[apply(df_names, 1, function(x) sd(x, na.rm = TRUE)!=0),]

clusters_annotated_unordered <- cbind(rownames(df_names_notnull),clusters$cluster)
colnames(clusters_annotated_unordered) <- c("Gene_ID","Cluster_raw")
write.table(clusters_annotated_unordered,"03-Capitella_teleta_clusters_annotation_raw.txt", 
            quote = F, sep = '\t', row.names = FALSE)

clusters_expression_unordered <- cbind(rownames(df_names_notnull),df_names_notnull,clusters$cluster)
colnames(clusters_expression_unordered) <- c("Gene_ID","blastula","gastrula","st4tt",
                                             "st5","st7","pre-competent","competent","Cluster_raw")
write.table(clusters_expression_unordered,"04-Capitella_teleta_clusters_expression_raw.txt", 
            quote = F, sep = '\t', row.names = FALSE)


#######################################################################
########################  2. Clustered heatmap ########################
#######################################################################

df_raw <- read.table("04-Capitella_teleta_clusters_expression_raw.txt", header = TRUE)
df_zscore <- data.frame(cbind(t(scale(t(df_raw[-c(1,9)]))),'cluster_raw'=df_raw$Cluster_raw))
df_zscore$cluster_raw <- factor(df_zscore$cluster_raw, 
                                levels = c("12","5","8","2","11","7","6","9","4","3","10","1"))

df_sorted <- df_zscore[order(df_zscore$cluster_raw),]

heatmap_color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

ComplexHeatmap::Heatmap(as.matrix(df_sorted[-c(8)]),
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        show_row_names = FALSE,
                        col = heatmap_color)

#re-write clusters
# equivalence is the following (old:new)
# 1:12, 2:4, 3:10, 4:9, 5:2, 6:7, 7:6, 8:3, 9:8, 10:11, 11:5. 12:1
df_reordered <- df_raw
df_reordered$Cluster_corrected <- c(rep(0,nrow(df_reordered)))
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 1] <- 12
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 2] <- 4
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 3] <- 10
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 4] <- 9
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 5] <- 2
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 6] <- 7
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 7] <- 6
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 8] <- 3
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 9] <- 8
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 10] <- 11
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 11] <- 5
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 12] <- 1

df_reordered_clean <- df_reordered[-c(9)]

write.table(df_reordered_clean[c(1,9)],"06-Capitella_teleta_clusters_annotation_corrected.txt", 
            quote = F, sep = '\t', row.names = FALSE)
write.table(df_reordered_clean,"07-Capitella_teleta_clusters_expression_corrected.txt", 
            quote = F, sep = '\t', row.names = FALSE)


#######################################################################
########################  3. LOESS smoothing ##########################
#######################################################################


library(ggplot2)
library(ggpubr)
library(gplots)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(matrixStats)
library(data.table)
library(mgcv)
library(scales)
library(viridis)

palette <- viridis_pal(option = "C", direction = 1)(30)
palette2 <- gsub('.{2}$', '', palette)
write.csv(rbind(palette2), file = 'palette.txt', row.names = FALSE)

palette <- viridis_pal(option = "C", direction = 1)(12)
time <- c(1:7)

write.table(c(palette), "palette.txt", sep = ',', quote = TRUE)

df_raw <- read.table("07-Capitella_teleta_clusters_expression_corrected.txt", header = T)

# For cluster 1
df_cluster1_raw <- df_raw[df_raw$Cluster_corrected == 12,-c(1,9)] # change 1 for each cluster and run all lines
mat_cluster1_raw <- as.matrix(df_cluster1_raw) 
mat_cluster1_norm <- t(scale(t(mat_cluster1_raw))) 
df_cluster1_norm <- data.frame(mat_cluster1_norm)
df_cluster1_clean <- data.frame(cbind(time,t(df_cluster1_norm))) 
df_cluster1_clean_long <- gather(df_cluster1_clean, gene_ID, expression, -time)

summary_cluster1 <- data.frame(time=df_cluster1_clean$time, 
                               n=tapply(df_cluster1_clean_long$expression, 
                                        df_cluster1_clean_long$time, length), 
                               mean=tapply(df_cluster1_clean_long$expression, 
                                           df_cluster1_clean_long$time, mean))

# Plot data
cluster12_plot <- ggplot(summary_cluster1, aes(x=time, y=mean)) + # change cluster1_plot to clusterX_plot for each cluster
  geom_line(data = df_cluster1_clean_long, aes(x=time, y=expression, group=gene_ID), 
            color="gray") +
  stat_smooth(colour = palette[12], # change 1 to match the number of the cluster, so that each cluster has its own colour of the continuous palette
              fill = palette[12], # change 1 to match the number of the cluster, so that each cluster has its own colour of the continuous palette
              size = 2, alpha = 0.5) +  
  geom_hline(yintercept = 0, linetype="dashed", color = "black") + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(y = "Normalised gene expression (z-score)", x = "developmental stage") +
  scale_y_continuous(breaks = c(-3:3), limits = c(-8,8)) +
  scale_x_continuous(breaks = c(1:7), expand = c(0,0), limits = c(0.8,7.2),
                     labels = c("blastula", "gastrula","st4tt", "st5", 
                                "st7", "pre-competent", "competent")) +
  theme(plot.title = element_text(color="black",size=14,hjust=0.5),
        axis.title.x = element_text(color="black",size=14),
        axis.title.y = element_text(color="black",size=14),
        axis.text.x = element_text(color="black",size=12,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=12,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"))


allfigures <- ggarrange(cluster1_plot, cluster2_plot, cluster3_plot, cluster4_plot,
                        cluster5_plot, cluster6_plot, cluster7_plot, cluster8_plot,
                        cluster9_plot, cluster10_plot, cluster11_plot, cluster12_plot,
                        ncol = 3, nrow = 4)


allfigures +
  ggsave("<name-of-file>.pdf", width = 25, height = 35, units = "cm")



#######################################################################
################  4. GO terms enrichment analysis #####################
#######################################################################


## PART 1. Generating the GO terms

library(ggplot2)
library(ggpubr)
library(topGO)

geneID2GO <- readMappings(file = "09-Capitella_teleta_geneID2GO_topGO_onlyGOgenes.txt")
geneUniverse <- names(geneID2GO)

df_clusters <- read.table("06-Capitella_teleta_clusters_annotation_corrected.txt", header = T)
df_clusters <- df_clusters[-c(1),] #remove post1 as it does not have GO terms annotated

cluster_1 <- df_clusters[df_clusters$Cluster_corrected == 1,]
cluster_2 <- df_clusters[df_clusters$Cluster_corrected == 2,]
cluster_3 <- df_clusters[df_clusters$Cluster_corrected == 3,]
cluster_4 <- df_clusters[df_clusters$Cluster_corrected == 4,]
cluster_5 <- df_clusters[df_clusters$Cluster_corrected == 5,]
cluster_6 <- df_clusters[df_clusters$Cluster_corrected == 6,]
cluster_7 <- df_clusters[df_clusters$Cluster_corrected == 7,]
cluster_8 <- df_clusters[df_clusters$Cluster_corrected == 8,]
cluster_9 <- df_clusters[df_clusters$Cluster_corrected == 9,]
cluster_10 <- df_clusters[df_clusters$Cluster_corrected == 10,]
cluster_11 <- df_clusters[df_clusters$Cluster_corrected == 11,]
cluster_12 <- df_clusters[df_clusters$Cluster_corrected == 12,]

cluster_1_names <- as.character(cluster_1$Gene_ID)
cluster_1_list_for_GO <- factor(as.integer(geneUniverse %in% cluster_1_names))
names(cluster_1_list_for_GO) <- geneUniverse

cluster_2_names <- as.character(cluster_2$Gene_ID)
cluster_2_list_for_GO <- factor(as.integer(geneUniverse %in% cluster_2_names))
names(cluster_2_list_for_GO) <- geneUniverse

cluster_3_names <- as.character(cluster_3$Gene_ID)
cluster_3_list_for_GO <- factor(as.integer(geneUniverse %in% cluster_3_names))
names(cluster_3_list_for_GO) <- geneUniverse

cluster_4_names <- as.character(cluster_4$Gene_ID)
cluster_4_list_for_GO <- factor(as.integer(geneUniverse %in% cluster_4_names))
names(cluster_4_list_for_GO) <- geneUniverse

cluster_5_names <- as.character(cluster_5$Gene_ID)
cluster_5_list_for_GO <- factor(as.integer(geneUniverse %in% cluster_5_names))
names(cluster_5_list_for_GO) <- geneUniverse

cluster_6_names <- as.character(cluster_6$Gene_ID)
cluster_6_list_for_GO <- factor(as.integer(geneUniverse %in% cluster_6_names))
names(cluster_6_list_for_GO) <- geneUniverse

cluster_7_names <- as.character(cluster_7$Gene_ID)
cluster_7_list_for_GO <- factor(as.integer(geneUniverse %in% cluster_7_names))
names(cluster_7_list_for_GO) <- geneUniverse

cluster_8_names <- as.character(cluster_8$Gene_ID)
cluster_8_list_for_GO <- factor(as.integer(geneUniverse %in% cluster_8_names))
names(cluster_8_list_for_GO) <- geneUniverse

cluster_9_names <- as.character(cluster_9$Gene_ID)
cluster_9_list_for_GO <- factor(as.integer(geneUniverse %in% cluster_9_names))
names(cluster_9_list_for_GO) <- geneUniverse

cluster_10_names <- as.character(cluster_10$Gene_ID)
cluster_10_list_for_GO <- factor(as.integer(geneUniverse %in% cluster_10_names))
names(cluster_10_list_for_GO) <- geneUniverse

cluster_11_names <- as.character(cluster_11$Gene_ID)
cluster_11_list_for_GO <- factor(as.integer(geneUniverse %in% cluster_11_names))
names(cluster_11_list_for_GO) <- geneUniverse

cluster_12_names <- as.character(cluster_12$Gene_ID)
cluster_12_list_for_GO <- factor(as.integer(geneUniverse %in% cluster_12_names))
names(cluster_12_list_for_GO) <- geneUniverse


cluster_1_GOdata_BP <- new("topGOdata", description="cluster_1_program",
                           ontology="BP", allGenes=cluster_1_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) 
cluster_1_resultFisher_BP <- runTest(cluster_1_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_cluster_1 <- GenTable(cluster_1_GOdata_BP, classicFisher = cluster_1_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30)

cluster_2_GOdata_BP <- new("topGOdata", description="cluster_2_program",
                           ontology="BP", allGenes=cluster_2_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) 
cluster_2_resultFisher_BP <- runTest(cluster_2_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_cluster_2 <- GenTable(cluster_2_GOdata_BP, classicFisher = cluster_2_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30) 

cluster_3_GOdata_BP <- new("topGOdata", description="cluster_3_program",
                           ontology="BP", allGenes=cluster_3_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) 
cluster_3_resultFisher_BP <- runTest(cluster_3_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_cluster_3 <- GenTable(cluster_3_GOdata_BP, classicFisher = cluster_3_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30)

cluster_4_GOdata_BP <- new("topGOdata", description="cluster_4_program",
                           ontology="BP", allGenes=cluster_4_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO)
cluster_4_resultFisher_BP <- runTest(cluster_4_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_cluster_4 <- GenTable(cluster_4_GOdata_BP, classicFisher = cluster_4_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30) 

cluster_5_GOdata_BP <- new("topGOdata", description="cluster_5_program",
                           ontology="BP", allGenes=cluster_5_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) 
cluster_5_resultFisher_BP <- runTest(cluster_5_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_cluster_5 <- GenTable(cluster_5_GOdata_BP, classicFisher = cluster_5_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30)

cluster_6_GOdata_BP <- new("topGOdata", description="cluster_6_program",
                           ontology="BP", allGenes=cluster_6_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) 
cluster_6_resultFisher_BP <- runTest(cluster_6_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_cluster_6 <- GenTable(cluster_6_GOdata_BP, classicFisher = cluster_6_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30) 

cluster_7_GOdata_BP <- new("topGOdata", description="cluster_7_program",
                           ontology="BP", allGenes=cluster_7_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) 
cluster_7_resultFisher_BP <- runTest(cluster_7_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_cluster_7 <- GenTable(cluster_7_GOdata_BP, classicFisher = cluster_7_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30)

cluster_8_GOdata_BP <- new("topGOdata", description="cluster_8_program",
                           ontology="BP", allGenes=cluster_8_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) 
cluster_8_resultFisher_BP <- runTest(cluster_8_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_cluster_8 <- GenTable(cluster_8_GOdata_BP, classicFisher = cluster_8_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30) 


cluster_9_GOdata_BP <- new("topGOdata", description="cluster_9_program",
                           ontology="BP", allGenes=cluster_9_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) 
cluster_9_resultFisher_BP <- runTest(cluster_9_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_cluster_9 <- GenTable(cluster_9_GOdata_BP, classicFisher = cluster_9_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30)

cluster_10_GOdata_BP <- new("topGOdata", description="cluster_10_program",
                           ontology="BP", allGenes=cluster_10_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) 
cluster_10_resultFisher_BP <- runTest(cluster_10_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_cluster_10 <- GenTable(cluster_10_GOdata_BP, classicFisher = cluster_10_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30) 

cluster_11_GOdata_BP <- new("topGOdata", description="cluster_11_program",
                           ontology="BP", allGenes=cluster_11_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) 
cluster_11_resultFisher_BP <- runTest(cluster_11_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_cluster_11 <- GenTable(cluster_11_GOdata_BP, classicFisher = cluster_11_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30)

cluster_12_GOdata_BP <- new("topGOdata", description="cluster_12_program",
                           ontology="BP", allGenes=cluster_12_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) 
cluster_12_resultFisher_BP <- runTest(cluster_12_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_cluster_12 <- GenTable(cluster_12_GOdata_BP, classicFisher = cluster_12_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30)

write.table(results_cluster_1, "10-Cluster_1_GO_terms.txt", quote=FALSE, sep='\t', row.names = F) 
write.table(results_cluster_2, "10-Cluster_2_GO_terms.txt", quote=FALSE, sep='\t', row.names = F)
write.table(results_cluster_3, "10-Cluster_3_GO_terms.txt", quote=FALSE, sep='\t', row.names = F) 
write.table(results_cluster_4, "10-Cluster_4_GO_terms.txt", quote=FALSE, sep='\t', row.names = F) 
write.table(results_cluster_5, "10-Cluster_5_GO_terms.txt", quote=FALSE, sep='\t', row.names = F) 
write.table(results_cluster_6, "10-Cluster_6_GO_terms.txt", quote=FALSE, sep='\t', row.names = F)
write.table(results_cluster_7, "10-Cluster_7_GO_terms.txt", quote=FALSE, sep='\t', row.names = F) 
write.table(results_cluster_8, "10-Cluster_8_GO_terms.txt", quote=FALSE, sep='\t', row.names = F) 
write.table(results_cluster_9, "10-Cluster_9_GO_terms.txt", quote=FALSE, sep='\t', row.names = F) 
write.table(results_cluster_10, "10-Cluster_10_GO_terms.txt", quote=FALSE, sep='\t', row.names = F)
write.table(results_cluster_11, "10-Cluster_11_GO_terms.txt", quote=FALSE, sep='\t', row.names = F) 
write.table(results_cluster_12, "10-Cluster_12_GO_terms.txt", quote=FALSE, sep='\t', row.names = F) 

results_clusters <- rbind(cbind(Cluster="cluster 1", results_cluster_1),
                          cbind(Cluster="cluster 2", results_cluster_2),
                          cbind(Cluster="cluster 3", results_cluster_3),
                          cbind(Cluster="cluster 4", results_cluster_4),
                          cbind(Cluster="cluster 5", results_cluster_5),
                          cbind(Cluster="cluster 6", results_cluster_6),
                          cbind(Cluster="cluster 7", results_cluster_7),
                          cbind(Cluster="cluster 8", results_cluster_8),
                          cbind(Cluster="cluster 9", results_cluster_9),
                          cbind(Cluster="cluster 10", results_cluster_10),
                          cbind(Cluster="cluster 11", results_cluster_11),
                          cbind(Cluster="cluster 12", results_cluster_12))

write.table(results_clusters, "10-All_clusters_GO_terms.txt", quote=FALSE, sep='\t', row.names = F)


## PART 2. Plotting bar plots of -log10(p-values) of GO terms for each cluster

library(topGO)
library(ggplot2)
library(ggpubr)
library(cowplot)

df_clusters <- read.table("10-All_clusters_GO_terms.txt", header = TRUE, sep = '\t')
df_clusters$classicFisher[df_clusters$classicFisher == '< 1e-30'] <- 1e-30

results_cluster_1 <- df_clusters[df_clusters$Cluster == 'cluster 1',-c(1)]
results_cluster_2 <- df_clusters[df_clusters$Cluster == 'cluster 2',-c(1)]
results_cluster_3 <- df_clusters[df_clusters$Cluster == 'cluster 3',-c(1)]
results_cluster_4 <- df_clusters[df_clusters$Cluster == 'cluster 4',-c(1)]
results_cluster_5 <- df_clusters[df_clusters$Cluster == 'cluster 5',-c(1)]
results_cluster_6 <- df_clusters[df_clusters$Cluster == 'cluster 6',-c(1)]
results_cluster_7 <- df_clusters[df_clusters$Cluster == 'cluster 7',-c(1)]
results_cluster_8 <- df_clusters[df_clusters$Cluster == 'cluster 8',-c(1)]
results_cluster_9 <- df_clusters[df_clusters$Cluster == 'cluster 9',-c(1)]
results_cluster_10 <- df_clusters[df_clusters$Cluster == 'cluster 10',-c(1)]
results_cluster_11 <- df_clusters[df_clusters$Cluster == 'cluster 11',-c(1)]
results_cluster_12 <- df_clusters[df_clusters$Cluster == 'cluster 12',-c(1)]

# Cluster 1
goEnrichment <- results_cluster_1
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

cluster_1_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Cluster 1") +
  scale_y_continuous(limits=c(0,30),breaks=round(seq(0,30, by = 2), 1)) +
  theme_classic() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, vjust=0.5),
    axis.title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

# Cluster 2
goEnrichment <- results_cluster_2
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

cluster_2_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Cluster 2") +
  scale_y_continuous(limits=c(0,30),breaks=round(seq(0,30, by = 2), 1)) +
  theme_classic() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, vjust=0.5),
    axis.title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

# Cluster 3
goEnrichment <- results_cluster_3
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

cluster_3_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Cluster 3") +
  scale_y_continuous(limits=c(0,30),breaks=round(seq(0,30, by = 2), 1)) +
  theme_classic() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, vjust=0.5),
    axis.title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

# Cluster 4
goEnrichment <- results_cluster_4
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

cluster_4_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Cluster 4") +
  scale_y_continuous(limits=c(0,30),breaks=round(seq(0,30, by = 2), 1)) +
  theme_classic() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, vjust=0.5),
    axis.title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

# Cluster 5
goEnrichment <- results_cluster_5
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

cluster_5_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Cluster 5") +
  scale_y_continuous(limits=c(0,30),breaks=round(seq(0,30, by = 2), 1)) +
  theme_classic() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, vjust=0.5),
    axis.title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

# Cluster 6
goEnrichment <- results_cluster_6
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

cluster_6_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Cluster 6") +
  scale_y_continuous(limits=c(0,30),breaks=round(seq(0,30, by = 2), 1)) +
  theme_classic() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, vjust=0.5),
    axis.title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

# Cluster 7
goEnrichment <- results_cluster_7
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

cluster_7_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Cluster 7") +
  scale_y_continuous(limits=c(0,30),breaks=round(seq(0,30, by = 2), 1)) +
  theme_classic() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, vjust=0.5),
    axis.title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

# Cluster 8
goEnrichment <- results_cluster_8
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

cluster_8_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Cluster 8") +
  scale_y_continuous(limits=c(0,30),breaks=round(seq(0,30, by = 2), 1)) +
  theme_classic() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, vjust=0.5),
    axis.title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

# Cluster 9
goEnrichment <- results_cluster_9
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

cluster_9_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Cluster 9") +
  scale_y_continuous(limits=c(0,30),breaks=round(seq(0,30, by = 2), 1)) +
  theme_classic() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, vjust=0.5),
    axis.title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

# Cluster 10
goEnrichment <- results_cluster_10
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

cluster_10_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Cluster 10") +
  scale_y_continuous(limits=c(0,30),breaks=round(seq(0,30, by = 2), 1)) +
  theme_classic() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, vjust=0.5),
    axis.title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

# Cluster 11
goEnrichment <- results_cluster_11
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

cluster_11_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Cluster 11") +
  scale_y_continuous(limits=c(0,30),breaks=round(seq(0,30, by = 2), 1)) +
  theme_classic() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, vjust=0.5),
    axis.title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

# Cluster 12
goEnrichment <- results_cluster_12
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

cluster_12_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Cluster 12") +
  scale_y_continuous(limits=c(0,30),breaks=round(seq(0,30, by = 2), 1)) +
  theme_classic() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, vjust=0.5),
    axis.title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

# Plot them all
plot_grid(cluster_1_plot, cluster_2_plot, cluster_3_plot, cluster_4_plot,
          cluster_5_plot, cluster_6_plot, cluster_7_plot, cluster_8_plot,
          cluster_9_plot, cluster_10_plot, cluster_11_plot, cluster_12_plot,
          ncol = 2, align = "v") +
  ggsave("<name-of-file>.pdf", height = 60, width = 40, units = "cm")


## PART 3. k-means clustering of GO terms and summary visualisation of
## enriched GO terms in clusters

# Necessary libraries
library(topGO)
library(ggplot2)
library(ggpubr)
library(simplifyEnrichment)
library(ComplexHeatmap)
library(tidyr)
library(RColorBrewer)

Capitella_GO_clusters <- read.table("12-Capitella_teleta_GO_terms.txt", header=TRUE, sep ='\t')
colnames(Capitella_GO_clusters) <- c("Cluster", "GO.ID", "Term", "Annotated", "Significant", 
                                     "Expected", "classicFisher")

Owenia_GO_clusters <- read.table("13-Owenia_fusiformis_GO_terms.txt", header=TRUE, sep =',')
colnames(Owenia_GO_clusters) <- c("Cluster", "GO.ID", "Term", "Annotated", "Significant", 
                                  "Expected", "classicFisher")

go_id <- rbind(Owenia_GO_clusters, Capitella_GO_clusters)
go_id$GO.ID <- as.character(go_id$GO.ID)
mat <- GO_similarity(go_id$GO.ID, ont = "BP")
set.seed(300) # ensures reproducibility
GO_cluster_annotation_raw <- simplifyGO(mat, method = "kmeans")

GO_cluster_annotation_clean <- GO_cluster_annotation_raw[2:ncol(GO_cluster_annotation_raw)]
rownames(GO_cluster_annotation_clean) <- t(GO_cluster_annotation_raw[1]) 

Capitella_GO_clusters$Cluster <- Capitella_GO_clusters$Cluster+c(12)
GO_clusters <- rbind(Owenia_GO_clusters, Capitella_GO_clusters)
GO_clusters$classicFisher[GO_clusters$classicFisher == '< 1e-30'] <- 1e-30
GO_clusters$kmeans_cluster <- GO_cluster_annotation_clean[GO_clusters$GO.ID,2] 
GO_clusters <- GO_clusters[order(GO_clusters$kmeans_cluster, 
                                 GO_clusters$classicFisher),] 

heatmap_data <- GO_clusters
heatmap_data$classicFisher <- as.numeric(heatmap_data$classicFisher)
heatmap_data <- heatmap_data[,c("Cluster","GO.ID","classicFisher")] 

heatmap_data_long <- data.frame(t(spread(heatmap_data, "GO.ID","classicFisher")))
heatmap_data_long <- heatmap_data_long[2:nrow(heatmap_data_long),]
heatmap_data_long <- -log(heatmap_data_long,10)
colnames(heatmap_data_long) <- c(1:24)

heatmap_data_long$kmeans_cluster <- GO_cluster_annotation_clean[rownames(heatmap_data_long),2]
heatmap_data_long$averageFisher <- rowMeans(heatmap_data_long[,c(1:24)], na.rm = TRUE)
heatmap_data_long <- heatmap_data_long[order(heatmap_data_long$kmeans_cluster, 
                                             -heatmap_data_long$averageFisher),]

palette <- colorRampPalette(brewer.pal(n = 7, name = "Reds"))(200)
palette[1] <- rgb(255,255,255, maxColorValue = 255)

mat_heatmap <- data.matrix(heatmap_data_long[-c(25,26)])
mat_heatmap[is.na(mat_heatmap)] <- 0 # NAs get converted to 0s, which in turn will be shown in the heatmap in white
ComplexHeatmap::Heatmap(mat_heatmap, 
                        row_split = heatmap_data_long$kmeans_cluster,
                        border = TRUE,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        show_row_names = FALSE,
                        col = palette)

# I fusioned clusters 6 and 7 into a single one
# In order: transcription (1), amino acid and lipid metabolism (2), glycoconjugate metabolism (3), 
# complex metabolism (4), RNA processing (5), Ion homeostasis (6,7), Motility (8), Intracellular transport (9),
# ion transport (10), cell differentiation, morphogenesis, and organogenesis (11), cellular structures assembly (12),
# response to stimuli (13), cell communication and signalling (14), OXPHOS metabolism (15), Immune system development (16)
