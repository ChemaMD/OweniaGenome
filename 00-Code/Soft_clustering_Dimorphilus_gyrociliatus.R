library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(tidyr)
library(ggalluvial)
library(RColorBrewer)
library(corrplot)
library(factoextra)

##############################################
## PART A. Quality check and basic analyses ##
##############################################


# 1. Correlation scatterplot
df_TPM <- read.table("01-Dimorphilus_gyrociliatus_RNAseq_TPM_replicates.txt", header = TRUE, sep='\t')

early_scatterplot <- ggplot(df_TPM, aes(x=early_development_rep1, y=early_development_rep2)) +
  geom_point(size = 0.6, color=rgb(86,59,69, maxColorValue = 255)) +
  scale_x_log10(limits = c(1,1e5), breaks = c(1,10,100,1000,10000,100000), expand = c(0,0)) +
  scale_y_log10(limits = c(1,1e5), breaks = c(1,10,100,1000,10000,100000), expand = c(0,0)) +
  stat_cor(method = "pearson", label.x = 3, label.y = 30) +
  annotate(x=1e1, y=1e4, 
           label=paste("R = ", round(cor(df_TPM$early_development_rep1, df_TPM$early_development_rep2),2)), 
           geom="text", size=5) +
  theme_classic() +
  theme(axis.title.x = element_text(color="black",size=10),
        axis.title.y = element_text(color="black",size=10),
        axis.text.x = element_text(color="black", size=10,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=10,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"))

late_scatterplot <- ggplot(df_TPM, aes(x=late_development_rep1, y=late_development_rep2)) +
  geom_point(size = 0.6, color=rgb(210,163,89, maxColorValue = 255)) +
  scale_x_log10(limits = c(1,1e5), breaks = c(1,10,100,1000,10000,100000), expand = c(0,0)) +
  scale_y_log10(limits = c(1,1e5), breaks = c(1,10,100,1000,10000,100000), expand = c(0,0)) +
  stat_cor(method = "pearson", label.x = 3, label.y = 30) +
  annotate(x=1e1, y=1e4, 
           label=paste("R = ", round(cor(df_TPM$late_development_rep1, df_TPM$late_development_rep2),2)), 
           geom="text", size=5) +
  theme_classic() +
  theme(axis.title.x = element_text(color="black",size=10),
        axis.title.y = element_text(color="black",size=10),
        axis.text.x = element_text(color="black", size=10,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=10,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"))

hatchling_scatterplot <- ggplot(df_TPM, aes(x=hatchling_rep1, y=hatchling_rep2)) +
  geom_point(size = 0.6, color=rgb(198,97,58, maxColorValue = 255)) +
  scale_x_log10(limits = c(1,1e5), breaks = c(1,10,100,1000,10000,100000), expand = c(0,0)) +
  scale_y_log10(limits = c(1,1e5), breaks = c(1,10,100,1000,10000,100000), expand = c(0,0)) +
  stat_cor(method = "pearson", label.x = 3, label.y = 30) +
  annotate(x=1e1, y=1e4, 
           label=paste("R = ", round(cor(df_TPM$hatchling_rep1, df_TPM$hatchling_rep2),2)), 
           geom="text", size=5) +
  theme_classic() +
  theme(axis.title.x = element_text(color="black",size=10),
        axis.title.y = element_text(color="black",size=10),
        axis.text.x = element_text(color="black", size=10,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=10,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"))

adult_scatterplot <- ggplot(df_TPM, aes(x=female_adult_rep1, y=female_adult_rep2)) +
  geom_point(size = 0.6, color=rgb(11,108,160, maxColorValue = 255)) +
  scale_x_log10(limits = c(1,1e5), breaks = c(1,10,100,1000,10000,100000), expand = c(0,0)) +
  scale_y_log10(limits = c(1,1e5), breaks = c(1,10,100,1000,10000,100000), expand = c(0,0)) +
  stat_cor(method = "pearson", label.x = 3, label.y = 30) +
  annotate(x=1e1, y=1e4, 
           label=paste("R = ", round(cor(df_TPM$female_adult_rep1, df_TPM$female_adult_rep2),2)), 
           geom="text", size=5) +
  theme_classic() +
  theme(axis.title.x = element_text(color="black",size=10),
        axis.title.y = element_text(color="black",size=10),
        axis.text.x = element_text(color="black", size=10,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=10,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"))

ggarrange(early_scatterplot + rremove("y.title") + rremove("x.title"), 
          late_scatterplot + rremove("y.title") + rremove("x.title"), 
          hatchling_scatterplot + rremove("y.title") + rremove("x.title"),
          adult_scatterplot + rremove("y.title") + rremove("x.title"),
          ncol = 2, nrow = 2) +
  ggsave("Dimorphilus_correlation_scatterplots.pdf", width = 21, height = 20, units = "cm")


# 2. Cumulative frequency plots
early_cumfreq <- ggplot(df_TPM) + 
  stat_ecdf(aes(x= early_development_rep1,y=1-..y..),geom = "point", size=0.2,
            show.legend = TRUE,
            color=rgb(86,59,69, maxColorValue = 255)) +
  stat_ecdf(aes(x= early_development_rep2,y=1-..y..),geom = "point", size=0.2,
            show.legend = TRUE,
            color=rgb(86,59,69, maxColorValue = 255), alpha = 0.4) +
  scale_x_log10(limits = c(1e-3,1e5), breaks = c(1e-3,1e-2,1e-1,1,10,100,1000,10000,100000), expand = c(0,0)) +
  theme_classic() +
  labs(y = "cumulative frequency", x = "gene expression (TPM)")

late_cumfreq <- ggplot(df_TPM) + 
  stat_ecdf(aes(x= late_development_rep1,y=1-..y..),geom = "point", size=0.2,
            show.legend = TRUE,
            color=rgb(210,163,89, maxColorValue = 255)) +
  stat_ecdf(aes(x= late_development_rep2,y=1-..y..),geom = "point", size=0.2,
            show.legend = TRUE,
            color=rgb(210,163,89, maxColorValue = 255), alpha = 0.4) +
  scale_x_log10(limits = c(1e-3,1e5), breaks = c(1e-3,1e-2,1e-1,1,10,100,1000,10000,100000), expand = c(0,0)) +
  theme_classic() +
  labs(y = "cumulative frequency", x = "gene expression (TPM)")

hatchling_cumfreq <- ggplot(df_TPM) + 
  stat_ecdf(aes(x= hatchling_rep1,y=1-..y..),geom = "point", size=0.2,
            show.legend = TRUE,
            color=rgb(198,97,58, maxColorValue = 255)) +
  stat_ecdf(aes(x= hatchling_rep2,y=1-..y..),geom = "point", size=0.2,
            show.legend = TRUE,
            color=rgb(198,97,58, maxColorValue = 255), alpha = 0.4) +
  scale_x_log10(limits = c(1e-3,1e5), breaks = c(1e-3,1e-2,1e-1,1,10,100,1000,10000,100000), expand = c(0,0)) +
  theme_classic() +
  labs(y = "cumulative frequency", x = "gene expression (TPM)")

adult_cumfreq <- ggplot(df_TPM) + 
  stat_ecdf(aes(x= female_adult_rep1,y=1-..y..),geom = "point", size=0.2,
            show.legend = TRUE,
            color=rgb(11,108,160, maxColorValue = 255)) +
  stat_ecdf(aes(x= female_adult_rep2,y=1-..y..),geom = "point", size=0.2,
            show.legend = TRUE,
            color=rgb(11,108,160, maxColorValue = 255), alpha = 0.4) +
  scale_x_log10(limits = c(1e-3,1e5), breaks = c(1e-3,1e-2,1e-1,1,10,100,1000,10000,100000), expand = c(0,0)) +
  theme_classic() +
  labs(y = "cumulative frequency", x = "gene expression (TPM)")

ggarrange(early_cumfreq + rremove("y.title") + rremove("x.title"), 
          late_cumfreq + rremove("y.title") + rremove("x.title"), 
          hatchling_cumfreq + rremove("y.title") + rremove("x.title"),
          adult_cumfreq + rremove("y.title") + rremove("x.title"),
          ncol = 2, nrow = 2)  +
  ggsave("Dimorphilus_cumulative_frequency_plots.pdf", width = 22, height = 20, units = "cm")


# 3. Correlation plots (TPM and DESeq2)
df_TPM <- read.table("01-Dimorphilus_gyrociliatus_RNAseq_TPM_replicates.txt", sep='\t', header=TRUE)
df_norm <- read.table("02-Dimorphilus_gyrociliatus_RNAseq_DESeq2_replicates.txt", sep='\t', header=TRUE)

palet1 = colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(100)
palet2 = colorRampPalette(rev(brewer.pal(n = 7, name = "Greys")))(100)
palet = cbind(palet2,palet1)

cc = cor(df_TPM, method = "pearson")
corrplot(cc, tl.col = "black", method = "color", col=palet, cl.lim=c(0,1))

cc = cor(df_norm, method = "pearson")
corrplot(cc, tl.col = "black", method = "color", col=palet, cl.lim=c(0,1))


# 4. PCA analysis (DESeq2)
df_pca <- df_norm[rowSums(df_norm) != 0,]
pca_results_raw <- prcomp(t(df_pca), scale = TRUE)

fviz_pca_ind(pca_results_raw, 
             repel = TRUE) +
  theme_classic() +
  ggsave("Dimorphilus_PCA.pdf", width = 20, height = 20, units = "cm")


# 5. Ridgeline plots
df_TPM <- read.table("01-Dimorphilus_gyrociliatus_RNAseq_TPM_average.txt", header = TRUE, sep='\t')
df_average_log <- log(df_TPM,2)
df_average_log[!is.finite(data.matrix(df_average_log))] <- NA
df_average_log_tidy <- gather(df_average_log, "stage", "expression")
df_average_log_tidy$stage <- factor(df_average_log_tidy$stage,
                                    levels = c("early_development", "late_development",
                                               "hatchling", "female_adult"), 
                                    ordered = TRUE)

ggplot(df_average_log_tidy, aes(x = expression, y = stage)) + 
  geom_density_ridges(aes(fill = stage), rel_min_height = 0.01) +
  scale_x_continuous(limits = c(-12,17), breaks = c(-10,-5,0,5,10,15),
                     expand = c(0,0)) +
  scale_fill_manual(values = c("#563B45","#D2A359","#C6613A",
                               "#0B6CA0","#032C65","#507810",
                               "#204858")) +
  theme_classic() +
  labs(x = "gene expression (log2(TPM))") +
  geom_vline(xintercept=1, linetype="dashed", 
             color = rgb(255,0,0, maxColorValue = 255)) +
  geom_vline(xintercept=0, linetype="dashed", 
             color = rgb(0,0,255, maxColorValue = 255)) +
  theme(axis.title.x = element_text(color="black",size=10),
        axis.title.y = element_text(color="black",size=10),
        axis.text.x = element_text(color="black", size=10,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=10,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"),
        legend.position = "none") +
  ggsave("Dimorphilus_ridgeline_plots.pdf", height = 15, width = 18, units = "cm")


# 6. Alluvial diagram
alluvial_raw <- df_TPM
alluvial_raw[alluvial_raw <= 2] <- 0
alluvial_raw[alluvial_raw > 2] <- 1
alluvial_curated <- as.data.frame(table(alluvial_raw))
is_alluvia_form(alluvial_curated, axes = 1:4, silent = TRUE)

group <- c(1:16)
alluvial_long <- alluvial_curated
alluvial_tidy <- gather(alluvial_long, "stage", "expression", -Freq)
alluvial_tidy <- cbind(alluvial_tidy, group)
alluvial_tidy$stage <- factor(alluvial_tidy$stage,
                                    levels = c("early_development", "late_development",
                                               "hatchling", "female_adult"),
                                    ordered = TRUE)

ggplot(alluvial_tidy,
       aes(x = stage, stratum = expression, alluvium = group,
           y = Freq,
           fill = expression, label = expression)) +
  geom_flow(alpha = .4) +
  geom_stratum(alpha = .6) +
  geom_text(stat = "stratum", size = 3) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y = "# of genes") +
  scale_x_discrete(expand = c(0.05, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black",size=14),
        axis.text.x = element_text(color="black",size=12,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=12,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length.y = unit(-0.15, "cm"),
        axis.ticks.length.x = unit(0, "cm")) +
  ggsave("Dimorphilus_alluvial_plots.pdf", height = 20, width = 15, units = "cm")



##############################################
###### PART B. Soft k-means clustering #######
##############################################

library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(corrplot)
library(RColorBrewer)
library(Mfuzz)
library(ComplexHeatmap)

# 1. Fuzzy soft k-means clustering
df_average <- read.table("02-Dimorphilus_gyrociliatus_RNAseq_DESeq2_average.txt", sep='\t', header=TRUE)
mat_average <- as.matrix(df_average)
eset_rna <- new('ExpressionSet', exprs = mat_average)
# 17,388 genes
clean_eset_rna <- filter.std(eset_rna,min.std=0) 
# 200 genes with null expression at all stages: remaining 17,188
norm_eset_rna <- standardise(clean_eset_rna)
fuzzifier <- round(mestimate(norm_eset_rna),2)
# fuzzifier: 1.43
centroids <- Dmin(norm_eset_rna, fuzzifier, crange = seq(4,20,1), repeats=10, visu=TRUE)
# optimal centroid with 9 clusters
cluster_number <- 9
# Run the two following lines together to ensure reproducibility
set.seed(123)
clusters <- mfuzz(norm_eset_rna, c = 9, m = fuzzifier)
clusters$size 
# 1897 1665 1847 1864 2021 1634 2682 1615 1963

df_names <- df_average
df_names_notnull <- df_names[apply(df_names, 1, function(x) sd(x, na.rm = TRUE)!=0),]

clusters_annotated_unordered <- cbind(rownames(df_names_notnull),clusters$cluster)
colnames(clusters_annotated_unordered) <- c("Gene_ID","Cluster_raw")
write.table(clusters_annotated_unordered,"03-Dimorphilus_gyrociliatus_clusters_annotation_raw.txt", 
            quote = F, sep = '\t', row.names = FALSE)

clusters_expression_unordered <- cbind(rownames(df_names_notnull),df_names_notnull,clusters$cluster)
colnames(clusters_expression_unordered) <- c("Gene_ID","early_development","late_development",
                                             "hatchling","female_adult","Cluster_raw")
write.table(clusters_expression_unordered,"04-Capitella_teleta_clusters_expression_raw.txt", 
            quote = F, sep = '\t', row.names = FALSE)


#

df_raw <- read.table("04-Dimorphilus_gyrociliatus_clusters_expression_raw.txt", header = TRUE)
df_zscore <- data.frame(cbind(t(scale(t(df_raw[-c(1,6)]))),'cluster_raw'=df_raw$Cluster_raw))
df_zscore$cluster_raw <- factor(df_zscore$cluster_raw, 
                                levels = c("3","7","9","6","5","1","4","2","8"))

df_sorted <- df_zscore[order(df_zscore$cluster_raw),]

heatmap_color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

ComplexHeatmap::Heatmap(as.matrix(df_sorted[-c(5)]),
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        show_row_names = FALSE,
                        col = heatmap_color)

#re-write clusters
# equivalence is the following (old:new)
# 1:6, 2:8, 3:1, 4:7, 5:5, 6:4, 7:2, 8:9, 9:3
df_reordered <- df_raw
df_reordered$Cluster_corrected <- c(rep(0,nrow(df_reordered)))
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 1] <- 6
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 2] <- 8
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 3] <- 1
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 4] <- 7
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 5] <- 5
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 6] <- 4
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 7] <- 2
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 8] <- 9
df_reordered$Cluster_corrected[df_reordered$Cluster_raw == 9] <- 3


df_reordered_clean <- df_reordered[-c(6)]

write.table(df_reordered_clean[c(1,6)],"06-Dimorphilus_gyrociliatus_clusters_annotation_corrected.txt", 
            quote = F, sep = '\t', row.names = FALSE)
write.table(df_reordered_clean,"07-Dimorphilus_gyrociliatus_clusters_expression_corrected.txt", 
            quote = F, sep = '\t', row.names = FALSE)


# 2. Loess smothing

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

palette <- viridis_pal(option = "C", direction = 1)(9)
time <- c(1:4)

df_raw <- read.table("07-Dimorphilus_gyrociliatus_clusters_expression_corrected.txt", header = T)

# For cluster 1
df_cluster1_raw <- df_raw[df_raw$Cluster_corrected == 9,-c(1,6)] # change 1 for each cluster and run all lines
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
cluster9_plot <- ggplot(summary_cluster1, aes(x=time, y=mean)) + # change cluster1_plot to clusterX_plot for each cluster
  geom_line(data = df_cluster1_clean_long, aes(x=time, y=expression, group=gene_ID), 
            color="gray") +
  stat_smooth(colour = palette[9], # change 1 to match the number of the cluster, so that each cluster has its own colour of the continuous palette
              fill = palette[9], # change 1 to match the number of the cluster, so that each cluster has its own colour of the continuous palette
              size = 2, alpha = 0.5) +  
  geom_hline(yintercept = 0, linetype="dashed", color = "black") + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(y = "Normalised gene expression (z-score)", x = "developmental stage") +
  scale_y_continuous(breaks = c(-3:3), limits = c(-3,3)) +
  scale_x_continuous(breaks = c(1:4), expand = c(0,0), limits = c(0.8,4.2),
                     labels = c("early_development","late_development",
                                "hatchling","female_adult")) +
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
                        cluster9_plot,
                        ncol = 3, nrow = 3)

allfigures +
  ggsave("08-Dimorphilus_gyrociliatus_RNAseq_clusters_loess.pdf", width = 25, height = 27, units = "cm")



## 3. GO terms calculation for RNAseq clusters
library(ggplot2)
library(ggpubr)
library(topGO)

geneID2GO <- readMappings(file = "09-Dimorphilus_gyrociliatus_geneID2GO_topGO_onlyGOgenes.txt")
geneUniverse <- names(geneID2GO)

df_clusters <- read.table("06-Dimorphilus_gyrociliatus_clusters_annotation_corrected.txt", header = T)

cluster_1 <- df_clusters[df_clusters$Cluster_corrected == 1,]
cluster_2 <- df_clusters[df_clusters$Cluster_corrected == 2,]
cluster_3 <- df_clusters[df_clusters$Cluster_corrected == 3,]
cluster_4 <- df_clusters[df_clusters$Cluster_corrected == 4,]
cluster_5 <- df_clusters[df_clusters$Cluster_corrected == 5,]
cluster_6 <- df_clusters[df_clusters$Cluster_corrected == 6,]
cluster_7 <- df_clusters[df_clusters$Cluster_corrected == 7,]
cluster_8 <- df_clusters[df_clusters$Cluster_corrected == 8,]
cluster_9 <- df_clusters[df_clusters$Cluster_corrected == 9,]

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


write.table(results_cluster_1, "10-Cluster_1_GO_terms.txt", quote=FALSE, sep='\t', row.names = F) 
write.table(results_cluster_2, "10-Cluster_2_GO_terms.txt", quote=FALSE, sep='\t', row.names = F)
write.table(results_cluster_3, "10-Cluster_3_GO_terms.txt", quote=FALSE, sep='\t', row.names = F) 
write.table(results_cluster_4, "10-Cluster_4_GO_terms.txt", quote=FALSE, sep='\t', row.names = F) 
write.table(results_cluster_5, "10-Cluster_5_GO_terms.txt", quote=FALSE, sep='\t', row.names = F) 
write.table(results_cluster_6, "10-Cluster_6_GO_terms.txt", quote=FALSE, sep='\t', row.names = F)
write.table(results_cluster_7, "10-Cluster_7_GO_terms.txt", quote=FALSE, sep='\t', row.names = F) 
write.table(results_cluster_8, "10-Cluster_8_GO_terms.txt", quote=FALSE, sep='\t', row.names = F) 
write.table(results_cluster_9, "10-Cluster_9_GO_terms.txt", quote=FALSE, sep='\t', row.names = F) 


results_clusters <- rbind(cbind(Cluster="cluster 1", results_cluster_1),
                          cbind(Cluster="cluster 2", results_cluster_2),
                          cbind(Cluster="cluster 3", results_cluster_3),
                          cbind(Cluster="cluster 4", results_cluster_4),
                          cbind(Cluster="cluster 5", results_cluster_5),
                          cbind(Cluster="cluster 6", results_cluster_6),
                          cbind(Cluster="cluster 7", results_cluster_7),
                          cbind(Cluster="cluster 8", results_cluster_8),
                          cbind(Cluster="cluster 9", results_cluster_9))

write.table(results_clusters, "10-All_clusters_GO_terms.txt", quote=FALSE, sep='\t', row.names = F)




# 3. Plotting bar plots of -log10(p-values) of GO terms for each cluster
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


## 4. k-means clustering of GO terms and summary visualisation of enriched GO terms in clusters
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

Dimorphilus_GO_clusters <- read.table("14-Dimorphilus_gyrociliatus_GO_terms.txt", header=TRUE, sep ='\t')
colnames(Owenia_GO_clusters) <- c("Cluster", "GO.ID", "Term", "Annotated", "Significant", 
                                  "Expected", "classicFisher")

go_id <- rbind(Owenia_GO_clusters, Capitella_GO_clusters, Dimorphilus_GO_clusters)
go_id$GO.ID <- as.character(go_id$GO.ID)
mat <- GO_similarity(go_id$GO.ID, ont = "BP")

set.seed(22) # ensures reproducibility
GO_cluster_annotation_raw <- simplifyGO(mat, method = "kmeans")

GO_cluster_annotation_clean <- GO_cluster_annotation_raw[2:ncol(GO_cluster_annotation_raw)]
rownames(GO_cluster_annotation_clean) <- t(GO_cluster_annotation_raw[1]) 
Capitella_GO_clusters$Cluster <- Capitella_GO_clusters$Cluster+c(12)
Dimorphilus_GO_clusters$Cluster <- Dimorphilus_GO_clusters$Cluster+c(24)

GO_clusters <- rbind(Owenia_GO_clusters, Capitella_GO_clusters, Dimorphilus_GO_clusters)
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
colnames(heatmap_data_long) <- c(1:33)

heatmap_data_long$kmeans_cluster <- GO_cluster_annotation_clean[rownames(heatmap_data_long),2]
heatmap_data_long$averageFisher <- rowMeans(heatmap_data_long[,c(1:33)], na.rm = TRUE)
heatmap_data_long <- heatmap_data_long[order(heatmap_data_long$kmeans_cluster, 
                                             -heatmap_data_long$averageFisher),]

#manually fix 4 GO terms that were not clustered:
# GO:0070838 metal ion transport -> cluster 1 (cellular transport) -> now replaced by GO:0030001
# GO:0055114 oxidation-reduction process -> DISCARDED (it's obsolote because it's a molecular function, not a biological process)
# GO:0072511 inorganic cation transmembrane transport -> cluster 1 (cellular transport)
# GO:0017144: xenobiotic metabolic process -> cluster 15 (complex metabolism)

heatmap_data_long$kmeans_cluster[rownames(heatmap_data_long) == 'GO:0070838'] <- 1
heatmap_data_long$kmeans_cluster[rownames(heatmap_data_long) == 'GO:0072511'] <- 1
heatmap_data_long$kmeans_cluster[rownames(heatmap_data_long) == 'GO:0017144'] <- 15

to_remove <- c("GO:0055114")
heatmap_data_final <- heatmap_data_long[!(row.names(heatmap_data_long) %in% to_remove),]

palette <- colorRampPalette(brewer.pal(n = 7, name = "Reds"))(200)
palette[1] <- rgb(255,255,255, maxColorValue = 255)

mat_heatmap <- data.matrix(heatmap_data_final[-c(34,35)])
mat_heatmap[is.na(mat_heatmap)] <- 0 # NAs get converted to 0s, which in turn will be shown in the heatmap in white
ComplexHeatmap::Heatmap(mat_heatmap, 
                        row_split = heatmap_data_final$kmeans_cluster,
                        border = TRUE,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        show_row_names = FALSE,
                        col = palette)

# I fusioned clusters 6 and 7 into a single one because they are almost indistinguishable from each other
# In order: cellular transport (1), motility and movement (2), ion homeostasis (3), response to stimuli (4), 
# cell signaling (5), cell differentiation, morpho- and organogenesis (6,7), cell structures assembly (8),
# cytoskeleton dynamics (9), OXPHOS metabolism (10), immune system development (11), RNA biosynthesis (12),
# amino acid and lipid metabolism (13), glycoconjugates metabolism (14), complex metabolism (15), RNA processing (16)
