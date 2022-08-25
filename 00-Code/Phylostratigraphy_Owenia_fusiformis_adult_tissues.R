library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(RColorBrewer)
library(viridis)
library(patchwork)
library(preprocessCore)
library(pheatmap)

# Import all data
raw_data <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/24-Phylostratigraphy_Owenia_RNAseq_adult_tissues/01-Raw_data_Owenia_fusiformis.txt", 
                       as.is=TRUE, header=TRUE)
raw_data_clean <- raw_data[-c(1,2)]
rownames(raw_data_clean) <- t(raw_data[c(1)])


# Import anterior intersection genes
anterior_genes <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/24-Phylostratigraphy_Owenia_RNAseq_adult_tissues/02-Anterior_intersection_genes.txt",
                             sep = '\t', header = T)
rownames(anterior_genes) <- t(anterior_genes[c(1)])
colnames(anterior_genes) <- c("Gene_ID","blastulae", "gastrulae", "elongation",
                              "early_larvae", "mature_larvae", "competent_larvae","juvenile")


# Import posterior intersection genes
posterior_genes <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/24-Phylostratigraphy_Owenia_RNAseq_adult_tissues/03-Posterior_intersection_genes.txt",
                              sep = '\t', header = T)
rownames(posterior_genes) <- t(posterior_genes[c(1)])
colnames(posterior_genes) <- c("Gene_ID","blastulae", "gastrulae", "elongation",
                               "early_larvae", "mature_larvae", "competent_larvae","juvenile")


# Quantile normalisation
norm <- data.frame(cbind(raw_data_clean[1],normalize.quantiles(data.matrix(raw_data_clean[,2:ncol(raw_data_clean)]))))
norm$ofus.origin[norm$ofus.origin == "Ofus-2"] <- "Ofus"
norm$ofus.origin[norm$ofus.origin == "Nephrozoa" | norm$ofus.origin == "Bilateria"] <- "Nephrozoa_Bilateria"

colnames(norm) <- c("origin","blastulae", "gastrulae", "elongation",
                    "early_larvae", "mature_larvae", "competent_larvae","juvenile")
norm_tidy <- gather(norm, stage, expr, -origin)
summary(norm)


# Subset anterior and posterior intersection genes from normalised set
`%!in%` <- Negate(`%in%`)
norm_anterior <- norm[rownames(norm) %in% rownames(anterior_genes),]
norm_posterior <- norm[rownames(norm) %in% rownames(posterior_genes),]
norm_other <- norm[rownames(norm) %!in% rownames(anterior_genes) & rownames(norm) %!in% rownames(posterior_genes),]


# Get 75% quantile of expression for these gene sets:

## 1. For anterior genes
norm_75_quant_anterior <- data.frame(tapply(norm_anterior$blastulae, norm_anterior$origin, quantile, p = 0.75, na.rm = TRUE),
                                     tapply(norm_anterior$gastrulae, norm_anterior$origin, quantile, p = 0.75, na.rm = TRUE),
                                     tapply(norm_anterior$elongation, norm_anterior$origin, quantile, p = 0.75, na.rm = TRUE),
                                     tapply(norm_anterior$early_larvae, norm_anterior$origin, quantile, p = 0.75, na.rm = TRUE),
                                     tapply(norm_anterior$mature_larvae, norm_anterior$origin, quantile, p = 0.75, na.rm = TRUE),
                                     tapply(norm_anterior$competent_larvae, norm_anterior$origin, quantile, p = 0.75, na.rm = TRUE),
                                     tapply(norm_anterior$juvenile, norm_anterior$origin, quantile, p = 0.75, na.rm = TRUE))
colnames(norm_75_quant_anterior) <- c("blastulae", "gastrulae", "elongation",
                                      "early_larvae", "mature_larvae", 
                                      "competent_larvae","juvenile")
rownames(norm_75_quant_anterior) <- c("Annelida","Lophotrochozoa","Metazoa",
                                      "Nephrozoa_Bilateria","Owenia_fusiformis","Other","Protostomia")
norm_75_quant_anterior <- norm_75_quant_anterior[c(3,4,7,2,1,5,6),]


## 2. For posterior genes
norm_75_quant_posterior <- data.frame(tapply(norm_posterior$blastulae, norm_posterior$origin, quantile, p = 0.75, na.rm = TRUE),
                                      tapply(norm_posterior$gastrulae, norm_posterior$origin, quantile, p = 0.75, na.rm = TRUE),
                                      tapply(norm_posterior$elongation, norm_posterior$origin, quantile, p = 0.75, na.rm = TRUE),
                                      tapply(norm_posterior$early_larvae, norm_posterior$origin, quantile, p = 0.75, na.rm = TRUE),
                                      tapply(norm_posterior$mature_larvae, norm_posterior$origin, quantile, p = 0.75, na.rm = TRUE),
                                      tapply(norm_posterior$competent_larvae, norm_posterior$origin, quantile, p = 0.75, na.rm = TRUE),
                                      tapply(norm_posterior$juvenile, norm_posterior$origin, quantile, p = 0.75, na.rm = TRUE))
colnames(norm_75_quant_posterior) <- c("blastulae", "gastrulae", "elongation",
                                       "early_larvae", "mature_larvae", 
                                       "competent_larvae","juvenile")
rownames(norm_75_quant_posterior) <- c("Annelida","Lophotrochozoa","Metazoa","Nephrozoa_Bilateria","Owenia_fusiformis","Other","Protostomia")
norm_75_quant_posterior <- norm_75_quant_posterior[c(3,4,7,2,1,5,6),]


# Plot contribution and dynamics heatmaps
heatmap_color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

#Heatmap 1: Scaled by row, anterior
pheatmap(norm_75_quant_anterior, 
         scale = "row",
         cluster_rows=F,
         cluster_cols=F, 
         show_rownames=T, 
         fontsize= 12,
         col=heatmap_color, 
         border_color=NA, 
         cellwidth = 20, 
         cellheight = 20)

#Heatmap 2: Scaled by row, posterior
pheatmap(norm_75_quant_posterior, 
         scale = "row",
         cluster_rows=F,
         cluster_cols=F, 
         show_rownames=T, 
         fontsize= 12,
         col=heatmap_color, 
         border_color=NA, 
         cellwidth = 20, 
         cellheight = 20)

#Heatmap 3: Scaled by column, anterior
pheatmap(norm_75_quant_anterior, 
         scale = "column",
         cluster_rows=F,
         cluster_cols=F, 
         show_rownames=T, 
         fontsize= 12,
         col=heatmap_color, 
         border_color=NA, 
         cellwidth = 20, 
         cellheight = 20)

#Heatmap 4: Scaled by column, posterior
pheatmap(norm_75_quant_posterior, 
         scale = "column",
         cluster_rows=F,
         cluster_cols=F, 
         show_rownames=T, 
         fontsize= 12,
         col=heatmap_color, 
         border_color=NA, 
         cellwidth = 20, 
         cellheight = 20)


# Plot pie charts
pie_chart_anterior <- data.frame(table(norm_anterior$origin))
pie_chart_posterior <- data.frame(table(norm_posterior$origin))
pie_chart_all <- data.frame(table(norm$origin))

pie_chart_anterior$percentage <- round(pie_chart_anterior$Freq/sum(pie_chart_anterior$Freq)*100,1)
pie_chart_posterior$percentage <- round(pie_chart_posterior$Freq/sum(pie_chart_posterior$Freq)*100,1)
pie_chart_all$percentage <- round(pie_chart_all$Freq/sum(pie_chart_all$Freq)*100,1)

colnames(pie_chart_anterior) <- c("origin", "count", "percentage")
colnames(pie_chart_posterior) <- c("origin", "count", "percentage")
colnames(pie_chart_all) <- c("origin", "count", "percentage")

pie_chart_anterior$percentage <- as.numeric(pie_chart_anterior$percentage)
pie_chart_posterior$percentage <- as.numeric(pie_chart_posterior$percentage)
pie_chart_all$percentage <- as.numeric(pie_chart_all$percentage)

pie_chart_anterior$origin <- c("Annelida","Lophotrochozoa","Metazoa",
                               "Nephrozoa_Bilateria","Owenia_fusiformis","Other","Protostomia")
pie_chart_posterior$origin <- c("Annelida","Lophotrochozoa","Metazoa",
                               "Nephrozoa_Bilateria","Owenia_fusiformis","Other","Protostomia")
pie_chart_all$origin <- c("Annelida","Lophotrochozoa","Metazoa",
                          "Nephrozoa_Bilateria","Owenia_fusiformis","Other","Protostomia")

pie_chart_anterior$origin <- factor(pie_chart_anterior$origin, levels = c("Metazoa","Nephrozoa_Bilateria","Protostomia","Lophotrochozoa",
                                                                          "Annelida", "Owenia_fusiformis","Other"))
pie_chart_posterior$origin <- factor(pie_chart_posterior$origin, levels = c("Metazoa","Nephrozoa_Bilateria","Protostomia","Lophotrochozoa",
                                                                            "Annelida", "Owenia_fusiformis","Other"))
pie_chart_all$origin <- factor(pie_chart_all$origin, levels = c("Metazoa","Nephrozoa_Bilateria","Protostomia","Lophotrochozoa",
                                                                "Annelida", "Owenia_fusiformis","Other"))

anterior_to_plot <- pie_chart_anterior %>% 
  mutate(end = 2 * pi * cumsum(pie_chart_anterior$percentage)/sum(pie_chart_anterior$percentage),
         start = lag(end, default = 0),
         middle = 0.5 * (start + end),
         hjust = ifelse(middle > pi, 1, 0),
         vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))

posterior_to_plot <- pie_chart_posterior %>% 
  mutate(end = 2 * pi * cumsum(pie_chart_posterior$percentage)/sum(pie_chart_posterior$percentage),
         start = lag(end, default = 0),
         middle = 0.5 * (start + end),
         hjust = ifelse(middle > pi, 1, 0),
         vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))

colours <- c("#F94144","#F3722C","#F8961E","#F9C74F","#90BE6D","#43AA8B","#577590")

ggplot(anterior_to_plot) + 
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                   start = start - pi/4, end = end - pi/4, fill = origin)) +
  geom_text(aes(x = 1.05 * sin(middle - pi/4) , y = 1.05 * cos(middle - pi/4), label = paste0(origin,": ",percentage,"% (",count,")"),
                hjust = hjust, vjust = vjust)) +
  coord_fixed() +
  scale_fill_manual(values = colours) +
  scale_x_continuous(limits = c(-2,2),
                     name = "", breaks = NULL, labels = NULL) +
  scale_y_continuous(limits = c(-1.5, 1.5),
                     name = "", breaks = NULL, labels = NULL) +
  theme_void() +
  theme(legend.position = "none")

ggplot(posterior_to_plot) + 
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                   start = start - pi/4, end = end - pi/4, fill = origin)) +
  geom_text(aes(x = 1.05 * sin(middle - pi/4) , y = 1.05 * cos(middle - pi/4), label = paste0(origin,": ",percentage,"% (",count,")"),
                hjust = hjust, vjust = vjust)) +
  coord_fixed() +
  scale_x_continuous(limits = c(-2,2),
                     name = "", breaks = NULL, labels = NULL) +
  scale_y_continuous(limits = c(-1.5, 1.5),
                     name = "", breaks = NULL, labels = NULL) +
  scale_fill_manual(values = colours) +
  theme_void() +
  theme(legend.position = "none")


# Barplot of gene sets vs. whole genome and random subset for comparison
random_subset <- norm[sample(nrow(norm),1000),] # subset of 1,000 genes
pie_chart_random <- data.frame(table(random_subset$origin))
pie_chart_random$percentage <- round(pie_chart_random$Freq/sum(pie_chart_random$Freq)*100,1)
colnames(pie_chart_random) <- c("origin", "count", "percentage")
pie_chart_random$percentage <- as.numeric(pie_chart_random$percentage)
pie_chart_random$origin <- c("Annelida","Lophotrochozoa","Metazoa",
                             "Nephrozoa_Bilateria","Owenia_fusiformis","Other","Protostomia")
pie_chart_random$origin <- factor(pie_chart_random$origin, levels = c("Metazoa","Nephrozoa_Bilateria","Protostomia","Lophotrochozoa",
                                                                      "Annelida", "Owenia_fusiformis","Other"))

barplot <- cbind(pie_chart_all[1],
                 pie_chart_all[3], 
                 pie_chart_anterior[3],
                 pie_chart_posterior[3],
                 pie_chart_random[3])
colnames(barplot) <- c("origin","whole_genome","anterior","posterior","random")
barplot$origin <- factor(barplot$origin, levels = c("Metazoa","Nephrozoa_Bilateria","Protostomia","Lophotrochozoa",
                                                    "Annelida", "Owenia_fusiformis","Other"))

barplot_tidy <- gather(barplot, gene_set, percentage, -origin)
barplot_tidy$gene_set <- factor(barplot_tidy$gene_set, levels = c("anterior", "posterior", "whole_genome","random"))

ggplot(barplot_tidy, aes(fill = origin, y = percentage, x = gene_set)) + 
  geom_bar(position="stack", stat="identity") +
  scale_x_discrete(labels = c("Anterior", "Posterior", "Whole genome","Random")) +
  labs(y = "% of transcription factors", x = NA) +
  scale_fill_manual(values = colours) +
  theme_classic() +
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
        axis.ticks.length = unit(-0.15, "cm"))