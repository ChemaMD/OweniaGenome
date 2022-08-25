library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ComplexHeatmap)
library(viridis)

species_palette <- viridis(3)

# Import 1-to-1 orthologues between Owenia and Capitella/Dimorphilus
ofus2ctel <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/03-Fran_plots_expression/01-Ofus2Ctel_orthologues.txt", header=T, sep ='\t')
ofus2dgyr <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/03-Fran_plots_expression/01-Ofus2Dgyr_orthologues.txt", header=T, sep ='\t')
rownames(ofus2ctel) <- ofus2ctel$Ofus
rownames(ofus2dgyr) <- ofus2dgyr$Ofus

# Import RNA-seq data
ofus_deseq2 <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/03-Fran_plots_expression/04-Owenia_fusiformis_RNAseq_DESeq2_average.txt",
                          sep = ' ', header = TRUE)
ofus_deseq2 <- ofus_deseq2[-c(2:8)]
ctel_deseq2 <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/03-Fran_plots_expression/04-Capitella_teleta_RNAseq_DESeq2_average.txt",
                          row.names = 1, sep = '\t', header = TRUE)
ctel_deseq2$Gene_ID <- rownames(ctel_deseq2)
ctel_deseq2 <- ctel_deseq2[c(8,1,2,3,4,5,6,7)]
rownames(ctel_deseq2) <- c(1:nrow(ctel_deseq2))
dgyr_deseq2 <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/03-Fran_plots_expression/04-Dimorphilus_gyrociliatus_RNAseq_DESeq2_average.txt",
                          sep = '\t', header = TRUE)
dgyr_deseq2$Gene_ID <- rownames(dgyr_deseq2)
dgyr_deseq2 <- dgyr_deseq2[c(5,1,2,3,4)]
rownames(dgyr_deseq2) <- c(1:nrow(dgyr_deseq2))

# Import TF classification for Owenia
tfclass <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/16-Transcription_factors_Owenia_Capitella_Dimorphilus/06-Owenia_fusiformis_TF_list_classified.txt",
                      sep = '\t', header = F)
rownames(tfclass) <- tfclass$V1
tfclass <- tfclass[-c(1)]
colnames(tfclass) <- c("PFAM_ID", "PFAM_annotation")



####################################################################################
### 1. Displaced genes between Capitella late, Owenia late and Dimorphilus early ###
####################################################################################

# Import list of genes to plot for all species
ofus_displaced <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/03-Fran_plots_expression/03-Common_displaced_Dimorphilus_short.txt", header=T, sep ='\t')
all_displaced <- cbind(ofus_displaced[c(2,1)], ofus2ctel[ofus_displaced$transcript_id,"Ctel"], ofus2dgyr[ofus_displaced$transcript_id,"Dgyr"])
colnames(all_displaced) <- c("gene_symbol", "ofus_id", "ctel_id", "dgyr_id")
all_displaced$pfam_annotation <- tfclass[all_displaced$ofus_id, "PFAM_annotation"]
all_displaced <- all_displaced[order(all_displaced$pfam_annotation),]

write.table(all_displaced, "09-Displaced_genes_supp_table.txt", sep ='\t', quote = F, row.names = F)

# Expression dynamics of the TFs we have found in deseq2 levels (from 0 to max and z-scored)

rescale_custom <- function(x) (x/(max(x))) # define function to make 0 to max scale

ofus_displaced_deseq2 <- ofus_deseq2[ofus_deseq2$Gene_ID %in% all_displaced$ofus_id,]
rownames(ofus_displaced_deseq2) <- ofus_displaced_deseq2$Gene_ID
ctel_displaced_deseq2 <- ctel_deseq2[ctel_deseq2$Gene_ID %in% all_displaced$ctel_id,]
rownames(ctel_displaced_deseq2) <- ctel_displaced_deseq2$Gene_ID
dgyr_displaced_deseq2 <- dgyr_deseq2[dgyr_deseq2$Gene_ID %in% all_displaced$dgyr_id,]
rownames(dgyr_displaced_deseq2) <- dgyr_displaced_deseq2$Gene_ID

ofus_normalised <- t(apply(ofus_displaced_deseq2[-c(1)], 1, rescale_custom))
ctel_normalised <- t(apply(ctel_displaced_deseq2[-c(1)], 1, rescale_custom))
dgyr_normalised <- t(apply(dgyr_displaced_deseq2[-c(1)], 1, rescale_custom))

ofus_zscore <- t(scale(t(ofus_displaced_deseq2[-c(1)])))
ctel_zscore <- t(scale(t(ctel_displaced_deseq2[-c(1)])))
dgyr_zscore <- t(scale(t(dgyr_displaced_deseq2[-c(1)])))

ofus_order <- all_displaced$ofus_id
ctel_order <- all_displaced$ctel_id
dgyr_order <- all_displaced$dgyr_id

ofus_normalised_sorted <- ofus_normalised[match(ofus_order, rownames(ofus_normalised)),]
ctel_normalised_sorted <- ctel_normalised[match(ctel_order, rownames(ctel_normalised)),]
dgyr_normalised_sorted <- dgyr_normalised[match(dgyr_order, rownames(dgyr_normalised)),]

ofus_zscore_sorted <- ofus_zscore[match(ofus_order, rownames(ofus_zscore)),]
ctel_zscore_sorted <- ctel_zscore[match(ctel_order, rownames(ctel_zscore)),]
dgyr_zscore_sorted <- dgyr_zscore[match(dgyr_order, rownames(dgyr_zscore)),]

#1/2 rdbu (substitute of reds for 0 to 1)
heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color <- heatmap_color[50:100]

ofus_0_to_max <- ComplexHeatmap::Heatmap(ofus_normalised_sorted,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        col = heatmap_color,
                        row_labels = all_displaced$gene_symbol,
                        left_annotation = rowAnnotation(foo = anno_text(tfclass[all_displaced$ofus_id,"PFAM_annotation"])),
                        column_labels = c("blastula", "gastrula", "elongation",
                                          "early larva", "mitraria larva",
                                          "competent larva", "juvenile"),
                        heatmap_legend_param = list(color_bar = "continuous"))

ctel_0_to_max <- ComplexHeatmap::Heatmap(ctel_normalised_sorted,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        col = heatmap_color,
                        column_labels = c("64 cells", "gastrula", "stage 4tt larva",
                                          "stage 5 larva", "stage 7 larva",
                                          "pre-competent larva", "competent larva"),
                        row_labels = all_displaced$gene_symbol,
                        heatmap_legend_param = list(color_bar = "continuous"))

dgyr_0_to_max <- ComplexHeatmap::Heatmap(dgyr_normalised_sorted,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        col = heatmap_color,
                        row_labels = all_displaced$gene_symbol,
                        column_labels = c("early development","late development",
                                          "hatchling", "female adult"),
                        heatmap_legend_param = list(color_bar = "continuous"))

ofus_0_to_max + ctel_0_to_max + dgyr_0_to_max



#full rdbu
heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

ofus_neg_to_pos <- ComplexHeatmap::Heatmap(ofus_zscore_sorted,
                                         cluster_columns = FALSE,
                                         cluster_rows = FALSE,
                                         col = heatmap_color,
                                         row_labels = all_displaced$gene_symbol,
                                         left_annotation = rowAnnotation(foo = anno_text(tfclass[all_displaced$ofus_id,"PFAM_annotation"])),
                                         column_labels = c("blastula", "gastrula", "elongation",
                                                           "early larva", "mitraria larva",
                                                           "competent larva", "juvenile"),
                                         heatmap_legend_param = list(color_bar = "continuous"))

ctel_neg_to_pos <- ComplexHeatmap::Heatmap(ctel_zscore_sorted,
                                         cluster_columns = FALSE,
                                         cluster_rows = FALSE,
                                         col = heatmap_color,
                                         column_labels = c("64 cells", "gastrula", "stage 4tt larva",
                                                           "stage 5 larva", "stage 7 larva",
                                                           "pre-competent larva", "competent larva"),
                                         row_labels = all_displaced$gene_symbol,
                                         heatmap_legend_param = list(color_bar = "continuous"))

dgyr_neg_to_pos <- ComplexHeatmap::Heatmap(dgyr_zscore_sorted,
                                         cluster_columns = FALSE,
                                         cluster_rows = FALSE,
                                         col = heatmap_color,
                                         row_labels = all_displaced$gene_symbol,
                                         column_labels = c("early development","late development",
                                                           "hatchling", "female adult"),
                                         heatmap_legend_param = list(color_bar = "continuous"))

ofus_neg_to_pos + ctel_neg_to_pos + dgyr_neg_to_pos

#alternative_viz
heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

ofus_0_to_max <- ComplexHeatmap::Heatmap(ofus_normalised_sorted,
                                         cluster_columns = FALSE,
                                         cluster_rows = FALSE,
                                         col = heatmap_color,
                                         row_labels = all_displaced$gene_symbol,
                                         left_annotation = rowAnnotation(foo = anno_text(tfclass[all_displaced$ofus_id,"PFAM_annotation"])),
                                         column_labels = c("blastula", "gastrula", "elongation",
                                                           "early larva", "mitraria larva",
                                                           "competent larva", "juvenile"),
                                         heatmap_legend_param = list(color_bar = "continuous"))

ctel_0_to_max <- ComplexHeatmap::Heatmap(ctel_normalised_sorted,
                                         cluster_columns = FALSE,
                                         cluster_rows = FALSE,
                                         col = heatmap_color,
                                         column_labels = c("64 cells", "gastrula", "stage 4tt larva",
                                                           "stage 5 larva", "stage 7 larva",
                                                           "pre-competent larva", "competent larva"),
                                         row_labels = all_displaced$gene_symbol,
                                         heatmap_legend_param = list(color_bar = "continuous"))

dgyr_0_to_max <- ComplexHeatmap::Heatmap(dgyr_normalised_sorted,
                                         cluster_columns = FALSE,
                                         cluster_rows = FALSE,
                                         col = heatmap_color,
                                         row_labels = all_displaced$gene_symbol,
                                         column_labels = c("early development","late development",
                                                           "hatchling", "female adult"),
                                         heatmap_legend_param = list(color_bar = "continuous"))

ofus_0_to_max + ctel_0_to_max + dgyr_0_to_max


# Relative expression lineplots

species_palette <- viridis(3)
time <- c(1:7)
time2 <- c(1:4)


## Owenia
ofus_clean <- data.frame(cbind(time,t(ofus_normalised_sorted)))
ofus_clean_long <- gather(ofus_clean, gene_ID, expression, -time)

ofus_clean_long %>%
  ggplot(aes(x=time, y=expression)) +
  geom_smooth(method ='loess', size = 1, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "DESeq2 normalised expression", x = "developmental stage") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
  scale_x_continuous(expand = c(0,0), limits = c(0.8,7.2), breaks = c(1:7),
                     labels = c("blastula", "gastrula", "elongation",
                                "early larva", "mitraria larva",
                                "competent larva", "juvenile")) +
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


## Capitella
ctel_clean <- data.frame(cbind(time,t(ctel_normalised_sorted)))
ctel_clean_long <- gather(ctel_clean, gene_ID, expression, -time)

ctel_clean_long %>%
  ggplot(aes(x=time, y=expression)) +
  geom_smooth(method ='loess', size = 1, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "DESeq2 normalised expression", x = "developmental stage") +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
  scale_x_continuous(expand = c(0,0), limits = c(0.8,7.2), breaks = c(1:7),
                     labels = c("64 cells", "gastrula", "stage 4tt larva", 
                                "stage 5 larva", "stage 7 larva", 
                                "pre-competent larva", "competent larva")) +
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

## Dimorphilus
dgyr_clean <- data.frame(cbind(time2,t(dgyr_normalised_sorted)))
dgyr_clean_long <- gather(dgyr_clean, gene_ID, expression, -time2)

dgyr_clean_long %>%
  ggplot(aes(x=time2, y=expression)) +
  geom_smooth(method ='loess', size = 1, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "DESeq2 normalised expression", x = "developmental stage") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
  scale_x_continuous(expand = c(0,0), limits = c(0.8,4.2), breaks = c(1:4),
                     labels = c("early development","late development",
                                "hatchling","female adult")) +
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



##################################################################
### 2. Displaced genes between Owenia late and Capitella early ###
##################################################################

# Import list of genes to plot for all species
ofus_displaced <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/03-Fran_plots_expression/07-Owenia_late_Capitella_early_short.txt", header=T, sep ='\t')
all_displaced <- cbind(ofus_displaced[c(2,1)], ofus2ctel[ofus_displaced$transcript_id,"Ctel"], ofus2dgyr[ofus_displaced$transcript_id,"Dgyr"])
colnames(all_displaced) <- c("gene_symbol", "ofus_id", "ctel_id", "dgyr_id")
all_displaced$pfam_annotation <- tfclass[all_displaced$ofus_id, "PFAM_annotation"]
all_displaced <- all_displaced[order(all_displaced$pfam_annotation),]

write.table(all_displaced[-c(4)], "10-Displaced_TFs_Capi_early_Owenia_late_supp_table.txt",
            sep = '\t', quote = F, row.names = F)

# Expression dynamics of the TFs we have found in deseq2 levels (from 0 to max and z-scored)

rescale_custom <- function(x) (x/(max(x))) # define function to make 0 to max scale

ofus_displaced_deseq2 <- ofus_deseq2[ofus_deseq2$Gene_ID %in% all_displaced$ofus_id,]
rownames(ofus_displaced_deseq2) <- ofus_displaced_deseq2$Gene_ID
ctel_displaced_deseq2 <- ctel_deseq2[ctel_deseq2$Gene_ID %in% all_displaced$ctel_id,]
rownames(ctel_displaced_deseq2) <- ctel_displaced_deseq2$Gene_ID
dgyr_displaced_deseq2 <- dgyr_deseq2[dgyr_deseq2$Gene_ID %in% all_displaced$dgyr_id,]
rownames(dgyr_displaced_deseq2) <- dgyr_displaced_deseq2$Gene_ID

ofus_normalised <- t(apply(ofus_displaced_deseq2[-c(1)], 1, rescale_custom))
ctel_normalised <- t(apply(ctel_displaced_deseq2[-c(1)], 1, rescale_custom))
dgyr_normalised <- t(apply(dgyr_displaced_deseq2[-c(1)], 1, rescale_custom))

ofus_zscore <- t(scale(t(ofus_displaced_deseq2[-c(1)])))
ctel_zscore <- t(scale(t(ctel_displaced_deseq2[-c(1)])))
dgyr_zscore <- t(scale(t(dgyr_displaced_deseq2[-c(1)])))

ofus_order <- all_displaced$ofus_id
ctel_order <- all_displaced$ctel_id
dgyr_order <- all_displaced$dgyr_id

ofus_normalised_sorted <- ofus_normalised[match(ofus_order, rownames(ofus_normalised)),]
ctel_normalised_sorted <- ctel_normalised[match(ctel_order, rownames(ctel_normalised)),]
dgyr_normalised_sorted <- dgyr_normalised[match(dgyr_order, rownames(dgyr_normalised)),]

ofus_zscore_sorted <- ofus_zscore[match(ofus_order, rownames(ofus_zscore)),]
ctel_zscore_sorted <- ctel_zscore[match(ctel_order, rownames(ctel_zscore)),]
dgyr_zscore_sorted <- dgyr_zscore[match(dgyr_order, rownames(dgyr_zscore)),]

#1/2 rdbu (substitute of reds for 0 to 1)
heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color <- heatmap_color[50:100]

ofus_0_to_max <- ComplexHeatmap::Heatmap(ofus_normalised_sorted,
                                         cluster_columns = FALSE,
                                         cluster_rows = FALSE,
                                         col = heatmap_color,
                                         row_labels = all_displaced$gene_symbol,
                                         left_annotation = rowAnnotation(foo = anno_text(tfclass[all_displaced$ofus_id,"PFAM_annotation"])),
                                         column_labels = c("blastula", "gastrula", "elongation",
                                                           "early larva", "mitraria larva",
                                                           "competent larva", "juvenile"),
                                         heatmap_legend_param = list(color_bar = "continuous"))

ctel_0_to_max <- ComplexHeatmap::Heatmap(ctel_normalised_sorted,
                                         cluster_columns = FALSE,
                                         cluster_rows = FALSE,
                                         col = heatmap_color,
                                         column_labels = c("64 cells", "gastrula", "stage 4tt larva",
                                                           "stage 5 larva", "stage 7 larva",
                                                           "pre-competent larva", "competent larva"),
                                         row_labels = all_displaced$gene_symbol,
                                         heatmap_legend_param = list(color_bar = "continuous"))

dgyr_0_to_max <- ComplexHeatmap::Heatmap(dgyr_normalised_sorted,
                                         cluster_columns = FALSE,
                                         cluster_rows = FALSE,
                                         col = heatmap_color,
                                         row_labels = all_displaced$gene_symbol,
                                         column_labels = c("early development","late development",
                                                           "hatchling", "female adult"),
                                         heatmap_legend_param = list(color_bar = "continuous"))

ofus_0_to_max + ctel_0_to_max + dgyr_0_to_max
ofus_0_to_max + ctel_0_to_max



#full rdbu
heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

ofus_neg_to_pos <- ComplexHeatmap::Heatmap(ofus_zscore_sorted,
                                           cluster_columns = FALSE,
                                           cluster_rows = FALSE,
                                           col = heatmap_color,
                                           row_labels = all_displaced$gene_symbol,
                                           left_annotation = rowAnnotation(foo = anno_text(tfclass[all_displaced$ofus_id,"PFAM_annotation"])),
                                           column_labels = c("blastula", "gastrula", "elongation",
                                                             "early larva", "mitraria larva",
                                                             "competent larva", "juvenile"),
                                           heatmap_legend_param = list(color_bar = "continuous"))

ctel_neg_to_pos <- ComplexHeatmap::Heatmap(ctel_zscore_sorted,
                                           cluster_columns = FALSE,
                                           cluster_rows = FALSE,
                                           col = heatmap_color,
                                           column_labels = c("64 cells", "gastrula", "stage 4tt larva",
                                                             "stage 5 larva", "stage 7 larva",
                                                             "pre-competent larva", "competent larva"),
                                           row_labels = all_displaced$gene_symbol,
                                           heatmap_legend_param = list(color_bar = "continuous"))

dgyr_neg_to_pos <- ComplexHeatmap::Heatmap(dgyr_zscore_sorted,
                                           cluster_columns = FALSE,
                                           cluster_rows = FALSE,
                                           col = heatmap_color,
                                           row_labels = all_displaced$gene_symbol,
                                           column_labels = c("early development","late development",
                                                             "hatchling", "female adult"),
                                           heatmap_legend_param = list(color_bar = "continuous"))

ofus_neg_to_pos + ctel_neg_to_pos + dgyr_neg_to_pos
ofus_neg_to_pos + ctel_neg_to_pos


#alternative_viz
heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

ofus_0_to_max <- ComplexHeatmap::Heatmap(ofus_normalised_sorted,
                                         cluster_columns = FALSE,
                                         cluster_rows = FALSE,
                                         col = heatmap_color,
                                         row_labels = all_displaced$gene_symbol,
                                         left_annotation = rowAnnotation(foo = anno_text(tfclass[all_displaced$ofus_id,"PFAM_annotation"])),
                                         column_labels = c("blastula", "gastrula", "elongation",
                                                           "early larva", "mitraria larva",
                                                           "competent larva", "juvenile"),
                                         heatmap_legend_param = list(color_bar = "continuous"))

ctel_0_to_max <- ComplexHeatmap::Heatmap(ctel_normalised_sorted,
                                         cluster_columns = FALSE,
                                         cluster_rows = FALSE,
                                         col = heatmap_color,
                                         column_labels = c("64 cells", "gastrula", "stage 4tt larva",
                                                           "stage 5 larva", "stage 7 larva",
                                                           "pre-competent larva", "competent larva"),
                                         row_labels = all_displaced$gene_symbol,
                                         heatmap_legend_param = list(color_bar = "continuous"))

dgyr_0_to_max <- ComplexHeatmap::Heatmap(dgyr_normalised_sorted,
                                         cluster_columns = FALSE,
                                         cluster_rows = FALSE,
                                         col = heatmap_color,
                                         row_labels = all_displaced$gene_symbol,
                                         column_labels = c("early development","late development",
                                                           "hatchling", "female adult"),
                                         heatmap_legend_param = list(color_bar = "continuous"))

ofus_0_to_max + ctel_0_to_max + dgyr_0_to_max
ofus_0_to_max + ctel_0_to_max



# Relative expression lineplots

species_palette <- viridis(3)
time <- c(1:7)
time2 <- c(1:4)

## Owenia
ofus_clean <- data.frame(cbind(time,t(ofus_normalised_sorted)))
ofus_clean_long <- gather(ofus_clean, gene_ID, expression, -time)

ofus_clean_long %>%
  ggplot(aes(x=time, y=expression)) +
  geom_smooth(method ='loess', size = 1, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "DESeq2 normalised expression", x = "developmental stage") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
  scale_x_continuous(expand = c(0,0), limits = c(0.8,7.2), breaks = c(1:7),
                     labels = c("blastula", "gastrula", "elongation",
                                "early larva", "mitraria larva", "competent larva", 
                                "juvenile")) +
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


## Capitella
ctel_clean <- data.frame(cbind(time,t(ctel_normalised_sorted)))
ctel_clean_long <- gather(ctel_clean, gene_ID, expression, -time)

ctel_clean_long %>%
  ggplot(aes(x=time, y=expression)) +
  geom_smooth(method ='loess', size = 1, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "DESeq2 normalised expression", x = "developmental stage") +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
  scale_x_continuous(expand = c(0,0), limits = c(0.8,7.2), breaks = c(1:7),
                     labels = c("64 cells", "gastrula", "stage 4tt larva", 
                                "stage 5 larva", "stage 7 larva", 
                                "pre-competent larva", "competent larva")) +
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

## Dimorphilus
dgyr_clean <- data.frame(cbind(time2,t(dgyr_normalised_sorted)))
dgyr_clean_long <- gather(dgyr_clean, gene_ID, expression, -time2)

dgyr_clean_long %>%
  ggplot(aes(x=time2, y=expression)) +
  geom_smooth(method ='loess', size = 1, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "DESeq2 normalised expression", x = "developmental stage") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
  scale_x_continuous(expand = c(0,0), limits = c(0.8,4.2), breaks = c(1:4),
                     labels = c("early development","late development",
                                "hatchling","female adult")) +
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










##################################################################
### 3. Displaced genes between Owenia early and Capitella late ###
##################################################################

# Import list of genes to plot for all species
ofus_displaced <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/03-Fran_plots_expression/12-Owenia_early_Capitella_late_short.txt", header=T, sep ='\t')
all_displaced <- cbind(ofus_displaced[c(2,1)], ofus2ctel[ofus_displaced$transcript_id,"Ctel"], ofus2dgyr[ofus_displaced$transcript_id,"Dgyr"])
colnames(all_displaced) <- c("gene_symbol", "ofus_id", "ctel_id", "dgyr_id")
all_displaced$pfam_annotation <- tfclass[all_displaced$ofus_id, "PFAM_annotation"]
all_displaced <- all_displaced[order(all_displaced$pfam_annotation),]

write.table(all_displaced[-c(4)], "10b-Displaced_TFs_Capi_late_Owenia_late_supp_table.txt",
            sep = '\t', quote = F, row.names = F)

# Expression dynamics of the TFs we have found in deseq2 levels (from 0 to max and z-scored)

rescale_custom <- function(x) (x/(max(x))) # define function to make 0 to max scale

ofus_displaced_deseq2 <- ofus_deseq2[ofus_deseq2$Gene_ID %in% all_displaced$ofus_id,]
rownames(ofus_displaced_deseq2) <- ofus_displaced_deseq2$Gene_ID
ctel_displaced_deseq2 <- ctel_deseq2[ctel_deseq2$Gene_ID %in% all_displaced$ctel_id,]
rownames(ctel_displaced_deseq2) <- ctel_displaced_deseq2$Gene_ID
dgyr_displaced_deseq2 <- dgyr_deseq2[dgyr_deseq2$Gene_ID %in% all_displaced$dgyr_id,]
rownames(dgyr_displaced_deseq2) <- dgyr_displaced_deseq2$Gene_ID

ofus_normalised <- t(apply(ofus_displaced_deseq2[-c(1)], 1, rescale_custom))
ctel_normalised <- t(apply(ctel_displaced_deseq2[-c(1)], 1, rescale_custom))
dgyr_normalised <- t(apply(dgyr_displaced_deseq2[-c(1)], 1, rescale_custom))

ofus_zscore <- t(scale(t(ofus_displaced_deseq2[-c(1)])))
ctel_zscore <- t(scale(t(ctel_displaced_deseq2[-c(1)])))
dgyr_zscore <- t(scale(t(dgyr_displaced_deseq2[-c(1)])))

ofus_order <- all_displaced$ofus_id
ctel_order <- all_displaced$ctel_id
dgyr_order <- all_displaced$dgyr_id

ofus_normalised_sorted <- ofus_normalised[match(ofus_order, rownames(ofus_normalised)),]
ctel_normalised_sorted <- ctel_normalised[match(ctel_order, rownames(ctel_normalised)),]
dgyr_normalised_sorted <- dgyr_normalised[match(dgyr_order, rownames(dgyr_normalised)),]

ofus_zscore_sorted <- ofus_zscore[match(ofus_order, rownames(ofus_zscore)),]
ctel_zscore_sorted <- ctel_zscore[match(ctel_order, rownames(ctel_zscore)),]
dgyr_zscore_sorted <- dgyr_zscore[match(dgyr_order, rownames(dgyr_zscore)),]

#1/2 rdbu (substitute of reds for 0 to 1)
heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color <- heatmap_color[50:100]

ofus_0_to_max <- ComplexHeatmap::Heatmap(ofus_normalised_sorted,
                                         cluster_columns = FALSE,
                                         cluster_rows = FALSE,
                                         col = heatmap_color,
                                         row_labels = all_displaced$gene_symbol,
                                         left_annotation = rowAnnotation(foo = anno_text(tfclass[all_displaced$ofus_id,"PFAM_annotation"])),
                                         column_labels = c("blastula", "gastrula", "elongation",
                                                           "early larva", "mitraria larva",
                                                           "competent larva", "juvenile"),
                                         heatmap_legend_param = list(color_bar = "continuous"))

ctel_0_to_max <- ComplexHeatmap::Heatmap(ctel_normalised_sorted,
                                         cluster_columns = FALSE,
                                         cluster_rows = FALSE,
                                         col = heatmap_color,
                                         column_labels = c("64 cells", "gastrula", "stage 4tt larva",
                                                           "stage 5 larva", "stage 7 larva",
                                                           "pre-competent larva", "competent larva"),
                                         row_labels = all_displaced$gene_symbol,
                                         heatmap_legend_param = list(color_bar = "continuous"))

dgyr_0_to_max <- ComplexHeatmap::Heatmap(dgyr_normalised_sorted,
                                         cluster_columns = FALSE,
                                         cluster_rows = FALSE,
                                         col = heatmap_color,
                                         row_labels = all_displaced$gene_symbol,
                                         column_labels = c("early development","late development",
                                                           "hatchling", "female adult"),
                                         heatmap_legend_param = list(color_bar = "continuous"))

ofus_0_to_max + ctel_0_to_max + dgyr_0_to_max
ofus_0_to_max + ctel_0_to_max



#full rdbu
heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

ofus_neg_to_pos <- ComplexHeatmap::Heatmap(ofus_zscore_sorted,
                                           cluster_columns = FALSE,
                                           cluster_rows = FALSE,
                                           col = heatmap_color,
                                           row_labels = all_displaced$gene_symbol,
                                           left_annotation = rowAnnotation(foo = anno_text(tfclass[all_displaced$ofus_id,"PFAM_annotation"])),
                                           column_labels = c("blastula", "gastrula", "elongation",
                                                             "early larva", "mitraria larva",
                                                             "competent larva", "juvenile"),
                                           heatmap_legend_param = list(color_bar = "continuous"))

ctel_neg_to_pos <- ComplexHeatmap::Heatmap(ctel_zscore_sorted,
                                           cluster_columns = FALSE,
                                           cluster_rows = FALSE,
                                           col = heatmap_color,
                                           column_labels = c("64 cells", "gastrula", "stage 4tt larva",
                                                             "stage 5 larva", "stage 7 larva",
                                                             "pre-competent larva", "competent larva"),
                                           row_labels = all_displaced$gene_symbol,
                                           heatmap_legend_param = list(color_bar = "continuous"))

dgyr_neg_to_pos <- ComplexHeatmap::Heatmap(dgyr_zscore_sorted,
                                           cluster_columns = FALSE,
                                           cluster_rows = FALSE,
                                           col = heatmap_color,
                                           row_labels = all_displaced$gene_symbol,
                                           column_labels = c("early development","late development",
                                                             "hatchling", "female adult"),
                                           heatmap_legend_param = list(color_bar = "continuous"))

ofus_neg_to_pos + ctel_neg_to_pos + dgyr_neg_to_pos
ofus_neg_to_pos + ctel_neg_to_pos


#alternative_viz
heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

ofus_0_to_max <- ComplexHeatmap::Heatmap(ofus_normalised_sorted,
                                         cluster_columns = FALSE,
                                         cluster_rows = FALSE,
                                         col = heatmap_color,
                                         row_labels = all_displaced$gene_symbol,
                                         left_annotation = rowAnnotation(foo = anno_text(tfclass[all_displaced$ofus_id,"PFAM_annotation"])),
                                         column_labels = c("blastula", "gastrula", "elongation",
                                                           "early larva", "mitraria larva",
                                                           "competent larva", "juvenile"),
                                         heatmap_legend_param = list(color_bar = "continuous"))

ctel_0_to_max <- ComplexHeatmap::Heatmap(ctel_normalised_sorted,
                                         cluster_columns = FALSE,
                                         cluster_rows = FALSE,
                                         col = heatmap_color,
                                         column_labels = c("64 cells", "gastrula", "stage 4tt larva",
                                                           "stage 5 larva", "stage 7 larva",
                                                           "pre-competent larva", "competent larva"),
                                         row_labels = all_displaced$gene_symbol,
                                         heatmap_legend_param = list(color_bar = "continuous"))

dgyr_0_to_max <- ComplexHeatmap::Heatmap(dgyr_normalised_sorted,
                                         cluster_columns = FALSE,
                                         cluster_rows = FALSE,
                                         col = heatmap_color,
                                         row_labels = all_displaced$gene_symbol,
                                         column_labels = c("early development","late development",
                                                           "hatchling", "female adult"),
                                         heatmap_legend_param = list(color_bar = "continuous"))

ofus_0_to_max + ctel_0_to_max + dgyr_0_to_max
ofus_0_to_max + ctel_0_to_max



# Relative expression lineplots

species_palette <- viridis(3)
time <- c(1:7)
time2 <- c(1:4)

## Owenia
ofus_clean <- data.frame(cbind(time,t(ofus_normalised_sorted)))
ofus_clean_long <- gather(ofus_clean, gene_ID, expression, -time)

ofus_clean_long %>%
  ggplot(aes(x=time, y=expression)) +
  geom_smooth(method ='loess', size = 1, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "DESeq2 normalised expression", x = "developmental stage") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
  scale_x_continuous(expand = c(0,0), limits = c(0.8,7.2), breaks = c(1:7),
                     labels = c("blastula", "gastrula", "elongation",
                                "early larva", "mitraria larva", "competent larva", 
                                "juvenile")) +
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


## Capitella
ctel_clean <- data.frame(cbind(time,t(ctel_normalised_sorted)))
ctel_clean_long <- gather(ctel_clean, gene_ID, expression, -time)

ctel_clean_long %>%
  ggplot(aes(x=time, y=expression)) +
  geom_smooth(method ='loess', size = 1, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "DESeq2 normalised expression", x = "developmental stage") +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
  scale_x_continuous(expand = c(0,0), limits = c(0.8,7.2), breaks = c(1:7),
                     labels = c("64 cells", "gastrula", "stage 4tt larva", 
                                "stage 5 larva", "stage 7 larva", 
                                "pre-competent larva", "competent larva")) +
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

## Dimorphilus
dgyr_clean <- data.frame(cbind(time2,t(dgyr_normalised_sorted)))
dgyr_clean_long <- gather(dgyr_clean, gene_ID, expression, -time2)

dgyr_clean_long %>%
  ggplot(aes(x=time2, y=expression)) +
  geom_smooth(method ='loess', size = 1, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "DESeq2 normalised expression", x = "developmental stage") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
  scale_x_continuous(expand = c(0,0), limits = c(0.8,4.2), breaks = c(1:4),
                     labels = c("early development","late development",
                                "hatchling","female adult")) +
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
