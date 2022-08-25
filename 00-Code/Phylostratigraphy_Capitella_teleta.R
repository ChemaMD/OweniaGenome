library(ggplot2)
library(ggpubr)
library(gplots)
library(RColorBrewer)
library(ggridges)
library(tidyr)
library(dplyr)
library(preprocessCore)
library(pheatmap)
library(ComplexHeatmap)

raw_data <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/05-Phylostratigraphy_with_original_species_Owenia_and_Capitella/02-Capitella_teleta/01-Raw_data_Capitella_teleta.txt", 
                       as.is=TRUE, header=TRUE)
rownames(raw_data) <- t(raw_data[c(1)])
raw_data_clean <- raw_data[-c(1)]
colnames(raw_data_clean)[2] <- "64_cells"


#QUANTILE NORMALISATION
#In here I did the normalisation for all the genes, I think this is the appropriate thing to do.
norm <- data.frame(cbind(raw_data_clean[1],normalize.quantiles(data.matrix(raw_data_clean[,2:ncol(raw_data_clean)]))))
norm$origin[norm$origin == "Ctel-2"] <- "Ctel"
norm$origin[norm$origin == "Nephrozoa" | norm$origin == "Bilateria"] <- "Nephrozoa_Bilateria"

colnames(norm) <- c("origin","64_cells","gastrula","st4tt_larva","st5_larva",
                    "st7_larva","precompetent_larva","competent_larva")
norm_tidy <- gather(norm, stage, expr, -origin)
summary(norm)

ctel_norm <- norm_tidy %>% filter(grepl('Ctel', origin))
metazoa_norm <- norm_tidy %>% filter(grepl('Metazoa', origin))
nephrozoa_bilateria_norm <- norm_tidy %>% filter(grepl('Nephrozoa_Bilateria', origin))
protostomia_norm <- norm_tidy %>% filter(grepl('Protostomia', origin))
spiralia_norm <- norm_tidy %>% filter(grepl('Lophotrochozoa', origin))
annelida_norm <- norm_tidy %>% filter(grepl('Annelida', origin))
other_norm <- norm_tidy %>% filter(grepl('others', origin))

ctel_norm$stage <- factor(ctel_norm$stage, levels = unique(ctel_norm$stage))
metazoa_norm$stage <- factor(metazoa_norm$stage, levels = unique(metazoa_norm$stage))
nephrozoa_bilateria_norm$stage <- factor(nephrozoa_bilateria_norm$stage, levels = unique(nephrozoa_bilateria_norm$stage))
protostomia_norm$stage <- factor(protostomia_norm$stage, levels = unique(protostomia_norm$stage))
spiralia_norm$stage <- factor(spiralia_norm$stage, levels = unique(spiralia_norm$stage))
annelida_norm$stage <- factor(annelida_norm$stage, levels = unique(annelida_norm$stage))
other_norm$stage <- factor(other_norm$stage, levels = unique(other_norm$stage))

norm_75_quant <- data.frame(tapply(norm$`64_cells`, norm$origin, quantile, p = 0.75, na.rm = TRUE),
                            tapply(norm$gastrula, norm$origin, quantile, p = 0.75, na.rm = TRUE),
                            tapply(norm$st4tt_larva, norm$origin, quantile, p = 0.75, na.rm = TRUE),
                            tapply(norm$st5_larva, norm$origin, quantile, p = 0.75, na.rm = TRUE),
                            tapply(norm$st7_larva, norm$origin, quantile, p = 0.75, na.rm = TRUE),
                            tapply(norm$precompetent_larva, norm$origin, quantile, p = 0.75, na.rm = TRUE),
                            tapply(norm$competent_larva, norm$origin, quantile, p = 0.75, na.rm = TRUE))
colnames(norm_75_quant) <- c("64_cells","gastrula","st4tt_larva","st5_larva",
                             "st7_larva","precompetent_larva","competent_larva")
rownames(norm_75_quant) <- c("Annelida","Capitella_teleta","Lophotrochozoa","Metazoa",
                             "Nephrozoa_Bilateria","Other","Protostomia")
norm_75_quant <- norm_75_quant[c(4,5,7,3,1,2,6),]


heatmap_color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

#Heatmap 1: Scaled by row
pheatmap(norm_75_quant, 
         scale = "row",
         cluster_rows=F,
         cluster_cols=F, 
         show_rownames=T, 
         fontsize= 12,
         col=heatmap_color, 
         border_color=NA, 
         cellwidth = 20, 
         cellheight = 20)

#Heatmap 2: Scaled by column
pheatmap(norm_75_quant, 
         scale = "column",
         cluster_rows=F,
         cluster_cols=F, 
         show_rownames=T, 
         fontsize= 12,
         col=heatmap_color, 
         border_color=NA, 
         cellwidth = 20, 
         cellheight = 20)

metazoa_plot <- ggplot(metazoa_norm, aes(x = stage, y = expr)) +
  geom_boxplot(aes(x = stage, y = expr, fill = stage), outlier.shape = NA,
               width = 0.7) +
  scale_fill_manual(values = c("#563B45","#D2A359","#C6613A",
                               "#0B6CA0","#032C65","#507810",
                               "#204858")) +
  coord_cartesian(ylim = c(0,3500)) +
  theme_classic() +
  labs(y = "quantile-normalised gene expression") +
  theme(axis.title.x = element_text(color="black",size=10),
        axis.title.y = element_text(color="black",size=10),
        axis.text.x = element_text(color="black",size=12,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=10,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"),
        legend.position = "none")

nephrozoa_bilateria_plot <- ggplot(nephrozoa_bilateria_norm, aes(x = stage, y = expr)) +
  geom_boxplot(aes(x = stage, y = expr, fill = stage), outlier.shape = NA,
               width = 0.7) +
  scale_fill_manual(values = c("#563B45","#D2A359","#C6613A",
                               "#0B6CA0","#032C65","#507810",
                               "#204858")) +
  coord_cartesian(ylim = c(0,1500)) +
  theme_classic() +
  labs(y = "quantile-normalised gene expression") +
  theme(axis.title.x = element_text(color="black",size=10),
        axis.title.y = element_text(color="black",size=10),
        axis.text.x = element_text(color="black",size=12,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=10,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"),
        legend.position = "none")

protostomia_plot <- ggplot(protostomia_norm, aes(x = stage, y = expr)) +
  geom_boxplot(aes(x = stage, y = expr, fill = stage), outlier.shape = NA,
               width = 0.7) +
  scale_fill_manual(values = c("#563B45","#D2A359","#C6613A",
                               "#0B6CA0","#032C65","#507810",
                               "#204858")) +
  coord_cartesian(ylim = c(0,700)) +
  theme_classic() +
  labs(y = "quantile-normalised gene expression") +
  theme(axis.title.x = element_text(color="black",size=10),
        axis.title.y = element_text(color="black",size=10),
        axis.text.x = element_text(color="black",size=12,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=10,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"),
        legend.position = "none")

spiralia_plot <- ggplot(spiralia_norm, aes(x = stage, y = expr)) +
  geom_boxplot(aes(x = stage, y = expr, fill = stage), outlier.shape = NA,
               width = 0.7) +
  scale_fill_manual(values = c("#563B45","#D2A359","#C6613A",
                               "#0B6CA0","#032C65","#507810",
                               "#204858")) +
  coord_cartesian(ylim = c(0,1100)) +
  theme_classic() +
  labs(y = "quantile-normalised gene expression") +
  theme(axis.title.x = element_text(color="black",size=10),
        axis.title.y = element_text(color="black",size=10),
        axis.text.x = element_text(color="black",size=12,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=10,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"),
        legend.position = "none")

annelida_plot <- ggplot(annelida_norm, aes(x = stage, y = expr)) +
  geom_boxplot(aes(x = stage, y = expr, fill = stage), outlier.shape = NA,
               width = 0.7) +
  scale_fill_manual(values = c("#563B45","#D2A359","#C6613A",
                               "#0B6CA0","#032C65","#507810",
                               "#204858")) +
  coord_cartesian(ylim = c(0,1000)) +
  theme_classic() +
  labs(y = "quantile-normalised gene expression") +
  theme(axis.title.x = element_text(color="black",size=10),
        axis.title.y = element_text(color="black",size=10),
        axis.text.x = element_text(color="black",size=12,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=10,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"),
        legend.position = "none")

capitella_teleta_plot <- ggplot(ctel_norm, aes(x = stage, y = expr)) +
  geom_boxplot(aes(x = stage, y = expr, fill = stage), outlier.shape = NA,
               width = 0.7) +
  scale_fill_manual(values = c("#563B45","#D2A359","#C6613A",
                               "#0B6CA0","#032C65","#507810",
                               "#204858")) +
  coord_cartesian(ylim = c(0,350)) +
  theme_classic() +
  labs(y = "quantile-normalised gene expression") +
  theme(axis.title.x = element_text(color="black",size=10),
        axis.title.y = element_text(color="black",size=10),
        axis.text.x = element_text(color="black",size=12,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=10,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"),
        legend.position = "none")


other_plot <- ggplot(other_norm, aes(x = stage, y = expr)) +
  geom_boxplot(aes(x = stage, y = expr, fill = stage), outlier.shape = NA,
               width = 0.7) +
  scale_fill_manual(values = c("#563B45","#D2A359","#C6613A",
                               "#0B6CA0","#032C65","#507810",
                               "#204858")) +
  coord_cartesian(ylim = c(0,600)) +
  theme_classic() +
  labs(y = "quantile-normalised gene expression") +
  theme(axis.title.x = element_text(color="black",size=10),
        axis.title.y = element_text(color="black",size=10),
        axis.text.x = element_text(color="black",size=12,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=10,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"),
        legend.position = "none")

#Add random selection of 2,000 genes as control
random <- norm[sample(nrow(norm),2000),]
random_norm <- gather(random, stage, expr, -origin)
random_norm$stage <- factor(random_norm$stage, levels = unique(random_norm$stage))

random_plot <- ggplot(random_norm, aes(x = stage, y = expr)) +
  geom_boxplot(aes(x = stage, y = expr, fill = stage), outlier.shape = NA,
               width = 0.7) +
  scale_fill_manual(values = c("#563B45","#D2A359","#C6613A",
                               "#0B6CA0","#032C65","#507810",
                               "#204858")) +
  coord_cartesian(ylim = c(0,1200)) +
  theme_classic() +
  labs(y = "quantile-normalised gene expression") +
  theme(axis.title.x = element_text(color="black",size=10),
        axis.title.y = element_text(color="black",size=10),
        axis.text.x = element_text(color="black",size=12,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=10,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"),
        legend.position = "none")

ggarrange(metazoa_plot, nephrozoa_bilateria_plot, protostomia_plot,
          spiralia_plot, annelida_plot, capitella_teleta_plot, 
          other_plot, random_plot,
          ncol = 4, nrow = 2)


# Phylostratigraphy analysis by RNA-seq mfuzz clusters
# Check if different lineages of genes are under- or overrepresented in the
# different gene clusters we determined by mfuzz clustering for RNA-seq
# using a Fisher's exact test (two-tailed)

raw_data <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/05-Phylostratigraphy_with_original_species_Owenia_and_Capitella/02-Capitella_teleta/01-Raw_data_Capitella_teleta.txt", 
                       as.is=TRUE, header=TRUE, row.names = 1)
lineages <- cbind(rownames(raw_data),raw_data[1])
colnames(lineages) <- c("gene_ID","origin")

clusters <- read.table("/Users/franciscomanuelmartinzamora/Dropbox/02-OweniaGenome/06-Revision/00-DATA/01-Capitella_RNAseq_core_analyses/06-Capitella_teleta_clusters_annotation_corrected.txt",
                       as.is = T, header = T, row.names = 1)

lineages_expr <- cbind(lineages[rownames(clusters),],clusters$Cluster_corrected)
colnames(lineages_expr) <- c("gene_ID","origin","cluster")
rownames(lineages_expr) <- c(1:nrow(lineages_expr))
lineages_expr$origin[lineages_expr$origin == "others"] <- "Other"
lineages_expr$origin[lineages_expr$origin == "Ctel-2"] <- "Capitella_teleta"
lineages$origin[lineages$origin == "Ctel"] <- "Capitella_teleta"
lineages_expr$origin[lineages_expr$origin == "Nephrozoa" | lineages_expr$origin == "Bilateria"] <- "Nephrozoa_Bilateria"


pvalues <- data.frame(Metazoa=numeric(), Nephrozoa_Bilateria=numeric(), Protostomia=numeric(),
                      Lophotrochozoa=numeric(), Annelida=numeric(), Capitella_teleta=numeric(), Other=numeric())
oddsratio <- data.frame(Metazoa=numeric(), Nephrozoa_Bilateria=numeric(), Protostomia=numeric(),
                        Lophotrochozoa=numeric(), Annelida=numeric(), Capitella_teleta=numeric(), Other=numeric())

for (i in c("Metazoa","Nephrozoa_Bilateria","Protostomia","Lophotrochozoa","Annelida","Capitella_teleta","Other")){
  for (j in c(1:12)){
    test <- data.frame("in_lineage" = c(length(which(lineages_expr$cluster == j & lineages_expr$origin == i)),
                                        length(which(lineages_expr$cluster != j & lineages_expr$origin == i))),
                       "not_in_lineage" = c(length(which(lineages_expr$cluster == j & lineages_expr$origin != i)),
                                            length(which(lineages_expr$cluster != j & lineages_expr$origin != i))),
                       row.names = c("in_cluster", "not_in_cluster"),
                       stringsAsFactors = FALSE)
    oddsratio[j,i] <- fisher.test(test)$estimate[1]
    if (fisher.test(test)$estimate[1] > 1){
      pvalues[j,i] <- fisher.test(test)$p.value
    } else {
      pvalues[j,i] <- -fisher.test(test)$p.value
    }
}
}

write.table(oddsratio, sep="\t", file="Odds_ratio.txt")

pvalues$Metazoa <- p.adjust(pvalues$Metazoa, method = "bonferroni")
pvalues$Nephrozoa_Bilateria <- p.adjust(pvalues$Nephrozoa_Bilateria, method = "bonferroni")
pvalues$Protostomia <- p.adjust(pvalues$Protostomia, method = "bonferroni")
pvalues$Lophotrochozoa <- p.adjust(pvalues$Lophotrochozoa, method = "bonferroni")
pvalues$Annelida <- p.adjust(pvalues$Annelida, method = "bonferroni")
pvalues$Capitella_teleta <- p.adjust(pvalues$Capitella_teleta, method = "bonferroni")
pvalues$Other <- p.adjust(pvalues$Other, method = "bonferroni")

pvalues[pvalues>0.05] <- 1
pvalues[pvalues<(-0.05)] <- 1
pvalues_log <- data.matrix(-log10(abs(pvalues)))

for (i in c(1:12)){
  for (j in c(1:7)){
    if (pvalues[i,j] < 0) {
      pvalues_log[i,j] <- -pvalues_log[i,j]
    }
  }
}

paletteLength <- 100
myBreaks <- c(seq(min(pvalues_log), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(pvalues_log)/paletteLength, max(pvalues_log), length.out=floor(paletteLength/2)))

myBreaks <- c(seq(-50, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(50/paletteLength, 50, length.out=floor(paletteLength/2)))

heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color[51] <- rgb(1,1,1)

pheatmap(t(pvalues_log),
         cluster_rows = F,
         cluster_cols = F,
         border_color = NA,
         color = heatmap_color,
         breaks = myBreaks,
         display_numbers = T)



#E. Bubbles proportional to gene size
library(grid)

count <- as.data.frame(table(norm$origin))
count_table <- norm_75_quant
count_table["Metazoa",] <- c(count$Freq[count$Var1 == "Metazoa"])
count_table["Nephrozoa_Bilateria",] <- c(count$Freq[count$Var1 == "Nephrozoa_Bilateria"])
count_table["Protostomia",] <- c(count$Freq[count$Var1 == "Protostomia"])
count_table["Lophotrochozoa",] <- c(count$Freq[count$Var1 == "Lophotrochozoa"])
count_table["Annelida",] <- c(count$Freq[count$Var1 == "Annelida"])
count_table["Capitella_teleta",] <- c(count$Freq[count$Var1 == "Ctel"])
count_table["Other",] <- c(count$Freq[count$Var1 == "others"])
count_table <- round(count_table)

write.table(count_table[1], "02-Capitella_teleta/07-Phylostrata_gene_count.txt", sep ='\t', quote = F, col.names = F)
