library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(corrplot)
library(RColorBrewer)
library(Mfuzz)
library(ComplexHeatmap)
library(viridis)
library(gplots)
library(tidyr)
library(matrixStats)
library(data.table)
library(mgcv)
library(scales)
library(preprocessCore)



#######################################################################
####################  1. Amphimedon queenslandica #####################
#######################################################################


# Input data
stages <- c("cleavage","brown_stage","spot_stage","ring_stage","precompetent_larva",
            "competent_larva","early_postsettlement","postsettlement","juvenile")
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/03-DESeq2_expression_data/")
df_DESeq2 <- read.table("Aque_DESeq2_average_Metazoa_and_TRGs.txt", header = TRUE)
rownames(df_DESeq2) <- df_DESeq2$Gene_ID
colnames(df_DESeq2) <- c("Gene_ID",stages)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/04-Metazoa_and_TRGs/")
df_origin <- read.table("Aque_Metazoa_and_TRGs_gene_ID.txt", header = F)
colnames(df_origin) <- c("origin","Gene_ID")
df_origin <- df_origin[!is.na(df_origin$Gene_ID),]
rownames(df_origin) <- df_origin$Gene_ID

# Subset genes with non-null expression and with metazoan/species-specific origin
df_notnull <- df_DESeq2[rowSums(df_DESeq2[,-c(1)])>0,]
df_notnull$origin <- NA
df_notnull[df_notnull$Gene_ID,'origin'] <- df_origin[df_notnull$Gene_ID,'origin']
df_notnull$origin[df_notnull$origin == 'Aque'] <- "Species-specific"
df_withorigin <- df_notnull[!is.na(df_notnull$origin),]

# Quantile normalise expression data
df_withorigin_qn <- cbind(df_withorigin[c(1,ncol(df_withorigin))],normalize.quantiles(data.matrix(df_withorigin[-c(1,ncol(df_withorigin))]))) # quantile normalisation of data
df_withorigin_qn$origin <- factor(df_withorigin_qn$origin, levels = c("Metazoa","Species-specific"))
colnames(df_withorigin_qn)[3:ncol(df_withorigin_qn)] <- stages

# Obtain 75% quantile for data distribution representation purposes in heatmaps and export
quantile_75 <- data.frame(matrix(ncol = length(stages), nrow = 2))
rownames(quantile_75) <- c("Metazoa","Species-specific")
colnames(quantile_75) <- stages

for (i in stages){
  quantile_75[i] <- round(tapply(df_withorigin_qn[,i], df_withorigin_qn$origin, quantile, p = 0.75, na.rm = TRUE),2)
}

quantile_75_clean <- cbind('origin'=rownames(quantile_75),quantile_75)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/08-quantile_normalised_heatmaps/")
write.table(quantile_75_clean, "01-Aque_quantile_75.txt", quote = F, sep = '\t', row.names = F)

# Transformations for boxplot
df_clean <- data.frame(cbind(stages,t(df_withorigin_qn[-c(1,2)])))
df_clean_long <- gather(df_clean, gene_ID, expression, -stages)
df_clean_long$origin <- rep(df_withorigin_qn$origin, each = length(stages))
df_clean_long$expression <- as.numeric(df_clean_long$expression)
df_clean_long$origin <- factor(df_clean_long$origin, levels = c("Metazoa","Species-specific"))
df_clean_long$stages <- factor(df_clean_long$stages, levels = stages)
df_clean_long_metazoa <- df_clean_long[df_clean_long$origin == 'Metazoa',]
df_clean_long_specific <- df_clean_long[df_clean_long$origin == 'Species-specific',]

# Boxplots
metazoa_boxplot <- ggplot(df_clean_long_metazoa, aes(x=stages, y=expression, fill=origin)) +
  geom_boxplot(outlier.shape = NA, fill = '#FF0000') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "quantile-normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,1800)) +
  scale_y_continuous(breaks = c(seq(0,1800,200))) +
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

specific_boxplot <- ggplot(df_clean_long_specific, aes(x=stages, y=expression, fill=origin)) +
  geom_boxplot(outlier.shape = NA, fill = '#147DF5') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "quantile-normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,180)) +
  scale_y_continuous(breaks = c(seq(0,180,20))) +
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

ggarrange(metazoa_boxplot, specific_boxplot)





#######################################################################
######################  2. Clytia hemisphaerica #######################
#######################################################################


# Input data
stages <- c("early_gastrula","early_planula","mid_planula","late_planula",
                      "primary_polyp","gastrozooid","gonozooid","stolon","baby_medusa")
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/03-DESeq2_expression_data/")
df_DESeq2 <- read.table("Chem_DESeq2_average_Metazoa_and_TRGs.txt", header = TRUE)
rownames(df_DESeq2) <- df_DESeq2$Gene_ID
colnames(df_DESeq2) <- c("Gene_ID",stages)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/04-Metazoa_and_TRGs/")
df_origin <- read.table("Chem_Metazoa_and_TRGs_gene_ID.txt", header = F)
colnames(df_origin) <- c("origin","Gene_ID")
df_origin <- df_origin[!is.na(df_origin$Gene_ID),]
rownames(df_origin) <- df_origin$Gene_ID

# Subset genes with non-null expression and with metazoan/species-specific origin
df_notnull <- df_DESeq2[rowSums(df_DESeq2[,-c(1)])>0,]
df_notnull$origin <- NA
df_notnull[df_notnull$Gene_ID,'origin'] <- df_origin[df_notnull$Gene_ID,'origin']
df_notnull$origin[df_notnull$origin == 'Chem'] <- "Species-specific"
df_withorigin <- df_notnull[!is.na(df_notnull$origin),]

# Quantile normalise expression data
df_withorigin_qn <- cbind(df_withorigin[c(1,ncol(df_withorigin))],normalize.quantiles(data.matrix(df_withorigin[-c(1,ncol(df_withorigin))]))) # quantile normalisation of data
df_withorigin_qn$origin <- factor(df_withorigin_qn$origin, levels = c("Metazoa","Species-specific"))
colnames(df_withorigin_qn)[3:ncol(df_withorigin_qn)] <- stages

# Obtain 75% quantile for data distribution representation purposes in heatmaps and export
quantile_75 <- data.frame(matrix(ncol = length(stages), nrow = 2))
rownames(quantile_75) <- c("Metazoa","Species-specific")
colnames(quantile_75) <- stages

for (i in stages){
  quantile_75[i] <- round(tapply(df_withorigin_qn[,i], df_withorigin_qn$origin, quantile, p = 0.75, na.rm = TRUE),2)
}

quantile_75_clean <- cbind('origin'=rownames(quantile_75),quantile_75)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/08-quantile_normalised_heatmaps/")
write.table(quantile_75_clean, "02-Chem_quantile_75.txt", quote = F, sep = '\t', row.names = F)

# Transformations for boxplot
df_clean <- data.frame(cbind(stages,t(df_withorigin_qn[-c(1,2)])))
df_clean_long <- gather(df_clean, gene_ID, expression, -stages)
df_clean_long$origin <- rep(df_withorigin_qn$origin, each = length(stages))
df_clean_long$expression <- as.numeric(df_clean_long$expression)
df_clean_long$origin <- factor(df_clean_long$origin, levels = c("Metazoa","Species-specific"))
df_clean_long$stages <- factor(df_clean_long$stages, levels = stages)
df_clean_long_metazoa <- df_clean_long[df_clean_long$origin == 'Metazoa',]
df_clean_long_specific <- df_clean_long[df_clean_long$origin == 'Species-specific',]

# Boxplots
metazoa_boxplot <- ggplot(df_clean_long_metazoa, aes(x=stages, y=expression, fill=origin)) +
  geom_boxplot(outlier.shape = NA, fill = '#FF0000') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "quantile-normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,6000)) +
  scale_y_continuous(breaks = c(seq(0,6000,500))) +
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

specific_boxplot <- ggplot(df_clean_long_specific, aes(x=stages, y=expression, fill=origin)) +
  geom_boxplot(outlier.shape = NA, fill = '#147DF5') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "quantile-normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,1200)) +
  scale_y_continuous(breaks = c(seq(0,1200,100))) +
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

ggarrange(metazoa_boxplot, specific_boxplot)




#######################################################################
#####################  3. Nematostella vectensis ######################
#######################################################################


# Input data
stages <- c("24 hpf","48 hpf","72 hpf","96 hpf","120 hpf","144 hpf",
            "168 hpf","192 hpf")
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/03-DESeq2_expression_data/")
df_DESeq2 <- read.table("Nvec_DESeq2_average_Metazoa_and_TRGs.txt", header = TRUE)
rownames(df_DESeq2) <- df_DESeq2$Gene_ID
colnames(df_DESeq2) <- c("Gene_ID",stages)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/04-Metazoa_and_TRGs/")
df_origin <- read.table("Nvec_Metazoa_and_TRGs_gene_ID.txt", header = F)
colnames(df_origin) <- c("origin","Gene_ID")
df_origin <- df_origin[!is.na(df_origin$Gene_ID),]
rownames(df_origin) <- df_origin$Gene_ID

# Subset genes with non-null expression and with metazoan/species-specific origin
df_notnull <- df_DESeq2[rowSums(df_DESeq2[,-c(1)])>0,]
df_notnull$origin <- NA
df_notnull[df_notnull$Gene_ID,'origin'] <- df_origin[df_notnull$Gene_ID,'origin']
df_notnull$origin[df_notnull$origin == 'Nvec'] <- "Species-specific"
df_withorigin <- df_notnull[!is.na(df_notnull$origin),]

# Quantile normalise expression data
df_withorigin_qn <- cbind(df_withorigin[c(1,ncol(df_withorigin))],normalize.quantiles(data.matrix(df_withorigin[-c(1,ncol(df_withorigin))]))) # quantile normalisation of data
df_withorigin_qn$origin <- factor(df_withorigin_qn$origin, levels = c("Metazoa","Species-specific"))
colnames(df_withorigin_qn)[3:ncol(df_withorigin_qn)] <- stages

# Obtain 75% quantile for data distribution representation purposes in heatmaps and export
quantile_75 <- data.frame(matrix(ncol = length(stages), nrow = 2))
rownames(quantile_75) <- c("Metazoa","Species-specific")
colnames(quantile_75) <- stages

for (i in stages){
  quantile_75[i] <- round(tapply(df_withorigin_qn[,i], df_withorigin_qn$origin, quantile, p = 0.75, na.rm = TRUE),2)
}

quantile_75_clean <- cbind('origin'=rownames(quantile_75),quantile_75)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/08-quantile_normalised_heatmaps/")
write.table(quantile_75_clean, "03-Nvec_quantile_75.txt", quote = F, sep = '\t', row.names = F)

# Transformations for boxplot
df_clean <- data.frame(cbind(stages,t(df_withorigin_qn[-c(1,2)])))
df_clean_long <- gather(df_clean, gene_ID, expression, -stages)
df_clean_long$origin <- rep(df_withorigin_qn$origin, each = length(stages))
df_clean_long$expression <- as.numeric(df_clean_long$expression)
df_clean_long$origin <- factor(df_clean_long$origin, levels = c("Metazoa","Species-specific"))
df_clean_long$stages <- factor(df_clean_long$stages, levels = stages)
df_clean_long_metazoa <- df_clean_long[df_clean_long$origin == 'Metazoa',]
df_clean_long_specific <- df_clean_long[df_clean_long$origin == 'Species-specific',]

# Boxplots
metazoa_boxplot <- ggplot(df_clean_long_metazoa, aes(x=stages, y=expression, fill=origin)) +
  geom_boxplot(outlier.shape = NA, fill = '#FF0000') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "quantile-normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,2000)) +
  scale_y_continuous(breaks = c(seq(0,2000,200))) +
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

specific_boxplot <- ggplot(df_clean_long_specific, aes(x=stages, y=expression, fill=origin)) +
  geom_boxplot(outlier.shape = NA, fill = '#147DF5') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "quantile-normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,240)) +
  scale_y_continuous(breaks = c(seq(0,240,20))) +
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

ggarrange(metazoa_boxplot, specific_boxplot)



#######################################################################
#################  4. Strongylocentrotus purpuratus ###################
#######################################################################


# Input data
stages <- c("hatched_blastula","mesenchyme_blastula","early_gastrula","mid_blastula",
            "late_gastrula","prism","late_prism","pluteus","four-arm_stage","vestibular_stage",
            "pentagonal_stage","tube_foot_stage","post-metamorphic","juvenile")
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/03-DESeq2_expression_data/")
df_DESeq2 <- read.table("Spur_DESeq2_average_Metazoa_and_TRGs.txt", header = TRUE)
rownames(df_DESeq2) <- df_DESeq2$Gene_ID
colnames(df_DESeq2) <- c("Gene_ID",stages)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/04-Metazoa_and_TRGs/")
df_origin <- read.table("Spur_Metazoa_and_TRGs_gene_ID.txt", header = F)
colnames(df_origin) <- c("origin","Gene_ID")
df_origin <- df_origin[!is.na(df_origin$Gene_ID),]
rownames(df_origin) <- df_origin$Gene_ID

# Subset genes with non-null expression and with metazoan/species-specific origin
df_notnull <- df_DESeq2[rowSums(df_DESeq2[,-c(1)])>0,]
df_notnull$origin <- NA
df_notnull[df_notnull$Gene_ID,'origin'] <- df_origin[df_notnull$Gene_ID,'origin']
df_notnull$origin[df_notnull$origin == 'Spur'] <- "Species-specific"
df_withorigin <- df_notnull[!is.na(df_notnull$origin),]

# Quantile normalise expression data
df_withorigin_qn <- cbind(df_withorigin[c(1,ncol(df_withorigin))],normalize.quantiles(data.matrix(df_withorigin[-c(1,ncol(df_withorigin))]))) # quantile normalisation of data
df_withorigin_qn$origin <- factor(df_withorigin_qn$origin, levels = c("Metazoa","Species-specific"))
colnames(df_withorigin_qn)[3:ncol(df_withorigin_qn)] <- stages

# Obtain 75% quantile for data distribution representation purposes in heatmaps and export
quantile_75 <- data.frame(matrix(ncol = length(stages), nrow = 2))
rownames(quantile_75) <- c("Metazoa","Species-specific")
colnames(quantile_75) <- stages

for (i in stages){
  quantile_75[i] <- round(tapply(df_withorigin_qn[,i], df_withorigin_qn$origin, quantile, p = 0.75, na.rm = TRUE),2)
}

quantile_75_clean <- cbind('origin'=rownames(quantile_75),quantile_75)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/08-quantile_normalised_heatmaps/")
write.table(quantile_75_clean, "04-Spur_quantile_75.txt", quote = F, sep = '\t', row.names = F)

# Transformations for boxplot
df_clean <- data.frame(cbind(stages,t(df_withorigin_qn[-c(1,2)])))
df_clean_long <- gather(df_clean, gene_ID, expression, -stages)
df_clean_long$origin <- rep(df_withorigin_qn$origin, each = length(stages))
df_clean_long$expression <- as.numeric(df_clean_long$expression)
df_clean_long$origin <- factor(df_clean_long$origin, levels = c("Metazoa","Species-specific"))
df_clean_long$stages <- factor(df_clean_long$stages, levels = stages)
df_clean_long_metazoa <- df_clean_long[df_clean_long$origin == 'Metazoa',]
df_clean_long_specific <- df_clean_long[df_clean_long$origin == 'Species-specific',]

# Boxplots
metazoa_boxplot <- ggplot(df_clean_long_metazoa, aes(x=stages, y=expression, fill=origin)) +
  geom_boxplot(outlier.shape = NA, fill = '#FF0000') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "quantile-normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,1200)) +
  scale_y_continuous(breaks = c(seq(0,1200,100))) +
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

specific_boxplot <- ggplot(df_clean_long_specific, aes(x=stages, y=expression, fill=origin)) +
  geom_boxplot(outlier.shape = NA, fill = '#147DF5') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "quantile-normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,240)) +
  scale_y_continuous(breaks = c(seq(0,240,20))) +
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

ggarrange(metazoa_boxplot, specific_boxplot)


#######################################################################
##################  5. Branchiostoma lanceolatum ######################
#######################################################################


# Input data
stages <- c("blastula","7_hpf","8_hpf","10_hpf","11_hpf","15_hpf","18_hpf",
            "21_hpf","24_hpf","27_hpf","36_hpf","50_hpf","60_hpf","pre-met.")
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/03-DESeq2_expression_data/")
df_DESeq2 <- read.table("Blan_DESeq2_average_Metazoa_and_TRGs.txt", header = TRUE)
rownames(df_DESeq2) <- df_DESeq2$Gene_ID
colnames(df_DESeq2) <- c("Gene_ID",stages)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/04-Metazoa_and_TRGs/")
df_origin <- read.table("Blan_Metazoa_and_TRGs_gene_ID.txt", header = F)
colnames(df_origin) <- c("origin","Gene_ID")
df_origin <- df_origin[!is.na(df_origin$Gene_ID),]
rownames(df_origin) <- df_origin$Gene_ID

# Subset genes with non-null expression and with metazoan/species-specific origin
df_notnull <- df_DESeq2[rowSums(df_DESeq2[,-c(1)])>0,]
df_notnull$origin <- NA
df_notnull[df_notnull$Gene_ID,'origin'] <- df_origin[df_notnull$Gene_ID,'origin']
df_notnull$origin[df_notnull$origin == 'Blan'] <- "Species-specific"
df_withorigin <- df_notnull[!is.na(df_notnull$origin),]

# Quantile normalise expression data
df_withorigin_qn <- cbind(df_withorigin[c(1,ncol(df_withorigin))],normalize.quantiles(data.matrix(df_withorigin[-c(1,ncol(df_withorigin))]))) # quantile normalisation of data
df_withorigin_qn$origin <- factor(df_withorigin_qn$origin, levels = c("Metazoa","Species-specific"))
colnames(df_withorigin_qn)[3:ncol(df_withorigin_qn)] <- stages

# Obtain 75% quantile for data distribution representation purposes in heatmaps and export
quantile_75 <- data.frame(matrix(ncol = length(stages), nrow = 2))
rownames(quantile_75) <- c("Metazoa","Species-specific")
colnames(quantile_75) <- stages

for (i in stages){
  quantile_75[i] <- round(tapply(df_withorigin_qn[,i], df_withorigin_qn$origin, quantile, p = 0.75, na.rm = TRUE),2)
}

quantile_75_clean <- cbind('origin'=rownames(quantile_75),quantile_75)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/08-quantile_normalised_heatmaps/")
write.table(quantile_75_clean, "05-Blan_quantile_75.txt", quote = F, sep = '\t', row.names = F)

# Transformations for boxplot
df_clean <- data.frame(cbind(stages,t(df_withorigin_qn[-c(1,2)])))
df_clean_long <- gather(df_clean, gene_ID, expression, -stages)
df_clean_long$origin <- rep(df_withorigin_qn$origin, each = length(stages))
df_clean_long$expression <- as.numeric(df_clean_long$expression)
df_clean_long$origin <- factor(df_clean_long$origin, levels = c("Metazoa","Species-specific"))
df_clean_long$stages <- factor(df_clean_long$stages, levels = stages)
df_clean_long_metazoa <- df_clean_long[df_clean_long$origin == 'Metazoa',]
df_clean_long_specific <- df_clean_long[df_clean_long$origin == 'Species-specific',]

# Boxplots
metazoa_boxplot <- ggplot(df_clean_long_metazoa, aes(x=stages, y=expression, fill=origin)) +
  geom_boxplot(outlier.shape = NA, fill = '#FF0000') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "quantile-normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,1500)) +
  scale_y_continuous(breaks = c(seq(0,1500,150))) +
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

specific_boxplot <- ggplot(df_clean_long_specific, aes(x=stages, y=expression, fill=origin)) +
  geom_boxplot(outlier.shape = NA, fill = '#147DF5') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "quantile-normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,100)) +
  scale_y_continuous(breaks = c(seq(0,100,10))) +
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

ggarrange(metazoa_boxplot, specific_boxplot)


#######################################################################
#########################  6. Danio rerio #############################
#######################################################################


# Input data
stages <- c("2_hpf","1000_cells","dome","shield","6_hpf","8_hpf","bud",
            "12_hpf","16_hpf","20_hpf","26_hpf","28_hpf","30_hpf","36_hpf",
            "2_dpf","3_dpf","5_dpf","7_dpf")
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/03-DESeq2_expression_data/")
df_DESeq2 <- read.table("Drer_DESeq2_average_Metazoa_and_TRGs.txt", header = TRUE)
rownames(df_DESeq2) <- df_DESeq2$Gene_ID
colnames(df_DESeq2) <- c("Gene_ID",stages)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/04-Metazoa_and_TRGs/")
df_origin <- read.table("Drer_Metazoa_and_TRGs_gene_ID.txt", header = F)
colnames(df_origin) <- c("origin","Gene_ID")
df_origin <- df_origin[!is.na(df_origin$Gene_ID),]
rownames(df_origin) <- df_origin$Gene_ID

# Subset genes with non-null expression and with metazoan/species-specific origin
df_notnull <- df_DESeq2[rowSums(df_DESeq2[,-c(1)])>0,]
df_notnull$origin <- NA
df_notnull[df_notnull$Gene_ID,'origin'] <- df_origin[df_notnull$Gene_ID,'origin']
df_notnull$origin[df_notnull$origin == 'Drer'] <- "Species-specific"
df_withorigin <- df_notnull[!is.na(df_notnull$origin),]

# Quantile normalise expression data
df_withorigin_qn <- cbind(df_withorigin[c(1,ncol(df_withorigin))],normalize.quantiles(data.matrix(df_withorigin[-c(1,ncol(df_withorigin))]))) # quantile normalisation of data
df_withorigin_qn$origin <- factor(df_withorigin_qn$origin, levels = c("Metazoa","Species-specific"))
colnames(df_withorigin_qn)[3:ncol(df_withorigin_qn)] <- stages

# Obtain 75% quantile for data distribution representation purposes in heatmaps and export
quantile_75 <- data.frame(matrix(ncol = length(stages), nrow = 2))
rownames(quantile_75) <- c("Metazoa","Species-specific")
colnames(quantile_75) <- stages

for (i in stages){
  quantile_75[i] <- round(tapply(df_withorigin_qn[,i], df_withorigin_qn$origin, quantile, p = 0.75, na.rm = TRUE),2)
}

quantile_75_clean <- cbind('origin'=rownames(quantile_75),quantile_75)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/08-quantile_normalised_heatmaps/")
write.table(quantile_75_clean, "06-Drer_quantile_75.txt", quote = F, sep = '\t', row.names = F)

# Transformations for boxplot
df_clean <- data.frame(cbind(stages,t(df_withorigin_qn[-c(1,2)])))
df_clean_long <- gather(df_clean, gene_ID, expression, -stages)
df_clean_long$origin <- rep(df_withorigin_qn$origin, each = length(stages))
df_clean_long$expression <- as.numeric(df_clean_long$expression)
df_clean_long$origin <- factor(df_clean_long$origin, levels = c("Metazoa","Species-specific"))
df_clean_long$stages <- factor(df_clean_long$stages, levels = stages)
df_clean_long_metazoa <- df_clean_long[df_clean_long$origin == 'Metazoa',]
df_clean_long_specific <- df_clean_long[df_clean_long$origin == 'Species-specific',]

# Boxplots
metazoa_boxplot <- ggplot(df_clean_long_metazoa, aes(x=stages, y=expression, fill=origin)) +
  geom_boxplot(outlier.shape = NA, fill = '#FF0000') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "quantile-normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,3000)) +
  scale_y_continuous(breaks = c(seq(0,3000,300))) +
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

specific_boxplot <- ggplot(df_clean_long_specific, aes(x=stages, y=expression, fill=origin)) +
  geom_boxplot(outlier.shape = NA, fill = '#147DF5') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "quantile-normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,300)) +
  scale_y_continuous(breaks = c(seq(0,300,25))) +
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

ggarrange(metazoa_boxplot, specific_boxplot)




#######################################################################
###################  7. Drosophila melanogaster #######################
#######################################################################


# Input data
stages <- c("0–2_hpf","2–4_hpf","4–8_hpf","8–12_hpf","12–16_hpf","16–20_hpf",
            "20–24_hpf","L1_larva","L2_larva","L3_larva_12h","L3_larva_puffstage_1–9",
            "white_prepupa","late_prepupa/early_pupa","pupa_2–3d","pupa_4d","adult_1d")
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/03-DESeq2_expression_data/")
df_DESeq2 <- read.table("Dmel_DESeq2_average_Metazoa_and_TRGs.txt", header = TRUE)
rownames(df_DESeq2) <- df_DESeq2$Gene_ID
colnames(df_DESeq2) <- c("Gene_ID",stages)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/04-Metazoa_and_TRGs/")
df_origin <- read.table("Dmel_Metazoa_and_TRGs_gene_ID.txt", header = F)
colnames(df_origin) <- c("origin","Gene_ID")
df_origin <- df_origin[!is.na(df_origin$Gene_ID),]
rownames(df_origin) <- df_origin$Gene_ID

# Subset genes with non-null expression and with metazoan/species-specific origin
df_notnull <- df_DESeq2[rowSums(df_DESeq2[,-c(1)])>0,]
df_notnull$origin <- NA
df_notnull[df_notnull$Gene_ID,'origin'] <- df_origin[df_notnull$Gene_ID,'origin']
df_notnull$origin[df_notnull$origin == 'Dmel'] <- "Species-specific"
df_withorigin <- df_notnull[!is.na(df_notnull$origin),]

# Quantile normalise expression data
df_withorigin_qn <- cbind(df_withorigin[c(1,ncol(df_withorigin))],normalize.quantiles(data.matrix(df_withorigin[-c(1,ncol(df_withorigin))]))) # quantile normalisation of data
df_withorigin_qn$origin <- factor(df_withorigin_qn$origin, levels = c("Metazoa","Species-specific"))
colnames(df_withorigin_qn)[3:ncol(df_withorigin_qn)] <- stages

# Obtain 75% quantile for data distribution representation purposes in heatmaps and export
quantile_75 <- data.frame(matrix(ncol = length(stages), nrow = 2))
rownames(quantile_75) <- c("Metazoa","Species-specific")
colnames(quantile_75) <- stages

for (i in stages){
  quantile_75[i] <- round(tapply(df_withorigin_qn[,i], df_withorigin_qn$origin, quantile, p = 0.75, na.rm = TRUE),2)
}

quantile_75_clean <- cbind('origin'=rownames(quantile_75),quantile_75)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/08-quantile_normalised_heatmaps/")
write.table(quantile_75_clean, "07-Dmel_quantile_75.txt", quote = F, sep = '\t', row.names = F)

# Transformations for boxplot
df_clean <- data.frame(cbind(stages,t(df_withorigin_qn[-c(1,2)])))
df_clean_long <- gather(df_clean, gene_ID, expression, -stages)
df_clean_long$origin <- rep(df_withorigin_qn$origin, each = length(stages))
df_clean_long$expression <- as.numeric(df_clean_long$expression)
df_clean_long$origin <- factor(df_clean_long$origin, levels = c("Metazoa","Species-specific"))
df_clean_long$stages <- factor(df_clean_long$stages, levels = stages)
df_clean_long_metazoa <- df_clean_long[df_clean_long$origin == 'Metazoa',]
df_clean_long_specific <- df_clean_long[df_clean_long$origin == 'Species-specific',]

# Boxplots
metazoa_boxplot <- ggplot(df_clean_long_metazoa, aes(x=stages, y=expression, fill=origin)) +
  geom_boxplot(outlier.shape = NA, fill = '#FF0000') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "quantile-normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,1300)) +
  scale_y_continuous(breaks = c(seq(0,1300,100))) +
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

specific_boxplot <- ggplot(df_clean_long_specific, aes(x=stages, y=expression, fill=origin)) +
  geom_boxplot(outlier.shape = NA, fill = '#147DF5') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "quantile-normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,250)) +
  scale_y_continuous(breaks = c(seq(0,250,25))) +
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

ggarrange(metazoa_boxplot, specific_boxplot)




#######################################################################
###################  8. Caenorhabditis elegans ########################
#######################################################################


# Input data
stages <- c("0–90_mpf","120–180_mpf","210–240_mpf","270–330_mpf","360–390_mpf",
            "420–480_mpf","510–540_mpf","570–600_mpf","L1_larva","L3_larva","L4_larva")
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/03-DESeq2_expression_data/")
df_DESeq2 <- read.table("Cele_DESeq2_average_Metazoa_and_TRGs.txt", header = TRUE)
rownames(df_DESeq2) <- df_DESeq2$Gene_ID
colnames(df_DESeq2) <- c("Gene_ID",stages)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/04-Metazoa_and_TRGs/")
df_origin <- read.table("Cele_Metazoa_and_TRGs_gene_ID.txt", header = F)
colnames(df_origin) <- c("origin","Gene_ID")
df_origin <- df_origin[!is.na(df_origin$Gene_ID),]
rownames(df_origin) <- df_origin$Gene_ID

# Subset genes with non-null expression and with metazoan/species-specific origin
df_notnull <- df_DESeq2[rowSums(df_DESeq2[,-c(1)])>0,]
df_notnull$origin <- NA
df_notnull[df_notnull$Gene_ID,'origin'] <- df_origin[df_notnull$Gene_ID,'origin']
df_notnull$origin[df_notnull$origin == 'Cele'] <- "Species-specific"
df_withorigin <- df_notnull[!is.na(df_notnull$origin),-c(ncol(df_notnull)-1)]

# Quantile normalise expression data
df_withorigin_qn <- cbind(df_withorigin[c(1,ncol(df_withorigin))],normalize.quantiles(data.matrix(df_withorigin[-c(1,ncol(df_withorigin))]))) # quantile normalisation of data
df_withorigin_qn$origin <- factor(df_withorigin_qn$origin, levels = c("Metazoa","Species-specific"))
colnames(df_withorigin_qn)[3:ncol(df_withorigin_qn)] <- stages

# Obtain 75% quantile for data distribution representation purposes in heatmaps and export
quantile_75 <- data.frame(matrix(ncol = length(stages), nrow = 2))
rownames(quantile_75) <- c("Metazoa","Species-specific")
colnames(quantile_75) <- stages

for (i in stages){
  quantile_75[i] <- round(tapply(df_withorigin_qn[,i], df_withorigin_qn$origin, quantile, p = 0.75, na.rm = TRUE),2)
}

quantile_75_clean <- cbind('origin'=rownames(quantile_75),quantile_75)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/08-quantile_normalised_heatmaps/")
write.table(quantile_75_clean, "08-Cele_quantile_75.txt", quote = F, sep = '\t', row.names = F)

# Transformations for boxplot
df_clean <- data.frame(cbind(stages,t(df_withorigin_qn[-c(1,2)])))
df_clean_long <- gather(df_clean, gene_ID, expression, -stages)
df_clean_long$origin <- rep(df_withorigin_qn$origin, each = length(stages))
df_clean_long$expression <- as.numeric(df_clean_long$expression)
df_clean_long$origin <- factor(df_clean_long$origin, levels = c("Metazoa","Species-specific"))
df_clean_long$stages <- factor(df_clean_long$stages, levels = stages)
df_clean_long_metazoa <- df_clean_long[df_clean_long$origin == 'Metazoa',]
df_clean_long_specific <- df_clean_long[df_clean_long$origin == 'Species-specific',]

# Boxplots
metazoa_boxplot <- ggplot(df_clean_long_metazoa, aes(x=stages, y=expression, fill=origin)) +
  geom_boxplot(outlier.shape = NA, fill = '#FF0000') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "quantile-normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,350)) +
  scale_y_continuous(breaks = c(seq(0,350,50))) +
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

specific_boxplot <- ggplot(df_clean_long_specific, aes(x=stages, y=expression, fill=origin)) +
  geom_boxplot(outlier.shape = NA, fill = '#147DF5') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "quantile-normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,40)) +
  scale_y_continuous(breaks = c(seq(0,40,5))) +
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

ggarrange(metazoa_boxplot, specific_boxplot)




#######################################################################
#####################  9. Crassostrea gigas ###########################
#######################################################################


# Input data
stages <- c("early_morula","morula","blastula","rotary_movement","free-swimming",
            "early_gastrula","gastrula","trochophore_1","trochophore_2","trochophore_3",
            "trochophore_4","trochophore_5","early_D-stage_1","early_D-stage_2",
            "D-stage_1","D-stage_2","D-stage_3","D-stage_4","D-stage_5","D-stage_6",
            "D-stage_7","early_umbo_1","early_umbo_2","umbo_1","umbo_2","umbo_3",
            "umbo_4","umbo_5","umbo_6","late_umbo_1","late_umbo_2","pediveliger_1",
            "pediveliger_2","spat","juvenile")
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/03-DESeq2_expression_data/")
df_DESeq2 <- read.table("Cgig_DESeq2_average_Metazoa_and_TRGs.txt", header = TRUE)
rownames(df_DESeq2) <- df_DESeq2$Gene_ID
colnames(df_DESeq2) <- c("Gene_ID",stages)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/04-Metazoa_and_TRGs/")
df_origin <- read.table("Cgig_Metazoa_and_TRGs_gene_ID.txt", header = F)
colnames(df_origin) <- c("origin","Gene_ID")
df_origin <- df_origin[!is.na(df_origin$Gene_ID),]
rownames(df_origin) <- df_origin$Gene_ID

# Subset genes with non-null expression and with metazoan/species-specific origin
df_notnull <- df_DESeq2[rowSums(df_DESeq2[,-c(1)])>0,]
df_notnull$origin <- NA
df_notnull[df_notnull$Gene_ID,'origin'] <- df_origin[df_notnull$Gene_ID,'origin']
df_notnull$origin[df_notnull$origin == 'Cgig'] <- "Species-specific"
df_withorigin <- df_notnull[!is.na(df_notnull$origin),]

# Quantile normalise expression data
df_withorigin_qn <- cbind(df_withorigin[c(1,ncol(df_withorigin))],normalize.quantiles(data.matrix(df_withorigin[-c(1,ncol(df_withorigin))]))) # quantile normalisation of data
df_withorigin_qn$origin <- factor(df_withorigin_qn$origin, levels = c("Metazoa","Species-specific"))
colnames(df_withorigin_qn)[3:ncol(df_withorigin_qn)] <- stages

# Obtain 75% quantile for data distribution representation purposes in heatmaps and export
quantile_75 <- data.frame(matrix(ncol = length(stages), nrow = 2))
rownames(quantile_75) <- c("Metazoa","Species-specific")
colnames(quantile_75) <- stages

for (i in stages){
  quantile_75[i] <- round(tapply(df_withorigin_qn[,i], df_withorigin_qn$origin, quantile, p = 0.75, na.rm = TRUE),2)
}

quantile_75_clean <- cbind('origin'=rownames(quantile_75),quantile_75)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/08-quantile_normalised_heatmaps/")
write.table(quantile_75_clean, "09-Cgig_quantile_75.txt", quote = F, sep = '\t', row.names = F)

# Transformations for boxplot
df_clean <- data.frame(cbind(stages,t(df_withorigin_qn[-c(1,2)])))
df_clean_long <- gather(df_clean, gene_ID, expression, -stages)
df_clean_long$origin <- rep(df_withorigin_qn$origin, each = length(stages))
df_clean_long$expression <- as.numeric(df_clean_long$expression)
df_clean_long$origin <- factor(df_clean_long$origin, levels = c("Metazoa","Species-specific"))
df_clean_long$stages <- factor(df_clean_long$stages, levels = stages)
df_clean_long_metazoa <- df_clean_long[df_clean_long$origin == 'Metazoa',]
df_clean_long_specific <- df_clean_long[df_clean_long$origin == 'Species-specific',]

# Boxplots
metazoa_boxplot <- ggplot(df_clean_long_metazoa, aes(x=stages, y=expression, fill=origin)) +
  geom_boxplot(outlier.shape = NA, fill = '#FF0000') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "quantile-normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,900)) +
  scale_y_continuous(breaks = c(seq(0,900,100))) +
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

specific_boxplot <- ggplot(df_clean_long_specific, aes(x=stages, y=expression, fill=origin)) +
  geom_boxplot(outlier.shape = NA, fill = '#147DF5') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "quantile-normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,50)) +
  scale_y_continuous(breaks = c(seq(0,50,5))) +
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

ggarrange(metazoa_boxplot, specific_boxplot)



#######################################################################
#####################  10. Capitella teleta ###########################
#######################################################################


# Input data
stages <- c("64_cells","gastrula","stage_4tt_larva","stage_5_larva",
            "stage_7_larva","pre-competent_larva","competent_larva")
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/03-DESeq2_expression_data/")
df_DESeq2 <- read.table("Ctel_DESeq2_average_Metazoa_and_TRGs.txt", header = TRUE)
df_DESeq2 <- df_DESeq2[-c(2:8)]
rownames(df_DESeq2) <- df_DESeq2$Gene_ID
colnames(df_DESeq2) <- c("Gene_ID",stages)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/04-Metazoa_and_TRGs/")
df_origin <- read.table("Ctel_Metazoa_and_TRGs_gene_ID.txt", header = F)
colnames(df_origin) <- c("origin","Gene_ID")
df_origin <- df_origin[!is.na(df_origin$Gene_ID),]
rownames(df_origin) <- df_origin$Gene_ID

# Subset genes with non-null expression and with metazoan/species-specific origin
df_notnull <- df_DESeq2[rowSums(df_DESeq2[,-c(1)])>0,]
df_notnull$origin <- NA
df_notnull[df_notnull$Gene_ID,'origin'] <- df_origin[df_notnull$Gene_ID,'origin']
df_notnull$origin[df_notnull$origin == 'Ctel'] <- "Species-specific"
df_withorigin <- df_notnull[!is.na(df_notnull$origin),]

# Quantile normalise expression data
df_withorigin_qn <- cbind(df_withorigin[c(1,ncol(df_withorigin))],normalize.quantiles(data.matrix(df_withorigin[-c(1,ncol(df_withorigin))]))) # quantile normalisation of data
df_withorigin_qn$origin <- factor(df_withorigin_qn$origin, levels = c("Metazoa","Species-specific"))
colnames(df_withorigin_qn)[3:ncol(df_withorigin_qn)] <- stages

# Obtain 75% quantile for data distribution representation purposes in heatmaps and export
quantile_75 <- data.frame(matrix(ncol = length(stages), nrow = 2))
rownames(quantile_75) <- c("Metazoa","Species-specific")
colnames(quantile_75) <- stages

for (i in stages){
  quantile_75[i] <- round(tapply(df_withorigin_qn[,i], df_withorigin_qn$origin, quantile, p = 0.75, na.rm = TRUE),2)
}

quantile_75_clean <- cbind('origin'=rownames(quantile_75),quantile_75)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/08-quantile_normalised_heatmaps/")
write.table(quantile_75_clean, "10-Ctel_quantile_75.txt", quote = F, sep = '\t', row.names = F)

# Transformations for boxplot
df_clean <- data.frame(cbind(stages,t(df_withorigin_qn[-c(1,2)])))
df_clean_long <- gather(df_clean, gene_ID, expression, -stages)
df_clean_long$origin <- rep(df_withorigin_qn$origin, each = length(stages))
df_clean_long$expression <- as.numeric(df_clean_long$expression)
df_clean_long$origin <- factor(df_clean_long$origin, levels = c("Metazoa","Species-specific"))
df_clean_long$stages <- factor(df_clean_long$stages, levels = stages)
df_clean_long_metazoa <- df_clean_long[df_clean_long$origin == 'Metazoa',]
df_clean_long_specific <- df_clean_long[df_clean_long$origin == 'Species-specific',]

# Boxplots
metazoa_boxplot <- ggplot(df_clean_long_metazoa, aes(x=stages, y=expression, fill=origin)) +
  geom_boxplot(outlier.shape = NA, fill = '#FF0000') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "quantile-normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,3600)) +
  scale_y_continuous(breaks = c(seq(0,3600,500))) +
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

specific_boxplot <- ggplot(df_clean_long_specific, aes(x=stages, y=expression, fill=origin)) +
  geom_boxplot(outlier.shape = NA, fill = '#147DF5') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "quantile-normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,100)) +
  scale_y_continuous(breaks = c(seq(0,100,10))) +
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

ggarrange(metazoa_boxplot, specific_boxplot)




#######################################################################
####################  11. Owenia fusiformis ###########################
#######################################################################


# Input data
stages <- c("blastula","gastrula","elongation","early_larva","mitraria_larva",
            "competent_larva","juvenile")
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/03-DESeq2_expression_data/")
df_DESeq2 <- read.table("Ofus_DESeq2_average_Metazoa_and_TRGs.txt", header = TRUE)
df_DESeq2 <- df_DESeq2[-c(2:8)]
rownames(df_DESeq2) <- df_DESeq2$Gene_ID
colnames(df_DESeq2) <- c("Gene_ID",stages)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/04-Metazoa_and_TRGs/")
df_origin <- read.table("Ofus_Metazoa_and_TRGs_gene_ID.txt", header = F)
colnames(df_origin) <- c("origin","Gene_ID")
df_origin <- df_origin[!is.na(df_origin$Gene_ID),]
rownames(df_origin) <- df_origin$Gene_ID

# Subset genes with non-null expression and with metazoan/species-specific origin
df_notnull <- df_DESeq2[rowSums(df_DESeq2[,-c(1)])>0,]
df_notnull$origin <- NA
df_notnull[df_notnull$Gene_ID,'origin'] <- df_origin[df_notnull$Gene_ID,'origin']
df_notnull$origin[df_notnull$origin == 'Ofus'] <- "Species-specific"
df_withorigin <- df_notnull[!is.na(df_notnull$origin),]

# Quantile normalise expression data
df_withorigin_qn <- cbind(df_withorigin[c(1,ncol(df_withorigin))],normalize.quantiles(data.matrix(df_withorigin[-c(1,ncol(df_withorigin))]))) # quantile normalisation of data
df_withorigin_qn$origin <- factor(df_withorigin_qn$origin, levels = c("Metazoa","Species-specific"))
colnames(df_withorigin_qn)[3:ncol(df_withorigin_qn)] <- stages

# Obtain 75% quantile for data distribution representation purposes in heatmaps and export
quantile_75 <- data.frame(matrix(ncol = length(stages), nrow = 2))
rownames(quantile_75) <- c("Metazoa","Species-specific")
colnames(quantile_75) <- stages

for (i in stages){
  quantile_75[i] <- round(tapply(df_withorigin_qn[,i], df_withorigin_qn$origin, quantile, p = 0.75, na.rm = TRUE),2)
}

quantile_75_clean <- cbind('origin'=rownames(quantile_75),quantile_75)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/10-Phylostratigraphy_JSD_species_Metazoan_TRGs/08-quantile_normalised_heatmaps/")
write.table(quantile_75_clean, "11-Ofus_quantile_75.txt", quote = F, sep = '\t', row.names = F)

# Transformations for boxplot
df_clean <- data.frame(cbind(stages,t(df_withorigin_qn[-c(1,2)])))
df_clean_long <- gather(df_clean, gene_ID, expression, -stages)
df_clean_long$origin <- rep(df_withorigin_qn$origin, each = length(stages))
df_clean_long$expression <- as.numeric(df_clean_long$expression)
df_clean_long$origin <- factor(df_clean_long$origin, levels = c("Metazoa","Species-specific"))
df_clean_long$stages <- factor(df_clean_long$stages, levels = stages)
df_clean_long_metazoa <- df_clean_long[df_clean_long$origin == 'Metazoa',]
df_clean_long_specific <- df_clean_long[df_clean_long$origin == 'Species-specific',]

# Boxplots
metazoa_boxplot <- ggplot(df_clean_long_metazoa, aes(x=stages, y=expression, fill=origin)) +
  geom_boxplot(outlier.shape = NA, fill = '#FF0000') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "quantile-normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,4500)) +
  scale_y_continuous(breaks = c(seq(0,4500,500))) +
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

specific_boxplot <- ggplot(df_clean_long_specific, aes(x=stages, y=expression, fill=origin)) +
  geom_boxplot(outlier.shape = NA, fill = '#147DF5') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "quantile-normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,500)) +
  scale_y_continuous(breaks = c(seq(0,500,50))) +
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

ggarrange(metazoa_boxplot, specific_boxplot)



