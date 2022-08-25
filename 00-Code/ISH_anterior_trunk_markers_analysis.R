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

#######################################################################
########################  1. Heatmap 0 to max #########################
#######################################################################

rescale_custom <- function(x) (x/(max(x)))

# Capitella teleta
raw <- read.table("06-Fran_Data_Capitella.txt", header = TRUE)
rownames(raw) <- raw$Transcript_ID

df_raw <- raw[,-c(1:10)]
df_norm <- t(apply(df_raw, 1, rescale_custom))
df_norm <- na.omit(df_norm)

heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color <- heatmap_color[50:100]

ComplexHeatmap::Heatmap(df_norm,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        col = heatmap_color,
                        heatmap_legend_param = list(color_bar = "continuous"),
                        row_labels = raw$Gene_symbol)

# Owenia fusiformis
raw <- read.table("05-Fran_Data_Owenia.txt", header = TRUE)
rownames(raw) <- raw$Transcript_ID

df_raw <- raw[,-c(1:10)]
df_norm <- t(apply(df_raw, 1, rescale_custom))
df_norm <- na.omit(df_norm)

heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color <- heatmap_color[50:100]

ComplexHeatmap::Heatmap(df_norm,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        col = heatmap_color,
                        heatmap_legend_param = list(color_bar = "continuous"),
                        row_labels = raw$Gene_symbol)



#######################################################################
#######################  2. DESeq2 expression #########################
#######################################################################

palette_raw <- viridis_pal(option = "C", direction = 1)(5)
palette <- palette_raw[c(1,3,5)]
time <- c(1:7)

# Capitella teleta
df_DESeq2 <- read.table("06-Fran_Data_Capitella.txt", header = TRUE)
rownames(df_DESeq2) <- df_DESeq2$Transcript_ID

df_norm <- df_DESeq2[,-c(1:10)]
df_clean <- data.frame(cbind(time,t(df_norm)))
df_clean_long <- gather(df_clean, gene_ID, expression, -time)
df_clean_long$family <- rep(df_DESeq2$Family, each = 7)

df_clean_long %>%
  ggplot(aes(x=time, y=expression, color=family)) +
  geom_line(data = df_clean_long, aes(x=time, y=expression, group = gene_ID), alpha = 0.5) +
  geom_smooth(method ='loess', size = 2, alpha = 0.2) +  
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "DESeq2 normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,250)) +
  scale_y_continuous(breaks = c(seq(0,250,50))) +
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


# Owenia fusiformis
df_DESeq2 <- read.table("05-Fran_Data_Owenia.txt", header = TRUE)
rownames(df_DESeq2) <- df_DESeq2$Transcript_ID

df_norm <- df_DESeq2[,-c(1:10)]
df_clean <- data.frame(cbind(time,t(df_norm)))
df_clean_long <- gather(df_clean, gene_ID, expression, -time)
df_clean_long$family <- rep(df_DESeq2$Family, each = 7)

df_clean_long %>%
  ggplot(aes(x=time, y=expression, color=family)) +
  geom_line(data = df_clean_long, aes(x=time, y=expression, group = gene_ID), alpha = 0.5) +
  geom_smooth(method ='loess', size = 2, alpha = 0.2) +  
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "DESeq2 normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,5000)) +
  scale_y_continuous(breaks = c(seq(0,5000,500))) +
  scale_x_continuous(breaks = c(1:7), expand = c(0,0), limits = c(0.8,7.2),
                     labels = c("blastula", "gastrula","elongation",
                                "early larva", "mitraria larva",
                                "competent larva","juvenile")) +
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


#######################################################################
######################  3. Normalised dynamics ########################
#######################################################################

palette_raw <- viridis_pal(option = "C", direction = 1)(5)
palette <- palette_raw[c(1,3,5)]
time <- c(1:7)
time2 <- c(1:4)

# Dimorphilus gyrociliatus
df_DESeq2 <- read.table("07-Fran_Data_Dimorphilus.txt", header = TRUE)
rownames(df_DESeq2) <- df_DESeq2$Transcript_ID

df_copy <- df_DESeq2[,-c(1:3)]
df_norm <- data.frame(t(scale(t(data.matrix(df_copy)))))
df_clean <- data.frame(cbind(time2,t(df_norm)))
df_clean_long <- gather(df_clean, gene_ID, expression, -time2)
df_clean_long$family <- rep(df_DESeq2$Family, each = 4)

df_clean_long %>%
  ggplot(aes(x=time2, y=expression, color=family)) +
  geom_smooth(aes(fill = family), method ='loess', size = 2, alpha = 0.1) +  
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "DESeq2 normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(-3,3)) +
  scale_y_continuous(breaks = c(seq(-3,3,1))) +
  scale_x_continuous(breaks = c(1:4), expand = c(0,0), limits = c(0.8,4.2),
                     labels = c("early_development","late_development",
                                "hatchling","female_Adult")) +
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

# Capitella teleta
df_DESeq2 <- read.table("06-Fran_Data_Capitella.txt", header = TRUE)
rownames(df_DESeq2) <- df_DESeq2$Transcript_ID

df_copy <- df_DESeq2[,-c(1:10)]
df_norm <- data.frame(t(scale(t(data.matrix(df_copy)))))
df_clean <- data.frame(cbind(time,t(df_norm)))
df_clean_long <- gather(df_clean, gene_ID, expression, -time)
df_clean_long$family <- rep(df_DESeq2$Family, each = 7)

df_clean_long %>%
  ggplot(aes(x=time, y=expression, color=family)) +
  geom_smooth(aes(fill = family), method ='loess', size = 2, alpha = 0.1) +  
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "DESeq2 normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(-3,3)) +
  scale_y_continuous(breaks = c(seq(-3,3,1))) +
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


# Owenia fusiformis
df_DESeq2 <- read.table("05-Fran_Data_Owenia.txt", header = TRUE)
rownames(df_DESeq2) <- df_DESeq2$Transcript_ID

df_copy <- df_DESeq2[,-c(1:10)]
df_norm <- data.frame(t(scale(t(data.matrix(df_copy)))))
df_clean <- data.frame(cbind(time,t(df_norm)))
df_clean_long <- gather(df_clean, gene_ID, expression, -time)
df_clean_long$family <- rep(df_DESeq2$Family, each = 7)

df_clean_long %>%
  ggplot(aes(x=time, y=expression, color=family)) +
  geom_smooth(aes(fill = family), method ='loess', size = 2, alpha = 0.1) +  
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "DESeq2 normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(-3,3)) +
  scale_y_continuous(breaks = c(seq(-3,3,1))) +
  scale_x_continuous(breaks = c(1:7), expand = c(0,0), limits = c(0.8,7.2),
                     labels = c("blastula", "gastrula","elongation",
                                "early larva", "mitraria larva",
                                "competent larva","juvenile")) +
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



#######################################################################
###########  4. Boxplots/violin plots of DESeq2 expression ############
#######################################################################


# Capitella teleta
stages <- c("blastula", "gastrula","st4tt", "st5", 
            "st7", "pre-competent", "competent")
df_DESeq2 <- read.table("06-Fran_Data_Capitella.txt", header = TRUE)
rownames(df_DESeq2) <- df_DESeq2$Transcript_ID

df_copy <- df_DESeq2[,-c(1:10)]
df_clean <- data.frame(cbind(stages,t(df_copy)))
df_clean_long <- gather(df_clean, gene_ID, expression, -stages)
df_clean_long$family <- rep(df_DESeq2$Family, each = 7)
df_clean_long$expression <- as.numeric(df_clean_long$expression)
df_clean_long$family <- factor(df_clean_long$family, levels = c("Anterior", "Posterior", "Trunk"))
df_clean_long$stages <- factor(df_clean_long$stages, levels = c("blastula", "gastrula","st4tt", "st5", 
                                                                "st7", "pre-competent", "competent"))


ggplot(df_clean_long, aes(x=stages, y=expression, fill=family)) +
  geom_boxplot() +  
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "DESeq2 normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,250)) +
  scale_y_continuous(breaks = c(seq(0,250,50))) +
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


# Owenia fusiformis
stages <-  c("blastula", "gastrula","elongation",
             "early larva", "mitraria larva",
             "competent larva","juvenile")
df_DESeq2 <- read.table("05-Fran_Data_Owenia.txt", header = TRUE)
rownames(df_DESeq2) <- df_DESeq2$Transcript_ID

df_copy <- df_DESeq2[,-c(1:10)]
df_clean <- data.frame(cbind(stages,t(df_copy)))
df_clean_long <- gather(df_clean, gene_ID, expression, -stages)
df_clean_long$family <- rep(df_DESeq2$Family, each = 7)
df_clean_long$expression <- as.numeric(df_clean_long$expression)
df_clean_long$family <- factor(df_clean_long$family, levels = c("Anterior", "Posterior", "Trunk"))
df_clean_long$stages <- factor(df_clean_long$stages, levels = c("blastula", "gastrula","elongation",
                                                                "early larva", "mitraria larva",
                                                                "competent larva","juvenile"))


ggplot(df_clean_long, aes(x=stages, y=expression, fill=family)) +
  geom_boxplot() +  
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "DESeq2 normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,5000)) +
  scale_y_continuous(breaks = c(seq(0,5000,500))) +
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



#######################################################################
##########  5. Boxplots/violin plots of normalised dynamics ###########
#######################################################################

# Dimorphilus gyrociliatus
stages <- c("early_development","late_development","hatchling",
            "female_adult")
df_DESeq2 <- read.table("07-Fran_Data_Dimorphilus.txt", header = TRUE)
rownames(df_DESeq2) <- df_DESeq2$Transcript_ID

df_copy <- df_DESeq2[,-c(1:3)]
df_norm <- data.frame(t(scale(t(data.matrix(df_copy)))))
df_clean <- data.frame(cbind(stages,t(df_norm)))
df_clean_long <- gather(df_clean, gene_ID, expression, -stages)
df_clean_long$family <- rep(df_DESeq2$Family, each = 4)
df_clean_long$expression <- as.numeric(df_clean_long$expression)
df_clean_long$family <- factor(df_clean_long$family, levels = c("Anterior", "Trunk"))
df_clean_long$stages <- factor(df_clean_long$stages, levels = stages)

ggplot(df_clean_long, aes(x=stages, y=expression, fill=family)) +
  geom_boxplot(outlier.shape = NA) +  
  facet_wrap(~family) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "DESeq2 normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(-3,3)) +
  scale_y_continuous(breaks = c(seq(-3,3,1))) +
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

# Capitella teleta
stages <- c("blastula", "gastrula","st4tt", "st5", 
            "st7", "pre-competent", "competent")
df_DESeq2 <- read.table("06-Fran_Data_Capitella.txt", header = TRUE)
rownames(df_DESeq2) <- df_DESeq2$Transcript_ID

df_copy <- df_DESeq2[,-c(1:10)]
df_norm <- data.frame(t(scale(t(data.matrix(df_copy)))))
df_clean <- data.frame(cbind(stages,t(df_norm)))
df_clean_long <- gather(df_clean, gene_ID, expression, -stages)
df_clean_long$family <- rep(df_DESeq2$Family, each = 7)
df_clean_long$expression <- as.numeric(df_clean_long$expression)
df_clean_long$family <- factor(df_clean_long$family, levels = c("Anterior", "Posterior", "Trunk"))
df_clean_long$stages <- factor(df_clean_long$stages, levels = c("blastula", "gastrula","st4tt", "st5", 
                                                                "st7", "pre-competent", "competent"))

ggplot(df_clean_long, aes(x=stages, y=expression, fill=family)) +
  geom_boxplot(outlier.shape = NA) +  
  facet_wrap(~family) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "DESeq2 normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(-3,3)) +
  scale_y_continuous(breaks = c(seq(-3,3,1))) +
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

# Owenia fusiformis
stages <-  c("blastula", "gastrula","elongation",
             "early larva", "mitraria larva",
             "competent larva","juvenile")
df_DESeq2 <- read.table("05-Fran_Data_Owenia.txt", header = TRUE)
rownames(df_DESeq2) <- df_DESeq2$Transcript_ID


df_copy <- df_DESeq2[,-c(1:10)]
df_norm <- data.frame(t(scale(t(data.matrix(df_copy)))))
df_clean <- data.frame(cbind(stages,t(df_norm)))
df_clean_long <- gather(df_clean, gene_ID, expression, -stages)
df_clean_long$family <- rep(df_DESeq2$Family, each = 7)
df_clean_long$expression <- as.numeric(df_clean_long$expression)
df_clean_long$family <- factor(df_clean_long$family, levels = c("Anterior", "Posterior", "Trunk"))
df_clean_long$stages <- factor(df_clean_long$stages, levels = c("blastula", "gastrula","elongation",
                                                                "early larva", "mitraria larva",
                                                                "competent larva","juvenile"))

ggplot(df_clean_long, aes(x=stages, y=expression, fill=family)) +
  geom_boxplot(outlier.shape = NA) +  
  theme_bw() +
  facet_wrap(~family) +
  theme(panel.grid = element_blank()) +
  labs(y = "DESeq2 normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(-3,3)) +
  scale_y_continuous(breaks = c(seq(-3,3,1))) +
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
