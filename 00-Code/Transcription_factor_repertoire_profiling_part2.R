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


gene_count <- c(31907,41219,17432) #ofus, #ctel, #dgyr, in that order
column_names <- c("pfam", "tfclass", "count", "geneID")
species_palette <- viridis(3)
pfam_tf_list <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/16-Transcription_factors_Owenia_Capitella_Dimorphilus/01-TF_PFAM_domains.txt", sep ='\t', header = T)
tfclass <- pfam_tf_list$Custom_abbreviation

ofus_tf_stats <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/16-Transcription_factors_Owenia_Capitella_Dimorphilus/07-Owenia_fusiformis_TF_stats.txt", sep='\t', header = F)
ctel_tf_stats <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/16-Transcription_factors_Owenia_Capitella_Dimorphilus/07-Capitella_teleta_TF_stats.txt", sep='\t', header = F)
dgyr_tf_stats <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/16-Transcription_factors_Owenia_Capitella_Dimorphilus/07-Dimorphilus_gyrociliatus_TF_stats.txt", sep='\t', header = F)

colnames(ofus_tf_stats) <- column_names
colnames(ctel_tf_stats) <- column_names
colnames(dgyr_tf_stats) <- column_names




# 1. Absolute number of transcription factor transcripts

combined_stats <- cbind(ofus_tf_stats[c(1:2)], 'Owenia_fusiformis' = ofus_tf_stats$count, 'Capitella_teleta' = ctel_tf_stats$count, 'Dimorphilus_gyrociliatus' = dgyr_tf_stats$count)
combined_stats <- combined_stats[match(tfclass, combined_stats$tfclass),]
combined_stats$tfclass <- factor(combined_stats$tfclass, levels = tfclass)
combined_stats_tidy <- gather(combined_stats, species, count, -pfam, -tfclass)
combined_stats_tidy$species <- factor(combined_stats_tidy$species, levels = c("Owenia_fusiformis", "Capitella_teleta", "Dimorphilus_gyrociliatus"))
combined_stats_tidy$tfclass <- factor(combined_stats_tidy$tfclass, levels = tfclass)

ggplot(combined_stats_tidy, aes(fill = species, y = count, x = species)) + 
  geom_bar(position="stack", stat="identity") +
  scale_x_discrete(labels = c("Owenia fusiformis", "Capitella teleta", "Dimorphilus gyrociliatus")) +
  labs(y = "# of transcription factors", x = NA) +
  scale_fill_viridis_d() +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
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
  
ggsave("09-Barplot_TF_number_absolute.pdf",
       width = 3, height = 6)

# 2. Percentage of genome of transcription factors

combined_stats_proportion <- combined_stats
combined_stats_proportion$Owenia_fusiformis <- combined_stats_proportion$Owenia_fusiformis/gene_count[1]*100
combined_stats_proportion$Capitella_teleta <- combined_stats_proportion$Capitella_teleta/gene_count[2]*100
combined_stats_proportion$Dimorphilus_gyrociliatus <- combined_stats_proportion$Dimorphilus_gyrociliatus/gene_count[3]*100

combined_stats_proportion_tidy <- gather(combined_stats_proportion, species, count, -pfam, -tfclass)
combined_stats_proportion_tidy$species <- factor(combined_stats_proportion_tidy$species, levels = c("Owenia_fusiformis", "Capitella_teleta", "Dimorphilus_gyrociliatus"))
combined_stats_proportion_tidy$tfclass <- factor(combined_stats_proportion_tidy$tfclass, levels = tfclass)

ggplot(combined_stats_proportion_tidy, aes(fill = species, y = count, x = species)) + 
  geom_bar(position="stack", stat="identity") +
  scale_x_discrete(labels = c("Owenia fusiformis", "Capitella teleta", "Dimorphilus gyrociliatus")) +
  labs(y = "% of transcription factors", x = NA) +
  scale_fill_viridis_d() +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
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

ggsave("10-Barplot_TF_number_relative.pdf",
       width = 3, height = 6)

# 3. Pie chart of transcription factors by class

combined_stats_relative <- combined_stats
combined_stats_relative$Owenia_fusiformis <- round(combined_stats_relative$Owenia_fusiformis/sum(combined_stats_relative$Owenia_fusiformis)*100,1)
combined_stats_relative$Capitella_teleta <- round(combined_stats_relative$Capitella_teleta/sum(combined_stats_relative$Capitella_teleta)*100,1)
combined_stats_relative$Dimorphilus_gyrociliatus <- round(combined_stats_relative$Dimorphilus_gyrociliatus/sum(combined_stats_relative$Dimorphilus_gyrociliatus)*100,1)

combined_stats_relative_tidy <- gather(combined_stats_relative, species, percentage, -pfam, -tfclass)
combined_stats_relative_tidy$count <- combined_stats_tidy$count
combined_stats_relative_tidy$species <- factor(combined_stats_relative_tidy$species, levels = c("Owenia_fusiformis", "Capitella_teleta", "Dimorphilus_gyrociliatus"))
combined_stats_relative_tidy$tfclass <- factor(combined_stats_relative_tidy$tfclass, levels = tfclass)

ofus_piechart <- combined_stats_relative_tidy[combined_stats_relative_tidy$species == 'Owenia_fusiformis',] %>% 
  mutate(end = 2 * pi * cumsum(combined_stats_relative_tidy[combined_stats_relative_tidy$species == 'Owenia_fusiformis',"percentage"])/sum(combined_stats_relative_tidy[combined_stats_relative_tidy$species == 'Owenia_fusiformis',"percentage"]),
         start = lag(end, default = 0),
         middle = 0.5 * (start + end),
         hjust = ifelse(middle > pi, 1, 0),
         vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))

ctel_piechart <- combined_stats_relative_tidy[combined_stats_relative_tidy$species == 'Capitella_teleta',] %>% 
  mutate(end = 2 * pi * cumsum(combined_stats_relative_tidy[combined_stats_relative_tidy$species == 'Capitella_teleta',"percentage"])/sum(combined_stats_relative_tidy[combined_stats_relative_tidy$species == 'Capitella_teleta',"percentage"]),
         start = lag(end, default = 0),
         middle = 0.5 * (start + end),
         hjust = ifelse(middle > pi, 1, 0),
         vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))

dgyr_piechart <- combined_stats_relative_tidy[combined_stats_relative_tidy$species == 'Dimorphilus_gyrociliatus',] %>% 
  mutate(end = 2 * pi * cumsum(combined_stats_relative_tidy[combined_stats_relative_tidy$species == 'Dimorphilus_gyrociliatus',"percentage"])/sum(combined_stats_relative_tidy[combined_stats_relative_tidy$species == 'Dimorphilus_gyrociliatus',"percentage"]),
         start = lag(end, default = 0),
         middle = 0.5 * (start + end),
         hjust = ifelse(middle > pi, 1, 0),
         vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))

ggplot(ofus_piechart) + 
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                   start = start, end = end, fill = tfclass)) +
  geom_text(aes(x = 1.05 * sin(middle), y = 1.05 * cos(middle), label = paste0(tfclass,": ",percentage,"% (",count,")"),
                hjust = hjust, vjust = vjust)) +
  coord_fixed() +
  scale_x_continuous(limits = c(-2,2),
                     name = "", breaks = NULL, labels = NULL) +
  scale_y_continuous(limits = c(-1.5, 1.5),
                     name = "", breaks = NULL, labels = NULL) +
  theme_void()

ggsave("11a-Piechart_Owenia_fusiformis_TF_class.pdf",
       width = 10.5, height = 6)

ggplot(ctel_piechart) + 
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                   start = start, end = end, fill = tfclass)) +
  geom_text(aes(x = 1.05 * sin(middle), y = 1.05 * cos(middle), label = paste0(tfclass,": ",percentage,"% (",count,")"),
                hjust = hjust, vjust = vjust)) +
  coord_fixed() +
  scale_x_continuous(limits = c(-2,2),
                     name = "", breaks = NULL, labels = NULL) +
  scale_y_continuous(limits = c(-1.5, 1.5),
                     name = "", breaks = NULL, labels = NULL) +
  theme_void()

ggsave("11b-Piechart_Capitella_teleta_TF_class.pdf",
       width = 10.5, height = 6)

ggplot(dgyr_piechart) + 
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                   start = start, end = end, fill = tfclass)) +
  geom_text(aes(x = 1.05 * sin(middle), y = 1.05 * cos(middle), label = paste0(tfclass,": ",percentage,"% (",count,")"),
                hjust = hjust, vjust = vjust)) +
  coord_fixed() +
  scale_x_continuous(limits = c(-2,2),
                     name = "", breaks = NULL, labels = NULL) +
  scale_y_continuous(limits = c(-1.5, 1.5),
                     name = "", breaks = NULL, labels = NULL) +
  theme_void()

ggsave("11c-Piechart_Dimorphilus_gyrociliatus_TF_class.pdf",
       width = 10.5, height = 6)

# 4. Line plot of # of expressed genes and TFs during development time courses 

ofus_timecourse <- read.table("08-Owenia_fusiformis_RNAseq_TPM_average.txt", sep='\t', header = T)
ctel_timecourse <- read.table("08-Capitella_teleta_RNAseq_TPM_average.txt", sep='\t', header = T)
dgyr_timecourse <- read.table("08-Dimorphilus_gyrociliatus_RNAseq_TPM_average.txt", sep='\t', header = T)
colnames(ctel_timecourse)[2] <- "64_cells"

ofus_tf_list <- read.table("05-Owenia_fusiformis_TF_list_gene_ID_only.txt", sep='\t', header=F)
ofus_tf_list <- ofus_tf_list$V1
ctel_tf_list <- read.table("05-Capitella_teleta_TF_list_gene_ID_only.txt", sep='\t', header=F)
ctel_tf_list <- ctel_tf_list$V1
dgyr_tf_list <- read.table("05-Dimorphilus_gyrociliatus_TF_list_gene_ID_only.txt", sep='\t', header=F)
dgyr_tf_list <- dgyr_tf_list$V1

ofus_tf_timecourse <- subset(ofus_timecourse, Gene_ID %in% ofus_tf_list)
ctel_tf_timecourse <- subset(ctel_timecourse, Gene_ID %in% ctel_tf_list)
dgyr_tf_timecourse <- subset(dgyr_timecourse, Gene_ID %in% dgyr_tf_list)

ofus_xaxis <- c(1:7)
ctel_xaxis <- c(1:7)
dgyr_xaxis <- c(1:4)


# 4a. Owenia fusiformis
ofus_expressed <- cbind(sum((ofus_timecourse$blastula > 2), na.rm = TRUE),
                        sum((ofus_timecourse$gastrula > 2), na.rm = TRUE),
                        sum((ofus_timecourse$elongation > 2), na.rm = TRUE),
                        sum((ofus_timecourse$early_larva > 2), na.rm = TRUE),
                        sum((ofus_timecourse$mitraria_larva > 2), na.rm = TRUE),
                        sum((ofus_timecourse$competent_larva > 2), na.rm = TRUE),
                        sum((ofus_timecourse$juvenile > 2), na.rm = TRUE))

ofus_tf_expressed <- cbind(sum((ofus_tf_timecourse$blastula > 2), na.rm = TRUE),
                           sum((ofus_tf_timecourse$gastrula > 2), na.rm = TRUE),
                           sum((ofus_tf_timecourse$elongation > 2), na.rm = TRUE),
                           sum((ofus_tf_timecourse$early_larva > 2), na.rm = TRUE),
                           sum((ofus_tf_timecourse$mitraria_larva > 2), na.rm = TRUE),
                           sum((ofus_tf_timecourse$competent_larva > 2), na.rm = TRUE),
                           sum((ofus_tf_timecourse$juvenile > 2), na.rm = TRUE))

ofus_lineplot <- data.frame(cbind(ofus_xaxis, t(ofus_expressed), t(ofus_tf_expressed)))

p1 <- ggplot() +
  geom_point(data=ofus_lineplot, aes(x=ofus_xaxis, y=V3), color = species_palette[1], size = 2, alpha = 0.6) +
  geom_line(data=ofus_lineplot, aes(x=ofus_xaxis, y=V3), color = species_palette[1], alpha = 0.6) +
  theme_classic() +
  coord_cartesian(ylim = c(500,1000)) +
  labs(y = "# of expressed TFs") +
  scale_x_continuous(breaks = c(1:7), limits = c(1,7),
                     expand = c(0,0.75),
                     labels = c("blastula", "gastrula", "elongation", 
                                "early larva", "mitraria larva",
                                "competent larva", "juvenile")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black",size=14),
        axis.text.x = element_text(color="black",size=12,angle=90,
                                   hjust=0.95,vjust=0.5,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=12,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"))
  
p2 <- ggplot() +
  geom_point(data=ofus_lineplot, aes(x=ofus_xaxis, y=V2), color = species_palette[1], size = 2) +
  geom_line(data=ofus_lineplot, aes(x=ofus_xaxis, y=V2), color = species_palette[1]) +
  theme_classic() +
  coord_cartesian(ylim = c(12000,22000)) + 
  guides(x = "none")+ # remove x line
  labs(y = "# of expressed genes") +
  scale_y_continuous(breaks = c(12000,14000,16000,18000,20000,22000), limits = c(0,22000)) +
  scale_x_continuous(expand = c(0,0.75)) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(color="black",size=14),
        axis.text.y = element_text(color="black",size=12,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"))
  
p2/p1
ggsave("12a-Lineplot_Owenia_fusiformis_expressed_genes_TF_absolute.pdf",
       width = 3, height = 6)

# 4b. Capitella teleta
ctel_expressed <- cbind(sum((ctel_timecourse$`64_cells` > 2), na.rm = TRUE),
                        sum((ctel_timecourse$gastrula > 2), na.rm = TRUE),
                        sum((ctel_timecourse$st4tt_larva > 2), na.rm = TRUE),
                        sum((ctel_timecourse$st5_larva > 2), na.rm = TRUE),
                        sum((ctel_timecourse$st7_larva > 2), na.rm = TRUE),
                        sum((ctel_timecourse$precompetent_larva > 2), na.rm = TRUE),
                        sum((ctel_timecourse$competent_larva > 2), na.rm = TRUE))

ctel_tf_expressed <- cbind(sum((ctel_tf_timecourse$`64_cells` > 2), na.rm = TRUE),
                           sum((ctel_tf_timecourse$gastrula > 2), na.rm = TRUE),
                           sum((ctel_tf_timecourse$st4tt_larva > 2), na.rm = TRUE),
                           sum((ctel_tf_timecourse$st5_larva > 2), na.rm = TRUE),
                           sum((ctel_tf_timecourse$st7_larva > 2), na.rm = TRUE),
                           sum((ctel_tf_timecourse$precompetent_larva > 2), na.rm = TRUE),
                           sum((ctel_tf_timecourse$competent_larva > 2), na.rm = TRUE))

ctel_lineplot <- data.frame(cbind(ctel_xaxis, t(ctel_expressed), t(ctel_tf_expressed)))

p1 <- ggplot() +
  geom_point(data=ctel_lineplot, aes(x=ctel_xaxis, y=V3), color = species_palette[2], size = 2, alpha = 0.6) +
  geom_line(data=ctel_lineplot, aes(x=ctel_xaxis, y=V3), color = species_palette[2], alpha = 0.6) +
  theme_classic() +
  coord_cartesian(ylim = c(300,800)) +
  labs(y = "# of expressed TFs") +
  scale_x_continuous(breaks = c(1:7), limits = c(1,7),
                     expand = c(0,0.75),
                     labels = c("64 cells", "gastrula", "stage 4tt larva", 
                                "stage 5 larva", "stage 7 larva",
                                "pre-competent larva", "competent larva")) +
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

p2 <- ggplot() +
  geom_point(data=ctel_lineplot, aes(x=ctel_xaxis, y=V2), color = species_palette[2], size = 2) +
  geom_line(data=ctel_lineplot, aes(x=ctel_xaxis, y=V2), color = species_palette[2]) +
  theme_classic() +
  coord_cartesian(ylim = c(12000,22000)) + 
  guides(x = "none")+ # remove x line
  labs(y = "# of expressed genes") +
  scale_y_continuous(breaks = c(12000,14000,16000,18000,20000,22000), limits = c(0,22000))  +
  scale_x_continuous(expand = c(0,0.75)) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(color="black",size=14),
        axis.text.y = element_text(color="black",size=12,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"))

p2/p1
ggsave("12b-Lineplot_Capitella_teleta_expressed_genes_TF_absolute.pdf",
       width = 3, height = 6)


# 4c. Dimorphilus gyrociliatus
dgyr_expressed <- cbind(sum((dgyr_timecourse$early_development > 2), na.rm = TRUE),
                        sum((dgyr_timecourse$late_development > 2), na.rm = TRUE),
                        sum((dgyr_timecourse$hatchling > 2), na.rm = TRUE),
                        sum((dgyr_timecourse$female_adult > 2), na.rm = TRUE))

dgyr_tf_expressed <- cbind(sum((dgyr_tf_timecourse$early_development > 2), na.rm = TRUE),
                           sum((dgyr_tf_timecourse$late_development > 2), na.rm = TRUE),
                           sum((dgyr_tf_timecourse$hatchling > 2), na.rm = TRUE),
                           sum((dgyr_tf_timecourse$female_adult > 2), na.rm = TRUE))

dgyr_lineplot <- data.frame(cbind(dgyr_xaxis, t(dgyr_expressed), t(dgyr_tf_expressed)))

p1 <- ggplot() +
  geom_point(data=dgyr_lineplot, aes(x=dgyr_xaxis, y=V3), color = species_palette[3], size = 2, alpha = 0.6) +
  geom_line(data=dgyr_lineplot, aes(x=dgyr_xaxis, y=V3), color = species_palette[3], alpha = 0.6) +
  theme_classic() +
  coord_cartesian(ylim = c(200,600)) +
  labs(y = "# of expressed TFs") +
  scale_x_continuous(breaks = c(1:4), limits = c(1,4),
                     expand = c(0,0.75),
                     labels = c("early development", "late development",
                                "hatchling", "female adult")) +
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

p2 <- ggplot() +
  geom_point(data=dgyr_lineplot, aes(x=dgyr_xaxis, y=V2), color = species_palette[3], size = 2) +
  geom_line(data=dgyr_lineplot, aes(x=dgyr_xaxis, y=V2), color = species_palette[3]) +
  theme_classic() +
  coord_cartesian(ylim = c(8000,16000)) + 
  guides(x = "none")+ # remove x line
  labs(y = "# of expressed genes") +
  scale_y_continuous(breaks = c(8000, 10000, 12000,14000,16000), limits = c(8000,16000))  +
  scale_x_continuous(expand = c(0,0.75)) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(color="black",size=14),
        axis.text.y = element_text(color="black",size=12,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"))

p2/p1
ggsave("12c-Lineplot_Dimorphilus_gyrociliatus_expressed_genes_TF_absolute.pdf",
       width = 3, height = 6)


# 5. Line plot of % of expressed genes and TFs during development time courses
ofus_lineplot_relative <- ofus_lineplot
ofus_lineplot_relative[c(2,3)] <- ofus_lineplot[c(2,3)]/gene_count[1]*100

ctel_lineplot_relative <- ctel_lineplot
ctel_lineplot_relative[c(2,3)] <- ctel_lineplot[c(2,3)]/gene_count[2]*100

dgyr_lineplot_relative <- dgyr_lineplot
dgyr_lineplot_relative[c(2,3)] <- dgyr_lineplot[c(2,3)]/gene_count[3]*100


# 5a. Owenia fusiformis
p1 <- ggplot() +
  geom_point(data=ofus_lineplot_relative, aes(x=ofus_xaxis, y=V3), color = species_palette[1], size = 2, alpha = 0.6) +
  geom_line(data=ofus_lineplot_relative, aes(x=ofus_xaxis, y=V3), color = species_palette[1], alpha = 0.6) +
  theme_classic() +
  coord_cartesian(ylim = c(0.5,3)) +
  labs(y = "% of expressed TFs") +
  scale_x_continuous(breaks = c(1:7), limits = c(1,7),
                     expand = c(0,0.75),
                     labels = c("blastula", "gastrula", "elongation", 
                                "early larva", "mitraria larva",
                                "competent larva", "juvenile")) +
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

p2 <- ggplot() +
  geom_point(data=ofus_lineplot_relative, aes(x=ofus_xaxis, y=V2), color = species_palette[1], size = 2) +
  geom_line(data=ofus_lineplot_relative, aes(x=ofus_xaxis, y=V2), color = species_palette[1]) +
  theme_classic() +
  coord_cartesian(ylim = c(30,90)) + 
  guides(x = "none")+ # remove x line
  labs(y = "% of expressed genes") +
  scale_y_continuous(breaks = c(30,40,50,60,70,80,90), limits = c(30,90))  +
  scale_x_continuous(expand = c(0,0.75)) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(color="black",size=14),
        axis.text.y = element_text(color="black",size=12,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"))

p2/p1
ggsave("13a-Lineplot_Owenia_fusiformis_expressed_genes_TF_relative.pdf",
       width = 3, height = 6)



# 5b. Capitella teleta
p1 <- ggplot() +
  geom_point(data=ctel_lineplot_relative, aes(x=ctel_xaxis, y=V3), color = species_palette[2], size = 2, alpha = 0.6) +
  geom_line(data=ctel_lineplot_relative, aes(x=ctel_xaxis, y=V3), color = species_palette[2], alpha = 0.6) +
  theme_classic() +
  coord_cartesian(ylim = c(0.5,3)) +
  labs(y = "% of expressed TFs") +
  scale_x_continuous(breaks = c(1:7), limits = c(1,7),
                     expand = c(0,0.75),
                     labels = c("64 cells", "gastrula", "stage 4tt larva", 
                                "stage 5 larva", "stage 7 larva",
                                "pre-competent larva", "competent larva")) +
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

p2 <- ggplot() +
  geom_point(data=ctel_lineplot_relative, aes(x=ctel_xaxis, y=V2), color = species_palette[2], size = 2) +
  geom_line(data=ctel_lineplot_relative, aes(x=ctel_xaxis, y=V2), color = species_palette[2]) +
  theme_classic() +
  coord_cartesian(ylim = c(30,90)) + 
  guides(x = "none")+ # remove x line
  labs(y = "% of expressed genes") +
  scale_y_continuous(breaks = c(30,40,50,60,70,80,90), limits = c(30,90))  +
  scale_x_continuous(expand = c(0,0.75)) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(color="black",size=14),
        axis.text.y = element_text(color="black",size=12,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"))

p2/p1
ggsave("13b-Lineplot_Capitella_teleta_expressed_genes_TF_relative.pdf",
       width = 3, height = 6)

# 5c. Dimorphilus gyrociliatus
p1 <- ggplot() +
  geom_point(data=dgyr_lineplot_relative, aes(x=dgyr_xaxis, y=V3), color = species_palette[3], size = 2, alpha = 0.6) +
  geom_line(data=dgyr_lineplot_relative, aes(x=dgyr_xaxis, y=V3), color = species_palette[3], alpha = 0.6) +
  theme_classic() +
  coord_cartesian(ylim = c(0.5,3)) +
  labs(y = "% of expressed TFs") +
  scale_x_continuous(breaks = c(1:4), limits = c(1,4),
                     expand = c(0,0.75),
                     labels = c("early development", "late development",
                                "hatchling", "female adult")) +
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

p2 <- ggplot() +
  geom_point(data=dgyr_lineplot_relative, aes(x=dgyr_xaxis, y=V2), color = species_palette[3], size = 2) +
  geom_line(data=dgyr_lineplot_relative, aes(x=dgyr_xaxis, y=V2), color = species_palette[3]) +
  theme_classic() +
  coord_cartesian(ylim = c(30,90)) + 
  guides(x = "none")+ # remove x line
  labs(y = "% of expressed genes") +
  scale_y_continuous(breaks = c(30,40,50,60,70,80,90), limits = c(30,90))  +
  scale_x_continuous(expand = c(0,0.75)) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(color="black",size=14),
        axis.text.y = element_text(color="black",size=12,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"))

p2/p1
ggsave("13c-Lineplot_Dimorphilus_gyrociliatus_expressed_genes_TF_relative.pdf",
       width = 3, height = 6)


# 6. Stacked barplots with % of TFs by class per species

ggplot(combined_stats_relative_tidy, aes(fill = tfclass, y = percentage, x = species)) + 
  geom_bar(position="stack", stat="identity", color = 'black', size = 0.35) +
  scale_x_discrete(labels = c("Owenia fusiformis", "Capitella teleta", "Dimorphilus gyrociliatus")) +
  labs(y = "% of transcription factors", x = NA) +
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

ggsave("14-Barplot_TF_class.pdf",
       width = 4.5, height = 7.2)


# 7. Barplot of % of TFs expressed per class during development time courses

ofus_tf_list_classified <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/16-Transcription_factors_Owenia_Capitella_Dimorphilus/06-Owenia_fusiformis_TF_list_classified.txt", sep ='\t')
ctel_tf_list_classified <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/16-Transcription_factors_Owenia_Capitella_Dimorphilus/06-Capitella_teleta_TF_list_classified.txt", sep ='\t')
dgyr_tf_list_classified <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/16-Transcription_factors_Owenia_Capitella_Dimorphilus/06-Dimorphilus_gyrociliatus_TF_list_classified.txt", sep ='\t')

colnames(ofus_tf_list_classified) <- c("Gene_ID", "pfam", "tfclass")
colnames(ctel_tf_list_classified) <- c("Gene_ID", "pfam", "tfclass")
colnames(dgyr_tf_list_classified) <- c("Gene_ID", "pfam", "tfclass")

# 7a. Owenia fusiformis
ofus_tf_timecourse_classified <- cbind(ofus_tf_timecourse, "tfclass" = ofus_tf_list_classified$tfclass)
ofus_tf_timecourse_classified[ofus_tf_timecourse_classified$blastula < 2,2] <- FALSE
ofus_tf_timecourse_classified[ofus_tf_timecourse_classified$gastrula < 2,3] <- FALSE
ofus_tf_timecourse_classified[ofus_tf_timecourse_classified$elongation < 2,4] <- FALSE
ofus_tf_timecourse_classified[ofus_tf_timecourse_classified$early_larva < 2,5] <- FALSE
ofus_tf_timecourse_classified[ofus_tf_timecourse_classified$mitraria_larva < 2,6] <- FALSE
ofus_tf_timecourse_classified[ofus_tf_timecourse_classified$competent_larva < 2,7] <- FALSE
ofus_tf_timecourse_classified[ofus_tf_timecourse_classified$juvenile < 2,8] <- FALSE
ofus_tf_timecourse_classified[ofus_tf_timecourse_classified$blastula > 2,2] <- TRUE
ofus_tf_timecourse_classified[ofus_tf_timecourse_classified$gastrula > 2,3] <- TRUE
ofus_tf_timecourse_classified[ofus_tf_timecourse_classified$elongation > 2,4] <- TRUE
ofus_tf_timecourse_classified[ofus_tf_timecourse_classified$early_larva > 2,5] <- TRUE
ofus_tf_timecourse_classified[ofus_tf_timecourse_classified$mitraria_larva > 2,6] <- TRUE
ofus_tf_timecourse_classified[ofus_tf_timecourse_classified$competent_larva > 2,7] <- TRUE
ofus_tf_timecourse_classified[ofus_tf_timecourse_classified$juvenile > 2,8] <- TRUE

ofus_tf_timecourse_classified_tidy <- gather(ofus_tf_timecourse_classified, stage, expressed_value, -Gene_ID, -tfclass)
ofus_tf_timecourse_classified_tidy_expressed <- ofus_tf_timecourse_classified_tidy[ofus_tf_timecourse_classified_tidy$expressed_value == 1,]

ofus_tf_timecourse_expressed <- data.frame(matrix(NA, nrow = 36, ncol = 7))
colnames(ofus_tf_timecourse_expressed) <- colnames(ofus_tf_timecourse_classified)[c(2:8)]
rownames(ofus_tf_timecourse_expressed) <- tfclass

for (stage in colnames(ofus_tf_timecourse_expressed)){
  for (class in rownames(ofus_tf_timecourse_expressed)){
    ofus_tf_timecourse_expressed[class,stage] <- sum(ofus_tf_timecourse_classified_tidy_expressed[ofus_tf_timecourse_classified_tidy_expressed$stage == stage & ofus_tf_timecourse_classified_tidy_expressed$tfclass == class,"expressed_value"])
  }
}

ofus_tf_timecourse_expressed$tfclass <- rownames(ofus_tf_timecourse_expressed)
ofus_tf_timecourse_expressed_tidy <- gather(ofus_tf_timecourse_expressed, stage, count, -tfclass)

ofus_tf_timecourse_expressed_percentage <- ofus_tf_timecourse_expressed

for (stage in colnames(ofus_tf_timecourse_expressed)[-c(ncol(ofus_tf_timecourse_expressed))]){
  for (class in rownames(ofus_tf_timecourse_expressed)){
    ofus_tf_timecourse_expressed_percentage[class,stage] <- ofus_tf_timecourse_expressed[class,stage]/sum(ofus_tf_timecourse_expressed[,stage])*100
  }
}
  
ofus_tf_timecourse_expressed_percentage$tfclass <- rownames(ofus_tf_timecourse_expressed_percentage)
ofus_tf_timecourse_expressed_percentage_tidy <- gather(ofus_tf_timecourse_expressed_percentage, stage, percentage, -tfclass)

ggplot(ofus_tf_timecourse_expressed_percentage_tidy, aes(fill = tfclass, y = percentage, x = stage)) + 
  geom_bar(position="stack", stat="identity", color = 'black', size = 0.35) +
  scale_x_discrete(labels = c("blastula", "gastrula", "elongation", 
                              "early larva", "mitraria larva",
                              "competent larva","juvenile")) +
  labs(y = "% of transcription factors expressed", x = NA) +
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

ggsave("15a-Barplot_Owenia_fusiformis_TF_class_expressed.pdf",
       width = 5.54, height = 6.86)


# 7b. Capitella teleta
ctel_tf_timecourse_classified <- cbind(ctel_tf_timecourse, "tfclass" = ctel_tf_list_classified$tfclass)
ctel_tf_timecourse_classified[ctel_tf_timecourse_classified$`64_cells` < 2,2] <- FALSE
ctel_tf_timecourse_classified[ctel_tf_timecourse_classified$gastrula < 2,3] <- FALSE
ctel_tf_timecourse_classified[ctel_tf_timecourse_classified$st4tt_larva < 2,4] <- FALSE
ctel_tf_timecourse_classified[ctel_tf_timecourse_classified$st5_larva < 2,5] <- FALSE
ctel_tf_timecourse_classified[ctel_tf_timecourse_classified$st7_larva < 2,6] <- FALSE
ctel_tf_timecourse_classified[ctel_tf_timecourse_classified$precompetent_larva < 2,7] <- FALSE
ctel_tf_timecourse_classified[ctel_tf_timecourse_classified$competent_larva < 2,8] <- FALSE
ctel_tf_timecourse_classified[ctel_tf_timecourse_classified$`64_cells` > 2,2] <- TRUE
ctel_tf_timecourse_classified[ctel_tf_timecourse_classified$gastrula > 2,3] <- TRUE
ctel_tf_timecourse_classified[ctel_tf_timecourse_classified$st4tt_larva > 2,4] <- TRUE
ctel_tf_timecourse_classified[ctel_tf_timecourse_classified$st5_larva > 2,5] <- TRUE
ctel_tf_timecourse_classified[ctel_tf_timecourse_classified$st7_larva > 2,6] <- TRUE
ctel_tf_timecourse_classified[ctel_tf_timecourse_classified$precompetent_larva > 2,7] <- TRUE
ctel_tf_timecourse_classified[ctel_tf_timecourse_classified$competent_larva > 2,8] <- TRUE

ctel_tf_timecourse_classified_tidy <- gather(ctel_tf_timecourse_classified, stage, expressed_value, -Gene_ID, -tfclass)
ctel_tf_timecourse_classified_tidy_expressed <- ctel_tf_timecourse_classified_tidy[ctel_tf_timecourse_classified_tidy$expressed_value == 1,]

ctel_tf_timecourse_expressed <- data.frame(matrix(NA, nrow = 36, ncol = 7))
colnames(ctel_tf_timecourse_expressed) <- colnames(ctel_tf_timecourse_classified)[c(2:8)]
rownames(ctel_tf_timecourse_expressed) <- tfclass

for (stage in colnames(ctel_tf_timecourse_expressed)){
  for (class in rownames(ctel_tf_timecourse_expressed)){
    ctel_tf_timecourse_expressed[class,stage] <- sum(ctel_tf_timecourse_classified_tidy_expressed[ctel_tf_timecourse_classified_tidy_expressed$stage == stage & ctel_tf_timecourse_classified_tidy_expressed$tfclass == class,"expressed_value"])
  }
}

ctel_tf_timecourse_expressed$tfclass <- rownames(ctel_tf_timecourse_expressed)
ctel_tf_timecourse_expressed_tidy <- gather(ctel_tf_timecourse_expressed, stage, count, -tfclass)

ctel_tf_timecourse_expressed_percentage <- ctel_tf_timecourse_expressed

for (stage in colnames(ctel_tf_timecourse_expressed)[-c(ncol(ctel_tf_timecourse_expressed))]){
  for (class in rownames(ctel_tf_timecourse_expressed)){
    ctel_tf_timecourse_expressed_percentage[class,stage] <- ctel_tf_timecourse_expressed[class,stage]/sum(ctel_tf_timecourse_expressed[,stage])*100
  }
}

ctel_tf_timecourse_expressed_percentage$tfclass <- rownames(ctel_tf_timecourse_expressed_percentage)
ctel_tf_timecourse_expressed_percentage_tidy <- gather(ctel_tf_timecourse_expressed_percentage, stage, percentage, -tfclass)

ggplot(ctel_tf_timecourse_expressed_percentage_tidy, aes(fill = tfclass, y = percentage, x = stage)) + 
  geom_bar(position="stack", stat="identity", color = 'black', size = 0.35) +
  scale_x_discrete(labels = c("64 cells", "gastrula", "stage 4tt larva", 
                              "stage 5 larva", "stage 7 larva",
                              "pre-competent larva", "competent larva")) +
  labs(y = "% of transcription factors expressed", x = NA) +
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

ggsave("15b-Barplot_Capitella_teleta_TF_class_expressed.pdf",
       width = 5.54, height = 6.86)


# 7c. Dimorphilus gyrociliatus
# DgyrP12862 was not captured in RNAseq developmental timecourse, so it was removed
common <- intersect(dgyr_tf_timecourse$Gene_ID, dgyr_tf_list_classified$Gene_ID)  
dgyr_tf_timecourse_classified <- cbind(dgyr_tf_timecourse, "tfclass" = dgyr_tf_list_classified[dgyr_tf_list_classified$Gene_ID %in% common,"tfclass"])

dgyr_tf_timecourse_classified[dgyr_tf_timecourse_classified$early_development < 2,2] <- FALSE
dgyr_tf_timecourse_classified[dgyr_tf_timecourse_classified$late_development < 2,3] <- FALSE
dgyr_tf_timecourse_classified[dgyr_tf_timecourse_classified$hatchling < 2,4] <- FALSE
dgyr_tf_timecourse_classified[dgyr_tf_timecourse_classified$female_adult < 2,5] <- FALSE
dgyr_tf_timecourse_classified[dgyr_tf_timecourse_classified$early_development > 2,2] <- TRUE
dgyr_tf_timecourse_classified[dgyr_tf_timecourse_classified$late_development > 2,3] <- TRUE
dgyr_tf_timecourse_classified[dgyr_tf_timecourse_classified$hatchling > 2,4] <- TRUE
dgyr_tf_timecourse_classified[dgyr_tf_timecourse_classified$female_adult > 2,5] <- TRUE

dgyr_tf_timecourse_classified_tidy <- gather(dgyr_tf_timecourse_classified, stage, expressed_value, -Gene_ID, -tfclass)
dgyr_tf_timecourse_classified_tidy_expressed <- dgyr_tf_timecourse_classified_tidy[dgyr_tf_timecourse_classified_tidy$expressed_value == 1,]

dgyr_tf_timecourse_expressed <- data.frame(matrix(NA, nrow = 36, ncol = 4))
colnames(dgyr_tf_timecourse_expressed) <- colnames(dgyr_tf_timecourse_classified)[c(2:5)]
rownames(dgyr_tf_timecourse_expressed) <- tfclass

for (stage in colnames(dgyr_tf_timecourse_expressed)){
  for (class in rownames(dgyr_tf_timecourse_expressed)){
    dgyr_tf_timecourse_expressed[class,stage] <- sum(dgyr_tf_timecourse_classified_tidy_expressed[dgyr_tf_timecourse_classified_tidy_expressed$stage == stage & dgyr_tf_timecourse_classified_tidy_expressed$tfclass == class,"expressed_value"])
  }
}

dgyr_tf_timecourse_expressed$tfclass <- rownames(dgyr_tf_timecourse_expressed)
dgyr_tf_timecourse_expressed_tidy <- gather(dgyr_tf_timecourse_expressed, stage, count, -tfclass)

dgyr_tf_timecourse_expressed_percentage <- dgyr_tf_timecourse_expressed

for (stage in colnames(dgyr_tf_timecourse_expressed)[-c(ncol(dgyr_tf_timecourse_expressed))]){
  for (class in rownames(dgyr_tf_timecourse_expressed)){
    dgyr_tf_timecourse_expressed_percentage[class,stage] <- dgyr_tf_timecourse_expressed[class,stage]/sum(dgyr_tf_timecourse_expressed[,stage])*100
  }
}

dgyr_tf_timecourse_expressed_percentage$tfclass <- rownames(dgyr_tf_timecourse_expressed_percentage)
dgyr_tf_timecourse_expressed_percentage_tidy <- gather(dgyr_tf_timecourse_expressed_percentage, stage, percentage, -tfclass)

ggplot(dgyr_tf_timecourse_expressed_percentage_tidy, aes(fill = tfclass, y = percentage, x = stage)) + 
  geom_bar(position="stack", stat="identity", color = 'black', size = 0.35) +
  scale_x_discrete(labels = c("early development", "late development", "hatchling","female adult")) +
  labs(y = "% of transcription factors expressed", x = NA) +
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

ggsave("15c-Barplot_Dimorphilus_gyrociliatus_TF_class_expressed.pdf",
       width = 5.54, height = 6.86)



# 8. Heatmaps of expression dynamics of TF classes
heatmap_color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

# 8a. Owenia fusiformis
stages <- c("blastula","gastrula","elongation","early_larva",
            "mitraria_larva","competent_larva","juvenile")
ofus_expression <- cbind(ofus_tf_timecourse, "tfclass" = ofus_tf_list_classified$tfclass)
rownames(ofus_expression) <- ofus_expression$Gene_ID

ofus_notnull <- ofus_expression[rowSums(ofus_expression[,-c(1,9)])>0,]
ofus_notnull_qn <- cbind(ofus_notnull[c(1,ncol(ofus_notnull))],normalize.quantiles(data.matrix(ofus_notnull[-c(1,ncol(ofus_notnull))]))) # quantile normalisation of data
ofus_notnull_qn$tfclass <- factor(ofus_notnull_qn$tfclass, levels = tfclass)

colnames(ofus_notnull_qn)[3:ncol(ofus_notnull_qn)] <- stages

ofus_quantile_75 <- data.frame(matrix(ncol = length(stages), nrow = length(tfclass)))
rownames(ofus_quantile_75) <- c(tfclass)
colnames(ofus_quantile_75) <- stages

for (i in stages){
  ofus_quantile_75[i] <- round(tapply(ofus_notnull_qn[,i], ofus_notnull_qn$tfclass, quantile, p = 0.75, na.rm = TRUE),2)
}

pheatmap(ofus_quantile_75, 
         scale = "row",
         cluster_rows=T,
         cluster_cols=F, 
         show_rownames=T, 
         fontsize= 12,
         col=heatmap_color, 
         border_color=NA, 
         cellwidth = 10, 
         cellheight = 10)

pheatmap(ofus_quantile_75, 
         scale = "row",
         cluster_rows=F,
         cluster_cols=F, 
         show_rownames=T, 
         fontsize= 12,
         col=heatmap_color, 
         border_color=NA, 
         cellwidth = 10, 
         cellheight = 10)


# 8b. Capitella teleta
stages <- c("64_cells", "gastrula", "st4tt_larva"," st5_larva", "st7_larva",
            "precompetent_larva", "competent_larva")
ctel_expression <- cbind(ctel_tf_timecourse, "tfclass" = ctel_tf_list_classified$tfclass)
rownames(ctel_expression) <- ctel_expression$Gene_ID

ctel_notnull <- ctel_expression[rowSums(ctel_expression[,-c(1,9)])>0,]
ctel_notnull_qn <- cbind(ctel_notnull[c(1,ncol(ctel_notnull))],normalize.quantiles(data.matrix(ctel_notnull[-c(1,ncol(ctel_notnull))]))) # quantile normalisation of data
ctel_notnull_qn$tfclass <- factor(ctel_notnull_qn$tfclass, levels = tfclass)

colnames(ctel_notnull_qn)[3:ncol(ctel_notnull_qn)] <- stages

ctel_quantile_75 <- data.frame(matrix(ncol = length(stages), nrow = length(tfclass)))
rownames(ctel_quantile_75) <- c(tfclass)
colnames(ctel_quantile_75) <- stages

for (i in stages){
  ctel_quantile_75[i] <- round(tapply(ctel_notnull_qn[,i], ctel_notnull_qn$tfclass, quantile, p = 0.75, na.rm = TRUE),2)
}

pheatmap(ctel_quantile_75, 
         scale = "row",
         cluster_rows=T,
         cluster_cols=F, 
         show_rownames=T, 
         fontsize= 12,
         col=heatmap_color, 
         border_color=NA, 
         cellwidth = 10, 
         cellheight = 10)

pheatmap(ctel_quantile_75, 
         scale = "row",
         cluster_rows=F,
         cluster_cols=F, 
         show_rownames=T, 
         fontsize= 12,
         col=heatmap_color, 
         border_color=NA, 
         cellwidth = 10, 
         cellheight = 10)



# 8c. Dimorphilus gyrociliatus
# DgyrP12862 was not captured in RNAseq developmental timecourse, so it was removed
stages <- c("early_development","late_development","hatchling","female_adult")
common <- intersect(dgyr_tf_timecourse$Gene_ID, dgyr_tf_list_classified$Gene_ID)  
dgyr_expression <- cbind(dgyr_tf_timecourse, "tfclass" = dgyr_tf_list_classified[dgyr_tf_list_classified$Gene_ID %in% common,"tfclass"])
rownames(dgyr_expression) <- dgyr_expression$Gene_ID

dgyr_notnull <- dgyr_expression[rowSums(dgyr_expression[,-c(1,6)])>0,]
dgyr_notnull_qn <- cbind(dgyr_notnull[c(1,ncol(dgyr_notnull))],normalize.quantiles(data.matrix(dgyr_notnull[-c(1,ncol(dgyr_notnull))]))) # quantile normalisation of data
dgyr_notnull_qn$tfclass <- factor(dgyr_notnull_qn$tfclass, levels = tfclass)

colnames(dgyr_notnull_qn)[3:ncol(dgyr_notnull_qn)] <- stages

dgyr_quantile_75 <- data.frame(matrix(ncol = length(stages), nrow = length(tfclass)))
rownames(dgyr_quantile_75) <- c(tfclass)
colnames(dgyr_quantile_75) <- stages

for (i in stages){
  dgyr_quantile_75[i] <- round(tapply(dgyr_notnull_qn[,i], dgyr_notnull_qn$tfclass, quantile, p = 0.75, na.rm = TRUE),2)
}

pheatmap(dgyr_quantile_75, 
         scale = "row",
         cluster_rows=T,
         cluster_cols=F, 
         show_rownames=T, 
         fontsize= 12,
         col=heatmap_color, 
         border_color=NA, 
         cellwidth = 10, 
         cellheight = 10)

pheatmap(dgyr_quantile_75, 
         scale = "row",
         cluster_rows=F,
         cluster_cols=F, 
         show_rownames=T, 
         fontsize= 12,
         col=heatmap_color, 
         border_color=NA, 
         cellwidth = 10, 
         cellheight = 10)


## EXPORT FINAL DATASETS
write.table(ofus_notnull_qn, "17a-Owenia_fusiformis_TF_expression_DESeq2_average.txt", sep ='\t', quote = F, row.names = F)
write.table(ctel_notnull_qn, "17b-Capitella_teleta_TF_expression_DESeq2_average.txt", sep ='\t', quote = F, row.names = F)
write.table(dgyr_notnull_qn, "17c-Dimorphilus_gyrociliatus_TF_expression_DESeq2_average.txt", sep ='\t', quote = F, row.names = F)



