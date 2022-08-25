library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(corrplot)
library(philentropy)
library(tidyr)
library(dplyr)
library(preprocessCore)
library(plyr)
library(scales)


###################################
#### A. Calculate gene-wise JSD ###
###################################


# 1. Owenia fusiformis vs. Capitella teleta
# Subset the stage 5 larva and mitraria larva stages

Ofus2Ctel_orth_tpm_qn <- read.csv("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/02-quantile_transformed_data/Ofus2Ctel_TPM_mean_quantile_transform.csv", header = T)
Ctel2Ofus_orth_tpm_qn <- read.csv("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/02-quantile_transformed_data/Ctel2Ofus_TPM_mean_quantile_transform.csv", header = T)
both_ID <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/01-one2one_orthologs/Ofus2Ctel.txt", header = T)
match_Ofus <- Ofus2Ctel_orth_tpm_qn[match(both_ID$Ofus,Ofus2Ctel_orth_tpm_qn$Gene_ID),]
match_Ctel <- Ctel2Ofus_orth_tpm_qn[match(both_ID$Ctel,Ctel2Ofus_orth_tpm_qn$Gene_ID),]
match_Ctel_st5 <- data.frame(cbind("Gene_ID" = match_Ctel$Gene_ID, "st5" = match_Ctel$Ctel_TPM_mean_st5))
match_Ofus_mitraria <- data.frame(cbind("Gene_ID" = match_Ofus$Gene_ID, "mitraria" = match_Ofus$Ofus_TPM_mean_27h))

Ofus2Ctel_genewise_JSD <- data.frame(cbind("Ctel" = match_Ctel$Gene_ID, "Ofus" = match_Ofus$Gene_ID, 
                                           "gene_wise_JSD" = replicate(nrow(match_Ctel),0), "x_variable" = replicate(nrow(match_Ctel),"A")))

for (i in 1:length(rownames(match_Ctel_st5))){
  match <- rbind(as.numeric(match_Ctel_st5[i,2]), as.numeric(match_Ofus_mitraria[i,2]))
  Ofus2Ctel_genewise_JSD[i,3] <- JSD(match)
}

Ofus2Ctel_genewise_JSD$gene_wise_JSD <- as.numeric(Ofus2Ctel_genewise_JSD$gene_wise_JSD)
write.table(Ofus2Ctel_genewise_JSD[-c(4)],
            "~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/03-results/01-Ofus_mitraria_Ctel_stage5.txt",
            quote = F, sep = '\t', row.names = F)


# 2. Owenia fusiformis vs. Crassostrea gigas
# Subset the umbo 1 and mitraria larva stages

Ofus2Cgig_orth_tpm_qn <- read.csv("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/02-quantile_transformed_data/Ofus2Cgig_TPM_mean_quantile_transform.csv", header = T)
Cgig2Ofus_orth_tpm_qn <- read.csv("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/02-quantile_transformed_data/Cgig2Ofus_TPM_mean_quantile_transform.csv", header = T)
both_ID <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/01-one2one_orthologs/Ofus2Cgig.txt", header = T)
match_Ofus <- Ofus2Cgig_orth_tpm_qn[match(both_ID$Ofus,Ofus2Cgig_orth_tpm_qn$Gene_ID),]
match_Cgig <- Cgig2Ofus_orth_tpm_qn[match(both_ID$Cgig,Cgig2Ofus_orth_tpm_qn$Gene_ID),]
match_Cgig_umbo1 <- data.frame(cbind("Gene_ID" = match_Cgig$Gene_ID, "umbo1" = match_Cgig$umbo_larva_1))
match_Ofus_mitraria <- data.frame(cbind("Gene_ID" = match_Ofus$Gene_ID, "mitraria" = match_Ofus$Ofus_TPM_mean_27h))

Ofus2Cgig_genewise_JSD <- data.frame(cbind("Cgig" = match_Cgig$Gene_ID, "Ofus" = match_Ofus$Gene_ID, 
                                           "gene_wise_JSD" = replicate(nrow(match_Cgig),0), "x_variable" = replicate(nrow(match_Cgig),"A")))

for (i in 1:length(rownames(match_Cgig_umbo1))){
  match <- rbind(as.numeric(match_Cgig_umbo1[i,2]), as.numeric(match_Ofus_mitraria[i,2]))
  Ofus2Cgig_genewise_JSD[i,3] <- JSD(match)
}

Ofus2Cgig_genewise_JSD$gene_wise_JSD <- as.numeric(Ofus2Cgig_genewise_JSD$gene_wise_JSD)
write.table(Ofus2Cgig_genewise_JSD[-c(4)],
            "~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/03-results/02-Ofus_mitraria_Cgig_umbo1.txt",
            quote = F, sep = '\t', row.names = F)


# 3. Owenia fusiformis vs. Strongylocentrotus purpuratus
# Subset the pluteus and mitraria larva stages

Ofus2Spur_orth_tpm_qn <- read.csv("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/02-quantile_transformed_data/Ofus2Spur_TPM_mean_quantile_transform.csv", header = T)
Spur2Ofus_orth_tpm_qn <- read.csv("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/02-quantile_transformed_data/Spur2Ofus_TPM_mean_quantile_transform.csv", header = T)
both_ID <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/01-one2one_orthologs/Ofus2Spur.txt", header = T)
match_Ofus <- Ofus2Spur_orth_tpm_qn[match(both_ID$Ofus,Ofus2Spur_orth_tpm_qn$Gene_ID),]
match_Spur <- Spur2Ofus_orth_tpm_qn[match(both_ID$Spur,Spur2Ofus_orth_tpm_qn$Gene_ID),]
match_Spur_pluteus <- data.frame(cbind("Gene_ID" = match_Spur$Gene_ID, "pluteus" = match_Spur$X72.h))
match_Ofus_mitraria <- data.frame(cbind("Gene_ID" = match_Ofus$Gene_ID, "mitraria" = match_Ofus$Ofus_TPM_mean_27h))

Ofus2Spur_genewise_JSD <- data.frame(cbind("Spur" = match_Spur$Gene_ID, "Ofus" = match_Ofus$Gene_ID, 
                                           "gene_wise_JSD" = replicate(nrow(match_Spur),0), "x_variable" = replicate(nrow(match_Spur),"A")))

for (i in 1:length(rownames(match_Spur_pluteus))){
  match <- rbind(as.numeric(match_Spur_pluteus[i,2]), as.numeric(match_Ofus_mitraria[i,2]))
  Ofus2Spur_genewise_JSD[i,3] <- JSD(match)
}

Ofus2Spur_genewise_JSD$gene_wise_JSD <- as.numeric(Ofus2Spur_genewise_JSD$gene_wise_JSD)
write.table(Ofus2Spur_genewise_JSD[-c(4)],
            "~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/03-results/03-Ofus_mitraria_Spur_pluteus.txt",
            quote = F, sep = '\t', row.names = F)


# 4. Owenia fusiformis vs. Nematostella vectensis
# Subset the 72 hpf (planula larva) and mitraria larva stages

Ofus2Nvec_orth_tpm_qn <- read.csv("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/02-quantile_transformed_data/Ofus2Nvec_TPM_mean_quantile_transform.csv", header = T)
Nvec2Ofus_orth_tpm_qn <- read.csv("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/02-quantile_transformed_data/Nvec2Ofus_TPM_mean_quantile_transform.csv", header = T)
both_ID <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/01-one2one_orthologs/Ofus2Nvec.txt", header = T)
match_Ofus <- Ofus2Nvec_orth_tpm_qn[match(both_ID$Ofus,Ofus2Nvec_orth_tpm_qn$Gene_ID),]
match_Nvec <- Nvec2Ofus_orth_tpm_qn[match(both_ID$Nvec,Nvec2Ofus_orth_tpm_qn$Gene_ID),]
match_Nvec_72hpf <- data.frame(cbind("Gene_ID" = match_Nvec$Gene_ID, "72_hpf" = match_Nvec$h72_mean_TPM))
match_Ofus_mitraria <- data.frame(cbind("Gene_ID" = match_Ofus$Gene_ID, "mitraria" = match_Ofus$Ofus_TPM_mean_27h))

Ofus2Nvec_genewise_JSD <- data.frame(cbind("Nvec" = match_Nvec$Gene_ID, "Ofus" = match_Ofus$Gene_ID, 
                                           "gene_wise_JSD" = replicate(nrow(match_Nvec),0), "x_variable" = replicate(nrow(match_Nvec),"A")))

for (i in 1:length(rownames(match_Nvec_72hpf))){
  match <- rbind(as.numeric(match_Nvec_72hpf[i,2]), as.numeric(match_Ofus_mitraria[i,2]))
  Ofus2Nvec_genewise_JSD[i,3] <- JSD(match)
}

Ofus2Nvec_genewise_JSD$gene_wise_JSD <- as.numeric(Ofus2Nvec_genewise_JSD$gene_wise_JSD)
write.table(Ofus2Nvec_genewise_JSD[-c(4)],
            "~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/03-results/04-Ofus_mitraria_Nvec_72hpf.txt",
            quote = F, sep = '\t', row.names = F)


# 5. Owenia fusiformis vs. Clytia hemisphaerica
# Subset the planula larva and mitraria larva stages

Ofus2Chem_orth_tpm_qn <- read.csv("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/02-quantile_transformed_data/Ofus2Chem_TPM_mean_quantile_transform.csv", header = T)
Chem2Ofus_orth_tpm_qn <- read.csv("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/02-quantile_transformed_data/Chem2Ofus_TPM_mean_quantile_transform.csv", header = T)
both_ID <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/01-one2one_orthologs/Ofus2Chem.txt", header = T)
match_Ofus <- Ofus2Chem_orth_tpm_qn[match(both_ID$Ofus,Ofus2Chem_orth_tpm_qn$Gene_ID),]
match_Chem <- Chem2Ofus_orth_tpm_qn[match(both_ID$Chem,Chem2Ofus_orth_tpm_qn$Gene_ID),]
match_Chem_48hpf <- data.frame(cbind("Gene_ID" = match_Chem$Gene_ID, "48_hpf" = match_Chem$Planula_48h_mean_TPM))
match_Ofus_mitraria <- data.frame(cbind("Gene_ID" = match_Ofus$Gene_ID, "mitraria" = match_Ofus$Ofus_TPM_mean_27h))

Ofus2Chem_genewise_JSD <- data.frame(cbind("Chem" = match_Chem$Gene_ID, "Ofus" = match_Ofus$Gene_ID, 
                                           "gene_wise_JSD" = replicate(nrow(match_Chem),0), "x_variable" = replicate(nrow(match_Chem),"A")))

for (i in 1:length(rownames(match_Chem_48hpf))){
  match <- rbind(as.numeric(match_Chem_48hpf[i,2]), as.numeric(match_Ofus_mitraria[i,2]))
  Ofus2Chem_genewise_JSD[i,3] <- JSD(match)
}

Ofus2Chem_genewise_JSD$gene_wise_JSD <- as.numeric(Ofus2Chem_genewise_JSD$gene_wise_JSD)
write.table(Ofus2Chem_genewise_JSD[-c(4)],
            "~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/03-results/05-Ofus_mitraria_Chem_planula.txt",
            quote = F, sep = '\t', row.names = F)






############################################################
#### B. Plot violin plots of gene-wise JSD distribution ####
############################################################

# 0. Colours and aesthetics

colors <- hue_pal()(10)
palette <- colors[c(1,2,7,8,9)]
ctel_color <- colors[1]
cgig_color <- colors[2]
spur_color <- colors[7]
nvec_color <- colors[8]
chem_color <- colors[9]


# 1. Re-import results and prepare data for plotting

columns <- c("Target","Subject","gene_wise_JSD")
ctel <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/03-results/01-Ofus_mitraria_Ctel_stage5.txt", 
                   sep = '\t', header = T)
cgig <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/03-results/02-Ofus_mitraria_Cgig_umbo1.txt", 
                   sep = '\t', header = T)
spur <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/03-results/03-Ofus_mitraria_Spur_pluteus.txt", 
                   sep = '\t', header = T)
nvec <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/03-results/04-Ofus_mitraria_Nvec_72hpf.txt", 
                   sep = '\t', header = T)
chem <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/03-results/05-Ofus_mitraria_Chem_planula.txt", 
                   sep = '\t', header = T)

colnames(ctel) <- columns
colnames(cgig) <- columns
colnames(spur) <- columns
colnames(nvec) <- columns
colnames(chem) <- columns

ctel$species <- "C. teleta"
cgig$species <- "C. gigas"
spur$species <- "S. purpuratus"
nvec$species <- "N. vectensis"
chem$species <- "C. hemisphaerica"

all <- rbind(ctel,cgig,spur,nvec,chem)


# 2. Estimate the threshold point, which we establish at 1/4 of the point with the 
# highest probability of the kernel density estimation. Then subset the genes that
# are below this threshold, and export them.

ctel_density <- density(ctel$gene_wise_JSD)
ctel_threshold <- ctel_density$x[which.max(ctel_density$y)]/4
ctel_gene_set <- ctel[ctel$gene_wise_JSD < ctel_threshold,]
write.table(ctel_gene_set,
            "~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/05-gene_sets/01-Gene_set_ctel_stage5.txt",
            quote = F, sep = '\t', row.names = F)

cgig_density <- density(cgig$gene_wise_JSD)
cgig_threshold <- cgig_density$x[which.max(cgig_density$y)]/4
cgig_gene_set <- cgig[cgig$gene_wise_JSD < cgig_threshold,]
write.table(cgig_gene_set,
            "~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/05-gene_sets/02-Gene_set_cgig_umbo1.txt",
            quote = F, sep = '\t', row.names = F)

spur_density <- density(spur$gene_wise_JSD)
spur_threshold <- spur_density$x[which.max(spur_density$y)]/4
spur_gene_set <- spur[spur$gene_wise_JSD < spur_threshold,]
write.table(spur_gene_set,
            "~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/05-gene_sets/03-Gene_set_spur_pluteus.txt",
            quote = F, sep = '\t', row.names = F)

nvec_density <- density(nvec$gene_wise_JSD)
nvec_threshold <- nvec_density$x[which.max(nvec_density$y)]/4
nvec_gene_set <- nvec[nvec$gene_wise_JSD < nvec_threshold,]
write.table(nvec_gene_set,
            "~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/05-gene_sets/04-Gene_set_nvec_72hpf.txt",
            quote = F, sep = '\t', row.names = F)

chem_density <- density(chem$gene_wise_JSD)
chem_threshold <- chem_density$x[which.max(chem_density$y)]/4
chem_gene_set <- chem[chem$gene_wise_JSD < chem_threshold,]
write.table(chem_gene_set,
            "~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/05-gene_sets/05-Gene_set_chem_planula.txt",
            quote = F, sep = '\t', row.names = F)


# 3. Export thresholds in a separate file
thresholds <- data.frame(cbind("species" = c("Ctel","Cgig","Spur","Nvec","Chem"),
                               "gene_wise_JSD_threshold" = c(ctel_threshold,cgig_threshold,spur_threshold,nvec_threshold,chem_threshold),
                               "gene_number" = c(length(ctel_gene_set$Target),length(cgig_gene_set$Target),
                                             length(spur_gene_set$Target),length(nvec_gene_set$Target),
                                             length(chem_gene_set$Target))))

write.table(thresholds,
            "~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/05-gene_sets/00-Gene_set_thresholds_all_species.txt",
            quote = F, sep = '\t', row.names = F)


# 4. Plot violin plots
# 4a. Full scale violin plot (0 to default)
ggplot(all, aes(x=species, y=gene_wise_JSD, fill = species)) + 
  geom_violin() +
  theme_classic() +
  scale_fill_manual(values = palette) +
  labs(y = "Gene-wise Jensen-Shannon divergence (JSD)", x = NA) +
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
ggsave("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/04-violin_plots/01-All_species_full_scale_violin.pdf",
       width = 4.5, height = 8)

# 4a. Full scale violin plot with boxplot (0 to default)
ggplot(all, aes(x=species, y=gene_wise_JSD, fill = species)) + 
  geom_violin() +
  theme_classic() +
  scale_fill_manual(values = palette) +
  labs(y = "Gene-wise Jensen-Shannon divergence (JSD)", x = NA) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2, outlier.shape = NA) +
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
ggsave("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/04-violin_plots/02-All_species_full_scale_violin_boxplot.pdf",
       width = 4.5, height = 8)

# 4c. Limited scale violin plot (0 to 0.04)
ggplot(all, aes(x=species, y=gene_wise_JSD, fill = species)) + 
  geom_violin() +
  theme_classic() +
  scale_fill_manual(values = palette) +
  labs(y = "Gene-wise Jensen-Shannon divergence (JSD)", x = NA) +
  coord_cartesian(ylim = c(0,0.04)) +
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
ggsave("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/04-violin_plots/03-All_species_0.04_scale_violin.pdf",
       width = 4.5, height = 8)

# 4d. Limited scale violin plot (0 to 0.04) with mean estimation and 
# standard deviation as confidence intervals 
ggplot(all, aes(x=species, y=gene_wise_JSD, fill = species)) + 
  geom_violin() +
  theme_classic() +
  scale_fill_manual(values = palette) +
  labs(y = "Gene-wise Jensen-Shannon divergence (JSD)", x = NA) +
  coord_cartesian(ylim = c(0,0.04)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "crossbar",
               colour = "black", width = 0.4) +
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
ggsave("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/04-violin_plots/04-All_species_0.04_scale_violin_mean_plus_minus_sd.pdf",
       width = 4.5, height = 8)


# 4e. Limited scale violin plot (0 to 0.04) with mean estimation and 
# standard deviation as confidence intervals and dashed line indicating
# the points with maximum probability density
ggplot(all, aes(x=species, y=gene_wise_JSD, fill = species)) + 
  geom_violin() +
  theme_classic() +
  scale_fill_manual(values = palette) +
  labs(y = "Gene-wise Jensen-Shannon divergence (JSD)", x = NA) +
  coord_cartesian(ylim = c(0,0.04)) +
  geom_hline(yintercept = ctel_threshold*4, linetype = "dashed", color = palette[1]) +
  geom_hline(yintercept = cgig_threshold*4, linetype = "dashed", color = palette[2]) +
  geom_hline(yintercept = spur_threshold*4, linetype = "dashed", color = palette[3]) +
  geom_hline(yintercept = nvec_threshold*4, linetype = "dashed", color = palette[4]) +
  geom_hline(yintercept = chem_threshold*4, linetype = "dashed", color = palette[5]) +
  stat_summary(fun.data = "mean_cl_boot", geom = "crossbar",
               colour = "black", width = 0.4) +
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
ggsave("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/04-violin_plots/05-All_species_0.04_scale_violin_mean_plus_minus_sd_max_prob.pdf",
       width = 4.5, height = 8)




#############################################
#### C. GO terms enrichment of gene sets ####
#############################################

# 1. Import GO universe

library(topGO)
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics")
geneID2GO <- readMappings(file = "Owenia_fusiformis_geneID2GO_topGO_onlyGOgenes.txt")
geneUniverse <- names(geneID2GO)


# 2. Re-import gene sets

ctel_gene_set <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/05-gene_sets/01-Gene_set_ctel_stage5.txt",
                            sep = '\t', header = T)
cgig_gene_set <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/05-gene_sets/02-Gene_set_cgig_umbo1.txt",
                            sep = '\t', header = T)
spur_gene_set <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/05-gene_sets/03-Gene_set_spur_pluteus.txt",
                            sep = '\t', header = T)
nvec_gene_set <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/05-gene_sets/04-Gene_set_nvec_72hpf.txt",
                            sep = '\t', header = T)
chem_gene_set <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/05-gene_sets/05-Gene_set_chem_planula.txt",
                            sep = '\t', header = T)


# 3. Prepare GO enrichment

ctel_gene_set_names <- as.character(ctel_gene_set$Subject)
ctel_gene_set_list_for_GO <- factor(as.integer(geneUniverse %in% ctel_gene_set_names))
names(ctel_gene_set_list_for_GO) <- geneUniverse

cgig_gene_set_names <- as.character(cgig_gene_set$Subject)
cgig_gene_set_list_for_GO <- factor(as.integer(geneUniverse %in% cgig_gene_set_names))
names(cgig_gene_set_list_for_GO) <- geneUniverse

spur_gene_set_names <- as.character(spur_gene_set$Subject)
spur_gene_set_list_for_GO <- factor(as.integer(geneUniverse %in% spur_gene_set_names))
names(spur_gene_set_list_for_GO) <- geneUniverse

nvec_gene_set_names <- as.character(nvec_gene_set$Subject)
nvec_gene_set_list_for_GO <- factor(as.integer(geneUniverse %in% nvec_gene_set_names))
names(nvec_gene_set_list_for_GO) <- geneUniverse

chem_gene_set_names <- as.character(chem_gene_set$Subject)
chem_gene_set_list_for_GO <- factor(as.integer(geneUniverse %in% chem_gene_set_names))
names(chem_gene_set_list_for_GO) <- geneUniverse


# 4. Run GO enrichment 

ctel_gene_set_GOdata_BP <- new("topGOdata", description="ctel_gene_set_program_BP",
                           ontology="BP", allGenes=ctel_gene_set_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) 
ctel_gene_set_resultFisher_BP <- runTest(ctel_gene_set_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_ctel_gene_set_BP <- GenTable(ctel_gene_set_GOdata_BP, classicFisher = ctel_gene_set_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30)

cgig_gene_set_GOdata_BP <- new("topGOdata", description="cgig_gene_set_program_BP",
                           ontology="BP", allGenes=cgig_gene_set_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) 
cgig_gene_set_resultFisher_BP <- runTest(cgig_gene_set_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_cgig_gene_set_BP <- GenTable(cgig_gene_set_GOdata_BP, classicFisher = cgig_gene_set_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30)

spur_gene_set_GOdata_BP <- new("topGOdata", description="spur_gene_set_program_BP",
                           ontology="BP", allGenes=spur_gene_set_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) 
spur_gene_set_resultFisher_BP <- runTest(spur_gene_set_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_spur_gene_set_BP <- GenTable(spur_gene_set_GOdata_BP, classicFisher = spur_gene_set_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30)

nvec_gene_set_GOdata_BP <- new("topGOdata", description="nvec_gene_set_program_BP",
                               ontology="BP", allGenes=nvec_gene_set_list_for_GO,
                               annot = annFUN.gene2GO, gene2GO = geneID2GO) 
nvec_gene_set_resultFisher_BP <- runTest(nvec_gene_set_GOdata_BP,
                                         algorithm="classic", statistic="fisher") 
results_nvec_gene_set_BP <- GenTable(nvec_gene_set_GOdata_BP, classicFisher = nvec_gene_set_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30)

chem_gene_set_GOdata_BP <- new("topGOdata", description="chem_gene_set_program_BP",
                               ontology="BP", allGenes=chem_gene_set_list_for_GO,
                               annot = annFUN.gene2GO, gene2GO = geneID2GO) 
chem_gene_set_resultFisher_BP <- runTest(chem_gene_set_GOdata_BP,
                                         algorithm="classic", statistic="fisher") 
results_chem_gene_set_BP <- GenTable(chem_gene_set_GOdata_BP, classicFisher = chem_gene_set_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30)


# 5. Export results

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/06-GO_terms_enrichment")

write.table(results_ctel_gene_set_BP, "01-Gene_set_ctel_stage5_GO_terms.txt", quote=FALSE, sep='\t', row.names = F)
write.table(results_cgig_gene_set_BP, "02-Gene_set_cgig_umbo1_GO_terms.txt", quote=FALSE, sep='\t', row.names = F)
write.table(results_spur_gene_set_BP, "03-Gene_set_spur_pluteus_GO_terms.txt", quote=FALSE, sep='\t', row.names = F) 
write.table(results_nvec_gene_set_BP, "04-Gene_set_nvec_72hpf_GO_terms.txt", quote=FALSE, sep='\t', row.names = F)
write.table(results_chem_gene_set_BP, "05-Gene_set_chem_planula_GO_terms.txt", quote=FALSE, sep='\t', row.names = F) 


# 6. Plot barplots of GO terms

results_ctel_gene_set_BP$classicFisher[results_ctel_gene_set_BP$classicFisher == '< 1e-30'] <- 1e-30
results_cgig_gene_set_BP$classicFisher[results_cgig_gene_set_BP$classicFisher == '< 1e-30'] <- 1e-30
results_spur_gene_set_BP$classicFisher[results_spur_gene_set_BP$classicFisher == '< 1e-30'] <- 1e-30
results_nvec_gene_set_BP$classicFisher[results_nvec_gene_set_BP$classicFisher == '< 1e-30'] <- 1e-30
results_chem_gene_set_BP$classicFisher[results_chem_gene_set_BP$classicFisher == '< 1e-30'] <- 1e-30

goEnrichment <- results_ctel_gene_set_BP
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

ctel_gene_set_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("C. teleta stage 5 larva -- O. fusiformis mitraria larva") +
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

goEnrichment <- results_cgig_gene_set_BP
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

cgig_gene_set_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("C. gigas umbo 1 -- O. fusiformis mitraria larva") +
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

goEnrichment <- results_spur_gene_set_BP
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

spur_gene_set_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("S. purpuratus pluteus -- O. fusiformis mitraria larva") +
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

goEnrichment <- results_nvec_gene_set_BP
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

nvec_gene_set_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("N. vectensis 72 hpf -- O. fusiformis mitraria larva") +
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

goEnrichment <- results_chem_gene_set_BP
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

chem_gene_set_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("C. hemisphaerica planula -- O. fusiformis mitraria larva") +
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

library(cowplot)
plot_grid(ctel_gene_set_plot, cgig_gene_set_plot, spur_gene_set_plot,
          nvec_gene_set_plot, chem_gene_set_plot,
          align = "v", ncol = 2)
ggsave("06-Gene_set_GO_terms_barplots.pdf", width = 14, height = 13)


# 7. Compare different GO terms enrichments

results_ctel_gene_set_BP <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/06-GO_terms_enrichment/01-Gene_set_ctel_stage5_GO_terms.txt",
                                       sep = '\t', header = T)
results_cgig_gene_set_BP <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/06-GO_terms_enrichment/01-Gene_set_cgig_umbo1_GO_terms.txt",
                                       sep = '\t', header = T)
results_spur_gene_set_BP <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/06-GO_terms_enrichment/01-Gene_set_spur_pluteus_GO_terms.txt",
                                       sep = '\t', header = T)
results_nvec_gene_set_BP <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/06-GO_terms_enrichment/01-Gene_set_nvec_72hpf_GO_terms.txt",
                                       sep = '\t', header = T)
results_chem_gene_set_BP <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/06-GO_terms_enrichment/01-Gene_set_chem_planula_GO_terms.txt",
                                       sep = '\t', header = T)


library(simplifyEnrichment)
library(ComplexHeatmap)
go_id <- rbind(results_ctel_gene_set_BP, results_cgig_gene_set_BP, results_spur_gene_set_BP,
               results_nvec_gene_set_BP, results_chem_gene_set_BP)
go_id$GO.ID <- as.character(go_id$GO.ID)
mat <- GO_similarity(go_id$GO.ID, ont = "BP")
# 51 non-redundant GO terms 

set.seed(2498) # ensures reproducibility
GO_cluster_annotation_raw <- simplifyGO(mat, method = "kmeans")

GO_cluster_annotation_clean <- GO_cluster_annotation_raw[2:ncol(GO_cluster_annotation_raw)]
rownames(GO_cluster_annotation_clean) <- t(GO_cluster_annotation_raw[1])

results_ctel_gene_set_BP$species_id <- "1"
results_cgig_gene_set_BP$species_id <- "2"
results_spur_gene_set_BP$species_id <- "3"
results_nvec_gene_set_BP$species_id <- "4"
results_chem_gene_set_BP$species_id <- "5"

GO_clusters <- rbind(results_ctel_gene_set_BP, results_cgig_gene_set_BP, results_spur_gene_set_BP,
                     results_nvec_gene_set_BP, results_chem_gene_set_BP)
GO_clusters$classicFisher[GO_clusters$classicFisher == '< 1e-30'] <- 1e-30
GO_clusters$classicFisher[GO_clusters$classicFisher == '<1e-30'] <- 1e-30
GO_clusters$kmeans_cluster <- GO_cluster_annotation_clean[GO_clusters$GO.ID,2] 
GO_clusters <- GO_clusters[order(GO_clusters$kmeans_cluster, 
                                 GO_clusters$classicFisher),] 

heatmap_data <- GO_clusters
heatmap_data$classicFisher <- as.numeric(heatmap_data$classicFisher)
heatmap_data <- heatmap_data[,c("species_id","GO.ID","classicFisher")] 

heatmap_data_long <- data.frame(t(spread(heatmap_data, "GO.ID","classicFisher")))
heatmap_data_long <- heatmap_data_long[2:nrow(heatmap_data_long),]
colnames(heatmap_data_long) <- c("Ctel","Cgig","Nvec","Spur","Chem")
heatmap_data_long$Ctel <- as.numeric(heatmap_data_long$Ctel)
heatmap_data_long$Cgig <- as.numeric(heatmap_data_long$Cgig)
heatmap_data_long$Nvec <- as.numeric(heatmap_data_long$Nvec)
heatmap_data_long$Spur <- as.numeric(heatmap_data_long$Spur)
heatmap_data_long$Chem <- as.numeric(heatmap_data_long$Chem)
heatmap_data_long <- -log(heatmap_data_long,10)


heatmap_data_long$kmeans_cluster <- GO_cluster_annotation_clean[rownames(heatmap_data_long),2]
heatmap_data_long$averageFisher <- rowMeans(heatmap_data_long[,c(1:5)], na.rm = TRUE)
heatmap_data_long <- heatmap_data_long[order(heatmap_data_long$kmeans_cluster,
                                             -heatmap_data_long$averageFisher),]

palette <- colorRampPalette(brewer.pal(n = 7, name = "Reds"))(200)
palette[1] <- rgb(255,255,255, maxColorValue = 255)

mat_heatmap <- data.matrix(heatmap_data_long[-c(6,7)])
mat_heatmap[is.na(mat_heatmap)] <- 0 # NAs get converted to 0s, which in turn will be shown in the heatmap in white
ComplexHeatmap::Heatmap(mat_heatmap,
                        border = TRUE,
                        row_split = heatmap_data_long$kmeans_cluster,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        show_row_names = TRUE,
                        row_labels = GO_cluster_annotation_clean[rownames(heatmap_data_long),"term"],
                        col = palette)



##########################################
#### D. Transcription factor analyses ####
##########################################

# 0. Colours and aesthetics

colors <- hue_pal()(10)
palette <- colors[c(1,2,7,8,9)]
ctel_color <- colors[1]
cgig_color <- colors[2]
spur_color <- colors[7]
nvec_color <- colors[8]
chem_color <- colors[9]

# 1. Re-import gene sets

ctel_gene_set <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/05-gene_sets/01-Gene_set_ctel_stage5.txt",
                            sep = '\t', header = T)
cgig_gene_set <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/05-gene_sets/02-Gene_set_cgig_umbo1.txt",
                            sep = '\t', header = T)
spur_gene_set <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/05-gene_sets/03-Gene_set_spur_pluteus.txt",
                            sep = '\t', header = T)
nvec_gene_set <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/05-gene_sets/04-Gene_set_nvec_72hpf.txt",
                            sep = '\t', header = T)
chem_gene_set <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/05-gene_sets/05-Gene_set_chem_planula.txt",
                            sep = '\t', header = T)


# 2. Import transcription factor list

ofus_tf_list <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/16-Transcription_factors_Owenia_Capitella_Dimorphilus/06-Owenia_fusiformis_TF_list_classified.txt", 
                           header = F, sep ='\t')
colnames(ofus_tf_list) <- c("Gene_ID", "PFAM_ID", "PFAM_annotation")


# 3. Subset TFs from gene sets and get TF annotations

ctel_TF_set <- data.frame(ctel_gene_set[ctel_gene_set$Subject %in% ofus_tf_list$Gene_ID,-c(3)])
ctel_TF_set$PFAM_ID <- ofus_tf_list[ofus_tf_list$Gene_ID %in% ctel_TF_set$Subject,"PFAM_ID"]
ctel_TF_set$PFAM_annotation <- ofus_tf_list[ofus_tf_list$Gene_ID %in% ctel_TF_set$Subject,"PFAM_annotation"]

cgig_TF_set <- data.frame(cgig_gene_set[cgig_gene_set$Subject %in% ofus_tf_list$Gene_ID,-c(3)])
cgig_TF_set$PFAM_ID <- ofus_tf_list[ofus_tf_list$Gene_ID %in% cgig_TF_set$Subject,"PFAM_ID"]
cgig_TF_set$PFAM_annotation <- ofus_tf_list[ofus_tf_list$Gene_ID %in% cgig_TF_set$Subject,"PFAM_annotation"]

spur_TF_set <- data.frame(spur_gene_set[spur_gene_set$Subject %in% ofus_tf_list$Gene_ID,-c(3)])
spur_TF_set$PFAM_ID <- ofus_tf_list[ofus_tf_list$Gene_ID %in% spur_TF_set$Subject,"PFAM_ID"]
spur_TF_set$PFAM_annotation <- ofus_tf_list[ofus_tf_list$Gene_ID %in% spur_TF_set$Subject,"PFAM_annotation"]

nvec_TF_set <- data.frame(nvec_gene_set[nvec_gene_set$Subject %in% ofus_tf_list$Gene_ID,-c(3)])
nvec_TF_set$PFAM_ID <- ofus_tf_list[ofus_tf_list$Gene_ID %in% nvec_TF_set$Subject,"PFAM_ID"]
nvec_TF_set$PFAM_annotation <- ofus_tf_list[ofus_tf_list$Gene_ID %in% nvec_TF_set$Subject,"PFAM_annotation"]

chem_TF_set <- data.frame(chem_gene_set[chem_gene_set$Subject %in% ofus_tf_list$Gene_ID,-c(3)])
chem_TF_set$PFAM_ID <- ofus_tf_list[ofus_tf_list$Gene_ID %in% chem_TF_set$Subject,"PFAM_ID"]
chem_TF_set$PFAM_annotation <- ofus_tf_list[ofus_tf_list$Gene_ID %in% chem_TF_set$Subject,"PFAM_annotation"]


# 4. Export the files

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/28-Gene_specific_metazoan_comparative_transcriptomics/05-gene_sets")
write.table(ctel_TF_set, "06-Gene_set_ctel_stage5_TF_only.txt", sep = '\t', quote = F, row.names = F)
write.table(cgig_TF_set, "07-Gene_set_cgig_umbo1_TF_only.txt", sep = '\t', quote = F, row.names = F)
write.table(spur_TF_set, "08-Gene_set_spur_pluteus_TF_only.txt", sep = '\t', quote = F, row.names = F)
write.table(nvec_TF_set, "09-Gene_set_nvec_72hpf_TF_only.txt", sep = '\t', quote = F, row.names = F)
write.table(chem_TF_set, "10-Gene_set_chem_planula_TF_only.txt", sep = '\t', quote = F, row.names = F)


# 5. Explore potential commonalities and dissimilarities?

library(ggvenn)

venn_data <- list("Cgig" = cgig_TF_set$Subject, "Spur" = spur_TF_set$Subject,
                  "Nvec" = nvec_TF_set$Subject, "Chem" = chem_TF_set$Subject)

ggvenn(venn_data,
       fill_color = c(palette[-c(1)]),
       stroke_size = 0.5, set_name_size = 4)
ggsave("11-Gene_set_TF_only_Venn_diagram.pdf", width = 5, height = 5)

venn_data2 <- list("Ctel"= ctel_TF_set$Subject,
                   "Cgig" = cgig_TF_set$Subject, "Spur" = spur_TF_set$Subject)

ggvenn(venn_data2,
       fill_color = c(palette[-c(4,5)]),
       stroke_size = 0.5, set_name_size = 4)
ggsave("12-Gene_set_TF_only_Venn_diagram_Bilateria_only.pdf", width = 5, height = 5)





