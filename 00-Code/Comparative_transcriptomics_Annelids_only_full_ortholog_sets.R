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


###############################
## Annelid only comparisons ###
###############################


# 1. Owenia fusiformis vs. Capitella teleta
# Number of orthologues = 7,651 orthologues
# 14 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/30-Annelid_comparative_transcriptomics_JSD/Quantile_transformation")
Ofus2Ctel_orth_tpm_qn <- read.csv("Ofus2Ctel_TPM_mean_quantile_transform.csv", header = T)
Ctel2Ofus_orth_tpm_qn <- read.csv("Ctel2Ofus_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/30-Annelid_comparative_transcriptomics_JSD/one2one")
both_ID <- read.table("Ofus2Ctel.txt", header = T)
match_Ofus <- Ofus2Ctel_orth_tpm_qn[match(both_ID$Ofus,Ofus2Ctel_orth_tpm_qn$Gene_ID),]
match_Ctel <- Ctel2Ofus_orth_tpm_qn[match(both_ID$Ctel,Ctel2Ofus_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/30-Annelid_comparative_transcriptomics_JSD/order")
write.table(match_Ofus,"order_Ofus2Ctel_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Ctel,"order_Ctel2Ofus_tpm_qn.txt", sep="\t", quote=F)
match_Ofus <- read.table("order_Ofus2Ctel_tpm_qn.txt", header = T)
match_Ctel <- read.table("order_Ctel2Ofus_tpm_qn.txt", header = T)

Ofus2Ctel_full_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Ctel)-1,rep=TRUE)))
colnames(Ofus2Ctel_full_JSD) <- colnames(Ofus2Ctel_orth_tpm_qn)[-c(1)]
rownames(Ofus2Ctel_full_JSD) <- colnames(Ctel2Ofus_orth_tpm_qn)[-c(1)]
Ofus2Ctel_mean_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Ctel)-1,rep=TRUE)))
colnames(Ofus2Ctel_mean_JSD) <- colnames(Ofus2Ctel_orth_tpm_qn)[-c(1)]
rownames(Ofus2Ctel_mean_JSD) <- colnames(Ctel2Ofus_orth_tpm_qn)[-c(1)]
Ofus2Ctel_sd_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Ctel)-1,rep=TRUE)))
colnames(Ofus2Ctel_sd_JSD) <- colnames(Ofus2Ctel_orth_tpm_qn)[-c(1)]
rownames(Ofus2Ctel_sd_JSD) <- colnames(Ctel2Ofus_orth_tpm_qn)[-c(1)]

for (i in colnames(Ofus2Ctel_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Ctel2Ofus_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Ofus[[i]],match_Ctel[[j]])
    Ofus2Ctel_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Ofus2Ctel_mean_JSD[j,i] <- mean(all_JSD)
    Ofus2Ctel_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/30-Annelid_comparative_transcriptomics_JSD/Results")
write.table(Ofus2Ctel_full_JSD, "Ofus2Ctel_full_set_JSD.txt", sep ='\t')
write.table(Ofus2Ctel_mean_JSD, "Ofus2Ctel_subset_JSD_mean.txt", sep ='\t')
write.table(Ofus2Ctel_sd_JSD, "Ofus2Ctel_subset_JSD_sd.txt", sep ='\t')


# 2. Owenia fusiformis vs. Dimorphilus gyrociliatus
# Number of orthologues = 6,346 orthologues
# 4 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/30-Annelid_comparative_transcriptomics_JSD/Quantile_transformation")
Ofus2Dgyr_orth_tpm_qn <- read.csv("Ofus2Dgyr_TPM_mean_quantile_transform.csv", header = T)
Dgyr2Ofus_orth_tpm_qn <- read.csv("Dgyr2Ofus_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/30-Annelid_comparative_transcriptomics_JSD/one2one")
both_ID <- read.table("Dgyr2Ofus.txt", header = T)
match_Ofus <- Ofus2Dgyr_orth_tpm_qn[match(both_ID$Ofus,Ofus2Dgyr_orth_tpm_qn$Gene_ID),]
match_Dgyr <- Dgyr2Ofus_orth_tpm_qn[match(both_ID$Dgyr,Dgyr2Ofus_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/30-Annelid_comparative_transcriptomics_JSD/order")
write.table(match_Ofus,"order_Ofus2Dgyr_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Dgyr,"order_Dgyr2Ofus_tpm_qn.txt", sep="\t", quote=F)
match_Ofus <- read.table("order_Ofus2Dgyr_tpm_qn.txt", header = T)
match_Dgyr <- read.table("order_Dgyr2Ofus_tpm_qn.txt", header = T)

Ofus2Dgyr_full_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Dgyr)-1,rep=TRUE)))
colnames(Ofus2Dgyr_full_JSD) <- colnames(Ofus2Dgyr_orth_tpm_qn)[-c(1)]
rownames(Ofus2Dgyr_full_JSD) <- colnames(Dgyr2Ofus_orth_tpm_qn)[-c(1)]
Ofus2Dgyr_mean_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Dgyr)-1,rep=TRUE)))
colnames(Ofus2Dgyr_mean_JSD) <- colnames(Ofus2Dgyr_orth_tpm_qn)[-c(1)]
rownames(Ofus2Dgyr_mean_JSD) <- colnames(Dgyr2Ofus_orth_tpm_qn)[-c(1)]
Ofus2Dgyr_sd_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Dgyr)-1,rep=TRUE)))
colnames(Ofus2Dgyr_sd_JSD) <- colnames(Ofus2Dgyr_orth_tpm_qn)[-c(1)]
rownames(Ofus2Dgyr_sd_JSD) <- colnames(Dgyr2Ofus_orth_tpm_qn)[-c(1)]

for (i in colnames(Ofus2Dgyr_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Dgyr2Ofus_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Ofus[[i]],match_Dgyr[[j]])
    Ofus2Dgyr_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Ofus2Dgyr_mean_JSD[j,i] <- mean(all_JSD)
    Ofus2Dgyr_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/30-Annelid_comparative_transcriptomics_JSD/Results")
write.table(Ofus2Dgyr_full_JSD, "Ofus2Dgyr_full_set_JSD.txt", sep ='\t')
write.table(Ofus2Dgyr_mean_JSD, "Ofus2Dgyr_subset_JSD_mean.txt", sep ='\t')
write.table(Ofus2Dgyr_sd_JSD, "Ofus2Dgyr_subset_JSD_sd.txt", sep ='\t')


# 3. Capitella teleta vs. Dimorphilus gyrociliatus
# Number of orthologues = 6,236 orthologues
# 4 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/30-Annelid_comparative_transcriptomics_JSD/Quantile_transformation")
Ctel2Dgyr_orth_tpm_qn <- read.csv("Ctel2Dgyr_TPM_mean_quantile_transform.csv", header = T)
Dgyr2Ctel_orth_tpm_qn <- read.csv("Dgyr2Ctel_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/30-Annelid_comparative_transcriptomics_JSD/one2one")
both_ID <- read.table("Dgyr2Ctel.txt", header = T)
match_Ctel <- Ctel2Dgyr_orth_tpm_qn[match(both_ID$Ctel,Ctel2Dgyr_orth_tpm_qn$Gene_ID),]
match_Dgyr <- Dgyr2Ctel_orth_tpm_qn[match(both_ID$Dgyr,Dgyr2Ctel_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/30-Annelid_comparative_transcriptomics_JSD/order")
write.table(match_Ctel,"order_Ctel2Dgyr_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Dgyr,"order_Dgyr2Ctel_tpm_qn.txt", sep="\t", quote=F)
match_Ctel <- read.table("order_Ctel2Dgyr_tpm_qn.txt", header = T)
match_Dgyr <- read.table("order_Dgyr2Ctel_tpm_qn.txt", header = T)

Ctel2Dgyr_full_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Dgyr)-1,rep=TRUE)))
colnames(Ctel2Dgyr_full_JSD) <- colnames(Ctel2Dgyr_orth_tpm_qn)[-c(1)]
rownames(Ctel2Dgyr_full_JSD) <- colnames(Dgyr2Ctel_orth_tpm_qn)[-c(1)]
Ctel2Dgyr_mean_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Dgyr)-1,rep=TRUE)))
colnames(Ctel2Dgyr_mean_JSD) <- colnames(Ctel2Dgyr_orth_tpm_qn)[-c(1)]
rownames(Ctel2Dgyr_mean_JSD) <- colnames(Dgyr2Ctel_orth_tpm_qn)[-c(1)]
Ctel2Dgyr_sd_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Dgyr)-1,rep=TRUE)))
colnames(Ctel2Dgyr_sd_JSD) <- colnames(Ctel2Dgyr_orth_tpm_qn)[-c(1)]
rownames(Ctel2Dgyr_sd_JSD) <- colnames(Dgyr2Ctel_orth_tpm_qn)[-c(1)]

for (i in colnames(Ctel2Dgyr_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Dgyr2Ctel_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Ctel[[i]],match_Dgyr[[j]])
    Ctel2Dgyr_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Ctel2Dgyr_mean_JSD[j,i] <- mean(all_JSD)
    Ctel2Dgyr_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/30-Annelid_comparative_transcriptomics_JSD/Results")
write.table(Ctel2Dgyr_full_JSD, "Ctel2Dgyr_full_set_JSD.txt", sep ='\t')
write.table(Ctel2Dgyr_mean_JSD, "Ctel2Dgyr_subset_JSD_mean.txt", sep ='\t')
write.table(Ctel2Dgyr_sd_JSD, "Ctel2Dgyr_subset_JSD_sd.txt", sep ='\t')



#######################################
## Annelid only 1-to-all comparisons ##
#######################################

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/30-Annelid_comparative_transcriptomics_JSD/Results")

Ofus2Ctel_mean <- read.table("Ofus2Ctel_subset_JSD_mean.txt", header = T)[-c(1:7),-c(1:7)]
Ofus2Dgyr_mean <- read.table("Ofus2Dgyr_subset_JSD_mean.txt", header = T)[,-c(1:7)]
Ctel2Dgyr_mean <- read.table("Ctel2Dgyr_subset_JSD_mean.txt", header = T)[,-c(1:7)]

Ofus2Ctel_sd <- read.table("Ofus2Ctel_subset_JSD_sd.txt", header = T)[-c(1:7),-c(1:7)]
Ofus2Dgyr_sd <- read.table("Ofus2Dgyr_subset_JSD_sd.txt", header = T)[,-c(1:7)]
Ctel2Dgyr_sd <- read.table("Ctel2Dgyr_subset_JSD_sd.txt", header = T)[,-c(1:7)]

Ofus2Ctel_transposed <- data.frame(t(Ofus2Ctel_mean))
Ofus2Dgyr_transposed <- data.frame(t(Ofus2Dgyr_mean))
Ctel2Dgyr_transposed <- data.frame(t(Ctel2Dgyr_mean))

min_stages <- data.frame(replicate(3,sample(0:1,ncol(Ofus2Dgyr_mean),rep=TRUE)))
colnames(min_stages) <- c("Capitella_teleta","Dimorphilus_gyrociliatus_Ofus","Dimorphilus_gyrociliatus_Ctel")
rownames(min_stages) <- c("blastula/64 cells","gastrula/gastrula","elongation/stage 4tt larva",
                          "early larva/stage 5 larva","mitraria larva/stage 7 larva",
                          "competent larva/pre-competent larva", "juvenile/competent larva")
min_stages$Capitella_teleta <- names(Ofus2Ctel_transposed)[apply(Ofus2Ctel_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Dimorphilus_gyrociliatus_Ofus <- names(Ofus2Dgyr_transposed)[apply(Ofus2Dgyr_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Dimorphilus_gyrociliatus_Ctel <- names(Ctel2Dgyr_transposed)[apply(Ctel2Dgyr_transposed, MARGIN = 1, FUN = which.min)]

minimum_distance_mean <- data.frame(replicate(3,sample(0:1,ncol(Ofus2Dgyr_mean),rep=TRUE)))
minimum_distance_sd <- data.frame(replicate(3,sample(0:1,ncol(Ofus2Dgyr_mean),rep=TRUE)))
colnames(minimum_distance_mean) <- c("Capitella_teleta","Dimorphilus_gyrociliatus_Ofus","Dimorphilus_gyrociliatus_Ctel")
rownames(minimum_distance_mean) <- c("blastula/64 cells","gastrula/gastrula","elongation/stage 4tt larva",
                                     "early larva/stage 5 larva","mitraria larva/stage 7 larva",
                                     "competent larva/pre-competent larva", "juvenile/competent larva")
colnames(minimum_distance_sd) <- c("Capitella_teleta","Dimorphilus_gyrociliatus_Ofus","Dimorphilus_gyrociliatus_Ctel")
rownames(minimum_distance_sd) <- c("blastula/64 cells","gastrula/gastrula","elongation/stage 4tt larva",
                                   "early larva/stage 5 larva","mitraria larva/stage 7 larva",
                                   "competent larva/pre-competent larva", "juvenile/competent larva")

for (i in c(1:ncol(Ofus2Dgyr_mean))){
  minimum_distance_mean[i,"Capitella_teleta"] <- Ofus2Ctel_mean[min_stages$Capitella_teleta[i],i]
  minimum_distance_mean[i,"Dimorphilus_gyrociliatus_Ofus"] <- Ofus2Dgyr_mean[min_stages$Dimorphilus_gyrociliatus_Ofus[i],i]
  minimum_distance_mean[i,"Dimorphilus_gyrociliatus_Ctel"] <- Ctel2Dgyr_mean[min_stages$Dimorphilus_gyrociliatus_Ctel[i],i]
  minimum_distance_sd[i,"Capitella_teleta"] <- Ofus2Ctel_sd[min_stages$Capitella_teleta[i],i]
  minimum_distance_sd[i,"Dimorphilus_gyrociliatus_Ofus"] <- Ofus2Dgyr_sd[min_stages$Dimorphilus_gyrociliatus_Ofus[i],i]
  minimum_distance_sd[i,"Dimorphilus_gyrociliatus_Ctel"] <- Ctel2Dgyr_sd[min_stages$Dimorphilus_gyrociliatus_Ctel[i],i]
}

minimum_distance_mean$stage <- c(1:ncol(Ofus2Dgyr_mean))
minimum_distance_sd$stage <- c(1:ncol(Ofus2Dgyr_mean))
minimum_distance_mean_tidy <- gather(minimum_distance_mean, "species", "JSD", -stage)
minimum_distance_sd_tidy <- gather(minimum_distance_sd, "species", "JSD", -stage)
minimum_distance_mean_tidy$upper <- minimum_distance_mean_tidy["JSD"] + minimum_distance_sd_tidy["JSD"]
minimum_distance_mean_tidy$lower <- minimum_distance_mean_tidy["JSD"] - minimum_distance_sd_tidy["JSD"]
colnames(minimum_distance_mean_tidy)[4:5] <- c("upper","lower")
minimum_distance_mean_tidy$upper <- unlist(minimum_distance_mean_tidy$upper)
minimum_distance_mean_tidy$lower <- unlist(minimum_distance_mean_tidy$lower)

minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Capitella_teleta'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Capitella_teleta'),3:5]-min(minimum_distance_mean$Capitella_teleta))/(max(minimum_distance_mean$Capitella_teleta)-min(minimum_distance_mean$Capitella_teleta))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Dimorphilus_gyrociliatus_Ofus'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Dimorphilus_gyrociliatus_Ofus'),3:5]-min(minimum_distance_mean$Dimorphilus_gyrociliatus_Ofus))/(max(minimum_distance_mean$Dimorphilus_gyrociliatus_Ofus)-min(minimum_distance_mean$Dimorphilus_gyrociliatus_Ofus))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Dimorphilus_gyrociliatus_Ctel'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Dimorphilus_gyrociliatus_Ctel'),3:5]-min(minimum_distance_mean$Dimorphilus_gyrociliatus_Ctel))/(max(minimum_distance_mean$Dimorphilus_gyrociliatus_Ctel)-min(minimum_distance_mean$Dimorphilus_gyrociliatus_Ctel))

final_dataset <- minimum_distance_mean_tidy
final_dataset$species <- factor(final_dataset$species,
                                levels = c("Capitella_teleta","Dimorphilus_gyrociliatus_Ofus","Dimorphilus_gyrociliatus_Ctel"))

ggplot(final_dataset, aes(x=stage, y=JSD, colour=species)) + 
  geom_line(show.legend = FALSE) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=factor(species)), 
              colour = NA, show.legend = FALSE,
              alpha = 0.5) +
  facet_wrap(~species, nrow=3, ncol=1) +
  scale_x_continuous(labels = c(1:ncol(Ofus2Dgyr_mean)), breaks = c(1:ncol(Ofus2Dgyr_mean))) +
  theme_classic() +
  labs(x = "O. fusiformis/C. teleta stage", y = "Relative gene expression divergence (JSD)")

write.table(min_stages, "minimum_stages.txt", sep='\t', quote = FALSE)
write.table(minimum_distance_mean, "minimum_stages_distance_mean.txt", sep='\t', quote = FALSE)
write.table(minimum_distance_sd, "minimum_stages_distance_sd.txt", sep='\t', quote = FALSE)
write.table(minimum_distance_mean_tidy, "minimum_stages_to_plot.txt", sep='\t', quote = FALSE)



#########################################
## Annelid only normalised comparisons ##
#########################################

# Normalise results by the number of orthologues to make 
# the 1-to-1 comparisons comparable: heatmaps of NORMALISED data

# Number of orthologues in order: Ofus2Ctel, Ofus2Dgyr, Ofus2Ctel
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/30-Annelid_comparative_transcriptomics_JSD/Results")
Orthologues <- c(7651,6346,6236)

Ofus2Ctel_mean <- read.table("Ofus2Ctel_subset_JSD_mean.txt", header = T)[-c(1:7),-c(1:7)]
Ofus2Dgyr_mean <- read.table("Ofus2Dgyr_subset_JSD_mean.txt", header = T)[-c(1:7)]
Ctel2Dgyr_mean <- read.table("Ctel2Dgyr_subset_JSD_mean.txt", header = T)[-c(1:7)]

Ofus2Ctel_mean_divided <- Ofus2Ctel_mean/Orthologues[1]
Ofus2Dgyr_mean_divided <- Ofus2Dgyr_mean/Orthologues[2]
Ctel2Dgyr_mean_divided <- Ctel2Dgyr_mean/Orthologues[3]

all_values_divid1 <- rbind(Ofus2Ctel_mean_divided,Ofus2Dgyr_mean_divided)
all_values_divid2 <- rbind(Ctel2Dgyr_mean_divided)

maxima <- vector(length=2)
minima <- vector(length=2)

minima[1] <- min(all_values_divid1)
minima[2] <- min(all_values_divid2)
maxima[1] <- max(all_values_divid1)
maxima[2] <- max(all_values_divid2)

minimum_divided <- min(minima)
maximum_divided <- max(maxima)

Ofus2Ctel_final_normalised <- (Ofus2Ctel_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Ofus2Dgyr_final_normalised <- (Ofus2Dgyr_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Ctel2Dgyr_final_normalised <- (Ctel2Dgyr_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/30-Annelid_comparative_transcriptomics_JSD/Results")

write.table(Ofus2Ctel_final_normalised, paste0("Ofus2Ctel_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Ofus2Dgyr_final_normalised, paste0("Ofus2Dgyr_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Ctel2Dgyr_final_normalised, paste0("Ctel2Dgyr_final_normalised",".txt"), sep='\t', quote = FALSE)


## plots

heatmap_color <- colorRampPalette(brewer.pal(n = 7, name = "RdGy"))(100)
paletteLength <- 100
myBreaks <- c(seq(1/paletteLength, 1, length.out=floor(paletteLength)))
heatmap_color <- colorRampPalette(brewer.pal(n = 7, name = "RdGy"))(100)

pheatmap(Ofus2Ctel_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks,
         labels_col = c("blastula","gastrula","elongation","early larva",
                        "mitraria larva", "competent larva", "juvenile"),
         labels_row = c("64 cells", "stage 4tt larva", "stage 5 larva",
                        "stage 7 larva", "pre-competent larva",
                        "competent larva"))

pheatmap(Ofus2Dgyr_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks,
         labels_col = c("blastula","gastrula","elongation","early larva",
                        "mitraria larva", "competent larva", "juvenile"),
         labels_row = c("early development","late development","hatchling",
                        "female adult"))

pheatmap(Ctel2Dgyr_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks,
         labels_col = c("64 cells", "stage 4tt larva", "stage 5 larva",
                        "stage 7 larva", "pre-competent larva",
                        "competent larva"),
         labels_row = c("early development","late development","hatchling",
                        "female adult"))


## plots raw

heatmap_color <- colorRampPalette(brewer.pal(n = 7, name = "RdGy"))(100)
paletteLength <- 100
myBreaks <- c(seq(1/paletteLength, 1, length.out=floor(paletteLength)))
heatmap_color <- colorRampPalette(brewer.pal(n = 7, name = "RdGy"))(100)

pheatmap(Ofus2Ctel_mean, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color,
         labels_col = c("blastula","gastrula","elongation","early larva",
                        "mitraria larva", "competent larva", "juvenile"),
         labels_row = c("64 cells", "stage 4tt larva", "stage 5 larva",
                        "stage 7 larva", "pre-competent larva",
                        "competent larva"))

pheatmap(Ofus2Dgyr_mean, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color,
         labels_col = c("blastula","gastrula","elongation","early larva",
                        "mitraria larva", "competent larva", "juvenile"),
         labels_row = c("early development","late development","hatchling",
                        "female adult"))

pheatmap(Ctel2Dgyr_mean, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color,
         labels_col = c("64 cells", "stage 4tt larva", "stage 5 larva",
                        "stage 7 larva", "pre-competent larva",
                        "competent larva"),
         labels_row = c("early development","late development","hatchling",
                        "female adult"))



