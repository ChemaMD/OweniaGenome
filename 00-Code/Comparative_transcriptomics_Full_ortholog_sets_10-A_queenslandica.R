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


##########################################
## Amphimedon queenslandica 1-to-1 comparisons ##
##########################################

# 8. Amphimedon queenslandica vs. Clytica hemisphaerica
# Number of orthologues = 4,691 orthologues
# 9 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/10-A_queenslandica/Quantile_transformation")
Aque2Chem_orth_tpm_qn <- read.csv("Aque2Chem_TPM_mean_quantile_transform.csv", header = T)
Chem2Aque_orth_tpm_qn <- read.csv("Chem2Aque_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/10-A_queenslandica/one2one")
both_ID <- read.table("Aque2Chem.txt", header = T)
match_Aque <- Aque2Chem_orth_tpm_qn[match(both_ID$Aque,Aque2Chem_orth_tpm_qn$Gene_ID),]
match_Aque <- na.omit(match_Aque)
match_Chem <- Chem2Aque_orth_tpm_qn[match(both_ID$Chem,Chem2Aque_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/10-A_queenslandica/order")
write.table(match_Aque,"order_Aque2Chem_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Chem,"order_Chem2Aque_tpm_qn.txt", sep="\t", quote=F)
match_Aque <- read.table("order_Aque2Chem_tpm_qn.txt", header = T)
match_Chem <- read.table("order_Chem2Aque_tpm_qn.txt", header = T)

Aque2Chem_full_JSD <- data.frame(replicate(ncol(match_Aque)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Aque2Chem_full_JSD) <- colnames(Aque2Chem_orth_tpm_qn)[-c(1)]
rownames(Aque2Chem_full_JSD) <- colnames(Chem2Aque_orth_tpm_qn)[-c(1)]
Aque2Chem_mean_JSD <- data.frame(replicate(ncol(match_Aque)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Aque2Chem_mean_JSD) <- colnames(Aque2Chem_orth_tpm_qn)[-c(1)]
rownames(Aque2Chem_mean_JSD) <- colnames(Chem2Aque_orth_tpm_qn)[-c(1)]
Aque2Chem_sd_JSD <- data.frame(replicate(ncol(match_Aque)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Aque2Chem_sd_JSD) <- colnames(Aque2Chem_orth_tpm_qn)[-c(1)]
rownames(Aque2Chem_sd_JSD) <- colnames(Chem2Aque_orth_tpm_qn)[-c(1)]

for (i in colnames(Aque2Chem_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Chem2Aque_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Aque[[i]],match_Chem[[j]])
    Aque2Chem_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Aque2Chem_mean_JSD[j,i] <- mean(all_JSD)
    Aque2Chem_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/10-A_queenslandica/Results")
write.table(Aque2Chem_full_JSD, "Aque2Chem_full_set_JSD.txt", sep ='\t')
write.table(Aque2Chem_mean_JSD, "Aque2Chem_subset_JSD_mean.txt", sep ='\t')
write.table(Aque2Chem_sd_JSD, "Aque2Chem_subset_JSD_sd.txt", sep ='\t')