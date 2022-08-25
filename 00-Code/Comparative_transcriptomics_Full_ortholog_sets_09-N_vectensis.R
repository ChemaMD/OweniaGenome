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
## Nematostella vectensis 1-to-1 comparisons ##
##########################################

# 7. Nematostella vectensis vs. Amphimedon queenslandica
# Number of orthologues = 3,962 orthologues
# 9 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/09-N_vectensis/Quantile_transformation")
Nvec2Aque_orth_tpm_qn <- read.csv("Nvec2Aque_TPM_mean_quantile_transform.csv", header = T)
Aque2Nvec_orth_tpm_qn <- read.csv("Aque2Nvec_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/09-N_vectensis/one2one")
both_ID <- read.table("Nvec2Aque.txt", header = T)
match_Nvec <- Nvec2Aque_orth_tpm_qn[match(both_ID$Nvec,Nvec2Aque_orth_tpm_qn$Gene_ID),]
match_Nvec <- na.omit(match_Nvec)
match_Aque <- Aque2Nvec_orth_tpm_qn[match(both_ID$Aque,Aque2Nvec_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/09-N_vectensis/order")
write.table(match_Nvec,"order_Nvec2Aque_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Aque,"order_Aque2Nvec_tpm_qn.txt", sep="\t", quote=F)
match_Nvec <- read.table("order_Nvec2Aque_tpm_qn.txt", header = T)
match_Aque <- read.table("order_Aque2Nvec_tpm_qn.txt", header = T)

Nvec2Aque_full_JSD <- data.frame(replicate(ncol(match_Nvec)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Nvec2Aque_full_JSD) <- colnames(Nvec2Aque_orth_tpm_qn)[-c(1)]
rownames(Nvec2Aque_full_JSD) <- colnames(Aque2Nvec_orth_tpm_qn)[-c(1)]
Nvec2Aque_mean_JSD <- data.frame(replicate(ncol(match_Nvec)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Nvec2Aque_mean_JSD) <- colnames(Nvec2Aque_orth_tpm_qn)[-c(1)]
rownames(Nvec2Aque_mean_JSD) <- colnames(Aque2Nvec_orth_tpm_qn)[-c(1)]
Nvec2Aque_sd_JSD <- data.frame(replicate(ncol(match_Nvec)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Nvec2Aque_sd_JSD) <- colnames(Nvec2Aque_orth_tpm_qn)[-c(1)]
rownames(Nvec2Aque_sd_JSD) <- colnames(Aque2Nvec_orth_tpm_qn)[-c(1)]

for (i in colnames(Nvec2Aque_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Aque2Nvec_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Nvec[[i]],match_Aque[[j]])
    Nvec2Aque_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Nvec2Aque_mean_JSD[j,i] <- mean(all_JSD)
    Nvec2Aque_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/09-N_vectensis/Results")
write.table(Nvec2Aque_full_JSD, "Nvec2Aque_full_set_JSD.txt", sep ='\t')
write.table(Nvec2Aque_mean_JSD, "Nvec2Aque_subset_JSD_mean.txt", sep ='\t')
write.table(Nvec2Aque_sd_JSD, "Nvec2Aque_subset_JSD_sd.txt", sep ='\t')


# 8. Nematostella vectensis vs. Clytica hemisphaerica
# Number of orthologues = 4,691 orthologues
# 9 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/09-N_vectensis/Quantile_transformation")
Nvec2Chem_orth_tpm_qn <- read.csv("Nvec2Chem_TPM_mean_quantile_transform.csv", header = T)
Chem2Nvec_orth_tpm_qn <- read.csv("Chem2Nvec_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/09-N_vectensis/one2one")
both_ID <- read.table("Nvec2Chem.txt", header = T)
match_Nvec <- Nvec2Chem_orth_tpm_qn[match(both_ID$Nvec,Nvec2Chem_orth_tpm_qn$Gene_ID),]
match_Nvec <- na.omit(match_Nvec)
match_Chem <- Chem2Nvec_orth_tpm_qn[match(both_ID$Chem,Chem2Nvec_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/09-N_vectensis/order")
write.table(match_Nvec,"order_Nvec2Chem_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Chem,"order_Chem2Nvec_tpm_qn.txt", sep="\t", quote=F)
match_Nvec <- read.table("order_Nvec2Chem_tpm_qn.txt", header = T)
match_Chem <- read.table("order_Chem2Nvec_tpm_qn.txt", header = T)

Nvec2Chem_full_JSD <- data.frame(replicate(ncol(match_Nvec)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Nvec2Chem_full_JSD) <- colnames(Nvec2Chem_orth_tpm_qn)[-c(1)]
rownames(Nvec2Chem_full_JSD) <- colnames(Chem2Nvec_orth_tpm_qn)[-c(1)]
Nvec2Chem_mean_JSD <- data.frame(replicate(ncol(match_Nvec)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Nvec2Chem_mean_JSD) <- colnames(Nvec2Chem_orth_tpm_qn)[-c(1)]
rownames(Nvec2Chem_mean_JSD) <- colnames(Chem2Nvec_orth_tpm_qn)[-c(1)]
Nvec2Chem_sd_JSD <- data.frame(replicate(ncol(match_Nvec)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Nvec2Chem_sd_JSD) <- colnames(Nvec2Chem_orth_tpm_qn)[-c(1)]
rownames(Nvec2Chem_sd_JSD) <- colnames(Chem2Nvec_orth_tpm_qn)[-c(1)]

for (i in colnames(Nvec2Chem_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Chem2Nvec_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Nvec[[i]],match_Chem[[j]])
    Nvec2Chem_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Nvec2Chem_mean_JSD[j,i] <- mean(all_JSD)
    Nvec2Chem_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/09-N_vectensis/Results")
write.table(Nvec2Chem_full_JSD, "Nvec2Chem_full_set_JSD.txt", sep ='\t')
write.table(Nvec2Chem_mean_JSD, "Nvec2Chem_subset_JSD_mean.txt", sep ='\t')
write.table(Nvec2Chem_sd_JSD, "Nvec2Chem_subset_JSD_sd.txt", sep ='\t')