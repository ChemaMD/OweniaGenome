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

`
##########################################
## Strongylocentrotus purpuratus 1-to-1 comparisons ##
##########################################

# 1. Strongylocentrotus purpuratus vs. Nematostella vectensis
# Number of orthologues = 5,254 orthologues (-4)
# 8 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/08-S_purpuratus/Quantile_transformation")
Spur2Nvec_orth_tpm_qn <- read.csv("Spur2Nvec_TPM_mean_quantile_transform.csv", header = T)
Nvec2Spur_orth_tpm_qn <- read.csv("Nvec2Spur_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/08-S_purpuratus/one2one")
both_ID <- read.table("Spur2Nvec.txt", header = T)
match_Spur <- Spur2Nvec_orth_tpm_qn[match(both_ID$Spur,Spur2Nvec_orth_tpm_qn$Gene_ID),]
match_Nvec <- Nvec2Spur_orth_tpm_qn[match(both_ID$Nvec,Nvec2Spur_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/08-S_purpuratus/order")
write.table(match_Spur,"order_Spur2Nvec_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Nvec,"order_Nvec2Spur_tpm_qn.txt", sep="\t", quote=F)
match_Spur <- read.table("order_Spur2Nvec_tpm_qn.txt", header = T)
match_Nvec <- read.table("order_Nvec2Spur_tpm_qn.txt", header = T)

Spur2Nvec_full_JSD <- data.frame(replicate(ncol(match_Spur)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Spur2Nvec_full_JSD) <- colnames(Spur2Nvec_orth_tpm_qn)[-c(1)]
rownames(Spur2Nvec_full_JSD) <- colnames(Nvec2Spur_orth_tpm_qn)[-c(1)]
Spur2Nvec_mean_JSD <- data.frame(replicate(ncol(match_Spur)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Spur2Nvec_mean_JSD) <- colnames(Spur2Nvec_orth_tpm_qn)[-c(1)]
rownames(Spur2Nvec_mean_JSD) <- colnames(Nvec2Spur_orth_tpm_qn)[-c(1)]
Spur2Nvec_sd_JSD <- data.frame(replicate(ncol(match_Spur)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Spur2Nvec_sd_JSD) <- colnames(Spur2Nvec_orth_tpm_qn)[-c(1)]
rownames(Spur2Nvec_sd_JSD) <- colnames(Nvec2Spur_orth_tpm_qn)[-c(1)]

for (i in colnames(Spur2Nvec_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Nvec2Spur_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Spur[[i]],match_Nvec[[j]])
    Spur2Nvec_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Spur2Nvec_mean_JSD[j,i] <- mean(all_JSD)
    Spur2Nvec_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/08-S_purpuratus/Results")
write.table(Spur2Nvec_full_JSD, "Spur2Nvec_full_set_JSD.txt", sep ='\t')
write.table(Spur2Nvec_mean_JSD, "Spur2Nvec_subset_JSD_mean.txt", sep ='\t')
write.table(Spur2Nvec_sd_JSD, "Spur2Nvec_subset_JSD_sd.txt", sep ='\t')


# 7. Strongylocentrotus purpuratus vs. Amphimedon queenslandica
# Number of orthologues = 3,962 orthologues
# 9 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/08-S_purpuratus/Quantile_transformation")
Spur2Aque_orth_tpm_qn <- read.csv("Spur2Aque_TPM_mean_quantile_transform.csv", header = T)
Aque2Spur_orth_tpm_qn <- read.csv("Aque2Spur_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/08-S_purpuratus/one2one")
both_ID <- read.table("Spur2Aque.txt", header = T)
match_Spur <- Spur2Aque_orth_tpm_qn[match(both_ID$Spur,Spur2Aque_orth_tpm_qn$Gene_ID),]
match_Spur <- na.omit(match_Spur)
match_Aque <- Aque2Spur_orth_tpm_qn[match(both_ID$Aque,Aque2Spur_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/08-S_purpuratus/order")
write.table(match_Spur,"order_Spur2Aque_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Aque,"order_Aque2Spur_tpm_qn.txt", sep="\t", quote=F)
match_Spur <- read.table("order_Spur2Aque_tpm_qn.txt", header = T)
match_Aque <- read.table("order_Aque2Spur_tpm_qn.txt", header = T)

Spur2Aque_full_JSD <- data.frame(replicate(ncol(match_Spur)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Spur2Aque_full_JSD) <- colnames(Spur2Aque_orth_tpm_qn)[-c(1)]
rownames(Spur2Aque_full_JSD) <- colnames(Aque2Spur_orth_tpm_qn)[-c(1)]
Spur2Aque_mean_JSD <- data.frame(replicate(ncol(match_Spur)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Spur2Aque_mean_JSD) <- colnames(Spur2Aque_orth_tpm_qn)[-c(1)]
rownames(Spur2Aque_mean_JSD) <- colnames(Aque2Spur_orth_tpm_qn)[-c(1)]
Spur2Aque_sd_JSD <- data.frame(replicate(ncol(match_Spur)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Spur2Aque_sd_JSD) <- colnames(Spur2Aque_orth_tpm_qn)[-c(1)]
rownames(Spur2Aque_sd_JSD) <- colnames(Aque2Spur_orth_tpm_qn)[-c(1)]

for (i in colnames(Spur2Aque_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Aque2Spur_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Spur[[i]],match_Aque[[j]])
    Spur2Aque_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Spur2Aque_mean_JSD[j,i] <- mean(all_JSD)
    Spur2Aque_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/08-S_purpuratus/Results")
write.table(Spur2Aque_full_JSD, "Spur2Aque_full_set_JSD.txt", sep ='\t')
write.table(Spur2Aque_mean_JSD, "Spur2Aque_subset_JSD_mean.txt", sep ='\t')
write.table(Spur2Aque_sd_JSD, "Spur2Aque_subset_JSD_sd.txt", sep ='\t')


# 8. Strongylocentrotus purpuratus vs. Clytica hemisphaerica
# Number of orthologues = 4,691 orthologues
# 9 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/08-S_purpuratus/Quantile_transformation")
Spur2Chem_orth_tpm_qn <- read.csv("Spur2Chem_TPM_mean_quantile_transform.csv", header = T)
Chem2Spur_orth_tpm_qn <- read.csv("Chem2Spur_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/08-S_purpuratus/one2one")
both_ID <- read.table("Spur2Chem.txt", header = T)
match_Spur <- Spur2Chem_orth_tpm_qn[match(both_ID$Spur,Spur2Chem_orth_tpm_qn$Gene_ID),]
match_Spur <- na.omit(match_Spur)
match_Chem <- Chem2Spur_orth_tpm_qn[match(both_ID$Chem,Chem2Spur_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/08-S_purpuratus/order")
write.table(match_Spur,"order_Spur2Chem_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Chem,"order_Chem2Spur_tpm_qn.txt", sep="\t", quote=F)
match_Spur <- read.table("order_Spur2Chem_tpm_qn.txt", header = T)
match_Chem <- read.table("order_Chem2Spur_tpm_qn.txt", header = T)

Spur2Chem_full_JSD <- data.frame(replicate(ncol(match_Spur)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Spur2Chem_full_JSD) <- colnames(Spur2Chem_orth_tpm_qn)[-c(1)]
rownames(Spur2Chem_full_JSD) <- colnames(Chem2Spur_orth_tpm_qn)[-c(1)]
Spur2Chem_mean_JSD <- data.frame(replicate(ncol(match_Spur)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Spur2Chem_mean_JSD) <- colnames(Spur2Chem_orth_tpm_qn)[-c(1)]
rownames(Spur2Chem_mean_JSD) <- colnames(Chem2Spur_orth_tpm_qn)[-c(1)]
Spur2Chem_sd_JSD <- data.frame(replicate(ncol(match_Spur)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Spur2Chem_sd_JSD) <- colnames(Spur2Chem_orth_tpm_qn)[-c(1)]
rownames(Spur2Chem_sd_JSD) <- colnames(Chem2Spur_orth_tpm_qn)[-c(1)]

for (i in colnames(Spur2Chem_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Chem2Spur_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Spur[[i]],match_Chem[[j]])
    Spur2Chem_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Spur2Chem_mean_JSD[j,i] <- mean(all_JSD)
    Spur2Chem_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/08-S_purpuratus/Results")
write.table(Spur2Chem_full_JSD, "Spur2Chem_full_set_JSD.txt", sep ='\t')
write.table(Spur2Chem_mean_JSD, "Spur2Chem_subset_JSD_mean.txt", sep ='\t')
write.table(Spur2Chem_sd_JSD, "Spur2Chem_subset_JSD_sd.txt", sep ='\t')