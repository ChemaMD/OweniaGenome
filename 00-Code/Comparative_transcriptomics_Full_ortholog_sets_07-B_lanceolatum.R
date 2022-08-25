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
## Branchiostoma lanceolatum 1-to-1 comparisons ##
##########################################

# 1. Branchiostoma lanceolatum vs. Nematostella vectensis
# Number of orthologues = 5,254 orthologues (-4)
# 8 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/07-B_lanceolatum/Quantile_transformation")
Blan2Nvec_orth_tpm_qn <- read.csv("Blan2Nvec_TPM_mean_quantile_transform.csv", header = T)
Nvec2Blan_orth_tpm_qn <- read.csv("Nvec2Blan_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/07-B_lanceolatum/one2one")
both_ID <- read.table("Blan2Nvec.txt", header = T)
match_Blan <- Blan2Nvec_orth_tpm_qn[match(both_ID$Blan,Blan2Nvec_orth_tpm_qn$Gene_ID),]
match_Nvec <- Nvec2Blan_orth_tpm_qn[match(both_ID$Nvec,Nvec2Blan_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/07-B_lanceolatum/order")
write.table(match_Blan,"order_Blan2Nvec_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Nvec,"order_Nvec2Blan_tpm_qn.txt", sep="\t", quote=F)
match_Blan <- read.table("order_Blan2Nvec_tpm_qn.txt", header = T)
match_Nvec <- read.table("order_Nvec2Blan_tpm_qn.txt", header = T)

Blan2Nvec_full_JSD <- data.frame(replicate(ncol(match_Blan)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Blan2Nvec_full_JSD) <- colnames(Blan2Nvec_orth_tpm_qn)[-c(1)]
rownames(Blan2Nvec_full_JSD) <- colnames(Nvec2Blan_orth_tpm_qn)[-c(1)]
Blan2Nvec_mean_JSD <- data.frame(replicate(ncol(match_Blan)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Blan2Nvec_mean_JSD) <- colnames(Blan2Nvec_orth_tpm_qn)[-c(1)]
rownames(Blan2Nvec_mean_JSD) <- colnames(Nvec2Blan_orth_tpm_qn)[-c(1)]
Blan2Nvec_sd_JSD <- data.frame(replicate(ncol(match_Blan)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Blan2Nvec_sd_JSD) <- colnames(Blan2Nvec_orth_tpm_qn)[-c(1)]
rownames(Blan2Nvec_sd_JSD) <- colnames(Nvec2Blan_orth_tpm_qn)[-c(1)]

for (i in colnames(Blan2Nvec_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Nvec2Blan_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Blan[[i]],match_Nvec[[j]])
    Blan2Nvec_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Blan2Nvec_mean_JSD[j,i] <- mean(all_JSD)
    Blan2Nvec_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/07-B_lanceolatum/Results")
write.table(Blan2Nvec_full_JSD, "Blan2Nvec_full_set_JSD.txt", sep ='\t')
write.table(Blan2Nvec_mean_JSD, "Blan2Nvec_subset_JSD_mean.txt", sep ='\t')
write.table(Blan2Nvec_sd_JSD, "Blan2Nvec_subset_JSD_sd.txt", sep ='\t')


# 2. Branchiostoma lanceolatum vs. Strongylocentrotus purpuratus
# Number of orthologues = 5,015 orthologues (+0)
# 14 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/07-B_lanceolatum/Quantile_transformation")
Blan2Spur_orth_tpm_qn <- read.csv("Blan2Spur_TPM_mean_quantile_transform.csv", header = T)
Spur2Blan_orth_tpm_qn <- read.csv("Spur2Blan_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/07-B_lanceolatum/one2one")
both_ID <- read.table("Blan2Spur.txt", header = T)
match_Blan <- Blan2Spur_orth_tpm_qn[match(both_ID$Blan,Blan2Spur_orth_tpm_qn$Gene_ID),]
match_Blan <- na.omit(match_Blan)
match_Spur <- Spur2Blan_orth_tpm_qn[match(both_ID$Spur,Spur2Blan_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/07-B_lanceolatum/order")
write.table(match_Blan,"order_Blan2Spur_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Spur,"order_Spur2Blan_tpm_qn.txt", sep="\t", quote=F)
match_Blan <- read.table("order_Blan2Spur_tpm_qn.txt", header = T)
match_Spur <- read.table("order_Spur2Blan_tpm_qn.txt", header = T)

Blan2Spur_full_JSD <- data.frame(replicate(ncol(match_Blan)-1,sample(0:1,ncol(match_Spur)-1,rep=TRUE)))
colnames(Blan2Spur_full_JSD) <- colnames(Blan2Spur_orth_tpm_qn)[-c(1)]
rownames(Blan2Spur_full_JSD) <- colnames(Spur2Blan_orth_tpm_qn)[-c(1)]
Blan2Spur_mean_JSD <- data.frame(replicate(ncol(match_Blan)-1,sample(0:1,ncol(match_Spur)-1,rep=TRUE)))
colnames(Blan2Spur_mean_JSD) <- colnames(Blan2Spur_orth_tpm_qn)[-c(1)]
rownames(Blan2Spur_mean_JSD) <- colnames(Spur2Blan_orth_tpm_qn)[-c(1)]
Blan2Spur_sd_JSD <- data.frame(replicate(ncol(match_Blan)-1,sample(0:1,ncol(match_Spur)-1,rep=TRUE)))
colnames(Blan2Spur_sd_JSD) <- colnames(Blan2Spur_orth_tpm_qn)[-c(1)]
rownames(Blan2Spur_sd_JSD) <- colnames(Spur2Blan_orth_tpm_qn)[-c(1)]

for (i in colnames(Blan2Spur_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Spur2Blan_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Blan[[i]],match_Spur[[j]])
    Blan2Spur_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Blan2Spur_mean_JSD[j,i] <- mean(all_JSD)
    Blan2Spur_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/07-B_lanceolatum/Results")
write.table(Blan2Spur_full_JSD, "Blan2Spur_full_set_JSD.txt", sep ='\t')
write.table(Blan2Spur_mean_JSD, "Blan2Spur_subset_JSD_mean.txt", sep ='\t')
write.table(Blan2Spur_sd_JSD, "Blan2Spur_subset_JSD_sd.txt", sep ='\t')


# 7. Branchiostoma lanceolatum vs. Amphimedon queenslandica
# Number of orthologues = 3,962 orthologues
# 9 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/07-B_lanceolatum/Quantile_transformation")
Blan2Aque_orth_tpm_qn <- read.csv("Blan2Aque_TPM_mean_quantile_transform.csv", header = T)
Aque2Blan_orth_tpm_qn <- read.csv("Aque2Blan_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/07-B_lanceolatum/one2one")
both_ID <- read.table("Blan2Aque.txt", header = T)
match_Blan <- Blan2Aque_orth_tpm_qn[match(both_ID$Blan,Blan2Aque_orth_tpm_qn$Gene_ID),]
match_Blan <- na.omit(match_Blan)
match_Aque <- Aque2Blan_orth_tpm_qn[match(both_ID$Aque,Aque2Blan_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/07-B_lanceolatum/order")
write.table(match_Blan,"order_Blan2Aque_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Aque,"order_Aque2Blan_tpm_qn.txt", sep="\t", quote=F)
match_Blan <- read.table("order_Blan2Aque_tpm_qn.txt", header = T)
match_Aque <- read.table("order_Aque2Blan_tpm_qn.txt", header = T)

Blan2Aque_full_JSD <- data.frame(replicate(ncol(match_Blan)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Blan2Aque_full_JSD) <- colnames(Blan2Aque_orth_tpm_qn)[-c(1)]
rownames(Blan2Aque_full_JSD) <- colnames(Aque2Blan_orth_tpm_qn)[-c(1)]
Blan2Aque_mean_JSD <- data.frame(replicate(ncol(match_Blan)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Blan2Aque_mean_JSD) <- colnames(Blan2Aque_orth_tpm_qn)[-c(1)]
rownames(Blan2Aque_mean_JSD) <- colnames(Aque2Blan_orth_tpm_qn)[-c(1)]
Blan2Aque_sd_JSD <- data.frame(replicate(ncol(match_Blan)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Blan2Aque_sd_JSD) <- colnames(Blan2Aque_orth_tpm_qn)[-c(1)]
rownames(Blan2Aque_sd_JSD) <- colnames(Aque2Blan_orth_tpm_qn)[-c(1)]

for (i in colnames(Blan2Aque_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Aque2Blan_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Blan[[i]],match_Aque[[j]])
    Blan2Aque_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Blan2Aque_mean_JSD[j,i] <- mean(all_JSD)
    Blan2Aque_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/07-B_lanceolatum/Results")
write.table(Blan2Aque_full_JSD, "Blan2Aque_full_set_JSD.txt", sep ='\t')
write.table(Blan2Aque_mean_JSD, "Blan2Aque_subset_JSD_mean.txt", sep ='\t')
write.table(Blan2Aque_sd_JSD, "Blan2Aque_subset_JSD_sd.txt", sep ='\t')


# 8. Branchiostoma lanceolatum vs. Clytica hemisphaerica
# Number of orthologues = 4,691 orthologues
# 9 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/07-B_lanceolatum/Quantile_transformation")
Blan2Chem_orth_tpm_qn <- read.csv("Blan2Chem_TPM_mean_quantile_transform.csv", header = T)
Chem2Blan_orth_tpm_qn <- read.csv("Chem2Blan_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/07-B_lanceolatum/one2one")
both_ID <- read.table("Blan2Chem.txt", header = T)
match_Blan <- Blan2Chem_orth_tpm_qn[match(both_ID$Blan,Blan2Chem_orth_tpm_qn$Gene_ID),]
match_Blan <- na.omit(match_Blan)
match_Chem <- Chem2Blan_orth_tpm_qn[match(both_ID$Chem,Chem2Blan_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/07-B_lanceolatum/order")
write.table(match_Blan,"order_Blan2Chem_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Chem,"order_Chem2Blan_tpm_qn.txt", sep="\t", quote=F)
match_Blan <- read.table("order_Blan2Chem_tpm_qn.txt", header = T)
match_Chem <- read.table("order_Chem2Blan_tpm_qn.txt", header = T)

Blan2Chem_full_JSD <- data.frame(replicate(ncol(match_Blan)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Blan2Chem_full_JSD) <- colnames(Blan2Chem_orth_tpm_qn)[-c(1)]
rownames(Blan2Chem_full_JSD) <- colnames(Chem2Blan_orth_tpm_qn)[-c(1)]
Blan2Chem_mean_JSD <- data.frame(replicate(ncol(match_Blan)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Blan2Chem_mean_JSD) <- colnames(Blan2Chem_orth_tpm_qn)[-c(1)]
rownames(Blan2Chem_mean_JSD) <- colnames(Chem2Blan_orth_tpm_qn)[-c(1)]
Blan2Chem_sd_JSD <- data.frame(replicate(ncol(match_Blan)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Blan2Chem_sd_JSD) <- colnames(Blan2Chem_orth_tpm_qn)[-c(1)]
rownames(Blan2Chem_sd_JSD) <- colnames(Chem2Blan_orth_tpm_qn)[-c(1)]

for (i in colnames(Blan2Chem_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Chem2Blan_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Blan[[i]],match_Chem[[j]])
    Blan2Chem_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Blan2Chem_mean_JSD[j,i] <- mean(all_JSD)
    Blan2Chem_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/07-B_lanceolatum/Results")
write.table(Blan2Chem_full_JSD, "Blan2Chem_full_set_JSD.txt", sep ='\t')
write.table(Blan2Chem_mean_JSD, "Blan2Chem_subset_JSD_mean.txt", sep ='\t')
write.table(Blan2Chem_sd_JSD, "Blan2Chem_subset_JSD_sd.txt", sep ='\t')