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
## Danio rerio 1-to-1 comparisons ##
##########################################

# 1. Danio rerio vs. Nematostella vectensis
# Number of orthologues = 5,254 orthologues (-4)
# 8 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/06-D_rerio/Quantile_transformation")
Drer2Nvec_orth_tpm_qn <- read.csv("Drer2Nvec_TPM_mean_quantile_transform.csv", header = T)
Nvec2Drer_orth_tpm_qn <- read.csv("Nvec2Drer_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/06-D_rerio/one2one")
both_ID <- read.table("Drer2Nvec.txt", header = T)
match_Drer <- Drer2Nvec_orth_tpm_qn[match(both_ID$Drer,Drer2Nvec_orth_tpm_qn$Gene_ID),]
match_Nvec <- Nvec2Drer_orth_tpm_qn[match(both_ID$Nvec,Nvec2Drer_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/06-D_rerio/order")
write.table(match_Drer,"order_Drer2Nvec_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Nvec,"order_Nvec2Drer_tpm_qn.txt", sep="\t", quote=F)
match_Drer <- read.table("order_Drer2Nvec_tpm_qn.txt", header = T)
match_Nvec <- read.table("order_Nvec2Drer_tpm_qn.txt", header = T)

Drer2Nvec_full_JSD <- data.frame(replicate(ncol(match_Drer)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Drer2Nvec_full_JSD) <- colnames(Drer2Nvec_orth_tpm_qn)[-c(1)]
rownames(Drer2Nvec_full_JSD) <- colnames(Nvec2Drer_orth_tpm_qn)[-c(1)]
Drer2Nvec_mean_JSD <- data.frame(replicate(ncol(match_Drer)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Drer2Nvec_mean_JSD) <- colnames(Drer2Nvec_orth_tpm_qn)[-c(1)]
rownames(Drer2Nvec_mean_JSD) <- colnames(Nvec2Drer_orth_tpm_qn)[-c(1)]
Drer2Nvec_sd_JSD <- data.frame(replicate(ncol(match_Drer)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Drer2Nvec_sd_JSD) <- colnames(Drer2Nvec_orth_tpm_qn)[-c(1)]
rownames(Drer2Nvec_sd_JSD) <- colnames(Nvec2Drer_orth_tpm_qn)[-c(1)]

for (i in colnames(Drer2Nvec_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Nvec2Drer_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Drer[[i]],match_Nvec[[j]])
    Drer2Nvec_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Drer2Nvec_mean_JSD[j,i] <- mean(all_JSD)
    Drer2Nvec_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/06-D_rerio/Results")
write.table(Drer2Nvec_full_JSD, "Drer2Nvec_full_set_JSD.txt", sep ='\t')
write.table(Drer2Nvec_mean_JSD, "Drer2Nvec_subset_JSD_mean.txt", sep ='\t')
write.table(Drer2Nvec_sd_JSD, "Drer2Nvec_subset_JSD_sd.txt", sep ='\t')


# 2. Danio rerio vs. Strongylocentrotus purpuratus
# Number of orthologues = 5,015 orthologues (+0)
# 14 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/06-D_rerio/Quantile_transformation")
Drer2Spur_orth_tpm_qn <- read.csv("Drer2Spur_TPM_mean_quantile_transform.csv", header = T)
Spur2Drer_orth_tpm_qn <- read.csv("Spur2Drer_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/06-D_rerio/one2one")
both_ID <- read.table("Drer2Spur.txt", header = T)
match_Drer <- Drer2Spur_orth_tpm_qn[match(both_ID$Drer,Drer2Spur_orth_tpm_qn$Gene_ID),]
match_Drer <- na.omit(match_Drer)
match_Spur <- Spur2Drer_orth_tpm_qn[match(both_ID$Spur,Spur2Drer_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/06-D_rerio/order")
write.table(match_Drer,"order_Drer2Spur_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Spur,"order_Spur2Drer_tpm_qn.txt", sep="\t", quote=F)
match_Drer <- read.table("order_Drer2Spur_tpm_qn.txt", header = T)
match_Spur <- read.table("order_Spur2Drer_tpm_qn.txt", header = T)

Drer2Spur_full_JSD <- data.frame(replicate(ncol(match_Drer)-1,sample(0:1,ncol(match_Spur)-1,rep=TRUE)))
colnames(Drer2Spur_full_JSD) <- colnames(Drer2Spur_orth_tpm_qn)[-c(1)]
rownames(Drer2Spur_full_JSD) <- colnames(Spur2Drer_orth_tpm_qn)[-c(1)]
Drer2Spur_mean_JSD <- data.frame(replicate(ncol(match_Drer)-1,sample(0:1,ncol(match_Spur)-1,rep=TRUE)))
colnames(Drer2Spur_mean_JSD) <- colnames(Drer2Spur_orth_tpm_qn)[-c(1)]
rownames(Drer2Spur_mean_JSD) <- colnames(Spur2Drer_orth_tpm_qn)[-c(1)]
Drer2Spur_sd_JSD <- data.frame(replicate(ncol(match_Drer)-1,sample(0:1,ncol(match_Spur)-1,rep=TRUE)))
colnames(Drer2Spur_sd_JSD) <- colnames(Drer2Spur_orth_tpm_qn)[-c(1)]
rownames(Drer2Spur_sd_JSD) <- colnames(Spur2Drer_orth_tpm_qn)[-c(1)]

for (i in colnames(Drer2Spur_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Spur2Drer_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Drer[[i]],match_Spur[[j]])
    Drer2Spur_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Drer2Spur_mean_JSD[j,i] <- mean(all_JSD)
    Drer2Spur_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/06-D_rerio/Results")
write.table(Drer2Spur_full_JSD, "Drer2Spur_full_set_JSD.txt", sep ='\t')
write.table(Drer2Spur_mean_JSD, "Drer2Spur_subset_JSD_mean.txt", sep ='\t')
write.table(Drer2Spur_sd_JSD, "Drer2Spur_subset_JSD_sd.txt", sep ='\t')


# 3. Danio rerio vs. Branchiostoma lanceolatum
# Number of orthologues = 6,673 orthologues (+39)
# 14 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/06-D_rerio/Quantile_transformation")
Drer2Blan_orth_tpm_qn <- read.csv("Drer2Blan_TPM_mean_quantile_transform.csv", header = T)
Blan2Drer_orth_tpm_qn <- read.csv("Blan2Drer_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/06-D_rerio/one2one")
both_ID <- read.table("Drer2Blan.txt", header = T)
match_Drer <- Drer2Blan_orth_tpm_qn[match(both_ID$Drer,Drer2Blan_orth_tpm_qn$Gene_ID),]
match_Drer <- na.omit(match_Drer)
match_Blan <- Blan2Drer_orth_tpm_qn[match(both_ID$Blan,Blan2Drer_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/06-D_rerio/order")
write.table(match_Drer,"order_Drer2Blan_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Blan,"order_Blan2Drer_tpm_qn.txt", sep="\t", quote=F)
match_Drer <- read.table("order_Drer2Blan_tpm_qn.txt", header = T)
match_Blan <- read.table("order_Blan2Drer_tpm_qn.txt", header = T)

Drer2Blan_full_JSD <- data.frame(replicate(ncol(match_Drer)-1,sample(0:1,ncol(match_Blan)-1,rep=TRUE)))
colnames(Drer2Blan_full_JSD) <- colnames(Drer2Blan_orth_tpm_qn)[-c(1)]
rownames(Drer2Blan_full_JSD) <- colnames(Blan2Drer_orth_tpm_qn)[-c(1)]
Drer2Blan_mean_JSD <- data.frame(replicate(ncol(match_Drer)-1,sample(0:1,ncol(match_Blan)-1,rep=TRUE)))
colnames(Drer2Blan_mean_JSD) <- colnames(Drer2Blan_orth_tpm_qn)[-c(1)]
rownames(Drer2Blan_mean_JSD) <- colnames(Blan2Drer_orth_tpm_qn)[-c(1)]
Drer2Blan_sd_JSD <- data.frame(replicate(ncol(match_Drer)-1,sample(0:1,ncol(match_Blan)-1,rep=TRUE)))
colnames(Drer2Blan_sd_JSD) <- colnames(Drer2Blan_orth_tpm_qn)[-c(1)]
rownames(Drer2Blan_sd_JSD) <- colnames(Blan2Drer_orth_tpm_qn)[-c(1)]

for (i in colnames(Drer2Blan_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Blan2Drer_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Drer[[i]],match_Blan[[j]])
    Drer2Blan_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Drer2Blan_mean_JSD[j,i] <- mean(all_JSD)
    Drer2Blan_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/06-D_rerio/Results")
write.table(Drer2Blan_full_JSD, "Drer2Blan_full_set_JSD.txt", sep ='\t')
write.table(Drer2Blan_mean_JSD, "Drer2Blan_subset_JSD_mean.txt", sep ='\t')
write.table(Drer2Blan_sd_JSD, "Drer2Blan_subset_JSD_sd.txt", sep ='\t')


# 7. Danio rerio vs. Amphimedon queenslandica
# Number of orthologues = 3,962 orthologues
# 9 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/06-D_rerio/Quantile_transformation")
Drer2Aque_orth_tpm_qn <- read.csv("Drer2Aque_TPM_mean_quantile_transform.csv", header = T)
Aque2Drer_orth_tpm_qn <- read.csv("Aque2Drer_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/06-D_rerio/one2one")
both_ID <- read.table("Drer2Aque.txt", header = T)
match_Drer <- Drer2Aque_orth_tpm_qn[match(both_ID$Drer,Drer2Aque_orth_tpm_qn$Gene_ID),]
match_Drer <- na.omit(match_Drer)
match_Aque <- Aque2Drer_orth_tpm_qn[match(both_ID$Aque,Aque2Drer_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/06-D_rerio/order")
write.table(match_Drer,"order_Drer2Aque_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Aque,"order_Aque2Drer_tpm_qn.txt", sep="\t", quote=F)
match_Drer <- read.table("order_Drer2Aque_tpm_qn.txt", header = T)
match_Aque <- read.table("order_Aque2Drer_tpm_qn.txt", header = T)

Drer2Aque_full_JSD <- data.frame(replicate(ncol(match_Drer)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Drer2Aque_full_JSD) <- colnames(Drer2Aque_orth_tpm_qn)[-c(1)]
rownames(Drer2Aque_full_JSD) <- colnames(Aque2Drer_orth_tpm_qn)[-c(1)]
Drer2Aque_mean_JSD <- data.frame(replicate(ncol(match_Drer)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Drer2Aque_mean_JSD) <- colnames(Drer2Aque_orth_tpm_qn)[-c(1)]
rownames(Drer2Aque_mean_JSD) <- colnames(Aque2Drer_orth_tpm_qn)[-c(1)]
Drer2Aque_sd_JSD <- data.frame(replicate(ncol(match_Drer)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Drer2Aque_sd_JSD) <- colnames(Drer2Aque_orth_tpm_qn)[-c(1)]
rownames(Drer2Aque_sd_JSD) <- colnames(Aque2Drer_orth_tpm_qn)[-c(1)]

for (i in colnames(Drer2Aque_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Aque2Drer_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Drer[[i]],match_Aque[[j]])
    Drer2Aque_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Drer2Aque_mean_JSD[j,i] <- mean(all_JSD)
    Drer2Aque_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/06-D_rerio/Results")
write.table(Drer2Aque_full_JSD, "Drer2Aque_full_set_JSD.txt", sep ='\t')
write.table(Drer2Aque_mean_JSD, "Drer2Aque_subset_JSD_mean.txt", sep ='\t')
write.table(Drer2Aque_sd_JSD, "Drer2Aque_subset_JSD_sd.txt", sep ='\t')


# 8. Danio rerio vs. Clytica hemisphaerica
# Number of orthologues = 4,691 orthologues
# 9 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/06-D_rerio/Quantile_transformation")
Drer2Chem_orth_tpm_qn <- read.csv("Drer2Chem_TPM_mean_quantile_transform.csv", header = T)
Chem2Drer_orth_tpm_qn <- read.csv("Chem2Drer_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/06-D_rerio/one2one")
both_ID <- read.table("Drer2Chem.txt", header = T)
match_Drer <- Drer2Chem_orth_tpm_qn[match(both_ID$Drer,Drer2Chem_orth_tpm_qn$Gene_ID),]
match_Drer <- na.omit(match_Drer)
match_Chem <- Chem2Drer_orth_tpm_qn[match(both_ID$Chem,Chem2Drer_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/06-D_rerio/order")
write.table(match_Drer,"order_Drer2Chem_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Chem,"order_Chem2Drer_tpm_qn.txt", sep="\t", quote=F)
match_Drer <- read.table("order_Drer2Chem_tpm_qn.txt", header = T)
match_Chem <- read.table("order_Chem2Drer_tpm_qn.txt", header = T)

Drer2Chem_full_JSD <- data.frame(replicate(ncol(match_Drer)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Drer2Chem_full_JSD) <- colnames(Drer2Chem_orth_tpm_qn)[-c(1)]
rownames(Drer2Chem_full_JSD) <- colnames(Chem2Drer_orth_tpm_qn)[-c(1)]
Drer2Chem_mean_JSD <- data.frame(replicate(ncol(match_Drer)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Drer2Chem_mean_JSD) <- colnames(Drer2Chem_orth_tpm_qn)[-c(1)]
rownames(Drer2Chem_mean_JSD) <- colnames(Chem2Drer_orth_tpm_qn)[-c(1)]
Drer2Chem_sd_JSD <- data.frame(replicate(ncol(match_Drer)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Drer2Chem_sd_JSD) <- colnames(Drer2Chem_orth_tpm_qn)[-c(1)]
rownames(Drer2Chem_sd_JSD) <- colnames(Chem2Drer_orth_tpm_qn)[-c(1)]

for (i in colnames(Drer2Chem_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Chem2Drer_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Drer[[i]],match_Chem[[j]])
    Drer2Chem_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Drer2Chem_mean_JSD[j,i] <- mean(all_JSD)
    Drer2Chem_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/06-D_rerio/Results")
write.table(Drer2Chem_full_JSD, "Drer2Chem_full_set_JSD.txt", sep ='\t')
write.table(Drer2Chem_mean_JSD, "Drer2Chem_subset_JSD_mean.txt", sep ='\t')
write.table(Drer2Chem_sd_JSD, "Drer2Chem_subset_JSD_sd.txt", sep ='\t')