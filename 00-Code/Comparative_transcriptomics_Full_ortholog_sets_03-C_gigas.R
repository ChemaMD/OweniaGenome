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
## Crassostrea gigas 1-to-1 comparisons ##
##########################################

# 1. Crassostrea gigas vs. Nematostella vectensis
# Number of orthologues = 5,254 orthologues (-4)
# 8 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/Quantile_transformation")
Cgig2Nvec_orth_tpm_qn <- read.csv("Cgig2Nvec_TPM_mean_quantile_transform.csv", header = T)
Nvec2Cgig_orth_tpm_qn <- read.csv("Nvec2Cgig_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/one2one")
both_ID <- read.table("Cgig2Nvec.txt", header = T)
match_Cgig <- Cgig2Nvec_orth_tpm_qn[match(both_ID$Cgig,Cgig2Nvec_orth_tpm_qn$Gene_ID),]
match_Nvec <- Nvec2Cgig_orth_tpm_qn[match(both_ID$Nvec,Nvec2Cgig_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/order")
write.table(match_Cgig,"order_Cgig2Nvec_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Nvec,"order_Nvec2Cgig_tpm_qn.txt", sep="\t", quote=F)
match_Cgig <- read.table("order_Cgig2Nvec_tpm_qn.txt", header = T)
match_Nvec <- read.table("order_Nvec2Cgig_tpm_qn.txt", header = T)

Cgig2Nvec_full_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Cgig2Nvec_full_JSD) <- colnames(Cgig2Nvec_orth_tpm_qn)[-c(1)]
rownames(Cgig2Nvec_full_JSD) <- colnames(Nvec2Cgig_orth_tpm_qn)[-c(1)]
Cgig2Nvec_mean_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Cgig2Nvec_mean_JSD) <- colnames(Cgig2Nvec_orth_tpm_qn)[-c(1)]
rownames(Cgig2Nvec_mean_JSD) <- colnames(Nvec2Cgig_orth_tpm_qn)[-c(1)]
Cgig2Nvec_sd_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Cgig2Nvec_sd_JSD) <- colnames(Cgig2Nvec_orth_tpm_qn)[-c(1)]
rownames(Cgig2Nvec_sd_JSD) <- colnames(Nvec2Cgig_orth_tpm_qn)[-c(1)]

for (i in colnames(Cgig2Nvec_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Nvec2Cgig_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cgig[[i]],match_Nvec[[j]])
    Cgig2Nvec_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cgig2Nvec_mean_JSD[j,i] <- mean(all_JSD)
    Cgig2Nvec_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/Results")
write.table(Cgig2Nvec_full_JSD, "Cgig2Nvec_full_set_JSD.txt", sep ='\t')
write.table(Cgig2Nvec_mean_JSD, "Cgig2Nvec_subset_JSD_mean.txt", sep ='\t')
write.table(Cgig2Nvec_sd_JSD, "Cgig2Nvec_subset_JSD_sd.txt", sep ='\t')


# 2. Crassostrea gigas vs. Strongylocentrotus purpuratus
# Number of orthologues = 5,015 orthologues (+0)
# 14 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/Quantile_transformation")
Cgig2Spur_orth_tpm_qn <- read.csv("Cgig2Spur_TPM_mean_quantile_transform.csv", header = T)
Spur2Cgig_orth_tpm_qn <- read.csv("Spur2Cgig_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/one2one")
both_ID <- read.table("Cgig2Spur.txt", header = T)
match_Cgig <- Cgig2Spur_orth_tpm_qn[match(both_ID$Cgig,Cgig2Spur_orth_tpm_qn$Gene_ID),]
match_Cgig <- na.omit(match_Cgig)
match_Spur <- Spur2Cgig_orth_tpm_qn[match(both_ID$Spur,Spur2Cgig_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/order")
write.table(match_Cgig,"order_Cgig2Spur_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Spur,"order_Spur2Cgig_tpm_qn.txt", sep="\t", quote=F)
match_Cgig <- read.table("order_Cgig2Spur_tpm_qn.txt", header = T)
match_Spur <- read.table("order_Spur2Cgig_tpm_qn.txt", header = T)

Cgig2Spur_full_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Spur)-1,rep=TRUE)))
colnames(Cgig2Spur_full_JSD) <- colnames(Cgig2Spur_orth_tpm_qn)[-c(1)]
rownames(Cgig2Spur_full_JSD) <- colnames(Spur2Cgig_orth_tpm_qn)[-c(1)]
Cgig2Spur_mean_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Spur)-1,rep=TRUE)))
colnames(Cgig2Spur_mean_JSD) <- colnames(Cgig2Spur_orth_tpm_qn)[-c(1)]
rownames(Cgig2Spur_mean_JSD) <- colnames(Spur2Cgig_orth_tpm_qn)[-c(1)]
Cgig2Spur_sd_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Spur)-1,rep=TRUE)))
colnames(Cgig2Spur_sd_JSD) <- colnames(Cgig2Spur_orth_tpm_qn)[-c(1)]
rownames(Cgig2Spur_sd_JSD) <- colnames(Spur2Cgig_orth_tpm_qn)[-c(1)]

for (i in colnames(Cgig2Spur_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Spur2Cgig_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cgig[[i]],match_Spur[[j]])
    Cgig2Spur_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cgig2Spur_mean_JSD[j,i] <- mean(all_JSD)
    Cgig2Spur_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/Results")
write.table(Cgig2Spur_full_JSD, "Cgig2Spur_full_set_JSD.txt", sep ='\t')
write.table(Cgig2Spur_mean_JSD, "Cgig2Spur_subset_JSD_mean.txt", sep ='\t')
write.table(Cgig2Spur_sd_JSD, "Cgig2Spur_subset_JSD_sd.txt", sep ='\t')


# 3. Crassostrea gigas vs. Branchiostoma lanceolatum
# Number of orthologues = 6,673 orthologues (+39)
# 14 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/Quantile_transformation")
Cgig2Blan_orth_tpm_qn <- read.csv("Cgig2Blan_TPM_mean_quantile_transform.csv", header = T)
Blan2Cgig_orth_tpm_qn <- read.csv("Blan2Cgig_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/one2one")
both_ID <- read.table("Cgig2Blan.txt", header = T)
match_Cgig <- Cgig2Blan_orth_tpm_qn[match(both_ID$Cgig,Cgig2Blan_orth_tpm_qn$Gene_ID),]
match_Cgig <- na.omit(match_Cgig)
match_Blan <- Blan2Cgig_orth_tpm_qn[match(both_ID$Blan,Blan2Cgig_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/order")
write.table(match_Cgig,"order_Cgig2Blan_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Blan,"order_Blan2Cgig_tpm_qn.txt", sep="\t", quote=F)
match_Cgig <- read.table("order_Cgig2Blan_tpm_qn.txt", header = T)
match_Blan <- read.table("order_Blan2Cgig_tpm_qn.txt", header = T)

Cgig2Blan_full_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Blan)-1,rep=TRUE)))
colnames(Cgig2Blan_full_JSD) <- colnames(Cgig2Blan_orth_tpm_qn)[-c(1)]
rownames(Cgig2Blan_full_JSD) <- colnames(Blan2Cgig_orth_tpm_qn)[-c(1)]
Cgig2Blan_mean_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Blan)-1,rep=TRUE)))
colnames(Cgig2Blan_mean_JSD) <- colnames(Cgig2Blan_orth_tpm_qn)[-c(1)]
rownames(Cgig2Blan_mean_JSD) <- colnames(Blan2Cgig_orth_tpm_qn)[-c(1)]
Cgig2Blan_sd_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Blan)-1,rep=TRUE)))
colnames(Cgig2Blan_sd_JSD) <- colnames(Cgig2Blan_orth_tpm_qn)[-c(1)]
rownames(Cgig2Blan_sd_JSD) <- colnames(Blan2Cgig_orth_tpm_qn)[-c(1)]

for (i in colnames(Cgig2Blan_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Blan2Cgig_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cgig[[i]],match_Blan[[j]])
    Cgig2Blan_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cgig2Blan_mean_JSD[j,i] <- mean(all_JSD)
    Cgig2Blan_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/Results")
write.table(Cgig2Blan_full_JSD, "Cgig2Blan_full_set_JSD.txt", sep ='\t')
write.table(Cgig2Blan_mean_JSD, "Cgig2Blan_subset_JSD_mean.txt", sep ='\t')
write.table(Cgig2Blan_sd_JSD, "Cgig2Blan_subset_JSD_sd.txt", sep ='\t')


# 4. Crassostrea gigas vs. Danio rerio
# Number of orthologues = 4,316 orthologues (+21)
# 7 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/Quantile_transformation")
Cgig2Drer_orth_tpm_qn <- read.csv("Cgig2Drer_TPM_mean_quantile_transform.csv", header = T)
Drer2Cgig_orth_tpm_qn <- read.csv("Drer2Cgig_TPM_mean_quantile_transform.csv", header = T)
#Drer2Cgig_orth_tpm_qn <- cbind(Drer2Cgig_orth_tpm_qn_preliminary[1],Drer2Cgig_orth_tpm_qn_preliminary[12],Drer2Cgig_orth_tpm_qn_preliminary[c(3:11)],Drer2Cgig_orth_tpm_qn_preliminary[2],Drer2Cgig_orth_tpm_qn_preliminary[c(13:19)])

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/one2one")
both_ID <- read.table("Cgig2Drer.txt", header = T)
match_Cgig <- Cgig2Drer_orth_tpm_qn[match(both_ID$Cgig,Cgig2Drer_orth_tpm_qn$Gene_ID),]
match_Cgig <- na.omit(match_Cgig)
match_Drer <- Drer2Cgig_orth_tpm_qn[match(both_ID$Drer,Drer2Cgig_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/order")
write.table(match_Cgig,"order_Cgig2Drer_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Drer,"order_Drer2Cgig_tpm_qn.txt", sep="\t", quote=F)
match_Cgig <- read.table("order_Cgig2Drer_tpm_qn.txt", header = T)
match_Drer <- read.table("order_Drer2Cgig_tpm_qn.txt", header = T)

Cgig2Drer_full_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Drer)-1,rep=TRUE)))
colnames(Cgig2Drer_full_JSD) <- colnames(Cgig2Drer_orth_tpm_qn)[-c(1)]
rownames(Cgig2Drer_full_JSD) <- colnames(Drer2Cgig_orth_tpm_qn)[-c(1)]
Cgig2Drer_mean_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Drer)-1,rep=TRUE)))
colnames(Cgig2Drer_mean_JSD) <- colnames(Cgig2Drer_orth_tpm_qn)[-c(1)]
rownames(Cgig2Drer_mean_JSD) <- colnames(Drer2Cgig_orth_tpm_qn)[-c(1)]
Cgig2Drer_sd_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Drer)-1,rep=TRUE)))
colnames(Cgig2Drer_sd_JSD) <- colnames(Cgig2Drer_orth_tpm_qn)[-c(1)]
rownames(Cgig2Drer_sd_JSD) <- colnames(Drer2Cgig_orth_tpm_qn)[-c(1)]

for (i in colnames(Cgig2Drer_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Drer2Cgig_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cgig[[i]],match_Drer[[j]])
    Cgig2Drer_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cgig2Drer_mean_JSD[j,i] <- mean(all_JSD)
    Cgig2Drer_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/Results")
write.table(Cgig2Drer_full_JSD, "Cgig2Drer_full_set_JSD.txt", sep ='\t')
write.table(Cgig2Drer_mean_JSD, "Cgig2Drer_subset_JSD_mean.txt", sep ='\t')
write.table(Cgig2Drer_sd_JSD, "Cgig2Drer_subset_JSD_sd.txt", sep ='\t')


# 5. Crassostrea gigas vs. Drosophila melanogaster
# Number of orthologues = 4,635 orthologues (+45)
# 18 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/Quantile_transformation")
Cgig2Dmel_orth_tpm_qn <- read.csv("Cgig2Dmel_TPM_mean_quantile_transform.csv", header = T)
Dmel2Cgig_orth_tpm_qn <- read.csv("Dmel2Cgig_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/one2one")
both_ID <- read.table("Cgig2Dmel.txt", header = T)
match_Cgig <- Cgig2Dmel_orth_tpm_qn[match(both_ID$Cgig,Cgig2Dmel_orth_tpm_qn$Gene_ID),]
match_Cgig <- na.omit(match_Cgig)
match_Dmel <- Dmel2Cgig_orth_tpm_qn[match(both_ID$Dmel,Dmel2Cgig_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/order")
write.table(match_Cgig,"order_Cgig2Dmel_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Dmel,"order_Dmel2Cgig_tpm_qn.txt", sep="\t", quote=F)
match_Cgig <- read.table("order_Cgig2Dmel_tpm_qn.txt", header = T)
match_Dmel <- read.table("order_Dmel2Cgig_tpm_qn.txt", header = T)

Cgig2Dmel_full_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Dmel)-1,rep=TRUE)))
colnames(Cgig2Dmel_full_JSD) <- colnames(Cgig2Dmel_orth_tpm_qn)[-c(1)]
rownames(Cgig2Dmel_full_JSD) <- colnames(Dmel2Cgig_orth_tpm_qn)[-c(1)]
Cgig2Dmel_mean_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Dmel)-1,rep=TRUE)))
colnames(Cgig2Dmel_mean_JSD) <- colnames(Cgig2Dmel_orth_tpm_qn)[-c(1)]
rownames(Cgig2Dmel_mean_JSD) <- colnames(Dmel2Cgig_orth_tpm_qn)[-c(1)]
Cgig2Dmel_sd_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Dmel)-1,rep=TRUE)))
colnames(Cgig2Dmel_sd_JSD) <- colnames(Cgig2Dmel_orth_tpm_qn)[-c(1)]
rownames(Cgig2Dmel_sd_JSD) <- colnames(Dmel2Cgig_orth_tpm_qn)[-c(1)]

for (i in colnames(Cgig2Dmel_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Dmel2Cgig_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cgig[[i]],match_Dmel[[j]])
    Cgig2Dmel_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cgig2Dmel_mean_JSD[j,i] <- mean(all_JSD)
    Cgig2Dmel_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/Results")
write.table(Cgig2Dmel_full_JSD, "Cgig2Dmel_full_set_JSD.txt", sep ='\t')
write.table(Cgig2Dmel_mean_JSD, "Cgig2Dmel_subset_JSD_mean.txt", sep ='\t')
write.table(Cgig2Dmel_sd_JSD, "Cgig2Dmel_subset_JSD_sd.txt", sep ='\t')


# 6. Crassostrea gigas vs. Caenorhabidits elegans
# Number of orthologues = 3,767 orthologues (-43)
# 13 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/Quantile_transformation")
Cgig2Cele_orth_tpm_qn <- read.csv("Cgig2Cele_TPM_mean_quantile_transform.csv", header = T)
Cele2Cgig_orth_tpm_qn <- read.csv("Cele2Cgig_TPM_mean_quantile_transform.csv", header = T)
Cele2Cgig_orth_tpm_qn$Gene_ID <- paste0("CELE_", Cele2Cgig_orth_tpm_qn$Gene_ID)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/one2one")
both_ID <- read.table("Cgig2Cele.txt", header = T)
match_Cgig <- Cgig2Cele_orth_tpm_qn[match(both_ID$Cgig,Cgig2Cele_orth_tpm_qn$Gene_ID),]
match_Cgig <- na.omit(match_Cgig)
match_Cele <- Cele2Cgig_orth_tpm_qn[match(both_ID$Cele,Cele2Cgig_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/order")
write.table(match_Cgig,"order_Cgig2Cele_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Cele,"order_Cele2Cgig_tpm_qn.txt", sep="\t", quote=F)
match_Cgig <- read.table("order_Cgig2Cele_tpm_qn.txt", header = T)
match_Cele <- read.table("order_Cele2Cgig_tpm_qn.txt", header = T)

Cgig2Cele_full_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Cele)-1,rep=TRUE)))
colnames(Cgig2Cele_full_JSD) <- colnames(Cgig2Cele_orth_tpm_qn)[-c(1)]
rownames(Cgig2Cele_full_JSD) <- colnames(Cele2Cgig_orth_tpm_qn)[-c(1)]
Cgig2Cele_mean_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Cele)-1,rep=TRUE)))
colnames(Cgig2Cele_mean_JSD) <- colnames(Cgig2Cele_orth_tpm_qn)[-c(1)]
rownames(Cgig2Cele_mean_JSD) <- colnames(Cele2Cgig_orth_tpm_qn)[-c(1)]
Cgig2Cele_sd_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Cele)-1,rep=TRUE)))
colnames(Cgig2Cele_sd_JSD) <- colnames(Cgig2Cele_orth_tpm_qn)[-c(1)]
rownames(Cgig2Cele_sd_JSD) <- colnames(Cele2Cgig_orth_tpm_qn)[-c(1)]

for (i in colnames(Cgig2Cele_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Cele2Cgig_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cgig[[i]],match_Cele[[j]])
    Cgig2Cele_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cgig2Cele_mean_JSD[j,i] <- mean(all_JSD)
    Cgig2Cele_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/Results")
write.table(Cgig2Cele_full_JSD, "Cgig2Cele_full_set_JSD.txt", sep ='\t')
write.table(Cgig2Cele_mean_JSD, "Cgig2Cele_subset_JSD_mean.txt", sep ='\t')
write.table(Cgig2Cele_sd_JSD, "Cgig2Cele_subset_JSD_sd.txt", sep ='\t')


# 7. Crassostrea gigas vs. Amphimedon queenslandica
# Number of orthologues = 3,962 orthologues
# 9 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/Quantile_transformation")
Cgig2Aque_orth_tpm_qn <- read.csv("Cgig2Aque_TPM_mean_quantile_transform.csv", header = T)
Aque2Cgig_orth_tpm_qn <- read.csv("Aque2Cgig_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/one2one")
both_ID <- read.table("Cgig2Aque.txt", header = T)
match_Cgig <- Cgig2Aque_orth_tpm_qn[match(both_ID$Cgig,Cgig2Aque_orth_tpm_qn$Gene_ID),]
match_Cgig <- na.omit(match_Cgig)
match_Aque <- Aque2Cgig_orth_tpm_qn[match(both_ID$Aque,Aque2Cgig_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/order")
write.table(match_Cgig,"order_Cgig2Aque_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Aque,"order_Aque2Cgig_tpm_qn.txt", sep="\t", quote=F)
match_Cgig <- read.table("order_Cgig2Aque_tpm_qn.txt", header = T)
match_Aque <- read.table("order_Aque2Cgig_tpm_qn.txt", header = T)

Cgig2Aque_full_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Cgig2Aque_full_JSD) <- colnames(Cgig2Aque_orth_tpm_qn)[-c(1)]
rownames(Cgig2Aque_full_JSD) <- colnames(Aque2Cgig_orth_tpm_qn)[-c(1)]
Cgig2Aque_mean_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Cgig2Aque_mean_JSD) <- colnames(Cgig2Aque_orth_tpm_qn)[-c(1)]
rownames(Cgig2Aque_mean_JSD) <- colnames(Aque2Cgig_orth_tpm_qn)[-c(1)]
Cgig2Aque_sd_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Cgig2Aque_sd_JSD) <- colnames(Cgig2Aque_orth_tpm_qn)[-c(1)]
rownames(Cgig2Aque_sd_JSD) <- colnames(Aque2Cgig_orth_tpm_qn)[-c(1)]

for (i in colnames(Cgig2Aque_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Aque2Cgig_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cgig[[i]],match_Aque[[j]])
    Cgig2Aque_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cgig2Aque_mean_JSD[j,i] <- mean(all_JSD)
    Cgig2Aque_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/Results")
write.table(Cgig2Aque_full_JSD, "Cgig2Aque_full_set_JSD.txt", sep ='\t')
write.table(Cgig2Aque_mean_JSD, "Cgig2Aque_subset_JSD_mean.txt", sep ='\t')
write.table(Cgig2Aque_sd_JSD, "Cgig2Aque_subset_JSD_sd.txt", sep ='\t')


# 8. Crassostrea gigas vs. Clytica hemisphaerica
# Number of orthologues = 4,691 orthologues
# 9 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/Quantile_transformation")
Cgig2Chem_orth_tpm_qn <- read.csv("Cgig2Chem_TPM_mean_quantile_transform.csv", header = T)
Chem2Cgig_orth_tpm_qn <- read.csv("Chem2Cgig_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/one2one")
both_ID <- read.table("Cgig2Chem.txt", header = T)
match_Cgig <- Cgig2Chem_orth_tpm_qn[match(both_ID$Cgig,Cgig2Chem_orth_tpm_qn$Gene_ID),]
match_Cgig <- na.omit(match_Cgig)
match_Chem <- Chem2Cgig_orth_tpm_qn[match(both_ID$Chem,Chem2Cgig_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/order")
write.table(match_Cgig,"order_Cgig2Chem_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Chem,"order_Chem2Cgig_tpm_qn.txt", sep="\t", quote=F)
match_Cgig <- read.table("order_Cgig2Chem_tpm_qn.txt", header = T)
match_Chem <- read.table("order_Chem2Cgig_tpm_qn.txt", header = T)

Cgig2Chem_full_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Cgig2Chem_full_JSD) <- colnames(Cgig2Chem_orth_tpm_qn)[-c(1)]
rownames(Cgig2Chem_full_JSD) <- colnames(Chem2Cgig_orth_tpm_qn)[-c(1)]
Cgig2Chem_mean_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Cgig2Chem_mean_JSD) <- colnames(Cgig2Chem_orth_tpm_qn)[-c(1)]
rownames(Cgig2Chem_mean_JSD) <- colnames(Chem2Cgig_orth_tpm_qn)[-c(1)]
Cgig2Chem_sd_JSD <- data.frame(replicate(ncol(match_Cgig)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Cgig2Chem_sd_JSD) <- colnames(Cgig2Chem_orth_tpm_qn)[-c(1)]
rownames(Cgig2Chem_sd_JSD) <- colnames(Chem2Cgig_orth_tpm_qn)[-c(1)]

for (i in colnames(Cgig2Chem_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Chem2Cgig_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cgig[[i]],match_Chem[[j]])
    Cgig2Chem_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cgig2Chem_mean_JSD[j,i] <- mean(all_JSD)
    Cgig2Chem_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/Results")
write.table(Cgig2Chem_full_JSD, "Cgig2Chem_full_set_JSD.txt", sep ='\t')
write.table(Cgig2Chem_mean_JSD, "Cgig2Chem_subset_JSD_mean.txt", sep ='\t')
write.table(Cgig2Chem_sd_JSD, "Cgig2Chem_subset_JSD_sd.txt", sep ='\t')


# 9. Re-import data from bootstrapped JSDs and plot
# non-comparable comparisons: heatmaps of RAW data

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/Results")

Cgig2Nvec_mean <- read.table("Cgig2Nvec_subset_JSD_mean.txt", header = T)
Cgig2Spur_mean <- read.table("Cgig2Spur_subset_JSD_mean.txt", header = T)
Cgig2Blan_mean <- read.table("Cgig2Blan_subset_JSD_mean.txt", header = T)
Cgig2Drer_mean <- read.table("Cgig2Drer_subset_JSD_mean.txt", header = T)
Cgig2Dmel_mean <- read.table("Cgig2Dmel_subset_JSD_mean.txt", header = T)
Cgig2Cele_mean <- read.table("Cgig2Cele_subset_JSD_mean.txt", header = T)
Cgig2Aque_mean <- read.table("Cgig2Aque_subset_JSD_mean.txt", header = T)
Cgig2Chem_mean <- read.table("Cgig2Chem_subset_JSD_mean.txt", header = T)

#allvalues <- rbind(Cgig2Nvec_mean,Cgig2Spur_mean,Cgig2Cgig_mean,Cgig2Blan_mean,Cgig2Drer_mean,Cgig2Dmel_mean,Cgig2Cele_mean,Cgig2Cgig_mean)
#minimum <- round_any(min(allvalues), 10, f = floor)
#maximum <- round_any(max(allvalues), 10, f = ceiling)

heatmap_color <- colorRampPalette(brewer.pal(n = 7, name = "RdGy"))(200)

h1 <- ComplexHeatmap::Heatmap(data.matrix(Cgig2Nvec_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Cgig2Nvec_mean)*unit(5,"mm"),
                              height=nrow(Cgig2Nvec_mean)*unit(5,"mm"),
                              col=heatmap_color,
                              column_labels = c("early morula","morula","blastula","rotary movement","free-swimming",
                                                "early gastrula","gastrula","trochophore 1","trochophore 2","trochophore 3",
                                                "trochophore 4","trochophore 5","early D-stage 1","early D-stage 2",
                                                "D-stage 1","D-stage 2","D-stage 3","D-stage 4","D-stage 5","D-stage 6",
                                                "D-stage 7","early umbo 1","early umbo 2","umbo 1","umbo 2","umbo 3",
                                                "umbo 4","umbo 5","umbo 6","late umbo 1","late umbo 2","pediveliger 1",
                                                "pediveliger 2","spat","juvenile"),
                              row_labels = c("24 hpf","48 hpf","72 hpf","96 hpf","120 hpf","144 hpf","168 hpf","192 hpf"),
                              heatmap_legend_param = list(title="JSD_raw"))

h2 <- ComplexHeatmap::Heatmap(data.matrix(Cgig2Spur_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Cgig2Spur_mean)*unit(5,"mm"),
                              height=nrow(Cgig2Spur_mean)*unit(5,"mm"),
                              col=heatmap_color,
                              column_labels = c("early morula","morula","blastula","rotary movement","free-swimming",
                                                "early gastrula","gastrula","trochophore 1","trochophore 2","trochophore 3",
                                                "trochophore 4","trochophore 5","early D-stage 1","early D-stage 2",
                                                "D-stage 1","D-stage 2","D-stage 3","D-stage 4","D-stage 5","D-stage 6",
                                                "D-stage 7","early umbo 1","early umbo 2","umbo 1","umbo 2","umbo 3",
                                                "umbo 4","umbo 5","umbo 6","late umbo 1","late umbo 2","pediveliger 1",
                                                "pediveliger 2","spat","juvenile"),
                              row_labels = c("hatched blastula","mesenchyme blastula","early gastrula","mid gastrula",
                                             "late gastrula","prism","late prism","pluteus","four-arm stage",
                                             "vestibular stage","pentagonal stage","tube foot stage",
                                             "post-metamorphic","juvenile"),
                              heatmap_legend_param = list(title="JSD_raw"))

h4 <- ComplexHeatmap::Heatmap(data.matrix(Cgig2Blan_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Cgig2Blan_mean)*unit(5,"mm"),
                              height=nrow(Cgig2Blan_mean)*unit(5,"mm"),
                              col=heatmap_color,
                              column_labels = c("early morula","morula","blastula","rotary movement","free-swimming",
                                                "early gastrula","gastrula","trochophore 1","trochophore 2","trochophore 3",
                                                "trochophore 4","trochophore 5","early D-stage 1","early D-stage 2",
                                                "D-stage 1","D-stage 2","D-stage 3","D-stage 4","D-stage 5","D-stage 6",
                                                "D-stage 7","early umbo 1","early umbo 2","umbo 1","umbo 2","umbo 3",
                                                "umbo 4","umbo 5","umbo 6","late umbo 1","late umbo 2","pediveliger 1",
                                                "pediveliger 2","spat","juvenile"),
                              row_labels = c("blastula","7 hpf","8 hpf","10 hpf","11 hpf","15 hpf","18 hpf",
                                             "21 hpf","24 hpf","27 hpf","36 hpf","50 hpf","60 hpf","pre-metamorphic"),
                              heatmap_legend_param = list(title="JSD_raw"))

h5 <- ComplexHeatmap::Heatmap(data.matrix(Cgig2Drer_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Cgig2Drer_mean)*unit(5,"mm"),
                              height=nrow(Cgig2Drer_mean)*unit(5,"mm"),
                              col=heatmap_color,
                              column_labels = c("early morula","morula","blastula","rotary movement","free-swimming",
                                                "early gastrula","gastrula","trochophore 1","trochophore 2","trochophore 3",
                                                "trochophore 4","trochophore 5","early D-stage 1","early D-stage 2",
                                                "D-stage 1","D-stage 2","D-stage 3","D-stage 4","D-stage 5","D-stage 6",
                                                "D-stage 7","early umbo 1","early umbo 2","umbo 1","umbo 2","umbo 3",
                                                "umbo 4","umbo 5","umbo 6","late umbo 1","late umbo 2","pediveliger 1",
                                                "pediveliger 2","spat","juvenile"),
                              row_labels = c("2 hpf", "1,000 cells", "Dome", "Shield", "6 hpf", "8 hpf",
                                             "Bud", "12 hpf", "16 hpf", "20 hpf", "26 hpf", "28 hpf", "30 hpf",
                                             "36 hpf", "2 dpf", "3 dpf", "5 dpf", "7 dpf"),
                              heatmap_legend_param = list(title="JSD_raw"))

h6 <- ComplexHeatmap::Heatmap(data.matrix(Cgig2Dmel_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Cgig2Dmel_mean)*unit(5,"mm"),
                              height=nrow(Cgig2Dmel_mean)*unit(5,"mm"),
                              col=heatmap_color,
                              column_labels = c("early morula","morula","blastula","rotary movement","free-swimming",
                                                "early gastrula","gastrula","trochophore 1","trochophore 2","trochophore 3",
                                                "trochophore 4","trochophore 5","early D-stage 1","early D-stage 2",
                                                "D-stage 1","D-stage 2","D-stage 3","D-stage 4","D-stage 5","D-stage 6",
                                                "D-stage 7","early umbo 1","early umbo 2","umbo 1","umbo 2","umbo 3",
                                                "umbo 4","umbo 5","umbo 6","late umbo 1","late umbo 2","pediveliger 1",
                                                "pediveliger 2","spat","juvenile"),
                              row_labels = c("0-2 hpf","2-4 hpf","4-8 hpf","8-12 hpf","12-16 hpf","16-20 hpf",
                                             "20-24 hpf","L1 larva","L2 larva","L3 larva 12 h","L3 larva puffstage 1-9",
                                             "white prepupa","late prepupa/early pupa","pupa 2-3 d","pupa 4 d","adult 1 d"),
                              heatmap_legend_param = list(title="JSD_raw"))

h7 <- ComplexHeatmap::Heatmap(data.matrix(Cgig2Cele_mean[-c(10,13),]),cluster_rows=F,cluster_columns=F,
                              width=ncol(Cgig2Cele_mean[-c(10,13),])*unit(5,"mm"),
                              height=nrow(Cgig2Cele_mean[-c(10,13),])*unit(5,"mm"),
                              col=heatmap_color,
                              column_labels = c("early morula","morula","blastula","rotary movement","free-swimming",
                                                "early gastrula","gastrula","trochophore 1","trochophore 2","trochophore 3",
                                                "trochophore 4","trochophore 5","early D-stage 1","early D-stage 2",
                                                "D-stage 1","D-stage 2","D-stage 3","D-stage 4","D-stage 5","D-stage 6",
                                                "D-stage 7","early umbo 1","early umbo 2","umbo 1","umbo 2","umbo 3",
                                                "umbo 4","umbo 5","umbo 6","late umbo 1","late umbo 2","pediveliger 1",
                                                "pediveliger 2","spat","juvenile"),
                              row_labels = c("0-90 mpf","120-180 mpf", "210-240 mpf", "270-330 mpf",
                                            "360-390 mpf","420-480 mpf","510-540 mpf","570-600 mpf","L1 larva",
                                            "L3 larva","L4 larva"),
                              heatmap_legend_param = list(title="JSD_raw"))

h9 <- ComplexHeatmap::Heatmap(data.matrix(Cgig2Aque_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Cgig2Aque_mean)*unit(5,"mm"),
                              height=nrow(Cgig2Aque_mean)*unit(5,"mm"),
                              col=heatmap_color,
                              column_labels = c("early morula","morula","blastula","rotary movement","free-swimming",
                                                "early gastrula","gastrula","trochophore 1","trochophore 2","trochophore 3",
                                                "trochophore 4","trochophore 5","early D-stage 1","early D-stage 2",
                                                "D-stage 1","D-stage 2","D-stage 3","D-stage 4","D-stage 5","D-stage 6",
                                                "D-stage 7","early umbo 1","early umbo 2","umbo 1","umbo 2","umbo 3",
                                                "umbo 4","umbo 5","umbo 6","late umbo 1","late umbo 2","pediveliger 1",
                                                "pediveliger 2","spat","juvenile"),
                              row_labels = c("cleavage","brown 24 hpf", "spot 48 hpf", "ring 72 hpf",
                                             "pre-competent larva", "competent larva", "early post-settlement",
                                             "post-settlement", "juvenile"),
                              heatmap_legend_param = list(title="JSD_raw"))

h10 <- ComplexHeatmap::Heatmap(data.matrix(Cgig2Chem_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Cgig2Chem_mean)*unit(5,"mm"),
                              height=nrow(Cgig2Chem_mean)*unit(5,"mm"),
                              col=heatmap_color,
                              column_labels = c("early morula","morula","blastula","rotary movement","free-swimming",
                                                "early gastrula","gastrula","trochophore 1","trochophore 2","trochophore 3",
                                                "trochophore 4","trochophore 5","early D-stage 1","early D-stage 2",
                                                "D-stage 1","D-stage 2","D-stage 3","D-stage 4","D-stage 5","D-stage 6",
                                                "D-stage 7","early umbo 1","early umbo 2","umbo 1","umbo 2","umbo 3",
                                                "umbo 4","umbo 5","umbo 6","late umbo 1","late umbo 2","pediveliger 1",
                                                "pediveliger 2","spat","juvenile"),
                              row_labels = c("early gastrula", "planula 24 hpf", "planula 48 hpf", "planula 72 hpf",
                                             "primary polyp", "gastrozooid", "gonozooid", "stolon", "baby medusa"),
                              heatmap_legend_param = list(title="JSD_raw"))

h1
h2
h4
h5
h6
h7
h9
h10