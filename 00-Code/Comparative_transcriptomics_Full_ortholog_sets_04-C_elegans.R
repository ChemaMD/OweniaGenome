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
## Caenorhabditis elegans 1-to-1 comparisons ##
##########################################

# 1. Caenorhabditis elegans vs. Nematostella vectensis
# Number of orthologues = 5,254 orthologues (-4)
# 8 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/Quantile_transformation")
Cele2Nvec_orth_tpm_qn <- read.csv("Cele2Nvec_TPM_mean_quantile_transform.csv", header = T)
Nvec2Cele_orth_tpm_qn <- read.csv("Nvec2Cele_TPM_mean_quantile_transform.csv", header = T)
Cele2Nvec_orth_tpm_qn$Gene_ID <- paste0("CELE_", Cele2Nvec_orth_tpm_qn$Gene_ID)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/one2one")
both_ID <- read.table("Cele2Nvec.txt", header = T)
match_Cele <- Cele2Nvec_orth_tpm_qn[match(both_ID$Cele,Cele2Nvec_orth_tpm_qn$Gene_ID),]
match_Nvec <- Nvec2Cele_orth_tpm_qn[match(both_ID$Nvec,Nvec2Cele_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/order")
write.table(match_Cele,"order_Cele2Nvec_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Nvec,"order_Nvec2Cele_tpm_qn.txt", sep="\t", quote=F)
match_Cele <- read.table("order_Cele2Nvec_tpm_qn.txt", header = T)
match_Nvec <- read.table("order_Nvec2Cele_tpm_qn.txt", header = T)

Cele2Nvec_full_JSD <- data.frame(replicate(ncol(match_Cele)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Cele2Nvec_full_JSD) <- colnames(Cele2Nvec_orth_tpm_qn)[-c(1)]
rownames(Cele2Nvec_full_JSD) <- colnames(Nvec2Cele_orth_tpm_qn)[-c(1)]
Cele2Nvec_mean_JSD <- data.frame(replicate(ncol(match_Cele)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Cele2Nvec_mean_JSD) <- colnames(Cele2Nvec_orth_tpm_qn)[-c(1)]
rownames(Cele2Nvec_mean_JSD) <- colnames(Nvec2Cele_orth_tpm_qn)[-c(1)]
Cele2Nvec_sd_JSD <- data.frame(replicate(ncol(match_Cele)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Cele2Nvec_sd_JSD) <- colnames(Cele2Nvec_orth_tpm_qn)[-c(1)]
rownames(Cele2Nvec_sd_JSD) <- colnames(Nvec2Cele_orth_tpm_qn)[-c(1)]

for (i in colnames(Cele2Nvec_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Nvec2Cele_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cele[[i]],match_Nvec[[j]])
    Cele2Nvec_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cele2Nvec_mean_JSD[j,i] <- mean(all_JSD)
    Cele2Nvec_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/Results")
write.table(Cele2Nvec_full_JSD, "Cele2Nvec_full_set_JSD.txt", sep ='\t')
write.table(Cele2Nvec_mean_JSD, "Cele2Nvec_subset_JSD_mean.txt", sep ='\t')
write.table(Cele2Nvec_sd_JSD, "Cele2Nvec_subset_JSD_sd.txt", sep ='\t')


# 2. Caenorhabditis elegans vs. Strongylocentrotus purpuratus
# Number of orthologues = 5,015 orthologues (+0)
# 14 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/Quantile_transformation")
Cele2Spur_orth_tpm_qn <- read.csv("Cele2Spur_TPM_mean_quantile_transform.csv", header = T)
Spur2Cele_orth_tpm_qn <- read.csv("Spur2Cele_TPM_mean_quantile_transform.csv", header = T)
Cele2Spur_orth_tpm_qn$Gene_ID <- paste0("CELE_", Cele2Spur_orth_tpm_qn$Gene_ID)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/one2one")
both_ID <- read.table("Cele2Spur.txt", header = T)
match_Cele <- Cele2Spur_orth_tpm_qn[match(both_ID$Cele,Cele2Spur_orth_tpm_qn$Gene_ID),]
match_Cele <- na.omit(match_Cele)
match_Spur <- Spur2Cele_orth_tpm_qn[match(both_ID$Spur,Spur2Cele_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/order")
write.table(match_Cele,"order_Cele2Spur_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Spur,"order_Spur2Cele_tpm_qn.txt", sep="\t", quote=F)
match_Cele <- read.table("order_Cele2Spur_tpm_qn.txt", header = T)
match_Spur <- read.table("order_Spur2Cele_tpm_qn.txt", header = T)

Cele2Spur_full_JSD <- data.frame(replicate(ncol(match_Cele)-1,sample(0:1,ncol(match_Spur)-1,rep=TRUE)))
colnames(Cele2Spur_full_JSD) <- colnames(Cele2Spur_orth_tpm_qn)[-c(1)]
rownames(Cele2Spur_full_JSD) <- colnames(Spur2Cele_orth_tpm_qn)[-c(1)]
Cele2Spur_mean_JSD <- data.frame(replicate(ncol(match_Cele)-1,sample(0:1,ncol(match_Spur)-1,rep=TRUE)))
colnames(Cele2Spur_mean_JSD) <- colnames(Cele2Spur_orth_tpm_qn)[-c(1)]
rownames(Cele2Spur_mean_JSD) <- colnames(Spur2Cele_orth_tpm_qn)[-c(1)]
Cele2Spur_sd_JSD <- data.frame(replicate(ncol(match_Cele)-1,sample(0:1,ncol(match_Spur)-1,rep=TRUE)))
colnames(Cele2Spur_sd_JSD) <- colnames(Cele2Spur_orth_tpm_qn)[-c(1)]
rownames(Cele2Spur_sd_JSD) <- colnames(Spur2Cele_orth_tpm_qn)[-c(1)]

for (i in colnames(Cele2Spur_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Spur2Cele_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cele[[i]],match_Spur[[j]])
    Cele2Spur_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cele2Spur_mean_JSD[j,i] <- mean(all_JSD)
    Cele2Spur_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/Results")
write.table(Cele2Spur_full_JSD, "Cele2Spur_full_set_JSD.txt", sep ='\t')
write.table(Cele2Spur_mean_JSD, "Cele2Spur_subset_JSD_mean.txt", sep ='\t')
write.table(Cele2Spur_sd_JSD, "Cele2Spur_subset_JSD_sd.txt", sep ='\t')


# 3. Caenorhabditis elegans vs. Branchiostoma lanceolatum
# Number of orthologues = 6,673 orthologues (+39)
# 14 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/Quantile_transformation")
Cele2Blan_orth_tpm_qn <- read.csv("Cele2Blan_TPM_mean_quantile_transform.csv", header = T)
Blan2Cele_orth_tpm_qn <- read.csv("Blan2Cele_TPM_mean_quantile_transform.csv", header = T)
Cele2Blan_orth_tpm_qn$Gene_ID <- paste0("CELE_", Cele2Blan_orth_tpm_qn$Gene_ID)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/one2one")
both_ID <- read.table("Cele2Blan.txt", header = T)
match_Cele <- Cele2Blan_orth_tpm_qn[match(both_ID$Cele,Cele2Blan_orth_tpm_qn$Gene_ID),]
match_Cele <- na.omit(match_Cele)
match_Blan <- Blan2Cele_orth_tpm_qn[match(both_ID$Blan,Blan2Cele_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/order")
write.table(match_Cele,"order_Cele2Blan_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Blan,"order_Blan2Cele_tpm_qn.txt", sep="\t", quote=F)
match_Cele <- read.table("order_Cele2Blan_tpm_qn.txt", header = T)
match_Blan <- read.table("order_Blan2Cele_tpm_qn.txt", header = T)

Cele2Blan_full_JSD <- data.frame(replicate(ncol(match_Cele)-1,sample(0:1,ncol(match_Blan)-1,rep=TRUE)))
colnames(Cele2Blan_full_JSD) <- colnames(Cele2Blan_orth_tpm_qn)[-c(1)]
rownames(Cele2Blan_full_JSD) <- colnames(Blan2Cele_orth_tpm_qn)[-c(1)]
Cele2Blan_mean_JSD <- data.frame(replicate(ncol(match_Cele)-1,sample(0:1,ncol(match_Blan)-1,rep=TRUE)))
colnames(Cele2Blan_mean_JSD) <- colnames(Cele2Blan_orth_tpm_qn)[-c(1)]
rownames(Cele2Blan_mean_JSD) <- colnames(Blan2Cele_orth_tpm_qn)[-c(1)]
Cele2Blan_sd_JSD <- data.frame(replicate(ncol(match_Cele)-1,sample(0:1,ncol(match_Blan)-1,rep=TRUE)))
colnames(Cele2Blan_sd_JSD) <- colnames(Cele2Blan_orth_tpm_qn)[-c(1)]
rownames(Cele2Blan_sd_JSD) <- colnames(Blan2Cele_orth_tpm_qn)[-c(1)]

for (i in colnames(Cele2Blan_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Blan2Cele_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cele[[i]],match_Blan[[j]])
    Cele2Blan_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cele2Blan_mean_JSD[j,i] <- mean(all_JSD)
    Cele2Blan_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/Results")
write.table(Cele2Blan_full_JSD, "Cele2Blan_full_set_JSD.txt", sep ='\t')
write.table(Cele2Blan_mean_JSD, "Cele2Blan_subset_JSD_mean.txt", sep ='\t')
write.table(Cele2Blan_sd_JSD, "Cele2Blan_subset_JSD_sd.txt", sep ='\t')


# 4. Caenorhabditis elegans vs. Danio rerio
# Number of orthologues = 4,316 orthologues (+21)
# 7 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/Quantile_transformation")
Cele2Drer_orth_tpm_qn <- read.csv("Cele2Drer_TPM_mean_quantile_transform.csv", header = T)
Drer2Cele_orth_tpm_qn <- read.csv("Drer2Cele_TPM_mean_quantile_transform.csv", header = T)
Cele2Drer_orth_tpm_qn$Gene_ID <- paste0("CELE_", Cele2Drer_orth_tpm_qn$Gene_ID)


setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/one2one")
both_ID <- read.table("Cele2Drer.txt", header = T)
match_Cele <- Cele2Drer_orth_tpm_qn[match(both_ID$Cele,Cele2Drer_orth_tpm_qn$Gene_ID),]
match_Cele <- na.omit(match_Cele)
match_Drer <- Drer2Cele_orth_tpm_qn[match(both_ID$Drer,Drer2Cele_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/order")
write.table(match_Cele,"order_Cele2Drer_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Drer,"order_Drer2Cele_tpm_qn.txt", sep="\t", quote=F)
match_Cele <- read.table("order_Cele2Drer_tpm_qn.txt", header = T)
match_Drer <- read.table("order_Drer2Cele_tpm_qn.txt", header = T)

Cele2Drer_full_JSD <- data.frame(replicate(ncol(match_Cele)-1,sample(0:1,ncol(match_Drer)-1,rep=TRUE)))
colnames(Cele2Drer_full_JSD) <- colnames(Cele2Drer_orth_tpm_qn)[-c(1)]
rownames(Cele2Drer_full_JSD) <- colnames(Drer2Cele_orth_tpm_qn)[-c(1)]
Cele2Drer_mean_JSD <- data.frame(replicate(ncol(match_Cele)-1,sample(0:1,ncol(match_Drer)-1,rep=TRUE)))
colnames(Cele2Drer_mean_JSD) <- colnames(Cele2Drer_orth_tpm_qn)[-c(1)]
rownames(Cele2Drer_mean_JSD) <- colnames(Drer2Cele_orth_tpm_qn)[-c(1)]
Cele2Drer_sd_JSD <- data.frame(replicate(ncol(match_Cele)-1,sample(0:1,ncol(match_Drer)-1,rep=TRUE)))
colnames(Cele2Drer_sd_JSD) <- colnames(Cele2Drer_orth_tpm_qn)[-c(1)]
rownames(Cele2Drer_sd_JSD) <- colnames(Drer2Cele_orth_tpm_qn)[-c(1)]

for (i in colnames(Cele2Drer_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Drer2Cele_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cele[[i]],match_Drer[[j]])
    Cele2Drer_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cele2Drer_mean_JSD[j,i] <- mean(all_JSD)
    Cele2Drer_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/Results")
write.table(Cele2Drer_full_JSD, "Cele2Drer_full_set_JSD.txt", sep ='\t')
write.table(Cele2Drer_mean_JSD, "Cele2Drer_subset_JSD_mean.txt", sep ='\t')
write.table(Cele2Drer_sd_JSD, "Cele2Drer_subset_JSD_sd.txt", sep ='\t')


# 5. Caenorhabditis elegans vs. Drosophila melanogaster
# Number of orthologues = 4,635 orthologues (+45)
# 18 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/Quantile_transformation")
Cele2Dmel_orth_tpm_qn <- read.csv("Cele2Dmel_TPM_mean_quantile_transform.csv", header = T)
Dmel2Cele_orth_tpm_qn <- read.csv("Dmel2Cele_TPM_mean_quantile_transform.csv", header = T)
Cele2Dmel_orth_tpm_qn$Gene_ID <- paste0("CELE_", Cele2Dmel_orth_tpm_qn$Gene_ID)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/one2one")
both_ID <- read.table("Cele2Dmel.txt", header = T)
match_Cele <- Cele2Dmel_orth_tpm_qn[match(both_ID$Cele,Cele2Dmel_orth_tpm_qn$Gene_ID),]
match_Cele <- na.omit(match_Cele)
match_Dmel <- Dmel2Cele_orth_tpm_qn[match(both_ID$Dmel,Dmel2Cele_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/order")
write.table(match_Cele,"order_Cele2Dmel_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Dmel,"order_Dmel2Cele_tpm_qn.txt", sep="\t", quote=F)
match_Cele <- read.table("order_Cele2Dmel_tpm_qn.txt", header = T)
match_Dmel <- read.table("order_Dmel2Cele_tpm_qn.txt", header = T)

Cele2Dmel_full_JSD <- data.frame(replicate(ncol(match_Cele)-1,sample(0:1,ncol(match_Dmel)-1,rep=TRUE)))
colnames(Cele2Dmel_full_JSD) <- colnames(Cele2Dmel_orth_tpm_qn)[-c(1)]
rownames(Cele2Dmel_full_JSD) <- colnames(Dmel2Cele_orth_tpm_qn)[-c(1)]
Cele2Dmel_mean_JSD <- data.frame(replicate(ncol(match_Cele)-1,sample(0:1,ncol(match_Dmel)-1,rep=TRUE)))
colnames(Cele2Dmel_mean_JSD) <- colnames(Cele2Dmel_orth_tpm_qn)[-c(1)]
rownames(Cele2Dmel_mean_JSD) <- colnames(Dmel2Cele_orth_tpm_qn)[-c(1)]
Cele2Dmel_sd_JSD <- data.frame(replicate(ncol(match_Cele)-1,sample(0:1,ncol(match_Dmel)-1,rep=TRUE)))
colnames(Cele2Dmel_sd_JSD) <- colnames(Cele2Dmel_orth_tpm_qn)[-c(1)]
rownames(Cele2Dmel_sd_JSD) <- colnames(Dmel2Cele_orth_tpm_qn)[-c(1)]

for (i in colnames(Cele2Dmel_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Dmel2Cele_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cele[[i]],match_Dmel[[j]])
    Cele2Dmel_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cele2Dmel_mean_JSD[j,i] <- mean(all_JSD)
    Cele2Dmel_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/Results")
write.table(Cele2Dmel_full_JSD, "Cele2Dmel_full_set_JSD.txt", sep ='\t')
write.table(Cele2Dmel_mean_JSD, "Cele2Dmel_subset_JSD_mean.txt", sep ='\t')
write.table(Cele2Dmel_sd_JSD, "Cele2Dmel_subset_JSD_sd.txt", sep ='\t')


# 7. Caenorhabditis elegans vs. Amphimedon queenslandica
# Number of orthologues = 3,962 orthologues
# 9 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/Quantile_transformation")
Cele2Aque_orth_tpm_qn <- read.csv("Cele2Aque_TPM_mean_quantile_transform.csv", header = T)
Aque2Cele_orth_tpm_qn <- read.csv("Aque2Cele_TPM_mean_quantile_transform.csv", header = T)
Cele2Aque_orth_tpm_qn$Gene_ID <- paste0("CELE_", Cele2Aque_orth_tpm_qn$Gene_ID)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/one2one")
both_ID <- read.table("Cele2Aque.txt", header = T)
match_Cele <- Cele2Aque_orth_tpm_qn[match(both_ID$Cele,Cele2Aque_orth_tpm_qn$Gene_ID),]
match_Cele <- na.omit(match_Cele)
match_Aque <- Aque2Cele_orth_tpm_qn[match(both_ID$Aque,Aque2Cele_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/order")
write.table(match_Cele,"order_Cele2Aque_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Aque,"order_Aque2Cele_tpm_qn.txt", sep="\t", quote=F)
match_Cele <- read.table("order_Cele2Aque_tpm_qn.txt", header = T)
match_Aque <- read.table("order_Aque2Cele_tpm_qn.txt", header = T)

Cele2Aque_full_JSD <- data.frame(replicate(ncol(match_Cele)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Cele2Aque_full_JSD) <- colnames(Cele2Aque_orth_tpm_qn)[-c(1)]
rownames(Cele2Aque_full_JSD) <- colnames(Aque2Cele_orth_tpm_qn)[-c(1)]
Cele2Aque_mean_JSD <- data.frame(replicate(ncol(match_Cele)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Cele2Aque_mean_JSD) <- colnames(Cele2Aque_orth_tpm_qn)[-c(1)]
rownames(Cele2Aque_mean_JSD) <- colnames(Aque2Cele_orth_tpm_qn)[-c(1)]
Cele2Aque_sd_JSD <- data.frame(replicate(ncol(match_Cele)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Cele2Aque_sd_JSD) <- colnames(Cele2Aque_orth_tpm_qn)[-c(1)]
rownames(Cele2Aque_sd_JSD) <- colnames(Aque2Cele_orth_tpm_qn)[-c(1)]

for (i in colnames(Cele2Aque_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Aque2Cele_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cele[[i]],match_Aque[[j]])
    Cele2Aque_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cele2Aque_mean_JSD[j,i] <- mean(all_JSD)
    Cele2Aque_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/Results")
write.table(Cele2Aque_full_JSD, "Cele2Aque_full_set_JSD.txt", sep ='\t')
write.table(Cele2Aque_mean_JSD, "Cele2Aque_subset_JSD_mean.txt", sep ='\t')
write.table(Cele2Aque_sd_JSD, "Cele2Aque_subset_JSD_sd.txt", sep ='\t')


# 8. Caenorhabditis elegans vs. Clytica hemisphaerica
# Number of orthologues = 4,691 orthologues
# 9 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/Quantile_transformation")
Cele2Chem_orth_tpm_qn <- read.csv("Cele2Chem_TPM_mean_quantile_transform.csv", header = T)
Chem2Cele_orth_tpm_qn <- read.csv("Chem2Cele_TPM_mean_quantile_transform.csv", header = T)
Cele2Chem_orth_tpm_qn$Gene_ID <- paste0("CELE_", Cele2Chem_orth_tpm_qn$Gene_ID)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/one2one")
both_ID <- read.table("Cele2Chem.txt", header = T)
match_Cele <- Cele2Chem_orth_tpm_qn[match(both_ID$Cele,Cele2Chem_orth_tpm_qn$Gene_ID),]
match_Cele <- na.omit(match_Cele)
match_Chem <- Chem2Cele_orth_tpm_qn[match(both_ID$Chem,Chem2Cele_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/order")
write.table(match_Cele,"order_Cele2Chem_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Chem,"order_Chem2Cele_tpm_qn.txt", sep="\t", quote=F)
match_Cele <- read.table("order_Cele2Chem_tpm_qn.txt", header = T)
match_Chem <- read.table("order_Chem2Cele_tpm_qn.txt", header = T)

Cele2Chem_full_JSD <- data.frame(replicate(ncol(match_Cele)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Cele2Chem_full_JSD) <- colnames(Cele2Chem_orth_tpm_qn)[-c(1)]
rownames(Cele2Chem_full_JSD) <- colnames(Chem2Cele_orth_tpm_qn)[-c(1)]
Cele2Chem_mean_JSD <- data.frame(replicate(ncol(match_Cele)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Cele2Chem_mean_JSD) <- colnames(Cele2Chem_orth_tpm_qn)[-c(1)]
rownames(Cele2Chem_mean_JSD) <- colnames(Chem2Cele_orth_tpm_qn)[-c(1)]
Cele2Chem_sd_JSD <- data.frame(replicate(ncol(match_Cele)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Cele2Chem_sd_JSD) <- colnames(Cele2Chem_orth_tpm_qn)[-c(1)]
rownames(Cele2Chem_sd_JSD) <- colnames(Chem2Cele_orth_tpm_qn)[-c(1)]

for (i in colnames(Cele2Chem_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Chem2Cele_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Cele[[i]],match_Chem[[j]])
    Cele2Chem_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Cele2Chem_mean_JSD[j,i] <- mean(all_JSD)
    Cele2Chem_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/Results")
write.table(Cele2Chem_full_JSD, "Cele2Chem_full_set_JSD.txt", sep ='\t')
write.table(Cele2Chem_mean_JSD, "Cele2Chem_subset_JSD_mean.txt", sep ='\t')
write.table(Cele2Chem_sd_JSD, "Cele2Chem_subset_JSD_sd.txt", sep ='\t')