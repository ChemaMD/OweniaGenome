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
## Owenia fusiformis 1-to-1 comparisons ##
##########################################

# 1. Owenia fusiformis vs. Nematostella vectensis
# Number of orthologues = 5,254 orthologues (-4)
# 8 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Quantile_transformation")
Ofus2Nvec_orth_tpm_qn <- read.csv("Ofus2Nvec_TPM_mean_quantile_transform.csv", header = T)
Nvec2Ofus_orth_tpm_qn <- read.csv("Nvec2Ofus_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/one2one")
both_ID <- read.table("Ofus2Nvec.txt", header = T)
match_Ofus <- Ofus2Nvec_orth_tpm_qn[match(both_ID$Ofus,Ofus2Nvec_orth_tpm_qn$Gene_ID),]
match_Nvec <- Nvec2Ofus_orth_tpm_qn[match(both_ID$Nvec,Nvec2Ofus_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/order")
write.table(match_Ofus,"order_Ofus2Nvec_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Nvec,"order_Nvec2Ofus_tpm_qn.txt", sep="\t", quote=F)
match_Ofus <- read.table("order_Ofus2Nvec_tpm_qn.txt", header = T)
match_Nvec <- read.table("order_Nvec2Ofus_tpm_qn.txt", header = T)

Ofus2Nvec_full_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Ofus2Nvec_full_JSD) <- colnames(Ofus2Nvec_orth_tpm_qn)[-c(1)]
rownames(Ofus2Nvec_full_JSD) <- colnames(Nvec2Ofus_orth_tpm_qn)[-c(1)]
Ofus2Nvec_mean_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Ofus2Nvec_mean_JSD) <- colnames(Ofus2Nvec_orth_tpm_qn)[-c(1)]
rownames(Ofus2Nvec_mean_JSD) <- colnames(Nvec2Ofus_orth_tpm_qn)[-c(1)]
Ofus2Nvec_sd_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Ofus2Nvec_sd_JSD) <- colnames(Ofus2Nvec_orth_tpm_qn)[-c(1)]
rownames(Ofus2Nvec_sd_JSD) <- colnames(Nvec2Ofus_orth_tpm_qn)[-c(1)]

for (i in colnames(Ofus2Nvec_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Nvec2Ofus_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Ofus[[i]],match_Nvec[[j]])
    Ofus2Nvec_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Ofus2Nvec_mean_JSD[j,i] <- mean(all_JSD)
    Ofus2Nvec_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Results")
write.table(Ofus2Nvec_full_JSD, "Ofus2Nvec_full_set_JSD.txt", sep ='\t')
write.table(Ofus2Nvec_mean_JSD, "Ofus2Nvec_subset_JSD_mean.txt", sep ='\t')
write.table(Ofus2Nvec_sd_JSD, "Ofus2Nvec_subset_JSD_sd.txt", sep ='\t')


# 2. Owenia fusiformis vs. Strongylocentrotus purpuratus
# Number of orthologues = 5,015 orthologues (+0)
# 14 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Quantile_transformation")
Ofus2Spur_orth_tpm_qn <- read.csv("Ofus2Spur_TPM_mean_quantile_transform.csv", header = T)
Spur2Ofus_orth_tpm_qn <- read.csv("Spur2Ofus_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/one2one")
both_ID <- read.table("Ofus2Spur.txt", header = T)
match_Ofus <- Ofus2Spur_orth_tpm_qn[match(both_ID$Ofus,Ofus2Spur_orth_tpm_qn$Gene_ID),]
match_Spur <- Spur2Ofus_orth_tpm_qn[match(both_ID$Spur,Spur2Ofus_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/order")
write.table(match_Ofus,"order_Ofus2Spur_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Spur,"order_Spur2Ofus_tpm_qn.txt", sep="\t", quote=F)
match_Ofus <- read.table("order_Ofus2Spur_tpm_qn.txt", header = T)
match_Spur <- read.table("order_Spur2Ofus_tpm_qn.txt", header = T)

Ofus2Spur_full_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Spur)-1,rep=TRUE)))
colnames(Ofus2Spur_full_JSD) <- colnames(Ofus2Spur_orth_tpm_qn)[-c(1)]
rownames(Ofus2Spur_full_JSD) <- colnames(Spur2Ofus_orth_tpm_qn)[-c(1)]
Ofus2Spur_mean_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Spur)-1,rep=TRUE)))
colnames(Ofus2Spur_mean_JSD) <- colnames(Ofus2Spur_orth_tpm_qn)[-c(1)]
rownames(Ofus2Spur_mean_JSD) <- colnames(Spur2Ofus_orth_tpm_qn)[-c(1)]
Ofus2Spur_sd_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Spur)-1,rep=TRUE)))
colnames(Ofus2Spur_sd_JSD) <- colnames(Ofus2Spur_orth_tpm_qn)[-c(1)]
rownames(Ofus2Spur_sd_JSD) <- colnames(Spur2Ofus_orth_tpm_qn)[-c(1)]

for (i in colnames(Ofus2Spur_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Spur2Ofus_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Ofus[[i]],match_Spur[[j]])
    Ofus2Spur_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Ofus2Spur_mean_JSD[j,i] <- mean(all_JSD)
    Ofus2Spur_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Results")
write.table(Ofus2Spur_full_JSD, "Ofus2Spur_full_set_JSD.txt", sep ='\t')
write.table(Ofus2Spur_mean_JSD, "Ofus2Spur_subset_JSD_mean.txt", sep ='\t')
write.table(Ofus2Spur_sd_JSD, "Ofus2Spur_subset_JSD_sd.txt", sep ='\t')


# 3. Owenia fusiformis vs. Crassostrea gigas
# Number of orthologues = 6,737 orthologues (+62)
# 35 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Quantile_transformation")
Ofus2Cgig_orth_tpm_qn <- read.csv("Ofus2Cgig_TPM_mean_quantile_transform.csv", header = T)
Cgig2Ofus_orth_tpm_qn <- read.csv("Cgig2Ofus_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/one2one")
both_ID <- read.table("Ofus2Cgig.txt", header = T)
match_Ofus <- Ofus2Cgig_orth_tpm_qn[match(both_ID$Ofus,Ofus2Cgig_orth_tpm_qn$Gene_ID),]
match_Cgig <- Cgig2Ofus_orth_tpm_qn[match(both_ID$Cgig,Cgig2Ofus_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/order")
write.table(match_Ofus,"order_Ofus2Cgig_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Cgig,"order_Cgig2Ofus_tpm_qn.txt", sep="\t", quote=F)
match_Ofus <- read.table("order_Ofus2Cgig_tpm_qn.txt", header = T)
match_Cgig <- read.table("order_Cgig2Ofus_tpm_qn.txt", header = T)

Ofus2Cgig_full_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Cgig)-1,rep=TRUE)))
colnames(Ofus2Cgig_full_JSD) <- colnames(Ofus2Cgig_orth_tpm_qn)[-c(1)]
rownames(Ofus2Cgig_full_JSD) <- colnames(Cgig2Ofus_orth_tpm_qn)[-c(1)]
Ofus2Cgig_mean_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Cgig)-1,rep=TRUE)))
colnames(Ofus2Cgig_mean_JSD) <- colnames(Ofus2Cgig_orth_tpm_qn)[-c(1)]
rownames(Ofus2Cgig_mean_JSD) <- colnames(Cgig2Ofus_orth_tpm_qn)[-c(1)]
Ofus2Cgig_sd_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Cgig)-1,rep=TRUE)))
colnames(Ofus2Cgig_sd_JSD) <- colnames(Ofus2Cgig_orth_tpm_qn)[-c(1)]
rownames(Ofus2Cgig_sd_JSD) <- colnames(Cgig2Ofus_orth_tpm_qn)[-c(1)]

for (i in colnames(Ofus2Cgig_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Cgig2Ofus_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Ofus[[i]],match_Cgig[[j]])
    Ofus2Cgig_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Ofus2Cgig_mean_JSD[j,i] <- mean(all_JSD)
    Ofus2Cgig_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Results")
write.table(Ofus2Cgig_full_JSD, "Ofus2Cgig_full_set_JSD.txt", sep ='\t')
write.table(Ofus2Cgig_mean_JSD, "Ofus2Cgig_subset_JSD_mean.txt", sep ='\t')
write.table(Ofus2Cgig_sd_JSD, "Ofus2Cgig_subset_JSD_sd.txt", sep ='\t')


# 4. Owenia fusiformis vs. Branchiostoma lanceolatum
# Number of orthologues = 6,673 orthologues (+39)
# 14 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Quantile_transformation")
Ofus2Blan_orth_tpm_qn <- read.csv("Ofus2Blan_TPM_mean_quantile_transform.csv", header = T)
Blan2Ofus_orth_tpm_qn <- read.csv("Blan2Ofus_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/one2one")
both_ID <- read.table("Ofus2Blan.txt", header = T)
match_Ofus <- Ofus2Blan_orth_tpm_qn[match(both_ID$Ofus,Ofus2Blan_orth_tpm_qn$Gene_ID),]
match_Blan <- Blan2Ofus_orth_tpm_qn[match(both_ID$Blan,Blan2Ofus_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/order")
write.table(match_Ofus,"order_Ofus2Blan_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Blan,"order_Blan2Ofus_tpm_qn.txt", sep="\t", quote=F)
match_Ofus <- read.table("order_Ofus2Blan_tpm_qn.txt", header = T)
match_Blan <- read.table("order_Blan2Ofus_tpm_qn.txt", header = T)

Ofus2Blan_full_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Blan)-1,rep=TRUE)))
colnames(Ofus2Blan_full_JSD) <- colnames(Ofus2Blan_orth_tpm_qn)[-c(1)]
rownames(Ofus2Blan_full_JSD) <- colnames(Blan2Ofus_orth_tpm_qn)[-c(1)]
Ofus2Blan_mean_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Blan)-1,rep=TRUE)))
colnames(Ofus2Blan_mean_JSD) <- colnames(Ofus2Blan_orth_tpm_qn)[-c(1)]
rownames(Ofus2Blan_mean_JSD) <- colnames(Blan2Ofus_orth_tpm_qn)[-c(1)]
Ofus2Blan_sd_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Blan)-1,rep=TRUE)))
colnames(Ofus2Blan_sd_JSD) <- colnames(Ofus2Blan_orth_tpm_qn)[-c(1)]
rownames(Ofus2Blan_sd_JSD) <- colnames(Blan2Ofus_orth_tpm_qn)[-c(1)]

for (i in colnames(Ofus2Blan_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Blan2Ofus_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Ofus[[i]],match_Blan[[j]])
    Ofus2Blan_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Ofus2Blan_mean_JSD[j,i] <- mean(all_JSD)
    Ofus2Blan_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Results")
write.table(Ofus2Blan_full_JSD, "Ofus2Blan_full_set_JSD.txt", sep ='\t')
write.table(Ofus2Blan_mean_JSD, "Ofus2Blan_subset_JSD_mean.txt", sep ='\t')
write.table(Ofus2Blan_sd_JSD, "Ofus2Blan_subset_JSD_sd.txt", sep ='\t')


# 5. Owenia fusiformis vs. Danio rerio
# Number of orthologues = 4,316 orthologues (+21)
# 7 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Quantile_transformation")
Ofus2Drer_orth_tpm_qn <- read.csv("Ofus2Drer_TPM_mean_quantile_transform.csv", header = T)
Drer2Ofus_orth_tpm_qn <- read.csv("Drer2Ofus_TPM_mean_quantile_transform.csv", header = T)
#Drer2Ofus_orth_tpm_qn <- cbind(Drer2Ofus_orth_tpm_qn_preliminary[1],Drer2Ofus_orth_tpm_qn_preliminary[12],rer2Ofus_orth_tpm_qn_preliminary[c(3:11)],Drer2Ofus_orth_tpm_qn_preliminary[2],Drer2Ofus_orth_tpm_qn_preliminary[c(13:19)]) ## not needed anymore

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/one2one")
both_ID <- read.table("Ofus2Drer.txt", header = T)
match_Ofus <- Ofus2Drer_orth_tpm_qn[match(both_ID$Ofus,Ofus2Drer_orth_tpm_qn$Gene_ID),]
match_Drer <- Drer2Ofus_orth_tpm_qn[match(both_ID$Drer,Drer2Ofus_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/order")
write.table(match_Ofus,"order_Ofus2Drer_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Drer,"order_Drer2Ofus_tpm_qn.txt", sep="\t", quote=F)
match_Ofus <- read.table("order_Ofus2Drer_tpm_qn.txt", header = T)
match_Drer <- read.table("order_Drer2Ofus_tpm_qn.txt", header = T)

Ofus2Drer_full_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Drer)-1,rep=TRUE)))
colnames(Ofus2Drer_full_JSD) <- colnames(Ofus2Drer_orth_tpm_qn)[-c(1)]
rownames(Ofus2Drer_full_JSD) <- colnames(Drer2Ofus_orth_tpm_qn)[-c(1)]
Ofus2Drer_mean_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Drer)-1,rep=TRUE)))
colnames(Ofus2Drer_mean_JSD) <- colnames(Ofus2Drer_orth_tpm_qn)[-c(1)]
rownames(Ofus2Drer_mean_JSD) <- colnames(Drer2Ofus_orth_tpm_qn)[-c(1)]
Ofus2Drer_sd_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Drer)-1,rep=TRUE)))
colnames(Ofus2Drer_sd_JSD) <- colnames(Ofus2Drer_orth_tpm_qn)[-c(1)]
rownames(Ofus2Drer_sd_JSD) <- colnames(Drer2Ofus_orth_tpm_qn)[-c(1)]

for (i in colnames(Ofus2Drer_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Drer2Ofus_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Ofus[[i]],match_Drer[[j]])
    Ofus2Drer_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Ofus2Drer_mean_JSD[j,i] <- mean(all_JSD)
    Ofus2Drer_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Results")
write.table(Ofus2Drer_full_JSD, "Ofus2Drer_full_set_JSD.txt", sep ='\t')
write.table(Ofus2Drer_mean_JSD, "Ofus2Drer_subset_JSD_mean.txt", sep ='\t')
write.table(Ofus2Drer_sd_JSD, "Ofus2Drer_subset_JSD_sd.txt", sep ='\t')


# 6. Owenia fusiformis vs. Drosophila melanogaster
# Number of orthologues = 4,635 orthologues (+45)
# 18 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Quantile_transformation")
Ofus2Dmel_orth_tpm_qn <- read.csv("Ofus2Dmel_TPM_mean_quantile_transform.csv", header = T)
Dmel2Ofus_orth_tpm_qn <- read.csv("Dmel2Ofus_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/one2one")
both_ID <- read.table("Ofus2Dmel.txt", header = T)
match_Ofus <- Ofus2Dmel_orth_tpm_qn[match(both_ID$Ofus,Ofus2Dmel_orth_tpm_qn$Gene_ID),]
match_Dmel <- Dmel2Ofus_orth_tpm_qn[match(both_ID$Dmel,Dmel2Ofus_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/order")
write.table(match_Ofus,"order_Ofus2Dmel_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Dmel,"order_Dmel2Ofus_tpm_qn.txt", sep="\t", quote=F)
match_Ofus <- read.table("order_Ofus2Dmel_tpm_qn.txt", header = T)
match_Dmel <- read.table("order_Dmel2Ofus_tpm_qn.txt", header = T)

Ofus2Dmel_full_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Dmel)-1,rep=TRUE)))
colnames(Ofus2Dmel_full_JSD) <- colnames(Ofus2Dmel_orth_tpm_qn)[-c(1)]
rownames(Ofus2Dmel_full_JSD) <- colnames(Dmel2Ofus_orth_tpm_qn)[-c(1)]
Ofus2Dmel_mean_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Dmel)-1,rep=TRUE)))
colnames(Ofus2Dmel_mean_JSD) <- colnames(Ofus2Dmel_orth_tpm_qn)[-c(1)]
rownames(Ofus2Dmel_mean_JSD) <- colnames(Dmel2Ofus_orth_tpm_qn)[-c(1)]
Ofus2Dmel_sd_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Dmel)-1,rep=TRUE)))
colnames(Ofus2Dmel_sd_JSD) <- colnames(Ofus2Dmel_orth_tpm_qn)[-c(1)]
rownames(Ofus2Dmel_sd_JSD) <- colnames(Dmel2Ofus_orth_tpm_qn)[-c(1)]

for (i in colnames(Ofus2Dmel_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Dmel2Ofus_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Ofus[[i]],match_Dmel[[j]])
    Ofus2Dmel_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Ofus2Dmel_mean_JSD[j,i] <- mean(all_JSD)
    Ofus2Dmel_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Results")
write.table(Ofus2Dmel_full_JSD, "Ofus2Dmel_full_set_JSD.txt", sep ='\t')
write.table(Ofus2Dmel_mean_JSD, "Ofus2Dmel_subset_JSD_mean.txt", sep ='\t')
write.table(Ofus2Dmel_sd_JSD, "Ofus2Dmel_subset_JSD_sd.txt", sep ='\t')


# 7. Owenia fusiformis vs. Caenorhabidits elegans
# Number of orthologues = 3,767 orthologues (-43)
# 13 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Quantile_transformation")
Ofus2Cele_orth_tpm_qn <- read.csv("Ofus2Cele_TPM_mean_quantile_transform.csv", header = T)
Cele2Ofus_orth_tpm_qn <- read.csv("Cele2Ofus_TPM_mean_quantile_transform.csv", header = T)
Cele2Ofus_orth_tpm_qn$Gene_ID <- paste0("CELE_", Cele2Ofus_orth_tpm_qn$Gene_ID)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/one2one")
both_ID <- read.table("Ofus2Cele.txt", header = T)
match_Ofus <- Ofus2Cele_orth_tpm_qn[match(both_ID$Ofus,Ofus2Cele_orth_tpm_qn$Gene_ID),]
match_Cele <- Cele2Ofus_orth_tpm_qn[match(both_ID$Cele,Cele2Ofus_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/order")
write.table(match_Ofus,"order_Ofus2Cele_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Cele,"order_Cele2Ofus_tpm_qn.txt", sep="\t", quote=F)
match_Ofus <- read.table("order_Ofus2Cele_tpm_qn.txt", header = T)
match_Cele <- read.table("order_Cele2Ofus_tpm_qn.txt", header = T)

Ofus2Cele_full_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Cele)-1,rep=TRUE)))
colnames(Ofus2Cele_full_JSD) <- colnames(Ofus2Cele_orth_tpm_qn)[-c(1)]
rownames(Ofus2Cele_full_JSD) <- colnames(Cele2Ofus_orth_tpm_qn)[-c(1)]
Ofus2Cele_mean_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Cele)-1,rep=TRUE)))
colnames(Ofus2Cele_mean_JSD) <- colnames(Ofus2Cele_orth_tpm_qn)[-c(1)]
rownames(Ofus2Cele_mean_JSD) <- colnames(Cele2Ofus_orth_tpm_qn)[-c(1)]
Ofus2Cele_sd_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Cele)-1,rep=TRUE)))
colnames(Ofus2Cele_sd_JSD) <- colnames(Ofus2Cele_orth_tpm_qn)[-c(1)]
rownames(Ofus2Cele_sd_JSD) <- colnames(Cele2Ofus_orth_tpm_qn)[-c(1)]

for (i in colnames(Ofus2Cele_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Cele2Ofus_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Ofus[[i]],match_Cele[[j]])
    Ofus2Cele_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Ofus2Cele_mean_JSD[j,i] <- mean(all_JSD)
    Ofus2Cele_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Results")
write.table(Ofus2Cele_full_JSD, "Ofus2Cele_full_set_JSD.txt", sep ='\t')
write.table(Ofus2Cele_mean_JSD, "Ofus2Cele_subset_JSD_mean.txt", sep ='\t')
write.table(Ofus2Cele_sd_JSD, "Ofus2Cele_subset_JSD_sd.txt", sep ='\t')


# 8. Owenia fusiformis vs. Capitella teleta
# Number of orthologues = 7,651 orthologues
# 14 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Quantile_transformation")
Ofus2Ctel_orth_tpm_qn <- read.csv("Ofus2Ctel_TPM_mean_quantile_transform.csv", header = T)
Ctel2Ofus_orth_tpm_qn <- read.csv("Ctel2Ofus_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/one2one")
both_ID <- read.table("Ofus2Ctel.txt", header = T)
match_Ofus <- Ofus2Ctel_orth_tpm_qn[match(both_ID$Ofus,Ofus2Ctel_orth_tpm_qn$Gene_ID),]
match_Ctel <- Ctel2Ofus_orth_tpm_qn[match(both_ID$Ctel,Ctel2Ofus_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/order")
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

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Results")
write.table(Ofus2Ctel_full_JSD, "Ofus2Ctel_full_set_JSD.txt", sep ='\t')
write.table(Ofus2Ctel_mean_JSD, "Ofus2Ctel_subset_JSD_mean.txt", sep ='\t')
write.table(Ofus2Ctel_sd_JSD, "Ofus2Ctel_subset_JSD_sd.txt", sep ='\t')


# 9. Owenia fusiformis vs. Amphimedon queenslandica
# Number of orthologues = 3,962 orthologues
# 9 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Quantile_transformation")
Ofus2Aque_orth_tpm_qn <- read.csv("Ofus2Aque_TPM_mean_quantile_transform.csv", header = T)
Aque2Ofus_orth_tpm_qn <- read.csv("Aque2Ofus_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/one2one")
both_ID <- read.table("Ofus2Aque.txt", header = T)
match_Ofus <- Ofus2Aque_orth_tpm_qn[match(both_ID$Ofus,Ofus2Aque_orth_tpm_qn$Gene_ID),]
match_Aque <- Aque2Ofus_orth_tpm_qn[match(both_ID$Aque,Aque2Ofus_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/order")
write.table(match_Ofus,"order_Ofus2Aque_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Aque,"order_Aque2Ofus_tpm_qn.txt", sep="\t", quote=F)
match_Ofus <- read.table("order_Ofus2Aque_tpm_qn.txt", header = T)
match_Aque <- read.table("order_Aque2Ofus_tpm_qn.txt", header = T)

Ofus2Aque_full_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Ofus2Aque_full_JSD) <- colnames(Ofus2Aque_orth_tpm_qn)[-c(1)]
rownames(Ofus2Aque_full_JSD) <- colnames(Aque2Ofus_orth_tpm_qn)[-c(1)]
Ofus2Aque_mean_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Ofus2Aque_mean_JSD) <- colnames(Ofus2Aque_orth_tpm_qn)[-c(1)]
rownames(Ofus2Aque_mean_JSD) <- colnames(Aque2Ofus_orth_tpm_qn)[-c(1)]
Ofus2Aque_sd_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Ofus2Aque_sd_JSD) <- colnames(Ofus2Aque_orth_tpm_qn)[-c(1)]
rownames(Ofus2Aque_sd_JSD) <- colnames(Aque2Ofus_orth_tpm_qn)[-c(1)]

for (i in colnames(Ofus2Aque_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Aque2Ofus_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Ofus[[i]],match_Aque[[j]])
    Ofus2Aque_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Ofus2Aque_mean_JSD[j,i] <- mean(all_JSD)
    Ofus2Aque_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Results")
write.table(Ofus2Aque_full_JSD, "Ofus2Aque_full_set_JSD.txt", sep ='\t')
write.table(Ofus2Aque_mean_JSD, "Ofus2Aque_subset_JSD_mean.txt", sep ='\t')
write.table(Ofus2Aque_sd_JSD, "Ofus2Aque_subset_JSD_sd.txt", sep ='\t')


# 10. Owenia fusiformis vs. Clytica hemisphaerica
# Number of orthologues = 4,691 orthologues
# 9 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Quantile_transformation")
Ofus2Chem_orth_tpm_qn <- read.csv("Ofus2Chem_TPM_mean_quantile_transform.csv", header = T)
Chem2Ofus_orth_tpm_qn <- read.csv("Chem2Ofus_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/one2one")
both_ID <- read.table("Ofus2Chem.txt", header = T)
match_Ofus <- Ofus2Chem_orth_tpm_qn[match(both_ID$Ofus,Ofus2Chem_orth_tpm_qn$Gene_ID),]
match_Chem <- Chem2Ofus_orth_tpm_qn[match(both_ID$Chem,Chem2Ofus_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/order")
write.table(match_Ofus,"order_Ofus2Chem_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Chem,"order_Chem2Ofus_tpm_qn.txt", sep="\t", quote=F)
match_Ofus <- read.table("order_Ofus2Chem_tpm_qn.txt", header = T)
match_Chem <- read.table("order_Chem2Ofus_tpm_qn.txt", header = T)

Ofus2Chem_full_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Ofus2Chem_full_JSD) <- colnames(Ofus2Chem_orth_tpm_qn)[-c(1)]
rownames(Ofus2Chem_full_JSD) <- colnames(Chem2Ofus_orth_tpm_qn)[-c(1)]
Ofus2Chem_mean_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Ofus2Chem_mean_JSD) <- colnames(Ofus2Chem_orth_tpm_qn)[-c(1)]
rownames(Ofus2Chem_mean_JSD) <- colnames(Chem2Ofus_orth_tpm_qn)[-c(1)]
Ofus2Chem_sd_JSD <- data.frame(replicate(ncol(match_Ofus)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Ofus2Chem_sd_JSD) <- colnames(Ofus2Chem_orth_tpm_qn)[-c(1)]
rownames(Ofus2Chem_sd_JSD) <- colnames(Chem2Ofus_orth_tpm_qn)[-c(1)]

for (i in colnames(Ofus2Chem_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Chem2Ofus_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Ofus[[i]],match_Chem[[j]])
    Ofus2Chem_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Ofus2Chem_mean_JSD[j,i] <- mean(all_JSD)
    Ofus2Chem_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Results")
write.table(Ofus2Chem_full_JSD, "Ofus2Chem_full_set_JSD.txt", sep ='\t')
write.table(Ofus2Chem_mean_JSD, "Ofus2Chem_subset_JSD_mean.txt", sep ='\t')
write.table(Ofus2Chem_sd_JSD, "Ofus2Chem_subset_JSD_sd.txt", sep ='\t')


# 11. Re-import data from bootstrapped JSDs and plot
# non-comparable comparisons: heatmaps of RAW data

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Results")

Ofus2Nvec_mean <- read.table("Ofus2Nvec_subset_JSD_mean.txt", header = T)
Ofus2Spur_mean <- read.table("Ofus2Spur_subset_JSD_mean.txt", header = T)
Ofus2Cgig_mean <- read.table("Ofus2Cgig_subset_JSD_mean.txt", header = T)
Ofus2Blan_mean <- read.table("Ofus2Blan_subset_JSD_mean.txt", header = T)
Ofus2Drer_mean <- read.table("Ofus2Drer_subset_JSD_mean.txt", header = T)
Ofus2Dmel_mean <- read.table("Ofus2Dmel_subset_JSD_mean.txt", header = T)
Ofus2Cele_mean <- read.table("Ofus2Cele_subset_JSD_mean.txt", header = T)
Ofus2Ctel_mean <- read.table("Ofus2Ctel_subset_JSD_mean.txt", header = T)
Ofus2Aque_mean <- read.table("Ofus2Aque_subset_JSD_mean.txt", header = T)
Ofus2Chem_mean <- read.table("Ofus2Chem_subset_JSD_mean.txt", header = T)

#allvalues <- rbind(Ofus2Nvec_mean,Ofus2Spur_mean,Ofus2Cgig_mean,Ofus2Blan_mean,Ofus2Drer_mean,Ofus2Dmel_mean,Ofus2Cele_mean,Ofus2Ctel_mean)
#minimum <- round_any(min(allvalues), 10, f = floor)
#maximum <- round_any(max(allvalues), 10, f = ceiling)

heatmap_color <- colorRampPalette(brewer.pal(n = 7, name = "RdGy"))(200)

h1 <- ComplexHeatmap::Heatmap(data.matrix(Ofus2Nvec_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Ofus2Nvec_mean)*unit(5,"mm"),
                              height=nrow(Ofus2Nvec_mean)*unit(5,"mm"),
                              col=heatmap_color,
                              column_labels = c("oocyte","1 cell","2 cells","4 cells","8 cells",
                                                "3 hpf","4 hpf","blastula","gastrula","elongation",
                                                "early larva","mitraria larva","competent larva","juvenile"),
                              row_labels = c("24 hpf","48 hpf","72 hpf","96 hpf","120 hpf","144 hpf","168 hpf","192 hpf"),
                              heatmap_legend_param = list(title="JSD_raw"))

h2 <- ComplexHeatmap::Heatmap(data.matrix(Ofus2Spur_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Ofus2Spur_mean)*unit(5,"mm"),
                              height=nrow(Ofus2Spur_mean)*unit(5,"mm"),
                              col=heatmap_color,
                              column_labels = c("oocyte","1 cell","2 cells","4 cells","8 cells",
                                                "3 hpf","4 hpf","blastula","gastrula","elongation",
                                                "early larva","mitraria larva","competent larva","juvenile"),
                              row_labels = c("hatched blastula","mesenchyme blastula","early gastrula","mid gastrula",
                                             "late gastrula","prism","late prism","pluteus","four-arm stage",
                                             "vestibular stage","pentagonal stage","tube foot stage",
                                             "post-metamorphic","juvenile"),
                              heatmap_legend_param = list(title="JSD_raw"))

h3 <- ComplexHeatmap::Heatmap(data.matrix(Ofus2Cgig_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Ofus2Cgig_mean)*unit(5,"mm"),
                              height=nrow(Ofus2Cgig_mean)*unit(5,"mm"),
                              col=heatmap_color,
                              column_labels = c("oocyte","1 cell","2 cells","4 cells","8 cells",
                                                "3 hpf","4 hpf","blastula","gastrula","elongation",
                                                "early larva","mitraria larva","competent larva","juvenile"),
                              row_labels = c("early morula","morula","blastula","rotary movement","free-swimming",
                                             "early gastrula","gastrula","trochophore 1","trochophore 2","trochophore 3",
                                             "trochophore 4","trochophore 5","early D-stage 1","early D-stage 2",
                                             "D-stage 1","D-stage 2","D-stage 3","D-stage 4","D-stage 5","D-stage 6",
                                             "D-stage 7","early umbo 1","early umbo 2","umbo 1","umbo 2","umbo 3",
                                             "umbo 4","umbo 5","umbo 6","late umbo 1","late umbo 2","pediveliger 1",
                                             "pediveliger 2","spat","juvenile"),
                              heatmap_legend_param = list(title="JSD_raw"))

h4 <- ComplexHeatmap::Heatmap(data.matrix(Ofus2Blan_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Ofus2Blan_mean)*unit(5,"mm"),
                              height=nrow(Ofus2Blan_mean)*unit(5,"mm"),
                              col=heatmap_color,
                              column_labels = c("oocyte","1 cell","2 cells","4 cells","8 cells",
                                                "3 hpf","4 hpf","blastula","gastrula","elongation",
                                                "early larva","mitraria larva","competent larva","juvenile"),
                              row_labels = c("blastula","7 hpf","8 hpf","10 hpf","11 hpf","15 hpf","18 hpf",
                                             "21 hpf","24 hpf","27 hpf","36 hpf","50 hpf","60 hpf","pre-metamorphic"),
                              heatmap_legend_param = list(title="JSD_raw"))

h5 <- ComplexHeatmap::Heatmap(data.matrix(Ofus2Drer_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Ofus2Drer_mean)*unit(5,"mm"),
                              height=nrow(Ofus2Drer_mean)*unit(5,"mm"),
                              col=heatmap_color,
                              column_labels = c("oocyte","1 cell","2 cells","4 cells","8 cells",
                                                "3 hpf","4 hpf","blastula","gastrula","elongation",
                                                "early larva","mitraria larva","competent larva","juvenile"),
                              row_labels = c("2 hpf", "1,000 cells", "Dome", "Shield", "6 hpf", "8 hpf",
                                             "Bud", "12 hpf", "16 hpf", "20 hpf", "26 hpf", "28 hpf", "30 hpf",
                                             "36 hpf", "2 dpf", "3 dpf", "5 dpf", "7 dpf"),
                              heatmap_legend_param = list(title="JSD_raw"))

h6 <- ComplexHeatmap::Heatmap(data.matrix(Ofus2Dmel_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Ofus2Dmel_mean)*unit(5,"mm"),
                              height=nrow(Ofus2Dmel_mean)*unit(5,"mm"),
                              col=heatmap_color,
                              column_labels = c("oocyte","1 cell","2 cells","4 cells","8 cells",
                                                "3 hpf","4 hpf","blastula","gastrula","elongation",
                                                "early larva","mitraria larva","competent larva","juvenile"),
                              row_labels = c("0-2 hpf","2-4 hpf","4-8 hpf","8-12 hpf","12-16 hpf","16-20 hpf",
                                             "20-24 hpf","L1 larva","L2 larva","L3 larva 12 h","L3 larva puffstage 1-9",
                                             "white prepupa","late prepupa/early pupa","pupa 2-3 d","pupa 4 d","adult 1 d"),
                              heatmap_legend_param = list(title="JSD_raw"))

h7 <- ComplexHeatmap::Heatmap(data.matrix(Ofus2Cele_mean[-c(10,13),]),cluster_rows=F,cluster_columns=F,
                              width=ncol(Ofus2Cele_mean[-c(10,13),])*unit(5,"mm"),
                              height=nrow(Ofus2Cele_mean[-c(10,13),])*unit(5,"mm"),
                              col=heatmap_color,
                              column_labels = c("oocyte","1 cell","2 cells","4 cells","8 cells",
                                                "3 hpf","4 hpf","blastula","gastrula","elongation",
                                                "early larva","mitraria larva","competent larva","juvenile"),
                              row_labels = c("0-90 mpf","120-180 mpf", "210-240 mpf", "270-330 mpf",
                                            "360-390 mpf","420-480 mpf","510-540 mpf","570-600 mpf","L1 larva",
                                            "L3 larva","L4 larva"),
                              heatmap_legend_param = list(title="JSD_raw"))

h8 <- ComplexHeatmap::Heatmap(data.matrix(Ofus2Ctel_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Ofus2Ctel_mean)*unit(5,"mm"),
                              height=nrow(Ofus2Ctel_mean)*unit(5,"mm"),
                              col=heatmap_color,
                              column_labels = c("oocyte","1 cell","2 cells","4 cells","8 cells",
                                                "3 hpf","4 hpf","blastula","gastrula","elongation",
                                                "early larva","mitraria larva","competent larva","juvenile"),
                              row_labels = c("oocyte","1 cell","2 cells","4 cells","8 cells","16 cells",
                                             "32 cells","64 cells","gastrula","stage 4tt larva","stage 5 larva",
                                             "stage 7 larva","pre-competent larva","competent larva"),
                              heatmap_legend_param = list(title="JSD_raw"))

h9 <- ComplexHeatmap::Heatmap(data.matrix(Ofus2Aque_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Ofus2Aque_mean)*unit(5,"mm"),
                              height=nrow(Ofus2Aque_mean)*unit(5,"mm"),
                              col=heatmap_color,
                              column_labels = c("oocyte","1 cell","2 cells","4 cells","8 cells",
                                                "3 hpf","4 hpf","blastula","gastrula","elongation",
                                                "early larva","mitraria larva","competent larva","juvenile"),
                              row_labels = c("cleavage","brown 24 hpf", "spot 48 hpf", "ring 72 hpf",
                                             "pre-competent larva", "competent larva", "early post-settlement",
                                             "post-settlement", "juvenile"),
                              heatmap_legend_param = list(title="JSD_raw"))

h10 <- ComplexHeatmap::Heatmap(data.matrix(Ofus2Chem_mean),cluster_rows=F,cluster_columns=F,
                              width=ncol(Ofus2Chem_mean)*unit(5,"mm"),
                              height=nrow(Ofus2Chem_mean)*unit(5,"mm"),
                              col=heatmap_color,
                              column_labels = c("oocyte","1 cell","2 cells","4 cells","8 cells",
                                                "3 hpf","4 hpf","blastula","gastrula","elongation",
                                                "early larva","mitraria larva","competent larva","juvenile"),
                              row_labels = c("early gastrula", "planula 24 hpf", "planula 48 hpf", "planula 72 hpf",
                                             "primary polyp", "gastrozooid", "gonozooid", "stolon", "baby medusa"),
                              heatmap_legend_param = list(title="JSD_raw"))

h1
h2
h3
h4
h5
h6
h7
h8
h9
h10

############################################
## Owenia fusiformis 1-to-all comparisons ##
############################################

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Results")

Ofus2Nvec_mean <- read.table("Ofus2Nvec_subset_JSD_mean.txt", header = T)[,-c(1:7)]
Ofus2Spur_mean <- read.table("Ofus2Spur_subset_JSD_mean.txt", header = T)[,-c(1:7)]
Ofus2Cgig_mean <- read.table("Ofus2Cgig_subset_JSD_mean.txt", header = T)[,-c(1:7)]
Ofus2Blan_mean <- read.table("Ofus2Blan_subset_JSD_mean.txt", header = T)[,-c(1:7)]
Ofus2Drer_mean <- read.table("Ofus2Drer_subset_JSD_mean.txt", header = T)[,-c(1:7)]
Ofus2Dmel_mean <- read.table("Ofus2Dmel_subset_JSD_mean.txt", header = T)[,-c(1:7)]
Ofus2Cele_mean <- read.table("Ofus2Cele_subset_JSD_mean.txt", header = T)[-c(12),-c(1:7)]
Ofus2Ctel_mean <- read.table("Ofus2Ctel_subset_JSD_mean.txt", header = T)[-c(1:7),-c(1:7)]
Ofus2Aque_mean <- read.table("Ofus2Aque_subset_JSD_mean.txt", header = T)[,-c(1:7)]
Ofus2Chem_mean <- read.table("Ofus2Chem_subset_JSD_mean.txt", header = T)[,-c(1:7)]

Ofus2Nvec_sd <- read.table("Ofus2Nvec_subset_JSD_sd.txt", header = T)[,-c(1:7)]
Ofus2Spur_sd <- read.table("Ofus2Spur_subset_JSD_sd.txt", header = T)[,-c(1:7)]
Ofus2Cgig_sd <- read.table("Ofus2Cgig_subset_JSD_sd.txt", header = T)[,-c(1:7)]
Ofus2Blan_sd <- read.table("Ofus2Blan_subset_JSD_sd.txt", header = T)[,-c(1:7)]
Ofus2Drer_sd <- read.table("Ofus2Drer_subset_JSD_sd.txt", header = T)[,-c(1:7)]
Ofus2Dmel_sd <- read.table("Ofus2Dmel_subset_JSD_sd.txt", header = T)[,-c(1:7)]
Ofus2Cele_sd <- read.table("Ofus2Cele_subset_JSD_sd.txt", header = T)[-c(12),-c(1:7)]
Ofus2Ctel_sd <- read.table("Ofus2Ctel_subset_JSD_sd.txt", header = T)[-c(1:7),-c(1:7)]
Ofus2Aque_sd <- read.table("Ofus2Aque_subset_JSD_sd.txt", header = T)[,-c(1:7)]
Ofus2Chem_sd <- read.table("Ofus2Chem_subset_JSD_sd.txt", header = T)[,-c(1:7)]

Ofus2Nvec_transposed <- data.frame(t(Ofus2Nvec_mean))
Ofus2Spur_transposed <- data.frame(t(Ofus2Spur_mean))
Ofus2Cgig_transposed <- data.frame(t(Ofus2Cgig_mean))
Ofus2Blan_transposed <- data.frame(t(Ofus2Blan_mean))
Ofus2Drer_transposed <- data.frame(t(Ofus2Drer_mean))
Ofus2Dmel_transposed <- data.frame(t(Ofus2Dmel_mean))
Ofus2Cele_transposed <- data.frame(t(Ofus2Cele_mean))
Ofus2Ctel_transposed <- data.frame(t(Ofus2Ctel_mean))
Ofus2Aque_transposed <- data.frame(t(Ofus2Aque_mean))
Ofus2Chem_transposed <- data.frame(t(Ofus2Chem_mean))

min_stages <- data.frame(replicate(10,sample(0:1,ncol(Ofus2Drer_mean),rep=TRUE)))
colnames(min_stages) <- c("Nematostella_vectensis","Strongylocentrotus_purpuratus",
                          "Crassostrea_gigas","Branchiostoma_lanceolatum",
                          "Danio_rerio", "Caenorhabditis_elegans",
                          "Drosophila_melanogaster","Capitella_teleta",
                          "Amphimedon_queenslandica", "Clytia_hemisphaerica")
rownames(min_stages) <- c("blastula","gastrula","elongation",
                          "early_larva","mitraria_larva","competent_larva","juvenile")
min_stages$Nematostella_vectensis <- names(Ofus2Nvec_transposed)[apply(Ofus2Nvec_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Strongylocentrotus_purpuratus <- names(Ofus2Spur_transposed)[apply(Ofus2Spur_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Crassostrea_gigas <- names(Ofus2Cgig_transposed)[apply(Ofus2Cgig_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Branchiostoma_lanceolatum <- names(Ofus2Blan_transposed)[apply(Ofus2Blan_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Danio_rerio <- names(Ofus2Drer_transposed)[apply(Ofus2Drer_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Drosophila_melanogaster <- names(Ofus2Dmel_transposed)[apply(Ofus2Dmel_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Caenorhabditis_elegans <- names(Ofus2Cele_transposed)[apply(Ofus2Cele_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Capitella_teleta <- names(Ofus2Ctel_transposed)[apply(Ofus2Ctel_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Amphimedon_queenslandica <- names(Ofus2Aque_transposed)[apply(Ofus2Aque_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Clytia_hemisphaerica <- names(Ofus2Chem_transposed)[apply(Ofus2Chem_transposed, MARGIN = 1, FUN = which.min)]

minimum_distance_mean <- data.frame(replicate(10,sample(0:1,ncol(Ofus2Drer_mean),rep=TRUE)))
minimum_distance_sd <- data.frame(replicate(10,sample(0:1,ncol(Ofus2Drer_mean),rep=TRUE)))
colnames(minimum_distance_mean) <- c("Nematostella_vectensis","Strongylocentrotus_purpuratus",
                                     "Crassostrea_gigas","Branchiostoma_lanceolatum",
                                     "Danio_rerio", "Caenorhabditis_elegans",
                                     "Drosophila_melanogaster","Capitella_teleta",
                                     "Amphimedon_queenslandica", "Clytia_hemisphaerica")
rownames(minimum_distance_mean) <- c("blastula","gastrula","elongation",
                                     "early_larva","mitraria_larva","competent_larva","juvenile")
colnames(minimum_distance_sd) <- c("Nematostella_vectensis","Strongylocentrotus_purpuratus",
                                   "Crassostrea_gigas","Branchiostoma_lanceolatum",
                                   "Danio_rerio", "Caenorhabditis_elegans",
                                   "Drosophila_melanogaster","Capitella_teleta",
                                   "Amphimedon_queenslandica", "Clytia_hemisphaerica")
rownames(minimum_distance_sd) <- c("blastula","gastrula","elongation",
                                   "early_larva","mitraria_larva","competent_larva","juvenile")

for (i in c(1:ncol(Ofus2Drer_mean))){
  minimum_distance_mean[i,"Nematostella_vectensis"] <- Ofus2Nvec_mean[min_stages$Nematostella_vectensis[i],i]
  minimum_distance_mean[i,"Strongylocentrotus_purpuratus"] <- Ofus2Spur_mean[min_stages$Strongylocentrotus_purpuratus[i],i]
  minimum_distance_mean[i,"Crassostrea_gigas"] <- Ofus2Cgig_mean[min_stages$Crassostrea_gigas[i],i]
  minimum_distance_mean[i,"Branchiostoma_lanceolatum"] <- Ofus2Blan_mean[min_stages$Branchiostoma_lanceolatum[i],i]
  minimum_distance_mean[i,"Danio_rerio"] <- Ofus2Drer_mean[min_stages$Danio_rerio[i],i]
  minimum_distance_mean[i,"Drosophila_melanogaster"] <- Ofus2Dmel_mean[min_stages$Drosophila_melanogaster[i],i]
  minimum_distance_mean[i,"Caenorhabditis_elegans"] <- Ofus2Cele_mean[min_stages$Caenorhabditis_elegans[i],i]
  minimum_distance_mean[i,"Capitella_teleta"] <- Ofus2Ctel_mean[min_stages$Capitella_teleta[i],i]
  minimum_distance_mean[i,"Amphimedon_queenslandica"] <- Ofus2Aque_mean[min_stages$Amphimedon_queenslandica[i],i]
  minimum_distance_mean[i,"Clytia_hemisphaerica"] <- Ofus2Chem_mean[min_stages$Clytia_hemisphaerica[i],i]
  minimum_distance_sd[i,"Nematostella_vectensis"] <- Ofus2Nvec_sd[min_stages$Nematostella_vectensis[i],i]
  minimum_distance_sd[i,"Strongylocentrotus_purpuratus"] <- Ofus2Spur_sd[min_stages$Strongylocentrotus_purpuratus[i],i]
  minimum_distance_sd[i,"Crassostrea_gigas"] <- Ofus2Cgig_sd[min_stages$Crassostrea_gigas[i],i]
  minimum_distance_sd[i,"Branchiostoma_lanceolatum"] <- Ofus2Blan_sd[min_stages$Branchiostoma_lanceolatum[i],i]
  minimum_distance_sd[i,"Danio_rerio"] <- Ofus2Drer_sd[min_stages$Danio_rerio[i],i]
  minimum_distance_sd[i,"Drosophila_melanogaster"] <- Ofus2Dmel_sd[min_stages$Drosophila_melanogaster[i],i]
  minimum_distance_sd[i,"Caenorhabditis_elegans"] <- Ofus2Cele_sd[min_stages$Caenorhabditis_elegans[i],i]
  minimum_distance_sd[i,"Capitella_teleta"] <- Ofus2Ctel_sd[min_stages$Capitella_teleta[i],i]
  minimum_distance_sd[i,"Amphimedon_queenslandica"] <- Ofus2Aque_sd[min_stages$Amphimedon_queenslandica[i],i]
  minimum_distance_sd[i,"Clytia_hemisphaerica"] <- Ofus2Chem_sd[min_stages$Clytia_hemisphaerica[i],i]
}

minimum_distance_mean$stage <- c(1:ncol(Ofus2Drer_mean))
minimum_distance_sd$stage <- c(1:ncol(Ofus2Drer_mean))
minimum_distance_mean_tidy <- gather(minimum_distance_mean, "species", "JSD", -stage)
minimum_distance_sd_tidy <- gather(minimum_distance_sd, "species", "JSD", -stage)
minimum_distance_mean_tidy$upper <- minimum_distance_mean_tidy["JSD"] + minimum_distance_sd_tidy["JSD"]
minimum_distance_mean_tidy$lower <- minimum_distance_mean_tidy["JSD"] - minimum_distance_sd_tidy["JSD"]
colnames(minimum_distance_mean_tidy)[4:5] <- c("upper","lower")
minimum_distance_mean_tidy$upper <- unlist(minimum_distance_mean_tidy$upper)
minimum_distance_mean_tidy$lower <- unlist(minimum_distance_mean_tidy$lower)

minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Nematostella_vectensis'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Nematostella_vectensis'),3:5]-min(minimum_distance_mean$Nematostella_vectensis))/(max(minimum_distance_mean$Nematostella_vectensis)-min(minimum_distance_mean$Nematostella_vectensis))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Strongylocentrotus_purpuratus'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Strongylocentrotus_purpuratus'),3:5]-min(minimum_distance_mean$Strongylocentrotus_purpuratus))/(max(minimum_distance_mean$Strongylocentrotus_purpuratus)-min(minimum_distance_mean$Strongylocentrotus_purpuratus))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Crassostrea_gigas'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Crassostrea_gigas'),3:5]-min(minimum_distance_mean$Crassostrea_gigas))/(max(minimum_distance_mean$Crassostrea_gigas)-min(minimum_distance_mean$Crassostrea_gigas))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Branchiostoma_lanceolatum'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Branchiostoma_lanceolatum'),3:5]-min(minimum_distance_mean$Branchiostoma_lanceolatum))/(max(minimum_distance_mean$Branchiostoma_lanceolatum)-min(minimum_distance_mean$Branchiostoma_lanceolatum))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Danio_rerio'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Danio_rerio'),3:5]-min(minimum_distance_mean$Danio_rerio))/(max(minimum_distance_mean$Danio_rerio)-min(minimum_distance_mean$Danio_rerio))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Caenorhabditis_elegans'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Caenorhabditis_elegans'),3:5]-min(minimum_distance_mean$Caenorhabditis_elegans))/(max(minimum_distance_mean$Caenorhabditis_elegans)-min(minimum_distance_mean$Caenorhabditis_elegans))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Drosophila_melanogaster'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Drosophila_melanogaster'),3:5]-min(minimum_distance_mean$Drosophila_melanogaster))/(max(minimum_distance_mean$Drosophila_melanogaster)-min(minimum_distance_mean$Drosophila_melanogaster))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Capitella_teleta'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Capitella_teleta'),3:5]-min(minimum_distance_mean$Capitella_teleta))/(max(minimum_distance_mean$Capitella_teleta)-min(minimum_distance_mean$Capitella_teleta))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Amphimedon_queenslandica'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Amphimedon_queenslandica'),3:5]-min(minimum_distance_mean$Amphimedon_queenslandica))/(max(minimum_distance_mean$Amphimedon_queenslandica)-min(minimum_distance_mean$Amphimedon_queenslandica))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Clytia_hemisphaerica'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Clytia_hemisphaerica'),3:5]-min(minimum_distance_mean$Clytia_hemisphaerica))/(max(minimum_distance_mean$Clytia_hemisphaerica)-min(minimum_distance_mean$Clytia_hemisphaerica))

final_dataset <- minimum_distance_mean_tidy
final_dataset$species <- factor(final_dataset$species,
                                levels = c("Capitella_teleta","Crassostrea_gigas",
                                           "Caenorhabditis_elegans","Drosophila_melanogaster",
                                           "Danio_rerio","Branchiostoma_lanceolatum",
                                           "Strongylocentrotus_purpuratus","Nematostella_vectensis",
                                           "Clytia_hemisphaerica","Amphimedon_queenslandica"))

ggplot(final_dataset, aes(x=stage, y=JSD, colour=species)) + 
  geom_line(show.legend = FALSE) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=factor(species)), 
              colour = NA, show.legend = FALSE,
              alpha = 0.5) +
  facet_wrap(~species, nrow=10, ncol=1) +
  scale_x_continuous(labels = c(1:ncol(Ofus2Drer_mean)), breaks = c(1:ncol(Ofus2Drer_mean))) +
  theme_classic() +
  labs(x = "O. fusiformis stage", y = "Normalised gene expression divergence (JSD)")

write.table(min_stages, "minimum_stages.txt", sep='\t', quote = FALSE)
write.table(minimum_distance_mean, "minimum_stages_distance_mean.txt", sep='\t', quote = FALSE)
write.table(minimum_distance_sd, "minimum_stages_distance_sd.txt", sep='\t', quote = FALSE)
write.table(minimum_distance_mean_tidy, "minimum_stages_to_plot.txt", sep='\t', quote = FALSE)
