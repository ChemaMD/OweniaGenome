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
## Capitella teleta 1-to-1 comparisons ##
##########################################

# 1. Capitella teleta vs. Nematostella vectensis
# Number of orthologues = 5,254 orthologues (-4)
# 8 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/Quantile_transformation")
Ctel2Nvec_orth_tpm_qn <- read.csv("Ctel2Nvec_TPM_mean_quantile_transform.csv", header = T)
Nvec2Ctel_orth_tpm_qn <- read.csv("Nvec2Ctel_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/one2one")
both_ID <- read.table("Ctel2Nvec.txt", header = T)
match_Ctel <- Ctel2Nvec_orth_tpm_qn[match(both_ID$Ctel,Ctel2Nvec_orth_tpm_qn$Gene_ID),]
match_Nvec <- Nvec2Ctel_orth_tpm_qn[match(both_ID$Nvec,Nvec2Ctel_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/order")
write.table(match_Ctel,"order_Ctel2Nvec_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Nvec,"order_Nvec2Ctel_tpm_qn.txt", sep="\t", quote=F)
match_Ctel <- read.table("order_Ctel2Nvec_tpm_qn.txt", header = T)
match_Nvec <- read.table("order_Nvec2Ctel_tpm_qn.txt", header = T)

Ctel2Nvec_full_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Ctel2Nvec_full_JSD) <- colnames(Ctel2Nvec_orth_tpm_qn)[-c(1)]
rownames(Ctel2Nvec_full_JSD) <- colnames(Nvec2Ctel_orth_tpm_qn)[-c(1)]
Ctel2Nvec_mean_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Ctel2Nvec_mean_JSD) <- colnames(Ctel2Nvec_orth_tpm_qn)[-c(1)]
rownames(Ctel2Nvec_mean_JSD) <- colnames(Nvec2Ctel_orth_tpm_qn)[-c(1)]
Ctel2Nvec_sd_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Ctel2Nvec_sd_JSD) <- colnames(Ctel2Nvec_orth_tpm_qn)[-c(1)]
rownames(Ctel2Nvec_sd_JSD) <- colnames(Nvec2Ctel_orth_tpm_qn)[-c(1)]

for (i in colnames(Ctel2Nvec_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Nvec2Ctel_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Ctel[[i]],match_Nvec[[j]])
    Ctel2Nvec_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Ctel2Nvec_mean_JSD[j,i] <- mean(all_JSD)
    Ctel2Nvec_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/Results")
write.table(Ctel2Nvec_full_JSD, "Ctel2Nvec_full_set_JSD.txt", sep ='\t')
write.table(Ctel2Nvec_mean_JSD, "Ctel2Nvec_subset_JSD_mean.txt", sep ='\t')
write.table(Ctel2Nvec_sd_JSD, "Ctel2Nvec_subset_JSD_sd.txt", sep ='\t')


# 2. Capitella teleta vs. Strongylocentrotus purpuratus
# Number of orthologues = 5,015 orthologues (+0)
# 14 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/Quantile_transformation")
Ctel2Spur_orth_tpm_qn <- read.csv("Ctel2Spur_TPM_mean_quantile_transform.csv", header = T)
Spur2Ctel_orth_tpm_qn <- read.csv("Spur2Ctel_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/one2one")
both_ID <- read.table("Ctel2Spur.txt", header = T)
match_Ctel <- Ctel2Spur_orth_tpm_qn[match(both_ID$Ctel,Ctel2Spur_orth_tpm_qn$Gene_ID),]
match_Spur <- Spur2Ctel_orth_tpm_qn[match(both_ID$Spur,Spur2Ctel_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/order")
write.table(match_Ctel,"order_Ctel2Spur_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Spur,"order_Spur2Ctel_tpm_qn.txt", sep="\t", quote=F)
match_Ctel <- read.table("order_Ctel2Spur_tpm_qn.txt", header = T)
match_Spur <- read.table("order_Spur2Ctel_tpm_qn.txt", header = T)

Ctel2Spur_full_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Spur)-1,rep=TRUE)))
colnames(Ctel2Spur_full_JSD) <- colnames(Ctel2Spur_orth_tpm_qn)[-c(1)]
rownames(Ctel2Spur_full_JSD) <- colnames(Spur2Ctel_orth_tpm_qn)[-c(1)]
Ctel2Spur_mean_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Spur)-1,rep=TRUE)))
colnames(Ctel2Spur_mean_JSD) <- colnames(Ctel2Spur_orth_tpm_qn)[-c(1)]
rownames(Ctel2Spur_mean_JSD) <- colnames(Spur2Ctel_orth_tpm_qn)[-c(1)]
Ctel2Spur_sd_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Spur)-1,rep=TRUE)))
colnames(Ctel2Spur_sd_JSD) <- colnames(Ctel2Spur_orth_tpm_qn)[-c(1)]
rownames(Ctel2Spur_sd_JSD) <- colnames(Spur2Ctel_orth_tpm_qn)[-c(1)]

for (i in colnames(Ctel2Spur_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Spur2Ctel_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Ctel[[i]],match_Spur[[j]])
    Ctel2Spur_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Ctel2Spur_mean_JSD[j,i] <- mean(all_JSD)
    Ctel2Spur_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/Results")
write.table(Ctel2Spur_full_JSD, "Ctel2Spur_full_set_JSD.txt", sep ='\t')
write.table(Ctel2Spur_mean_JSD, "Ctel2Spur_subset_JSD_mean.txt", sep ='\t')
write.table(Ctel2Spur_sd_JSD, "Ctel2Spur_subset_JSD_sd.txt", sep ='\t')


# 3. Capitella teleta vs. Crassostrea gigas
# Number of orthologues = 6,737 orthologues (+62)
# 35 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/Quantile_transformation")
Ctel2Cgig_orth_tpm_qn <- read.csv("Ctel2Cgig_TPM_mean_quantile_transform.csv", header = T)
Cgig2Ctel_orth_tpm_qn <- read.csv("Cgig2Ctel_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/one2one")
both_ID <- read.table("Ctel2Cgig.txt", header = T)
match_Ctel <- Ctel2Cgig_orth_tpm_qn[match(both_ID$Ctel,Ctel2Cgig_orth_tpm_qn$Gene_ID),]
match_Cgig <- Cgig2Ctel_orth_tpm_qn[match(both_ID$Cgig,Cgig2Ctel_orth_tpm_qn$Gene_ID),]
match_Cgig <- na.omit(match_Cgig)
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/order")
write.table(match_Ctel,"order_Ctel2Cgig_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Cgig,"order_Cgig2Ctel_tpm_qn.txt", sep="\t", quote=F)
match_Ctel <- read.table("order_Ctel2Cgig_tpm_qn.txt", header = T)
match_Cgig <- read.table("order_Cgig2Ctel_tpm_qn.txt", header = T)

Ctel2Cgig_full_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Cgig)-1,rep=TRUE)))
colnames(Ctel2Cgig_full_JSD) <- colnames(Ctel2Cgig_orth_tpm_qn)[-c(1)]
rownames(Ctel2Cgig_full_JSD) <- colnames(Cgig2Ctel_orth_tpm_qn)[-c(1)]
Ctel2Cgig_mean_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Cgig)-1,rep=TRUE)))
colnames(Ctel2Cgig_mean_JSD) <- colnames(Ctel2Cgig_orth_tpm_qn)[-c(1)]
rownames(Ctel2Cgig_mean_JSD) <- colnames(Cgig2Ctel_orth_tpm_qn)[-c(1)]
Ctel2Cgig_sd_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Cgig)-1,rep=TRUE)))
colnames(Ctel2Cgig_sd_JSD) <- colnames(Ctel2Cgig_orth_tpm_qn)[-c(1)]
rownames(Ctel2Cgig_sd_JSD) <- colnames(Cgig2Ctel_orth_tpm_qn)[-c(1)]

for (i in colnames(Ctel2Cgig_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Cgig2Ctel_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Ctel[[i]],match_Cgig[[j]])
    Ctel2Cgig_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Ctel2Cgig_mean_JSD[j,i] <- mean(all_JSD)
    Ctel2Cgig_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/Results")
write.table(Ctel2Cgig_full_JSD, "Ctel2Cgig_full_set_JSD.txt", sep ='\t')
write.table(Ctel2Cgig_mean_JSD, "Ctel2Cgig_subset_JSD_mean.txt", sep ='\t')
write.table(Ctel2Cgig_sd_JSD, "Ctel2Cgig_subset_JSD_sd.txt", sep ='\t')


# 4. Capitella teleta vs. Branchiostoma lanceolatum
# Number of orthologues = 6,673 orthologues (+39)
# 14 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/Quantile_transformation")
Ctel2Blan_orth_tpm_qn <- read.csv("Ctel2Blan_TPM_mean_quantile_transform.csv", header = T)
Blan2Ctel_orth_tpm_qn <- read.csv("Blan2Ctel_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/one2one")
both_ID <- read.table("Ctel2Blan.txt", header = T)
match_Ctel <- Ctel2Blan_orth_tpm_qn[match(both_ID$Ctel,Ctel2Blan_orth_tpm_qn$Gene_ID),]
match_Blan <- Blan2Ctel_orth_tpm_qn[match(both_ID$Blan,Blan2Ctel_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/order")
write.table(match_Ctel,"order_Ctel2Blan_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Blan,"order_Blan2Ctel_tpm_qn.txt", sep="\t", quote=F)
match_Ctel <- read.table("order_Ctel2Blan_tpm_qn.txt", header = T)
match_Blan <- read.table("order_Blan2Ctel_tpm_qn.txt", header = T)

Ctel2Blan_full_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Blan)-1,rep=TRUE)))
colnames(Ctel2Blan_full_JSD) <- colnames(Ctel2Blan_orth_tpm_qn)[-c(1)]
rownames(Ctel2Blan_full_JSD) <- colnames(Blan2Ctel_orth_tpm_qn)[-c(1)]
Ctel2Blan_mean_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Blan)-1,rep=TRUE)))
colnames(Ctel2Blan_mean_JSD) <- colnames(Ctel2Blan_orth_tpm_qn)[-c(1)]
rownames(Ctel2Blan_mean_JSD) <- colnames(Blan2Ctel_orth_tpm_qn)[-c(1)]
Ctel2Blan_sd_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Blan)-1,rep=TRUE)))
colnames(Ctel2Blan_sd_JSD) <- colnames(Ctel2Blan_orth_tpm_qn)[-c(1)]
rownames(Ctel2Blan_sd_JSD) <- colnames(Blan2Ctel_orth_tpm_qn)[-c(1)]

for (i in colnames(Ctel2Blan_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Blan2Ctel_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Ctel[[i]],match_Blan[[j]])
    Ctel2Blan_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Ctel2Blan_mean_JSD[j,i] <- mean(all_JSD)
    Ctel2Blan_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/Results")
write.table(Ctel2Blan_full_JSD, "Ctel2Blan_full_set_JSD.txt", sep ='\t')
write.table(Ctel2Blan_mean_JSD, "Ctel2Blan_subset_JSD_mean.txt", sep ='\t')
write.table(Ctel2Blan_sd_JSD, "Ctel2Blan_subset_JSD_sd.txt", sep ='\t')


# 5. Capitella teleta vs. Danio rerio
# Number of orthologues = 4,316 orthologues (+21)
# 7 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/Quantile_transformation")
Ctel2Drer_orth_tpm_qn <- read.csv("Ctel2Drer_TPM_mean_quantile_transform.csv", header = T)
Drer2Ctel_orth_tpm_qn <- read.csv("Drer2Ctel_TPM_mean_quantile_transform.csv", header = T)
#Drer2Ctel_orth_tpm_qn <- cbind(Drer2Ctel_orth_tpm_qn_preliminary[1],Drer2Ctel_orth_tpm_qn_preliminary[12],Drer2Ctel_orth_tpm_qn_preliminary[c(3:11)],Drer2Ctel_orth_tpm_qn_preliminary[2],Drer2Ctel_orth_tpm_qn_preliminary[c(13:19)])

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/one2one")
both_ID <- read.table("Ctel2Drer.txt", header = T)
match_Ctel <- Ctel2Drer_orth_tpm_qn[match(both_ID$Ctel,Ctel2Drer_orth_tpm_qn$Gene_ID),]
match_Drer <- Drer2Ctel_orth_tpm_qn[match(both_ID$Drer,Drer2Ctel_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/order")
write.table(match_Ctel,"order_Ctel2Drer_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Drer,"order_Drer2Ctel_tpm_qn.txt", sep="\t", quote=F)
match_Ctel <- read.table("order_Ctel2Drer_tpm_qn.txt", header = T)
match_Drer <- read.table("order_Drer2Ctel_tpm_qn.txt", header = T)

Ctel2Drer_full_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Drer)-1,rep=TRUE)))
colnames(Ctel2Drer_full_JSD) <- colnames(Ctel2Drer_orth_tpm_qn)[-c(1)]
rownames(Ctel2Drer_full_JSD) <- colnames(Drer2Ctel_orth_tpm_qn)[-c(1)]
Ctel2Drer_mean_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Drer)-1,rep=TRUE)))
colnames(Ctel2Drer_mean_JSD) <- colnames(Ctel2Drer_orth_tpm_qn)[-c(1)]
rownames(Ctel2Drer_mean_JSD) <- colnames(Drer2Ctel_orth_tpm_qn)[-c(1)]
Ctel2Drer_sd_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Drer)-1,rep=TRUE)))
colnames(Ctel2Drer_sd_JSD) <- colnames(Ctel2Drer_orth_tpm_qn)[-c(1)]
rownames(Ctel2Drer_sd_JSD) <- colnames(Drer2Ctel_orth_tpm_qn)[-c(1)]

for (i in colnames(Ctel2Drer_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Drer2Ctel_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Ctel[[i]],match_Drer[[j]])
    Ctel2Drer_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Ctel2Drer_mean_JSD[j,i] <- mean(all_JSD)
    Ctel2Drer_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/Results")
write.table(Ctel2Drer_full_JSD, "Ctel2Drer_full_set_JSD.txt", sep ='\t')
write.table(Ctel2Drer_mean_JSD, "Ctel2Drer_subset_JSD_mean.txt", sep ='\t')
write.table(Ctel2Drer_sd_JSD, "Ctel2Drer_subset_JSD_sd.txt", sep ='\t')


# 6. Capitella teleta vs. Drosophila melanogaster
# Number of orthologues = 4,635 orthologues (+45)
# 18 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/Quantile_transformation")
Ctel2Dmel_orth_tpm_qn <- read.csv("Ctel2Dmel_TPM_mean_quantile_transform.csv", header = T)
Dmel2Ctel_orth_tpm_qn <- read.csv("Dmel2Ctel_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/one2one")
both_ID <- read.table("Ctel2Dmel.txt", header = T)
match_Ctel <- Ctel2Dmel_orth_tpm_qn[match(both_ID$Ctel,Ctel2Dmel_orth_tpm_qn$Gene_ID),]
match_Dmel <- Dmel2Ctel_orth_tpm_qn[match(both_ID$Dmel,Dmel2Ctel_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/order")
write.table(match_Ctel,"order_Ctel2Dmel_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Dmel,"order_Dmel2Ctel_tpm_qn.txt", sep="\t", quote=F)
match_Ctel <- read.table("order_Ctel2Dmel_tpm_qn.txt", header = T)
match_Dmel <- read.table("order_Dmel2Ctel_tpm_qn.txt", header = T)

Ctel2Dmel_full_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Dmel)-1,rep=TRUE)))
colnames(Ctel2Dmel_full_JSD) <- colnames(Ctel2Dmel_orth_tpm_qn)[-c(1)]
rownames(Ctel2Dmel_full_JSD) <- colnames(Dmel2Ctel_orth_tpm_qn)[-c(1)]
Ctel2Dmel_mean_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Dmel)-1,rep=TRUE)))
colnames(Ctel2Dmel_mean_JSD) <- colnames(Ctel2Dmel_orth_tpm_qn)[-c(1)]
rownames(Ctel2Dmel_mean_JSD) <- colnames(Dmel2Ctel_orth_tpm_qn)[-c(1)]
Ctel2Dmel_sd_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Dmel)-1,rep=TRUE)))
colnames(Ctel2Dmel_sd_JSD) <- colnames(Ctel2Dmel_orth_tpm_qn)[-c(1)]
rownames(Ctel2Dmel_sd_JSD) <- colnames(Dmel2Ctel_orth_tpm_qn)[-c(1)]

for (i in colnames(Ctel2Dmel_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Dmel2Ctel_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Ctel[[i]],match_Dmel[[j]])
    Ctel2Dmel_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Ctel2Dmel_mean_JSD[j,i] <- mean(all_JSD)
    Ctel2Dmel_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/Results")
write.table(Ctel2Dmel_full_JSD, "Ctel2Dmel_full_set_JSD.txt", sep ='\t')
write.table(Ctel2Dmel_mean_JSD, "Ctel2Dmel_subset_JSD_mean.txt", sep ='\t')
write.table(Ctel2Dmel_sd_JSD, "Ctel2Dmel_subset_JSD_sd.txt", sep ='\t')


# 7. Capitella teleta vs. Caenorhabidits elegans
# Number of orthologues = 3,767 orthologues (-43)
# 13 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/Quantile_transformation")
Ctel2Cele_orth_tpm_qn <- read.csv("Ctel2Cele_TPM_mean_quantile_transform.csv", header = T)
Cele2Ctel_orth_tpm_qn <- read.csv("Cele2Ctel_TPM_mean_quantile_transform.csv", header = T)
Cele2Ctel_orth_tpm_qn$Gene_ID <- paste0("CELE_", Cele2Ctel_orth_tpm_qn$Gene_ID)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/one2one")
both_ID <- read.table("Ctel2Cele.txt", header = T)
match_Ctel <- Ctel2Cele_orth_tpm_qn[match(both_ID$Ctel,Ctel2Cele_orth_tpm_qn$Gene_ID),]
match_Cele <- Cele2Ctel_orth_tpm_qn[match(both_ID$Cele,Cele2Ctel_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/order")
write.table(match_Ctel,"order_Ctel2Cele_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Cele,"order_Cele2Ctel_tpm_qn.txt", sep="\t", quote=F)
match_Ctel <- read.table("order_Ctel2Cele_tpm_qn.txt", header = T)
match_Cele <- read.table("order_Cele2Ctel_tpm_qn.txt", header = T)

Ctel2Cele_full_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Cele)-1,rep=TRUE)))
colnames(Ctel2Cele_full_JSD) <- colnames(Ctel2Cele_orth_tpm_qn)[-c(1)]
rownames(Ctel2Cele_full_JSD) <- colnames(Cele2Ctel_orth_tpm_qn)[-c(1)]
Ctel2Cele_mean_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Cele)-1,rep=TRUE)))
colnames(Ctel2Cele_mean_JSD) <- colnames(Ctel2Cele_orth_tpm_qn)[-c(1)]
rownames(Ctel2Cele_mean_JSD) <- colnames(Cele2Ctel_orth_tpm_qn)[-c(1)]
Ctel2Cele_sd_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Cele)-1,rep=TRUE)))
colnames(Ctel2Cele_sd_JSD) <- colnames(Ctel2Cele_orth_tpm_qn)[-c(1)]
rownames(Ctel2Cele_sd_JSD) <- colnames(Cele2Ctel_orth_tpm_qn)[-c(1)]

for (i in colnames(Ctel2Cele_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Cele2Ctel_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Ctel[[i]],match_Cele[[j]])
    Ctel2Cele_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Ctel2Cele_mean_JSD[j,i] <- mean(all_JSD)
    Ctel2Cele_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/Results")
write.table(Ctel2Cele_full_JSD, "Ctel2Cele_full_set_JSD.txt", sep ='\t')
write.table(Ctel2Cele_mean_JSD, "Ctel2Cele_subset_JSD_mean.txt", sep ='\t')
write.table(Ctel2Cele_sd_JSD, "Ctel2Cele_subset_JSD_sd.txt", sep ='\t')


# 8. Capitella teleta vs. Amphimedon queenslandica
# Number of orthologues = 3,962 orthologues
# 9 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/Quantile_transformation")
Ctel2Aque_orth_tpm_qn <- read.csv("Ctel2Aque_TPM_mean_quantile_transform.csv", header = T)
Aque2Ctel_orth_tpm_qn <- read.csv("Aque2Ctel_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/one2one")
both_ID <- read.table("Ctel2Aque.txt", header = T)
match_Ctel <- Ctel2Aque_orth_tpm_qn[match(both_ID$Ctel,Ctel2Aque_orth_tpm_qn$Gene_ID),]
match_Aque <- Aque2Ctel_orth_tpm_qn[match(both_ID$Aque,Aque2Ctel_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/order")
write.table(match_Ctel,"order_Ctel2Aque_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Aque,"order_Aque2Ctel_tpm_qn.txt", sep="\t", quote=F)
match_Ctel <- read.table("order_Ctel2Aque_tpm_qn.txt", header = T)
match_Aque <- read.table("order_Aque2Ctel_tpm_qn.txt", header = T)

Ctel2Aque_full_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Ctel2Aque_full_JSD) <- colnames(Ctel2Aque_orth_tpm_qn)[-c(1)]
rownames(Ctel2Aque_full_JSD) <- colnames(Aque2Ctel_orth_tpm_qn)[-c(1)]
Ctel2Aque_mean_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Ctel2Aque_mean_JSD) <- colnames(Ctel2Aque_orth_tpm_qn)[-c(1)]
rownames(Ctel2Aque_mean_JSD) <- colnames(Aque2Ctel_orth_tpm_qn)[-c(1)]
Ctel2Aque_sd_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Ctel2Aque_sd_JSD) <- colnames(Ctel2Aque_orth_tpm_qn)[-c(1)]
rownames(Ctel2Aque_sd_JSD) <- colnames(Aque2Ctel_orth_tpm_qn)[-c(1)]

for (i in colnames(Ctel2Aque_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Aque2Ctel_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Ctel[[i]],match_Aque[[j]])
    Ctel2Aque_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Ctel2Aque_mean_JSD[j,i] <- mean(all_JSD)
    Ctel2Aque_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/Results")
write.table(Ctel2Aque_full_JSD, "Ctel2Aque_full_set_JSD.txt", sep ='\t')
write.table(Ctel2Aque_mean_JSD, "Ctel2Aque_subset_JSD_mean.txt", sep ='\t')
write.table(Ctel2Aque_sd_JSD, "Ctel2Aque_subset_JSD_sd.txt", sep ='\t')


# 9. Capitella teleta vs. Clytica hemisphaerica
# Number of orthologues = 4,691 orthologues
# 9 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/Quantile_transformation")
Ctel2Chem_orth_tpm_qn <- read.csv("Ctel2Chem_TPM_mean_quantile_transform.csv", header = T)
Chem2Ctel_orth_tpm_qn <- read.csv("Chem2Ctel_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/one2one")
both_ID <- read.table("Ctel2Chem.txt", header = T)
match_Ctel <- Ctel2Chem_orth_tpm_qn[match(both_ID$Ctel,Ctel2Chem_orth_tpm_qn$Gene_ID),]
match_Chem <- Chem2Ctel_orth_tpm_qn[match(both_ID$Chem,Chem2Ctel_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/order")
write.table(match_Ctel,"order_Ctel2Chem_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Chem,"order_Chem2Ctel_tpm_qn.txt", sep="\t", quote=F)
match_Ctel <- read.table("order_Ctel2Chem_tpm_qn.txt", header = T)
match_Chem <- read.table("order_Chem2Ctel_tpm_qn.txt", header = T)

Ctel2Chem_full_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Ctel2Chem_full_JSD) <- colnames(Ctel2Chem_orth_tpm_qn)[-c(1)]
rownames(Ctel2Chem_full_JSD) <- colnames(Chem2Ctel_orth_tpm_qn)[-c(1)]
Ctel2Chem_mean_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Ctel2Chem_mean_JSD) <- colnames(Ctel2Chem_orth_tpm_qn)[-c(1)]
rownames(Ctel2Chem_mean_JSD) <- colnames(Chem2Ctel_orth_tpm_qn)[-c(1)]
Ctel2Chem_sd_JSD <- data.frame(replicate(ncol(match_Ctel)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Ctel2Chem_sd_JSD) <- colnames(Ctel2Chem_orth_tpm_qn)[-c(1)]
rownames(Ctel2Chem_sd_JSD) <- colnames(Chem2Ctel_orth_tpm_qn)[-c(1)]

for (i in colnames(Ctel2Chem_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Chem2Ctel_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Ctel[[i]],match_Chem[[j]])
    Ctel2Chem_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Ctel2Chem_mean_JSD[j,i] <- mean(all_JSD)
    Ctel2Chem_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/Results")
write.table(Ctel2Chem_full_JSD, "Ctel2Chem_full_set_JSD.txt", sep ='\t')
write.table(Ctel2Chem_mean_JSD, "Ctel2Chem_subset_JSD_mean.txt", sep ='\t')
write.table(Ctel2Chem_sd_JSD, "Ctel2Chem_subset_JSD_sd.txt", sep ='\t')




############################################
## Capitella teleta 1-to-all comparisons ##
############################################

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/Results")

Ctel2Nvec_mean <- read.table("Ctel2Nvec_subset_JSD_mean.txt", header = T)[,-c(1:7)]
Ctel2Spur_mean <- read.table("Ctel2Spur_subset_JSD_mean.txt", header = T)[,-c(1:7)]
Ctel2Cgig_mean <- read.table("Ctel2Cgig_subset_JSD_mean.txt", header = T)[,-c(1:7)]
Ctel2Blan_mean <- read.table("Ctel2Blan_subset_JSD_mean.txt", header = T)[,-c(1:7)]
Ctel2Drer_mean <- read.table("Ctel2Drer_subset_JSD_mean.txt", header = T)[,-c(1:7)]
Ctel2Dmel_mean <- read.table("Ctel2Dmel_subset_JSD_mean.txt", header = T)[-c(17:18),-c(1:7)]
Ctel2Cele_mean <- read.table("Ctel2Cele_subset_JSD_mean.txt", header = T)[-c(10,13),-c(1:7)]
Ctel2Aque_mean <- read.table("Ctel2Aque_subset_JSD_mean.txt", header = T)[,-c(1:7)]
Ctel2Chem_mean <- read.table("Ctel2Chem_subset_JSD_mean.txt", header = T)[,-c(1:7)]

Ctel2Nvec_sd <- read.table("Ctel2Nvec_subset_JSD_sd.txt", header = T)[,-c(1:7)]
Ctel2Spur_sd <- read.table("Ctel2Spur_subset_JSD_sd.txt", header = T)[,-c(1:7)]
Ctel2Cgig_sd <- read.table("Ctel2Cgig_subset_JSD_sd.txt", header = T)[,-c(1:7)]
Ctel2Blan_sd <- read.table("Ctel2Blan_subset_JSD_sd.txt", header = T)[,-c(1:7)]
Ctel2Drer_sd <- read.table("Ctel2Drer_subset_JSD_sd.txt", header = T)[,-c(1:7)]
Ctel2Dmel_sd <- read.table("Ctel2Dmel_subset_JSD_sd.txt", header = T)[-c(17:18),-c(1:7)]
Ctel2Cele_sd <- read.table("Ctel2Cele_subset_JSD_sd.txt", header = T)[-c(10,13),-c(1:7)]
Ctel2Aque_sd <- read.table("Ctel2Aque_subset_JSD_sd.txt", header = T)[,-c(1:7)]
Ctel2Chem_sd <- read.table("Ctel2Chem_subset_JSD_sd.txt", header = T)[,-c(1:7)]

Ctel2Nvec_transposed <- data.frame(t(Ctel2Nvec_mean))
Ctel2Spur_transposed <- data.frame(t(Ctel2Spur_mean))
Ctel2Cgig_transposed <- data.frame(t(Ctel2Cgig_mean))
Ctel2Blan_transposed <- data.frame(t(Ctel2Blan_mean))
Ctel2Drer_transposed <- data.frame(t(Ctel2Drer_mean))
Ctel2Dmel_transposed <- data.frame(t(Ctel2Dmel_mean))
Ctel2Cele_transposed <- data.frame(t(Ctel2Cele_mean))
Ctel2Aque_transposed <- data.frame(t(Ctel2Aque_mean))
Ctel2Chem_transposed <- data.frame(t(Ctel2Chem_mean))

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Results")

Ctel2Ofus_transposed <- read.table("Ofus2Ctel_subset_JSD_mean.txt", header = T)[-c(1:7),]
Ctel2Ofus_mean <- data.frame(t(Ctel2Ofus_transposed))

Ofus2Ctel_sd <- read.table("Ofus2Ctel_subset_JSD_sd.txt", header = T)[-c(1:7),]
Ctel2Ofus_sd <- data.frame(t(Ofus2Ctel_sd))

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/Results")


min_stages <- data.frame(replicate(10,sample(0:1,ncol(Ctel2Drer_mean),rep=TRUE)))
colnames(min_stages) <- c("Nematostella_vectensis","Strongylocentrotus_purpuratus",
                          "Crassostrea_gigas","Branchiostoma_lanceolatum",
                          "Danio_rerio", "Caenorhabditis_elegans",
                          "Drosophila_melanogaster","Owenia_fusiformis",
                          "Amphimedon_queenslandica", "Clytia_hemisphaerica")
rownames(min_stages) <- c("64 cells","gastrula","stage 4tt larva","stage 5 larva",
                          "stage 7 larva","pre-competent larva","competent larva")
min_stages$Nematostella_vectensis <- names(Ctel2Nvec_transposed)[apply(Ctel2Nvec_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Strongylocentrotus_purpuratus <- names(Ctel2Spur_transposed)[apply(Ctel2Spur_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Crassostrea_gigas <- names(Ctel2Cgig_transposed)[apply(Ctel2Cgig_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Branchiostoma_lanceolatum <- names(Ctel2Blan_transposed)[apply(Ctel2Blan_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Danio_rerio <- names(Ctel2Drer_transposed)[apply(Ctel2Drer_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Drosophila_melanogaster <- names(Ctel2Dmel_transposed)[apply(Ctel2Dmel_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Caenorhabditis_elegans <- names(Ctel2Cele_transposed)[apply(Ctel2Cele_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Amphimedon_queenslandica <- names(Ctel2Aque_transposed)[apply(Ctel2Aque_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Clytia_hemisphaerica <- names(Ctel2Chem_transposed)[apply(Ctel2Chem_transposed, MARGIN = 1, FUN = which.min)]
min_stages$Owenia_fusiformis <- names(Ctel2Ofus_transposed)[apply(Ctel2Ofus_transposed, MARGIN = 1, FUN = which.min)]

minimum_distance_mean <- data.frame(replicate(10,sample(0:1,ncol(Ctel2Drer_mean),rep=TRUE)))
minimum_distance_sd <- data.frame(replicate(10,sample(0:1,ncol(Ctel2Drer_mean),rep=TRUE)))
colnames(minimum_distance_mean) <- c("Nematostella_vectensis","Strongylocentrotus_purpuratus",
                                     "Crassostrea_gigas","Branchiostoma_lanceolatum",
                                     "Danio_rerio", "Caenorhabditis_elegans",
                                     "Drosophila_melanogaster","Owenia_fusiformis",
                                     "Amphimedon_queenslandica", "Clytia_hemisphaerica")
rownames(minimum_distance_mean) <- c("64 cells","gastrula","stage 4tt larva","stage 5 larva",
                                     "stage 7 larva","pre-competent larva","competent larva")
colnames(minimum_distance_sd) <- c("Nematostella_vectensis","Strongylocentrotus_purpuratus",
                                   "Crassostrea_gigas","Branchiostoma_lanceolatum",
                                   "Danio_rerio", "Caenorhabditis_elegans",
                                   "Drosophila_melanogaster","Owenia_fusiformis",
                                   "Amphimedon_queenslandica", "Clytia_hemisphaerica")
rownames(minimum_distance_sd) <- c("64 cells","gastrula","stage 4tt larva","stage 5 larva",
                                   "stage 7 larva","pre-competent larva","competent larva")

for (i in c(1:ncol(Ctel2Drer_mean))){
  minimum_distance_mean[i,"Nematostella_vectensis"] <- Ctel2Nvec_mean[min_stages$Nematostella_vectensis[i],i]
  minimum_distance_mean[i,"Strongylocentrotus_purpuratus"] <- Ctel2Spur_mean[min_stages$Strongylocentrotus_purpuratus[i],i]
  minimum_distance_mean[i,"Crassostrea_gigas"] <- Ctel2Cgig_mean[min_stages$Crassostrea_gigas[i],i]
  minimum_distance_mean[i,"Branchiostoma_lanceolatum"] <- Ctel2Blan_mean[min_stages$Branchiostoma_lanceolatum[i],i]
  minimum_distance_mean[i,"Danio_rerio"] <- Ctel2Drer_mean[min_stages$Danio_rerio[i],i]
  minimum_distance_mean[i,"Drosophila_melanogaster"] <- Ctel2Dmel_mean[min_stages$Drosophila_melanogaster[i],i]
  minimum_distance_mean[i,"Caenorhabditis_elegans"] <- Ctel2Cele_mean[min_stages$Caenorhabditis_elegans[i],i]
  minimum_distance_mean[i,"Amphimedon_queenslandica"] <- Ctel2Aque_mean[min_stages$Amphimedon_queenslandica[i],i]
  minimum_distance_mean[i,"Clytia_hemisphaerica"] <- Ctel2Chem_mean[min_stages$Clytia_hemisphaerica[i],i]
  minimum_distance_mean[i,"Owenia_fusiformis"] <- Ctel2Ofus_mean[min_stages$Owenia_fusiformis[i],i]
  minimum_distance_sd[i,"Nematostella_vectensis"] <- Ctel2Nvec_sd[min_stages$Nematostella_vectensis[i],i]
  minimum_distance_sd[i,"Strongylocentrotus_purpuratus"] <- Ctel2Spur_sd[min_stages$Strongylocentrotus_purpuratus[i],i]
  minimum_distance_sd[i,"Crassostrea_gigas"] <- Ctel2Cgig_sd[min_stages$Crassostrea_gigas[i],i]
  minimum_distance_sd[i,"Branchiostoma_lanceolatum"] <- Ctel2Blan_sd[min_stages$Branchiostoma_lanceolatum[i],i]
  minimum_distance_sd[i,"Danio_rerio"] <- Ctel2Drer_sd[min_stages$Danio_rerio[i],i]
  minimum_distance_sd[i,"Drosophila_melanogaster"] <- Ctel2Dmel_sd[min_stages$Drosophila_melanogaster[i],i]
  minimum_distance_sd[i,"Caenorhabditis_elegans"] <- Ctel2Cele_sd[min_stages$Caenorhabditis_elegans[i],i]
  minimum_distance_sd[i,"Amphimedon_queenslandica"] <- Ctel2Aque_sd[min_stages$Amphimedon_queenslandica[i],i]
  minimum_distance_sd[i,"Clytia_hemisphaerica"] <- Ctel2Chem_sd[min_stages$Clytia_hemisphaerica[i],i]
  minimum_distance_sd[i,"Owenia_fusiformis"] <- Ctel2Ofus_sd[min_stages$Owenia_fusiformis[i],i]
}

minimum_distance_mean$stage <- c(1:ncol(Ctel2Drer_mean))
minimum_distance_sd$stage <- c(1:ncol(Ctel2Drer_mean))
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
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Amphimedon_queenslandica'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Amphimedon_queenslandica'),3:5]-min(minimum_distance_mean$Amphimedon_queenslandica))/(max(minimum_distance_mean$Amphimedon_queenslandica)-min(minimum_distance_mean$Amphimedon_queenslandica))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Clytia_hemisphaerica'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Clytia_hemisphaerica'),3:5]-min(minimum_distance_mean$Clytia_hemisphaerica))/(max(minimum_distance_mean$Clytia_hemisphaerica)-min(minimum_distance_mean$Clytia_hemisphaerica))
minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Owenia_fusiformis'),3:5] <- (minimum_distance_mean_tidy[which(minimum_distance_mean_tidy$species=='Owenia_fusiformis'),3:5]-min(minimum_distance_mean$Owenia_fusiformis))/(max(minimum_distance_mean$Owenia_fusiformis)-min(minimum_distance_mean$Owenia_fusiformis))

final_dataset <- minimum_distance_mean_tidy
final_dataset$species <- factor(final_dataset$species,
                                levels = c("Owenia_fusiformis","Crassostrea_gigas",
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
  scale_x_continuous(labels = c(1:ncol(Ctel2Drer_mean)), breaks = c(1:ncol(Ctel2Drer_mean))) +
  theme_classic() +
  labs(x = "C. teleta stage", y = "Normalised gene expression divergence (JSD)")

write.table(min_stages, "minimum_stages.txt", sep='\t', quote = FALSE)
write.table(minimum_distance_mean, "minimum_stages_distance_mean.txt", sep='\t', quote = FALSE)
write.table(minimum_distance_sd, "minimum_stages_distance_sd.txt", sep='\t', quote = FALSE)
write.table(minimum_distance_mean_tidy, "minimum_stages_to_plot.txt", sep='\t', quote = FALSE)

