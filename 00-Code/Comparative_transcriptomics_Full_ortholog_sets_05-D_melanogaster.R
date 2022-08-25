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
## Drosophila melanogaster 1-to-1 comparisons ##
##########################################

# 1. Drosophila melanogaster vs. Nematostella vectensis
# Number of orthologues = 5,254 orthologues (-4)
# 8 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/Quantile_transformation")
Dmel2Nvec_orth_tpm_qn <- read.csv("Dmel2Nvec_TPM_mean_quantile_transform.csv", header = T)
Nvec2Dmel_orth_tpm_qn <- read.csv("Nvec2Dmel_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/one2one")
both_ID <- read.table("Dmel2Nvec.txt", header = T)
match_Dmel <- Dmel2Nvec_orth_tpm_qn[match(both_ID$Dmel,Dmel2Nvec_orth_tpm_qn$Gene_ID),]
match_Nvec <- Nvec2Dmel_orth_tpm_qn[match(both_ID$Nvec,Nvec2Dmel_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/order")
write.table(match_Dmel,"order_Dmel2Nvec_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Nvec,"order_Nvec2Dmel_tpm_qn.txt", sep="\t", quote=F)
match_Dmel <- read.table("order_Dmel2Nvec_tpm_qn.txt", header = T)
match_Nvec <- read.table("order_Nvec2Dmel_tpm_qn.txt", header = T)

Dmel2Nvec_full_JSD <- data.frame(replicate(ncol(match_Dmel)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Dmel2Nvec_full_JSD) <- colnames(Dmel2Nvec_orth_tpm_qn)[-c(1)]
rownames(Dmel2Nvec_full_JSD) <- colnames(Nvec2Dmel_orth_tpm_qn)[-c(1)]
Dmel2Nvec_mean_JSD <- data.frame(replicate(ncol(match_Dmel)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Dmel2Nvec_mean_JSD) <- colnames(Dmel2Nvec_orth_tpm_qn)[-c(1)]
rownames(Dmel2Nvec_mean_JSD) <- colnames(Nvec2Dmel_orth_tpm_qn)[-c(1)]
Dmel2Nvec_sd_JSD <- data.frame(replicate(ncol(match_Dmel)-1,sample(0:1,ncol(match_Nvec)-1,rep=TRUE)))
colnames(Dmel2Nvec_sd_JSD) <- colnames(Dmel2Nvec_orth_tpm_qn)[-c(1)]
rownames(Dmel2Nvec_sd_JSD) <- colnames(Nvec2Dmel_orth_tpm_qn)[-c(1)]

for (i in colnames(Dmel2Nvec_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Nvec2Dmel_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Dmel[[i]],match_Nvec[[j]])
    Dmel2Nvec_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Dmel2Nvec_mean_JSD[j,i] <- mean(all_JSD)
    Dmel2Nvec_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/Results")
write.table(Dmel2Nvec_full_JSD, "Dmel2Nvec_full_set_JSD.txt", sep ='\t')
write.table(Dmel2Nvec_mean_JSD, "Dmel2Nvec_subset_JSD_mean.txt", sep ='\t')
write.table(Dmel2Nvec_sd_JSD, "Dmel2Nvec_subset_JSD_sd.txt", sep ='\t')


# 2. Drosophila melanogaster vs. Strongylocentrotus purpuratus
# Number of orthologues = 5,015 orthologues (+0)
# 14 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/Quantile_transformation")
Dmel2Spur_orth_tpm_qn <- read.csv("Dmel2Spur_TPM_mean_quantile_transform.csv", header = T)
Spur2Dmel_orth_tpm_qn <- read.csv("Spur2Dmel_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/one2one")
both_ID <- read.table("Dmel2Spur.txt", header = T)
match_Dmel <- Dmel2Spur_orth_tpm_qn[match(both_ID$Dmel,Dmel2Spur_orth_tpm_qn$Gene_ID),]
match_Dmel <- na.omit(match_Dmel)
match_Spur <- Spur2Dmel_orth_tpm_qn[match(both_ID$Spur,Spur2Dmel_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/order")
write.table(match_Dmel,"order_Dmel2Spur_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Spur,"order_Spur2Dmel_tpm_qn.txt", sep="\t", quote=F)
match_Dmel <- read.table("order_Dmel2Spur_tpm_qn.txt", header = T)
match_Spur <- read.table("order_Spur2Dmel_tpm_qn.txt", header = T)

Dmel2Spur_full_JSD <- data.frame(replicate(ncol(match_Dmel)-1,sample(0:1,ncol(match_Spur)-1,rep=TRUE)))
colnames(Dmel2Spur_full_JSD) <- colnames(Dmel2Spur_orth_tpm_qn)[-c(1)]
rownames(Dmel2Spur_full_JSD) <- colnames(Spur2Dmel_orth_tpm_qn)[-c(1)]
Dmel2Spur_mean_JSD <- data.frame(replicate(ncol(match_Dmel)-1,sample(0:1,ncol(match_Spur)-1,rep=TRUE)))
colnames(Dmel2Spur_mean_JSD) <- colnames(Dmel2Spur_orth_tpm_qn)[-c(1)]
rownames(Dmel2Spur_mean_JSD) <- colnames(Spur2Dmel_orth_tpm_qn)[-c(1)]
Dmel2Spur_sd_JSD <- data.frame(replicate(ncol(match_Dmel)-1,sample(0:1,ncol(match_Spur)-1,rep=TRUE)))
colnames(Dmel2Spur_sd_JSD) <- colnames(Dmel2Spur_orth_tpm_qn)[-c(1)]
rownames(Dmel2Spur_sd_JSD) <- colnames(Spur2Dmel_orth_tpm_qn)[-c(1)]

for (i in colnames(Dmel2Spur_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Spur2Dmel_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Dmel[[i]],match_Spur[[j]])
    Dmel2Spur_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Dmel2Spur_mean_JSD[j,i] <- mean(all_JSD)
    Dmel2Spur_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/Results")
write.table(Dmel2Spur_full_JSD, "Dmel2Spur_full_set_JSD.txt", sep ='\t')
write.table(Dmel2Spur_mean_JSD, "Dmel2Spur_subset_JSD_mean.txt", sep ='\t')
write.table(Dmel2Spur_sd_JSD, "Dmel2Spur_subset_JSD_sd.txt", sep ='\t')


# 3. Drosophila melanogaster vs. Branchiostoma lanceolatum
# Number of orthologues = 6,673 orthologues (+39)
# 14 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/Quantile_transformation")
Dmel2Blan_orth_tpm_qn <- read.csv("Dmel2Blan_TPM_mean_quantile_transform.csv", header = T)
Blan2Dmel_orth_tpm_qn <- read.csv("Blan2Dmel_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/one2one")
both_ID <- read.table("Dmel2Blan.txt", header = T)
match_Dmel <- Dmel2Blan_orth_tpm_qn[match(both_ID$Dmel,Dmel2Blan_orth_tpm_qn$Gene_ID),]
match_Dmel <- na.omit(match_Dmel)
match_Blan <- Blan2Dmel_orth_tpm_qn[match(both_ID$Blan,Blan2Dmel_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/order")
write.table(match_Dmel,"order_Dmel2Blan_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Blan,"order_Blan2Dmel_tpm_qn.txt", sep="\t", quote=F)
match_Dmel <- read.table("order_Dmel2Blan_tpm_qn.txt", header = T)
match_Blan <- read.table("order_Blan2Dmel_tpm_qn.txt", header = T)

Dmel2Blan_full_JSD <- data.frame(replicate(ncol(match_Dmel)-1,sample(0:1,ncol(match_Blan)-1,rep=TRUE)))
colnames(Dmel2Blan_full_JSD) <- colnames(Dmel2Blan_orth_tpm_qn)[-c(1)]
rownames(Dmel2Blan_full_JSD) <- colnames(Blan2Dmel_orth_tpm_qn)[-c(1)]
Dmel2Blan_mean_JSD <- data.frame(replicate(ncol(match_Dmel)-1,sample(0:1,ncol(match_Blan)-1,rep=TRUE)))
colnames(Dmel2Blan_mean_JSD) <- colnames(Dmel2Blan_orth_tpm_qn)[-c(1)]
rownames(Dmel2Blan_mean_JSD) <- colnames(Blan2Dmel_orth_tpm_qn)[-c(1)]
Dmel2Blan_sd_JSD <- data.frame(replicate(ncol(match_Dmel)-1,sample(0:1,ncol(match_Blan)-1,rep=TRUE)))
colnames(Dmel2Blan_sd_JSD) <- colnames(Dmel2Blan_orth_tpm_qn)[-c(1)]
rownames(Dmel2Blan_sd_JSD) <- colnames(Blan2Dmel_orth_tpm_qn)[-c(1)]

for (i in colnames(Dmel2Blan_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Blan2Dmel_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Dmel[[i]],match_Blan[[j]])
    Dmel2Blan_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Dmel2Blan_mean_JSD[j,i] <- mean(all_JSD)
    Dmel2Blan_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/Results")
write.table(Dmel2Blan_full_JSD, "Dmel2Blan_full_set_JSD.txt", sep ='\t')
write.table(Dmel2Blan_mean_JSD, "Dmel2Blan_subset_JSD_mean.txt", sep ='\t')
write.table(Dmel2Blan_sd_JSD, "Dmel2Blan_subset_JSD_sd.txt", sep ='\t')


# 4. Drosophila melanogaster vs. Danio rerio
# Number of orthologues = 4,316 orthologues (+21)
# 7 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/Quantile_transformation")
Dmel2Drer_orth_tpm_qn <- read.csv("Dmel2Drer_TPM_mean_quantile_transform.csv", header = T)
Drer2Dmel_orth_tpm_qn <- read.csv("Drer2Dmel_TPM_mean_quantile_transform.csv", header = T)
#Drer2Dmel_orth_tpm_qn <- cbind(Drer2Dmel_orth_tpm_qn_preliminary[1],Drer2Dmel_orth_tpm_qn_preliminary[12],Drer2Dmel_orth_tpm_qn_preliminary[c(3:11)],Drer2Dmel_orth_tpm_qn_preliminary[2],Drer2Dmel_orth_tpm_qn_preliminary[c(13:19)])


setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/one2one")
both_ID <- read.table("Dmel2Drer.txt", header = T)
match_Dmel <- Dmel2Drer_orth_tpm_qn[match(both_ID$Dmel,Dmel2Drer_orth_tpm_qn$Gene_ID),]
match_Dmel <- na.omit(match_Dmel)
match_Drer <- Drer2Dmel_orth_tpm_qn[match(both_ID$Drer,Drer2Dmel_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/order")
write.table(match_Dmel,"order_Dmel2Drer_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Drer,"order_Drer2Dmel_tpm_qn.txt", sep="\t", quote=F)
match_Dmel <- read.table("order_Dmel2Drer_tpm_qn.txt", header = T)
match_Drer <- read.table("order_Drer2Dmel_tpm_qn.txt", header = T)

Dmel2Drer_full_JSD <- data.frame(replicate(ncol(match_Dmel)-1,sample(0:1,ncol(match_Drer)-1,rep=TRUE)))
colnames(Dmel2Drer_full_JSD) <- colnames(Dmel2Drer_orth_tpm_qn)[-c(1)]
rownames(Dmel2Drer_full_JSD) <- colnames(Drer2Dmel_orth_tpm_qn)[-c(1)]
Dmel2Drer_mean_JSD <- data.frame(replicate(ncol(match_Dmel)-1,sample(0:1,ncol(match_Drer)-1,rep=TRUE)))
colnames(Dmel2Drer_mean_JSD) <- colnames(Dmel2Drer_orth_tpm_qn)[-c(1)]
rownames(Dmel2Drer_mean_JSD) <- colnames(Drer2Dmel_orth_tpm_qn)[-c(1)]
Dmel2Drer_sd_JSD <- data.frame(replicate(ncol(match_Dmel)-1,sample(0:1,ncol(match_Drer)-1,rep=TRUE)))
colnames(Dmel2Drer_sd_JSD) <- colnames(Dmel2Drer_orth_tpm_qn)[-c(1)]
rownames(Dmel2Drer_sd_JSD) <- colnames(Drer2Dmel_orth_tpm_qn)[-c(1)]

for (i in colnames(Dmel2Drer_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Drer2Dmel_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Dmel[[i]],match_Drer[[j]])
    Dmel2Drer_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Dmel2Drer_mean_JSD[j,i] <- mean(all_JSD)
    Dmel2Drer_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/Results")
write.table(Dmel2Drer_full_JSD, "Dmel2Drer_full_set_JSD.txt", sep ='\t')
write.table(Dmel2Drer_mean_JSD, "Dmel2Drer_subset_JSD_mean.txt", sep ='\t')
write.table(Dmel2Drer_sd_JSD, "Dmel2Drer_subset_JSD_sd.txt", sep ='\t')



# 7. Drosophila melanogaster vs. Amphimedon queenslandica
# Number of orthologues = 3,962 orthologues
# 9 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/Quantile_transformation")
Dmel2Aque_orth_tpm_qn <- read.csv("Dmel2Aque_TPM_mean_quantile_transform.csv", header = T)
Aque2Dmel_orth_tpm_qn <- read.csv("Aque2Dmel_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/one2one")
both_ID <- read.table("Dmel2Aque.txt", header = T)
match_Dmel <- Dmel2Aque_orth_tpm_qn[match(both_ID$Dmel,Dmel2Aque_orth_tpm_qn$Gene_ID),]
match_Dmel <- na.omit(match_Dmel)
match_Aque <- Aque2Dmel_orth_tpm_qn[match(both_ID$Aque,Aque2Dmel_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/order")
write.table(match_Dmel,"order_Dmel2Aque_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Aque,"order_Aque2Dmel_tpm_qn.txt", sep="\t", quote=F)
match_Dmel <- read.table("order_Dmel2Aque_tpm_qn.txt", header = T)
match_Aque <- read.table("order_Aque2Dmel_tpm_qn.txt", header = T)

Dmel2Aque_full_JSD <- data.frame(replicate(ncol(match_Dmel)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Dmel2Aque_full_JSD) <- colnames(Dmel2Aque_orth_tpm_qn)[-c(1)]
rownames(Dmel2Aque_full_JSD) <- colnames(Aque2Dmel_orth_tpm_qn)[-c(1)]
Dmel2Aque_mean_JSD <- data.frame(replicate(ncol(match_Dmel)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Dmel2Aque_mean_JSD) <- colnames(Dmel2Aque_orth_tpm_qn)[-c(1)]
rownames(Dmel2Aque_mean_JSD) <- colnames(Aque2Dmel_orth_tpm_qn)[-c(1)]
Dmel2Aque_sd_JSD <- data.frame(replicate(ncol(match_Dmel)-1,sample(0:1,ncol(match_Aque)-1,rep=TRUE)))
colnames(Dmel2Aque_sd_JSD) <- colnames(Dmel2Aque_orth_tpm_qn)[-c(1)]
rownames(Dmel2Aque_sd_JSD) <- colnames(Aque2Dmel_orth_tpm_qn)[-c(1)]

for (i in colnames(Dmel2Aque_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Aque2Dmel_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Dmel[[i]],match_Aque[[j]])
    Dmel2Aque_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Dmel2Aque_mean_JSD[j,i] <- mean(all_JSD)
    Dmel2Aque_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/Results")
write.table(Dmel2Aque_full_JSD, "Dmel2Aque_full_set_JSD.txt", sep ='\t')
write.table(Dmel2Aque_mean_JSD, "Dmel2Aque_subset_JSD_mean.txt", sep ='\t')
write.table(Dmel2Aque_sd_JSD, "Dmel2Aque_subset_JSD_sd.txt", sep ='\t')


# 8. Drosophila melanogaster vs. Clytica hemisphaerica
# Number of orthologues = 4,691 orthologues
# 9 stages

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/Quantile_transformation")
Dmel2Chem_orth_tpm_qn <- read.csv("Dmel2Chem_TPM_mean_quantile_transform.csv", header = T)
Chem2Dmel_orth_tpm_qn <- read.csv("Chem2Dmel_TPM_mean_quantile_transform.csv", header = T)

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/one2one")
both_ID <- read.table("Dmel2Chem.txt", header = T)
match_Dmel <- Dmel2Chem_orth_tpm_qn[match(both_ID$Dmel,Dmel2Chem_orth_tpm_qn$Gene_ID),]
match_Dmel <- na.omit(match_Dmel)
match_Chem <- Chem2Dmel_orth_tpm_qn[match(both_ID$Chem,Chem2Dmel_orth_tpm_qn$Gene_ID),]
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/order")
write.table(match_Dmel,"order_Dmel2Chem_tpm_qn.txt", sep="\t", quote=F)
write.table(match_Chem,"order_Chem2Dmel_tpm_qn.txt", sep="\t", quote=F)
match_Dmel <- read.table("order_Dmel2Chem_tpm_qn.txt", header = T)
match_Chem <- read.table("order_Chem2Dmel_tpm_qn.txt", header = T)

Dmel2Chem_full_JSD <- data.frame(replicate(ncol(match_Dmel)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Dmel2Chem_full_JSD) <- colnames(Dmel2Chem_orth_tpm_qn)[-c(1)]
rownames(Dmel2Chem_full_JSD) <- colnames(Chem2Dmel_orth_tpm_qn)[-c(1)]
Dmel2Chem_mean_JSD <- data.frame(replicate(ncol(match_Dmel)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Dmel2Chem_mean_JSD) <- colnames(Dmel2Chem_orth_tpm_qn)[-c(1)]
rownames(Dmel2Chem_mean_JSD) <- colnames(Chem2Dmel_orth_tpm_qn)[-c(1)]
Dmel2Chem_sd_JSD <- data.frame(replicate(ncol(match_Dmel)-1,sample(0:1,ncol(match_Chem)-1,rep=TRUE)))
colnames(Dmel2Chem_sd_JSD) <- colnames(Dmel2Chem_orth_tpm_qn)[-c(1)]
rownames(Dmel2Chem_sd_JSD) <- colnames(Chem2Dmel_orth_tpm_qn)[-c(1)]

for (i in colnames(Dmel2Chem_orth_tpm_qn)[-c(1)]){
  for (j in colnames(Chem2Dmel_orth_tpm_qn)[-c(1)]){
    all_JSD <- vector()
    match <- rbind(match_Dmel[[i]],match_Chem[[j]])
    Dmel2Chem_full_JSD[j,i] <- JSD(match)
    for (k in c(1:250)){
      sub_match <- match[,sample(ncol(match),replace=TRUE)]
      all_JSD[k] <- JSD(sub_match)
    }
    Dmel2Chem_mean_JSD[j,i] <- mean(all_JSD)
    Dmel2Chem_sd_JSD[j,i] <- sd(all_JSD)
  }
}

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/Results")
write.table(Dmel2Chem_full_JSD, "Dmel2Chem_full_set_JSD.txt", sep ='\t')
write.table(Dmel2Chem_mean_JSD, "Dmel2Chem_subset_JSD_mean.txt", sep ='\t')
write.table(Dmel2Chem_sd_JSD, "Dmel2Chem_subset_JSD_sd.txt", sep ='\t')