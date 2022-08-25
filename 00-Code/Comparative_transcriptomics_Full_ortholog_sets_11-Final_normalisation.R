library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(corrplot)
library(philentropy)
library(tidyr)
library(dplyr)
library(preprocessCore)
library(plyr)
library(circlize)
library(grid)
library(gridGraphics)
library(pheatmap)


#######################################
## All-to-all normalised comparisons ##
#######################################

# Normalise results by the number of orthologues to make 
# the 1-to-1 comparisons comparable: heatmaps of NORMALISED data

# 1. Owenia fusiformis vs. 10 species
# Number of orthologues in order: Nvec, Spur, Cgig, Blan, Drer, Dmel, Cele, Ctel, Aque, Chem
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/01-O_fusiformis/Results")
Ofus_orthologues <- c(5254,5015,6737,6673,4316,4635,3767,7651,3962,4691)

Ofus2Nvec_mean <- read.table("Ofus2Nvec_subset_JSD_mean.txt", header = T)[-c(1:7)]
Ofus2Spur_mean <- read.table("Ofus2Spur_subset_JSD_mean.txt", header = T)[-c(1:7)]
Ofus2Cgig_mean <- read.table("Ofus2Cgig_subset_JSD_mean.txt", header = T)[-c(1:7)]
Ofus2Blan_mean <- read.table("Ofus2Blan_subset_JSD_mean.txt", header = T)[-c(1:7)]
Ofus2Drer_mean <- read.table("Ofus2Drer_subset_JSD_mean.txt", header = T)[-c(1:7)]
Ofus2Dmel_mean <- read.table("Ofus2Dmel_subset_JSD_mean.txt", header = T)[,-c(1:7)]
Ofus2Cele_mean <- read.table("Ofus2Cele_subset_JSD_mean.txt", header = T)[-c(12),-c(1:7)]
Ofus2Ctel_mean <- read.table("Ofus2Ctel_subset_JSD_mean.txt", header = T)[-c(1:7),-c(1:7)]
Ofus2Aque_mean <- read.table("Ofus2Aque_subset_JSD_mean.txt", header = T)[-c(1:7)]
Ofus2Chem_mean <- read.table("Ofus2Chem_subset_JSD_mean.txt", header = T)[-c(1:7)]

Ofus2Nvec_mean_divided <- Ofus2Nvec_mean/Ofus_orthologues[1]
Ofus2Spur_mean_divided <- Ofus2Spur_mean/Ofus_orthologues[2]
Ofus2Cgig_mean_divided <- Ofus2Cgig_mean/Ofus_orthologues[3]
Ofus2Blan_mean_divided <- Ofus2Blan_mean/Ofus_orthologues[4]
Ofus2Drer_mean_divided <- Ofus2Drer_mean/Ofus_orthologues[5]
Ofus2Dmel_mean_divided <- Ofus2Dmel_mean/Ofus_orthologues[6]
Ofus2Cele_mean_divided <- Ofus2Cele_mean/Ofus_orthologues[7]
Ofus2Ctel_mean_divided <- Ofus2Ctel_mean/Ofus_orthologues[8]
Ofus2Aque_mean_divided <- Ofus2Aque_mean/Ofus_orthologues[9]
Ofus2Chem_mean_divided <- Ofus2Chem_mean/Ofus_orthologues[10]


# 2. Capitella teleta vs. 9 species
# Number of orthologues in order: Nvec, Spur, Cgig, Blan, Drer, Dmel, Cele, Aque, Chem
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/02-C_teleta/Results")
Ctel_orthologues <- c(4945,4651,6108,6073,3951,4263,3411,3703,4246)

Ctel2Nvec_mean <- read.table("Ctel2Nvec_subset_JSD_mean.txt", header = T)[-c(1:7)]
Ctel2Spur_mean <- read.table("Ctel2Spur_subset_JSD_mean.txt", header = T)[-c(1:7)]
Ctel2Cgig_mean <- read.table("Ctel2Cgig_subset_JSD_mean.txt", header = T)[-c(1:7)]
Ctel2Blan_mean <- read.table("Ctel2Blan_subset_JSD_mean.txt", header = T)[-c(1:7)]
Ctel2Drer_mean <- read.table("Ctel2Drer_subset_JSD_mean.txt", header = T)[-c(1:7)]
Ctel2Dmel_mean <- read.table("Ctel2Dmel_subset_JSD_mean.txt", header = T)[,-c(1:7)]
Ctel2Cele_mean <- read.table("Ctel2Cele_subset_JSD_mean.txt", header = T)[-c(12),-c(1:7)]
Ctel2Aque_mean <- read.table("Ctel2Aque_subset_JSD_mean.txt", header = T)[-c(1:7)]
Ctel2Chem_mean <- read.table("Ctel2Chem_subset_JSD_mean.txt", header = T)[-c(1:7)]

Ctel2Nvec_mean_divided <- Ctel2Nvec_mean/Ctel_orthologues[1]
Ctel2Spur_mean_divided <- Ctel2Spur_mean/Ctel_orthologues[2]
Ctel2Cgig_mean_divided <- Ctel2Cgig_mean/Ctel_orthologues[3]
Ctel2Blan_mean_divided <- Ctel2Blan_mean/Ctel_orthologues[4]
Ctel2Drer_mean_divided <- Ctel2Drer_mean/Ctel_orthologues[5]
Ctel2Dmel_mean_divided <- Ctel2Dmel_mean/Ctel_orthologues[6]
Ctel2Cele_mean_divided <- Ctel2Cele_mean/Ctel_orthologues[7]
Ctel2Aque_mean_divided <- Ctel2Aque_mean/Ctel_orthologues[8]
Ctel2Chem_mean_divided <- Ctel2Chem_mean/Ctel_orthologues[9]


# 3. Crassostrea gigas vs. 8 species
# Number of orthologues in order: Nvec, Spur, Blan, Drer, Dmel, Cele, Aque, Chem
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/03-C_gigas/Results")
Cgig_orthologues <- c(4334,4199,5388,3424,3792,3053,3224,3786)

Cgig2Nvec_mean <- read.table("Cgig2Nvec_subset_JSD_mean.txt", header = T)
Cgig2Spur_mean <- read.table("Cgig2Spur_subset_JSD_mean.txt", header = T)
Cgig2Blan_mean <- read.table("Cgig2Blan_subset_JSD_mean.txt", header = T)
Cgig2Drer_mean <- read.table("Cgig2Drer_subset_JSD_mean.txt", header = T)
Cgig2Dmel_mean <- read.table("Cgig2Dmel_subset_JSD_mean.txt", header = T)[,]
Cgig2Cele_mean <- read.table("Cgig2Cele_subset_JSD_mean.txt", header = T)[-c(12),]
Cgig2Aque_mean <- read.table("Cgig2Aque_subset_JSD_mean.txt", header = T)
Cgig2Chem_mean <- read.table("Cgig2Chem_subset_JSD_mean.txt", header = T)

Cgig2Nvec_mean_divided <- Cgig2Nvec_mean/Cgig_orthologues[1]
Cgig2Spur_mean_divided <- Cgig2Spur_mean/Cgig_orthologues[2]
Cgig2Blan_mean_divided <- Cgig2Blan_mean/Cgig_orthologues[3]
Cgig2Drer_mean_divided <- Cgig2Drer_mean/Cgig_orthologues[4]
Cgig2Dmel_mean_divided <- Cgig2Dmel_mean/Cgig_orthologues[5]
Cgig2Cele_mean_divided <- Cgig2Cele_mean/Cgig_orthologues[6]
Cgig2Aque_mean_divided <- Cgig2Aque_mean/Cgig_orthologues[7]
Cgig2Chem_mean_divided <- Cgig2Chem_mean/Cgig_orthologues[8]


# 4. Caenorhabditis elegans vs. 7 species
# Number of orthologues in order: Nvec, Spur, Blan, Drer, Dmel, Aque, Chem
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/04-C_elegans/Results")
Cele_orthologues <- c(2890,2640,3497,2486,3432,2631,2943)

Cele2Nvec_mean <- read.table("Cele2Nvec_subset_JSD_mean.txt", header = T)[-c(12)]
Cele2Spur_mean <- read.table("Cele2Spur_subset_JSD_mean.txt", header = T)[-c(12)]
Cele2Blan_mean <- read.table("Cele2Blan_subset_JSD_mean.txt", header = T)[-c(12)]
Cele2Drer_mean <- read.table("Cele2Drer_subset_JSD_mean.txt", header = T)[-c(12)]
Cele2Dmel_mean <- read.table("Cele2Dmel_subset_JSD_mean.txt", header = T)[,-c(12)]
Cele2Aque_mean <- read.table("Cele2Aque_subset_JSD_mean.txt", header = T)[-c(12)]
Cele2Chem_mean <- read.table("Cele2Chem_subset_JSD_mean.txt", header = T)[-c(12)]

Cele2Nvec_mean_divided <- Cele2Nvec_mean/Cele_orthologues[1]
Cele2Spur_mean_divided <- Cele2Spur_mean/Cele_orthologues[2]
Cele2Blan_mean_divided <- Cele2Blan_mean/Cele_orthologues[3]
Cele2Drer_mean_divided <- Cele2Drer_mean/Cele_orthologues[4]
Cele2Dmel_mean_divided <- Cele2Dmel_mean/Cele_orthologues[5]
Cele2Aque_mean_divided <- Cele2Aque_mean/Cele_orthologues[6]
Cele2Chem_mean_divided <- Cele2Chem_mean/Cele_orthologues[7]


# 5. Drosophila melanogaster vs. 6 species
# Number of orthologues in order: Nvec, Spur, Blan, Drer, Aque, Chem
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/05-D_melanogaster/Results")
Dmel_orthologues <- c(3473,3241,4197,3021,3003,3359)

Dmel2Nvec_mean <- read.table("Dmel2Nvec_subset_JSD_mean.txt", header = T)[]
Dmel2Spur_mean <- read.table("Dmel2Spur_subset_JSD_mean.txt", header = T)[]
Dmel2Blan_mean <- read.table("Dmel2Blan_subset_JSD_mean.txt", header = T)[]
Dmel2Drer_mean <- read.table("Dmel2Drer_subset_JSD_mean.txt", header = T)[]
Dmel2Aque_mean <- read.table("Dmel2Aque_subset_JSD_mean.txt", header = T)[]
Dmel2Chem_mean <- read.table("Dmel2Chem_subset_JSD_mean.txt", header = T)[]

Dmel2Nvec_mean_divided <- Dmel2Nvec_mean/Dmel_orthologues[1]
Dmel2Spur_mean_divided <- Dmel2Spur_mean/Dmel_orthologues[2]
Dmel2Blan_mean_divided <- Dmel2Blan_mean/Dmel_orthologues[3]
Dmel2Drer_mean_divided <- Dmel2Drer_mean/Dmel_orthologues[4]
Dmel2Aque_mean_divided <- Dmel2Aque_mean/Dmel_orthologues[5]
Dmel2Chem_mean_divided <- Dmel2Chem_mean/Dmel_orthologues[6]


# 6. Danio rerio vs. 5 species
# Number of orthologues in order: Nvec, Spur, Blan, Aque, Chem
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/06-D_rerio/Results")
Drer_orthologues <- c(3479,3162,4314,2816,3001)

Drer2Nvec_mean <- read.table("Drer2Nvec_subset_JSD_mean.txt", header = T)
Drer2Spur_mean <- read.table("Drer2Spur_subset_JSD_mean.txt", header = T)
Drer2Blan_mean <- read.table("Drer2Blan_subset_JSD_mean.txt", header = T)
Drer2Aque_mean <- read.table("Drer2Aque_subset_JSD_mean.txt", header = T)
Drer2Chem_mean <- read.table("Drer2Chem_subset_JSD_mean.txt", header = T)

Drer2Nvec_mean_divided <- Drer2Nvec_mean/Drer_orthologues[1]
Drer2Spur_mean_divided <- Drer2Spur_mean/Drer_orthologues[2]
Drer2Blan_mean_divided <- Drer2Blan_mean/Drer_orthologues[3]
Drer2Aque_mean_divided <- Drer2Aque_mean/Drer_orthologues[4]
Drer2Chem_mean_divided <- Drer2Chem_mean/Drer_orthologues[5]


# 7. Branchiostoma lanceolatum vs. 4 species
# Number of orthologues in order: Nvec, Spur, Aque, Chem
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/07-B_lanceolatum/Results")
Blan_orthologues <- c(5148,5047,3926,4385)

Blan2Nvec_mean <- read.table("Blan2Nvec_subset_JSD_mean.txt", header = T)
Blan2Spur_mean <- read.table("Blan2Spur_subset_JSD_mean.txt", header = T)
Blan2Aque_mean <- read.table("Blan2Aque_subset_JSD_mean.txt", header = T)
Blan2Chem_mean <- read.table("Blan2Chem_subset_JSD_mean.txt", header = T)

Blan2Nvec_mean_divided <- Blan2Nvec_mean/Blan_orthologues[1]
Blan2Spur_mean_divided <- Blan2Spur_mean/Blan_orthologues[2]
Blan2Aque_mean_divided <- Blan2Aque_mean/Blan_orthologues[3]
Blan2Chem_mean_divided <- Blan2Chem_mean/Blan_orthologues[4]


# 8. Strongylocentrotus purpuratus vs. 3 species
# Number of orthologues in order: Nvec, Aque, Chem
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/08-S_purpuratus/Results")
Spur_orthologues <- c(4045,3021,3340)

Spur2Nvec_mean <- read.table("Spur2Nvec_subset_JSD_mean.txt", header = T)
Spur2Aque_mean <- read.table("Spur2Aque_subset_JSD_mean.txt", header = T)
Spur2Chem_mean <- read.table("Spur2Chem_subset_JSD_mean.txt", header = T)

Spur2Nvec_mean_divided <- Spur2Nvec_mean/Spur_orthologues[1]
Spur2Aque_mean_divided <- Spur2Aque_mean/Spur_orthologues[2]
Spur2Chem_mean_divided <- Spur2Chem_mean/Spur_orthologues[3]


# 9. Nematostella vectensis vs. 2 species
# Number of orthologues in order: Aque, Chem
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/09-N_vectensis/Results")
Nvec_orthologues <- c(3461,4461)

Nvec2Aque_mean <- read.table("Nvec2Aque_subset_JSD_mean.txt", header = T)
Nvec2Chem_mean <- read.table("Nvec2Chem_subset_JSD_mean.txt", header = T)

Nvec2Aque_mean_divided <- Nvec2Aque_mean/Nvec_orthologues[1]
Nvec2Chem_mean_divided <- Nvec2Chem_mean/Nvec_orthologues[2]


# 10. Amphimedon queenslandica vs. 1 species
# Number of orthologues in order: Chem
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/10-A_queenslandica/Results")
Aque_orthologues <- c(3274)

Aque2Chem_mean <- read.table("Aque2Chem_subset_JSD_mean.txt", header = T)

Aque2Chem_mean_divided <- Aque2Chem_mean/Aque_orthologues[1]


# 11. Get minimum and maximum values and normalise
# (JSD-JSDmin)/(JSDmax-JSDmin)

allvalues_divid_1 <- rbind(Ofus2Nvec_mean_divided,Ofus2Spur_mean_divided,Ofus2Cgig_mean_divided,
                           Ofus2Blan_mean_divided,Ofus2Drer_mean_divided,Ofus2Dmel_mean_divided,
                           Ofus2Cele_mean_divided,Ofus2Ctel_mean_divided,Ofus2Aque_mean_divided,
                           Ofus2Chem_mean_divided)
allvalues_divid_2 <- rbind(Ctel2Nvec_mean_divided,Ctel2Spur_mean_divided,Ctel2Cgig_mean_divided,
                           Ctel2Blan_mean_divided,Ctel2Drer_mean_divided,Ctel2Dmel_mean_divided,
                           Ctel2Cele_mean_divided,Ctel2Aque_mean_divided,Ctel2Chem_mean_divided)
allvalues_divid_3 <- rbind(Cgig2Nvec_mean_divided,Cgig2Spur_mean_divided,Cgig2Blan_mean_divided,
                           Cgig2Drer_mean_divided,Cgig2Dmel_mean_divided,Cgig2Cele_mean_divided,
                           Cgig2Aque_mean_divided,Cgig2Chem_mean_divided)
allvalues_divid_4 <- rbind(Cele2Nvec_mean_divided,Cele2Spur_mean_divided,Cele2Blan_mean_divided,
                           Cele2Drer_mean_divided,Cele2Dmel_mean_divided,Cele2Aque_mean_divided,
                           Cele2Chem_mean_divided)
allvalues_divid_5 <- rbind(Dmel2Nvec_mean_divided,Dmel2Spur_mean_divided,Dmel2Blan_mean_divided,
                           Dmel2Drer_mean_divided,Dmel2Aque_mean_divided,Dmel2Chem_mean_divided)
allvalues_divid_6 <- rbind(Drer2Nvec_mean_divided,Drer2Spur_mean_divided,Drer2Blan_mean_divided,
                           Drer2Aque_mean_divided,Drer2Chem_mean_divided)
allvalues_divid_7 <- rbind(Blan2Nvec_mean_divided,Blan2Spur_mean_divided,Blan2Aque_mean_divided,
                           Blan2Chem_mean_divided)
allvalues_divid_8 <- rbind(Spur2Nvec_mean_divided,Spur2Aque_mean_divided,Spur2Chem_mean_divided)
allvalues_divid_9 <- rbind(Nvec2Aque_mean_divided,Nvec2Chem_mean_divided)
allvalues_divid_10 <- rbind(Aque2Chem_mean_divided)

maxima <- vector(length=10)
minima <- vector(length=10)

minima[1] <- min(allvalues_divid_1)
minima[2] <- min(allvalues_divid_2)
minima[3] <- min(allvalues_divid_3)
minima[4] <- min(allvalues_divid_4)
minima[5] <- min(allvalues_divid_5)
minima[6] <- min(allvalues_divid_6)
minima[7] <- min(allvalues_divid_7)
minima[8] <- min(allvalues_divid_8)
minima[9] <- min(allvalues_divid_9)
minima[10] <- min(allvalues_divid_10)

maxima[1] <- max(allvalues_divid_1)
maxima[2] <- max(allvalues_divid_2)
maxima[3] <- max(allvalues_divid_3)
maxima[4] <- max(allvalues_divid_4)
maxima[5] <- max(allvalues_divid_5)
maxima[6] <- max(allvalues_divid_6)
maxima[7] <- max(allvalues_divid_7)
maxima[8] <- max(allvalues_divid_8)
maxima[9] <- max(allvalues_divid_9)
maxima[10] <- max(allvalues_divid_10)

minimum_divided <- min(minima)
maximum_divided <- max(maxima)

Ofus2Nvec_final_normalised <- (Ofus2Nvec_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Ofus2Spur_final_normalised <- (Ofus2Spur_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Ofus2Cgig_final_normalised <- (Ofus2Cgig_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Ofus2Blan_final_normalised <- (Ofus2Blan_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Ofus2Drer_final_normalised <- (Ofus2Drer_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Ofus2Dmel_final_normalised <- (Ofus2Dmel_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Ofus2Cele_final_normalised <- (Ofus2Cele_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Ofus2Ctel_final_normalised <- (Ofus2Ctel_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Ofus2Aque_final_normalised <- (Ofus2Aque_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Ofus2Chem_final_normalised <- (Ofus2Chem_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)

Ctel2Nvec_final_normalised <- (Ctel2Nvec_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Ctel2Spur_final_normalised <- (Ctel2Spur_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Ctel2Cgig_final_normalised <- (Ctel2Cgig_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Ctel2Blan_final_normalised <- (Ctel2Blan_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Ctel2Drer_final_normalised <- (Ctel2Drer_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Ctel2Dmel_final_normalised <- (Ctel2Dmel_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Ctel2Cele_final_normalised <- (Ctel2Cele_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Ctel2Aque_final_normalised <- (Ctel2Aque_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Ctel2Chem_final_normalised <- (Ctel2Chem_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)

Cgig2Nvec_final_normalised <- (Cgig2Nvec_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cgig2Spur_final_normalised <- (Cgig2Spur_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cgig2Blan_final_normalised <- (Cgig2Blan_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cgig2Drer_final_normalised <- (Cgig2Drer_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cgig2Dmel_final_normalised <- (Cgig2Dmel_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cgig2Cele_final_normalised <- (Cgig2Cele_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cgig2Aque_final_normalised <- (Cgig2Aque_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cgig2Chem_final_normalised <- (Cgig2Chem_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)

Cele2Nvec_final_normalised <- (Cele2Nvec_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cele2Spur_final_normalised <- (Cele2Spur_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cele2Blan_final_normalised <- (Cele2Blan_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cele2Drer_final_normalised <- (Cele2Drer_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cele2Dmel_final_normalised <- (Cele2Dmel_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cele2Aque_final_normalised <- (Cele2Aque_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Cele2Chem_final_normalised <- (Cele2Chem_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)

Dmel2Nvec_final_normalised <- (Dmel2Nvec_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Dmel2Spur_final_normalised <- (Dmel2Spur_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Dmel2Blan_final_normalised <- (Dmel2Blan_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Dmel2Drer_final_normalised <- (Dmel2Drer_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Dmel2Aque_final_normalised <- (Dmel2Aque_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Dmel2Chem_final_normalised <- (Dmel2Chem_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)

Drer2Nvec_final_normalised <- (Drer2Nvec_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Drer2Spur_final_normalised <- (Drer2Spur_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Drer2Blan_final_normalised <- (Drer2Blan_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Drer2Aque_final_normalised <- (Drer2Aque_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Drer2Chem_final_normalised <- (Drer2Chem_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)

Blan2Nvec_final_normalised <- (Blan2Nvec_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Blan2Spur_final_normalised <- (Blan2Spur_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Blan2Aque_final_normalised <- (Blan2Aque_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Blan2Chem_final_normalised <- (Blan2Chem_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)

Spur2Nvec_final_normalised <- (Spur2Nvec_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Spur2Aque_final_normalised <- (Spur2Aque_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Spur2Chem_final_normalised <- (Spur2Chem_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)

Nvec2Aque_final_normalised <- (Nvec2Aque_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)
Nvec2Chem_final_normalised <- (Nvec2Chem_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)

Aque2Chem_final_normalised <- (Aque2Chem_mean_divided-minimum_divided)/(maximum_divided-minimum_divided)


# 10. Export final JSD datasets
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/03-Metazoan_comparative_transcriptomics/12-Final_results")

write.table(Ofus2Nvec_final_normalised, paste0("Ofus2Nvec_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Ofus2Spur_final_normalised, paste0("Ofus2Spur_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Ofus2Cgig_final_normalised, paste0("Ofus2Cgig_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Ofus2Blan_final_normalised, paste0("Ofus2Blan_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Ofus2Drer_final_normalised, paste0("Ofus2Drer_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Ofus2Dmel_final_normalised, paste0("Ofus2Dmel_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Ofus2Cele_final_normalised, paste0("Ofus2Cele_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Ofus2Ctel_final_normalised, paste0("Ofus2Ctel_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Ofus2Aque_final_normalised, paste0("Ofus2Aque_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Ofus2Chem_final_normalised, paste0("Ofus2Chem_final_normalised",".txt"), sep='\t', quote = FALSE)

write.table(Ctel2Nvec_final_normalised, paste0("Ctel2Nvec_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Ctel2Spur_final_normalised, paste0("Ctel2Spur_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Ctel2Cgig_final_normalised, paste0("Ctel2Cgig_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Ctel2Blan_final_normalised, paste0("Ctel2Blan_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Ctel2Drer_final_normalised, paste0("Ctel2Drer_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Ctel2Dmel_final_normalised, paste0("Ctel2Dmel_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Ctel2Cele_final_normalised, paste0("Ctel2Cele_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Ctel2Aque_final_normalised, paste0("Ctel2Aque_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Ctel2Chem_final_normalised, paste0("Ctel2Chem_final_normalised",".txt"), sep='\t', quote = FALSE)

write.table(Cgig2Nvec_final_normalised, paste0("Cgig2Nvec_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cgig2Spur_final_normalised, paste0("Cgig2Spur_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cgig2Blan_final_normalised, paste0("Cgig2Blan_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cgig2Drer_final_normalised, paste0("Cgig2Drer_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cgig2Dmel_final_normalised, paste0("Cgig2Dmel_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cgig2Cele_final_normalised, paste0("Cgig2Cele_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cgig2Aque_final_normalised, paste0("Cgig2Aque_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cgig2Chem_final_normalised, paste0("Cgig2Chem_final_normalised",".txt"), sep='\t', quote = FALSE)

write.table(Cele2Nvec_final_normalised, paste0("Cele2Nvec_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cele2Spur_final_normalised, paste0("Cele2Spur_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cele2Blan_final_normalised, paste0("Cele2Blan_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cele2Drer_final_normalised, paste0("Cele2Drer_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cele2Dmel_final_normalised, paste0("Cele2Dmel_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cele2Aque_final_normalised, paste0("Cele2Aque_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Cele2Chem_final_normalised, paste0("Cele2Chem_final_normalised",".txt"), sep='\t', quote = FALSE)

write.table(Dmel2Nvec_final_normalised, paste0("Dmel2Nvec_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Dmel2Spur_final_normalised, paste0("Dmel2Spur_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Dmel2Blan_final_normalised, paste0("Dmel2Blan_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Dmel2Drer_final_normalised, paste0("Dmel2Drer_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Dmel2Aque_final_normalised, paste0("Dmel2Aque_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Dmel2Chem_final_normalised, paste0("Dmel2Chem_final_normalised",".txt"), sep='\t', quote = FALSE)

write.table(Drer2Nvec_final_normalised, paste0("Drer2Nvec_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Drer2Spur_final_normalised, paste0("Drer2Spur_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Drer2Blan_final_normalised, paste0("Drer2Blan_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Drer2Aque_final_normalised, paste0("Drer2Aque_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Drer2Chem_final_normalised, paste0("Drer2Chem_final_normalised",".txt"), sep='\t', quote = FALSE)

write.table(Blan2Nvec_final_normalised, paste0("Blan2Nvec_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Blan2Spur_final_normalised, paste0("Blan2Spur_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Blan2Aque_final_normalised, paste0("Blan2Aque_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Blan2Chem_final_normalised, paste0("Blan2Chem_final_normalised",".txt"), sep='\t', quote = FALSE)

write.table(Spur2Nvec_final_normalised, paste0("Spur2Nvec_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Spur2Aque_final_normalised, paste0("Spur2Aque_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Spur2Chem_final_normalised, paste0("Spur2Chem_final_normalised",".txt"), sep='\t', quote = FALSE)

write.table(Nvec2Aque_final_normalised, paste0("Nvec2Aque_final_normalised",".txt"), sep='\t', quote = FALSE)
write.table(Nvec2Chem_final_normalised, paste0("Nvec2Chem_final_normalised",".txt"), sep='\t', quote = FALSE)

write.table(Aque2Chem_final_normalised, paste0("Aque2Chem_final_normalised",".txt"), sep='\t', quote = FALSE)




# 12. Plot all 36 heatmaps: PENDING
# Create 9x9 grid
heatmap_color <- colorRampPalette(brewer.pal(n = 7, name = "RdGy"))(100)


paletteLength <- 100
myBreaks <- c(seq(1/paletteLength, 1, length.out=floor(paletteLength)))
heatmap_color <- colorRampPalette(brewer.pal(n = 7, name = "RdGy"))(100)

# Column 1: Owenia fusiformis vs. 10 species
pheatmap(Ofus2Ctel_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Ofus2Cgig_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Ofus2Cele_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Ofus2Dmel_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Ofus2Drer_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Ofus2Blan_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Ofus2Spur_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Ofus2Nvec_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Ofus2Chem_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Ofus2Aque_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)


# Column 2: Capitella teleta vs. 9 species 
pheatmap(Ctel2Cgig_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Ctel2Cele_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Ctel2Dmel_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Ctel2Drer_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Ctel2Blan_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Ctel2Spur_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Ctel2Nvec_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Ctel2Chem_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Ctel2Aque_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)


# Column 3: Crassostrea gigas vs. 8 species
pheatmap(Cgig2Cele_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Cgig2Dmel_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Cgig2Drer_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Cgig2Blan_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Cgig2Spur_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Cgig2Nvec_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Cgig2Chem_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Cgig2Aque_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)


# Column 4: Caenorhabditis elegans vs. 7 species
pheatmap(Cele2Dmel_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Cele2Drer_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Cele2Blan_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Cele2Spur_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Cele2Nvec_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Cele2Chem_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Cele2Aque_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)


# Column 5: Drosophila melanogaster vs. 6 species
pheatmap(Dmel2Drer_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Dmel2Blan_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Dmel2Spur_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Dmel2Nvec_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Dmel2Chem_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Dmel2Aque_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)


# Column 6: Danio rerio vs. 5 species
pheatmap(Drer2Blan_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Drer2Spur_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Drer2Nvec_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Drer2Chem_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Drer2Aque_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)


# Column 7: Branchiostoma lanceolatum vs. 4 species
pheatmap(Blan2Spur_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Blan2Nvec_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Blan2Chem_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Blan2Aque_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)


# Column 8: Strongylocentrotus purpuratus vs. 3 species
pheatmap(Spur2Nvec_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Spur2Chem_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Spur2Aque_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)


# Column 9: Nematostella vectensis vs. 2 species
pheatmap(Nvec2Chem_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

pheatmap(Nvec2Aque_final_normalised, 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)


# Column 10: Clytia hemisphaerica vs. 1 species
pheatmap(t(Aque2Chem_final_normalised), 
         cluster_rows = F, cluster_cols = F, cellheight = 10, cellwidth = 10,
         border_color = NA, color = heatmap_color, breaks = myBreaks)

