library(ggplot2)
library(tidyr)
library(preprocessCore)
library(RColorBrewer)
library(pheatmap)


#############################
#### 1. Capitella teleta ####  
#############################

# Import order of TFclass entries
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/16-Transcription_factors_Owenia_Capitella_Dimorphilus/")
tfclass <- read.table("01-TF_PFAM_domains.txt", sep = '\t', header = T)[4]
colnames(tfclass) <- "tfclass"

# Import expression data with TFclass
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/19-Transcription_factors_RNAseq_clusters_enrichment/")
stages <- c("64_cells","gastrula","st4tt_larva","st5_larva","st7_larva",
            "precompetent_larva","competent_larva")
ctel_raw <- read.table("01b-Capitella_teleta_TF_expression_DESeq2_average.txt", sep ='\t', header=T, row.names = 1)
colnames(ctel_raw)[2] <- stages[1]


# Import cluster annotation from mfuzz
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/01-Capitella_RNAseq_core_analyses/")
ctel_clusters <- read.table("06-Capitella_teleta_clusters_annotation_corrected.txt",
                       sep ='\t', header = T, row.names = 1)
ctel_data <- cbind("gene_ID" = rownames(ctel_raw), ctel_raw[1], "cluster" = ctel_clusters[rownames(ctel_raw),"Cluster_corrected"])


# Some TFs belong to more than one family, so we had to remove them for the enrichment testing
# Other TFs were not assigned a cluster during the clustering, so they were removed too
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/19-Transcription_factors_RNAseq_clusters_enrichment/")
ctel_data_unique <- na.omit(ctel_data)
rownames(ctel_data_unique) <- c(1:nrow(ctel_data_unique))

pvalues <- data.frame(matrix(ncol = length(unique(ctel_data_unique$tfclass)), nrow = length(unique(ctel_data_unique$cluster))))
oddsratio <- data.frame(matrix(ncol = length(unique(ctel_data_unique$tfclass)), nrow = length(unique(ctel_data_unique$cluster))))

colnames(pvalues) <- tfclass$tfclass
colnames(oddsratio) <- tfclass$tfclass
rownames(pvalues) <- c(1:length(unique(ctel_data_unique$cluster)))
rownames(oddsratio) <- c(1:length(unique(ctel_data_unique$cluster)))


for (i in tfclass$tfclass){
  for (j in c(1:length(unique(ctel_data_unique$cluster)))){
    test <- data.frame("in_tfclass" = c(length(which(ctel_data_unique$cluster == j & ctel_data_unique$tfclass == i)),
                                        length(which(ctel_data_unique$cluster != j & ctel_data_unique$tfclass == i))),
                       "not_in_tfclass" = c(length(which(ctel_data_unique$cluster == j & ctel_data_unique$tfclass != i)),
                                            length(which(ctel_data_unique$cluster != j & ctel_data_unique$tfclass != i))),
                       row.names = c("in_cluster", "not_in_cluster"),
                       stringsAsFactors = FALSE)
    oddsratio[j,i] <- fisher.test(test)$estimate[1]
    if (fisher.test(test)$estimate[1] > 1){
      pvalues[j,i] <- fisher.test(test)$p.value
    } else {
      pvalues[j,i] <- -fisher.test(test)$p.value
    }
}
}

write.table(oddsratio, sep="\t", file="Odds_ratio_Capitella.txt")

for (i in tfclass$tfclass){
  pvalues[,i] <- p.adjust(pvalues[,i], method = "BH")
}

pvalues[pvalues>0.05] <- 1
pvalues[pvalues<(-0.05)] <- 1
pvalues_log <- data.matrix(-log10(abs(pvalues)))

for (i in c(1:length(unique(ctel_data_unique$cluster)))){
  for (j in c(1:length(unique(ctel_data_unique$tfclass)))){
    if (pvalues[i,j] < 0) {
      pvalues_log[i,j] <- -pvalues_log[i,j]
    }
  }
}

paletteLength <- 100
myBreaks <- c(seq(min(pvalues_log), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(pvalues_log)/paletteLength, max(pvalues_log), length.out=floor(paletteLength/2)))

myBreaks <- c(seq(-5, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(50/paletteLength, 5, length.out=floor(paletteLength/2)))

heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color[50] <- rgb(1,1,1)

pheatmap(t(pvalues_log),
         cluster_rows = F,
         cluster_cols = F,
         border_color = NA,
         color = heatmap_color,
         breaks = myBreaks,
         display_numbers = T)




#####################################
#### 2. Dimorphilus gyrociliatus ####  
#####################################

# Import order of TFclass entries
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/16-Transcription_factors_Owenia_Capitella_Dimorphilus/")
tfclass <- read.table("01-TF_PFAM_domains.txt", sep = '\t', header = T)[4]
colnames(tfclass) <- "tfclass"

# Import expression data with TFclass
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/19-Transcription_factors_RNAseq_clusters_enrichment/")
stages <- c("early_development","late_development","hatchling","female_adult")
dgyr_raw <- read.table("01c-Dimorphilus_gyrociliatus_TF_expression_DESeq2_average.txt", sep ='\t', header=T, row.names = 1)

# Import cluster annotation from mfuzz
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/07-Dimorphilus_RNAseq_core_analyses/")
dgyr_clusters <- read.table("06-Dimorphilus_gyrociliatus_clusters_annotation_corrected.txt",
                            sep ='\t', header = T, row.names = 1)
dgyr_data <- cbind("gene_ID" = rownames(dgyr_raw), dgyr_raw[1], "cluster" = dgyr_clusters[rownames(dgyr_raw),"Cluster_corrected"])


# Some TFs belong to more than one family, so we had to remove them for the enrichment testing
# Other TFs were not assigned a cluster during the clustering, so they were removed too
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/19-Transcription_factors_RNAseq_clusters_enrichment/")
dgyr_data_unique <- na.omit(dgyr_data)
rownames(dgyr_data_unique) <- c(1:nrow(dgyr_data_unique))

pvalues <- data.frame(matrix(ncol = length(unique(dgyr_data_unique$tfclass)), nrow = length(unique(dgyr_data_unique$cluster))))
oddsratio <- data.frame(matrix(ncol = length(unique(dgyr_data_unique$tfclass)), nrow = length(unique(dgyr_data_unique$cluster))))

colnames(pvalues) <- tfclass$tfclass
colnames(oddsratio) <- tfclass$tfclass
rownames(pvalues) <- c(1:length(unique(dgyr_data_unique$cluster)))
rownames(oddsratio) <- c(1:length(unique(dgyr_data_unique$cluster)))


for (i in tfclass$tfclass){
  for (j in c(1:length(unique(dgyr_data_unique$cluster)))){
    test <- data.frame("in_tfclass" = c(length(which(dgyr_data_unique$cluster == j & dgyr_data_unique$tfclass == i)),
                                        length(which(dgyr_data_unique$cluster != j & dgyr_data_unique$tfclass == i))),
                       "not_in_tfclass" = c(length(which(dgyr_data_unique$cluster == j & dgyr_data_unique$tfclass != i)),
                                            length(which(dgyr_data_unique$cluster != j & dgyr_data_unique$tfclass != i))),
                       row.names = c("in_cluster", "not_in_cluster"),
                       stringsAsFactors = FALSE)
    oddsratio[j,i] <- fisher.test(test)$estimate[1]
    if (fisher.test(test)$estimate[1] > 1){
      pvalues[j,i] <- fisher.test(test)$p.value
    } else {
      pvalues[j,i] <- -fisher.test(test)$p.value
    }
  }
}

write.table(oddsratio, sep="\t", file="Odds_ratio_Dimorphilus.txt")

for (i in tfclass$tfclass){
  pvalues[,i] <- p.adjust(pvalues[,i], method = "BH")
}

pvalues[pvalues>0.05] <- 1
pvalues[pvalues<(-0.05)] <- 1
pvalues_log <- data.matrix(-log10(abs(pvalues)))

for (i in c(1:length(unique(dgyr_data_unique$cluster)))){
  for (j in c(1:length(unique(dgyr_data_unique$tfclass)))){
    if (pvalues[i,j] < 0) {
      pvalues_log[i,j] <- -pvalues_log[i,j]
    }
  }
}

paletteLength <- 100
myBreaks <- c(seq(min(pvalues_log), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(pvalues_log)/paletteLength, max(pvalues_log), length.out=floor(paletteLength/2)))

myBreaks <- c(seq(-5, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(50/paletteLength, 5, length.out=floor(paletteLength/2)))

heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color[50] <- rgb(1,1,1)

pheatmap(t(pvalues_log),
         cluster_rows = F,
         cluster_cols = F,
         border_color = NA,
         color = heatmap_color,
         breaks = myBreaks,
         display_numbers = T)




##############################
#### 3. Owenia fusiformis ####  
##############################

# Import order of TFclass entries
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/16-Transcription_factors_Owenia_Capitella_Dimorphilus/")
tfclass <- read.table("01-TF_PFAM_domains.txt", sep = '\t', header = T)[4]
colnames(tfclass) <- "tfclass"

# Import expression data with TFclass
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/19-Transcription_factors_RNAseq_clusters_enrichment/")
stages <- c("blastula","gastrula","elongation","early_larva",
            "mitraria_larva","competent_larva","juvenile")
ofus_raw <- read.table("01a-Owenia_fusiformis_TF_expression_DESeq2_average.txt", sep ='\t', header=T, row.names = 1)

# Import cluster annotation from mfuzz
cluster1 <- cbind(read.table("~/Dropbox/02-OweniaGenome/01-Data/07-RNAseq/01-New_genome/time_course_genes_expression/mfuzz/new_C1_average.txt", 
                             as.is=TRUE, header=TRUE, row.names = 1),cluster=c(1))
cluster2 <- cbind(read.table("~/Dropbox/02-OweniaGenome/01-Data/07-RNAseq/01-New_genome/time_course_genes_expression/mfuzz/new_C2_average.txt", 
                             as.is=TRUE, header=TRUE, row.names = 1),cluster=c(2))
cluster3 <- cbind(read.table("~/Dropbox/02-OweniaGenome/01-Data/07-RNAseq/01-New_genome/time_course_genes_expression/mfuzz/new_C3_average.txt", 
                             as.is=TRUE, header=TRUE, row.names = 1),cluster=c(3))
cluster4 <- cbind(read.table("~/Dropbox/02-OweniaGenome/01-Data/07-RNAseq/01-New_genome/time_course_genes_expression/mfuzz/new_C4_average.txt", 
                             as.is=TRUE, header=TRUE, row.names = 1),cluster=c(4))
cluster5 <- cbind(read.table("~/Dropbox/02-OweniaGenome/01-Data/07-RNAseq/01-New_genome/time_course_genes_expression/mfuzz/new_C5_average.txt", 
                             as.is=TRUE, header=TRUE, row.names = 1),cluster=c(5))
cluster6 <- cbind(read.table("~/Dropbox/02-OweniaGenome/01-Data/07-RNAseq/01-New_genome/time_course_genes_expression/mfuzz/new_C6_average.txt", 
                             as.is=TRUE, header=TRUE, row.names = 1),cluster=c(6))
cluster7 <- cbind(read.table("~/Dropbox/02-OweniaGenome/01-Data/07-RNAseq/01-New_genome/time_course_genes_expression/mfuzz/new_C7_average.txt", 
                             as.is=TRUE, header=TRUE, row.names = 1),cluster=c(7))
cluster8 <- cbind(read.table("~/Dropbox/02-OweniaGenome/01-Data/07-RNAseq/01-New_genome/time_course_genes_expression/mfuzz/new_C8_average.txt", 
                             as.is=TRUE, header=TRUE, row.names = 1),cluster=c(8))
cluster9 <- cbind(read.table("~/Dropbox/02-OweniaGenome/01-Data/07-RNAseq/01-New_genome/time_course_genes_expression/mfuzz/new_C9_average.txt", 
                             as.is=TRUE, header=TRUE, row.names = 1),cluster=c(9))
cluster10 <- cbind(read.table("~/Dropbox/02-OweniaGenome/01-Data/07-RNAseq/01-New_genome/time_course_genes_expression/mfuzz/new_C10_average.txt", 
                              as.is=TRUE, header=TRUE, row.names = 1),cluster=c(10))
cluster11 <- cbind(read.table("~/Dropbox/02-OweniaGenome/01-Data/07-RNAseq/01-New_genome/time_course_genes_expression/mfuzz/new_C11_average.txt", 
                              as.is=TRUE, header=TRUE, row.names = 1),cluster=c(11))
cluster12 <- cbind(read.table("~/Dropbox/02-OweniaGenome/01-Data/07-RNAseq/01-New_genome/time_course_genes_expression/mfuzz/new_C12_average.txt", 
                              as.is=TRUE, header=TRUE, row.names = 1),cluster=c(12))
ofus_clusters <- rbind(cluster1,cluster2,cluster3,cluster4,cluster5,cluster6,
                       cluster7,cluster8,cluster9,cluster10,cluster11,cluster12)

ofus_data <- cbind("gene_ID" = rownames(ofus_raw), ofus_raw[1], "cluster" = ofus_clusters[rownames(ofus_raw),"cluster"])


# Some TFs belong to more than one family, so we had to remove them for the enrichment testing
# Other TFs were not assigned a cluster during the clustering, so they were removed too
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/19-Transcription_factors_RNAseq_clusters_enrichment/")
ofus_data_unique <- na.omit(ofus_data)
rownames(ofus_data_unique) <- c(1:nrow(ofus_data_unique))

pvalues <- data.frame(matrix(ncol = length(unique(ofus_data_unique$tfclass)), nrow = length(unique(ofus_data_unique$cluster))))
oddsratio <- data.frame(matrix(ncol = length(unique(ofus_data_unique$tfclass)), nrow = length(unique(ofus_data_unique$cluster))))

colnames(pvalues) <- tfclass$tfclass
colnames(oddsratio) <- tfclass$tfclass
rownames(pvalues) <- c(1:length(unique(ofus_data_unique$cluster)))
rownames(oddsratio) <- c(1:length(unique(ofus_data_unique$cluster)))


for (i in tfclass$tfclass){
  for (j in c(1:length(unique(ofus_data_unique$cluster)))){
    test <- data.frame("in_tfclass" = c(length(which(ofus_data_unique$cluster == j & ofus_data_unique$tfclass == i)),
                                        length(which(ofus_data_unique$cluster != j & ofus_data_unique$tfclass == i))),
                       "not_in_tfclass" = c(length(which(ofus_data_unique$cluster == j & ofus_data_unique$tfclass != i)),
                                            length(which(ofus_data_unique$cluster != j & ofus_data_unique$tfclass != i))),
                       row.names = c("in_cluster", "not_in_cluster"),
                       stringsAsFactors = FALSE)
    oddsratio[j,i] <- fisher.test(test)$estimate[1]
    if (fisher.test(test)$estimate[1] > 1){
      pvalues[j,i] <- fisher.test(test)$p.value
    } else {
      pvalues[j,i] <- -fisher.test(test)$p.value
    }
  }
}

write.table(oddsratio, sep="\t", file="Odds_ratio_Owenia.txt")

for (i in tfclass$tfclass){
  pvalues[,i] <- p.adjust(pvalues[,i], method = "benjamini-hochberg")
}

pvalues[pvalues>0.05] <- 1
pvalues[pvalues<(-0.05)] <- 1
pvalues_log <- data.matrix(-log10(abs(pvalues)))

for (i in c(1:length(unique(ofus_data_unique$cluster)))){
  for (j in c(1:length(unique(ofus_data_unique$tfclass)))){
    if (pvalues[i,j] < 0) {
      pvalues_log[i,j] <- -pvalues_log[i,j]
    }
  }
}

paletteLength <- 100
myBreaks <- c(seq(min(pvalues_log), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(pvalues_log)/paletteLength, max(pvalues_log), length.out=floor(paletteLength/2)))

myBreaks <- c(seq(-5, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(50/paletteLength, 5, length.out=floor(paletteLength/2)))

heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color[50] <- rgb(1,1,1)

pheatmap(t(pvalues_log),
         cluster_rows = F,
         cluster_cols = F,
         border_color = NA,
         color = heatmap_color,
         breaks = myBreaks,
         display_numbers = T)




############################################# 
#### 4. Owenia fusiformis new clustering ####  
############################################# 

setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/19-Transcription_factors_RNAseq_clusters_enrichment/")
ofus_raw <- read.table("01a-Owenia_fusiformis_TF_expression_DESeq2_average.txt", sep ='\t', header=T, row.names = 1)

ofus_clusters <- read.table("/Users/franciscomanuelmartinzamora/Dropbox/02-OweniaGenome/06-Revision/00-DATA/27-Owenia_RNAseq_core_analyses_new_clustering/06-Owenia_fusiformis_clusters_annotation_corrected.txt",
                            as.is = T, header = T, row.names = 1)

ofus_data <- cbind("gene_ID" = rownames(ofus_raw), ofus_raw[1], "cluster" = ofus_clusters[rownames(ofus_raw),"Cluster_corrected"])


# Some TFs belong to more than one family, so we had to remove them for the enrichment testing
# Other TFs were not assigned a cluster during the clustering, so they were removed too
setwd("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/19-Transcription_factors_RNAseq_clusters_enrichment/")
ofus_data_unique <- na.omit(ofus_data)
rownames(ofus_data_unique) <- c(1:nrow(ofus_data_unique))

pvalues <- data.frame(matrix(ncol = length(unique(ofus_data_unique$tfclass)), nrow = length(unique(ofus_data_unique$cluster))))
oddsratio <- data.frame(matrix(ncol = length(unique(ofus_data_unique$tfclass)), nrow = length(unique(ofus_data_unique$cluster))))

colnames(pvalues) <- tfclass$tfclass
colnames(oddsratio) <- tfclass$tfclass
rownames(pvalues) <- c(1:length(unique(ofus_data_unique$cluster)))
rownames(oddsratio) <- c(1:length(unique(ofus_data_unique$cluster)))


for (i in tfclass$tfclass){
  for (j in c(1:length(unique(ofus_data_unique$cluster)))){
    test <- data.frame("in_tfclass" = c(length(which(ofus_data_unique$cluster == j & ofus_data_unique$tfclass == i)),
                                        length(which(ofus_data_unique$cluster != j & ofus_data_unique$tfclass == i))),
                       "not_in_tfclass" = c(length(which(ofus_data_unique$cluster == j & ofus_data_unique$tfclass != i)),
                                            length(which(ofus_data_unique$cluster != j & ofus_data_unique$tfclass != i))),
                       row.names = c("in_cluster", "not_in_cluster"),
                       stringsAsFactors = FALSE)
    oddsratio[j,i] <- fisher.test(test)$estimate[1]
    if (fisher.test(test)$estimate[1] > 1){
      pvalues[j,i] <- fisher.test(test)$p.value
    } else {
      pvalues[j,i] <- -fisher.test(test)$p.value
    }
  }
}

write.table(oddsratio, sep="\t", file="Odds_ratio_Owenia_new_clustering.txt")

for (i in tfclass$tfclass){
  pvalues[,i] <- p.adjust(pvalues[,i], method = "BH")
}

pvalues[pvalues>0.05] <- 1
pvalues[pvalues<(-0.05)] <- 1
pvalues_log <- data.matrix(-log10(abs(pvalues)))

for (i in c(1:length(unique(ofus_data_unique$cluster)))){
  for (j in c(1:length(unique(ofus_data_unique$tfclass)))){
    if (pvalues[i,j] < 0) {
      pvalues_log[i,j] <- -pvalues_log[i,j]
    }
  }
}

paletteLength <- 100
myBreaks <- c(seq(min(pvalues_log), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(pvalues_log)/paletteLength, max(pvalues_log), length.out=floor(paletteLength/2)))

myBreaks <- c(seq(-5, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(50/paletteLength, 5, length.out=floor(paletteLength/2)))

heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color[50] <- rgb(1,1,1)

pheatmap(t(pvalues_log),
         cluster_rows = F,
         cluster_cols = F,
         border_color = NA,
         color = heatmap_color,
         breaks = myBreaks,
         display_numbers = T)
