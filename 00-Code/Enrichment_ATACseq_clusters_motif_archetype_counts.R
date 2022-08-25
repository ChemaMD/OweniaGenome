library(ggplot2)
library(tidyr)
library(preprocessCore)
library(RColorBrewer)
library(pheatmap)

#################################
## 0. Prepare data for testing ##
#################################

mock_ctel <- read.table("02-Ctel_peakID_motifCount.txt", sep = '\t', header = T)
mock_ofus <- read.table("02-Ofus_peakID_motifCount.txt", sep = '\t', header = T)

motif_list_ctel <- data.frame("motif_name" = sort(colnames(mock_ctel)))
motif_list_ofus <- data.frame("motif_name" = sort(colnames(mock_ofus)))

motif_list_ofus_reduced <- data.frame("motif_name" = t(motif_list_ofus)[-c(27:45,51,53:65,67:96)])
motif_list_ctel_reduced <- data.frame("motif_name" = t(motif_list_ctel)[-c(27:45,52:62,64:92)])

#############################
#### 1. Capitella teleta ####  
#############################

ctel_clusters <- read.table("01-Capitella_teleta_peak_clusters.txt", sep ='\t', header = T)
rownames(ctel_clusters) <- ctel_clusters$Peak_ID

ctel_motif_count <- read.table("02-Ctel_peakID_motifCount.txt", sep ='\t', header = T)
ctel_motif_count_labelled <- cbind(ctel_motif_count, ctel_clusters[rownames(ctel_motif_count),"Cluster"])
colnames(ctel_motif_count_labelled)[93] <- "Cluster"
ctel_motif_count_final <- aggregate(ctel_motif_count_labelled[-c(93)], list(ctel_motif_count_labelled$Cluster), FUN=sum) 
colnames(ctel_motif_count_final)[1] <- "Cluster"

ctel_rowsummed <- data.frame("cluster_motif_total" = rowSums(ctel_motif_count_final[-c(1)], na.rm = T))

ctel_pvalues <- data.frame(matrix(ncol = length(motif_list_ctel$motif_name), nrow = length(ctel_motif_count_final$Cluster)))
ctel_oddsratio <- data.frame(matrix(ncol = length(motif_list_ctel$motif_name), nrow = length(ctel_motif_count_final$Cluster)))

colnames(ctel_pvalues) <- motif_list_ctel$motif_name
colnames(ctel_oddsratio) <- motif_list_ctel$motif_name
rownames(ctel_pvalues) <- ctel_motif_count_final$Cluster
rownames(ctel_pvalues) <- ctel_motif_count_final$Cluster

for (i in motif_list_ctel$motif_name){
  for (j in ctel_motif_count_final$Cluster){
    test <- data.frame("motif_occurences" = c(ctel_motif_count_final[j,i],
                                              ctel_rowsummed[j,] - ctel_motif_count_final[j,i]),
                       "other_motif_occurences" = c(sum(ctel_motif_count_final[,i]) - ctel_motif_count_final[j,i],
                                                    sum(ctel_rowsummed) - (sum(ctel_motif_count_final[,i]) - ctel_motif_count_final[j,i]) - (ctel_rowsummed[j,] - ctel_motif_count_final[j,i])),
                       row.names = c("in_cluster", "not_in_cluster"),
                       stringsAsFactors = FALSE)
    ctel_oddsratio[j,i] <- fisher.test(test)$estimate[1]
    if (fisher.test(test)$estimate[1] > 1){
      ctel_pvalues[j,i] <- fisher.test(test)$p.value
    } else {
      ctel_pvalues[j,i] <- -fisher.test(test)$p.value
    }
  }
}

write.table(ctel_oddsratio, sep="\t", file="03-Odds_ratio_Capitella.txt", quote = F)

# Subet motifs of interest
ctel_pvalues <- ctel_pvalues[colnames(ctel_pvalues) %in% motif_list_ctel_reduced$motif_name]


for (i in motif_list_ctel_reduced$motif_name){
  ctel_pvalues[,i] <- p.adjust(ctel_pvalues[,i], method = "bonferroni")
}

ctel_pvalues[ctel_pvalues>0.05] <- 1
ctel_pvalues[ctel_pvalues<(-0.05)] <- 1
ctel_pvalues_log <- data.matrix(-log10(abs(ctel_pvalues)))

for (i in c(1:length(unique(ctel_motif_count_final$Cluster)))){
  for (j in c(1:length(unique(motif_list_ctel_reduced$motif_name)))){
    if (ctel_pvalues[i,j] < 0) {
      ctel_pvalues_log[i,j] <- -ctel_pvalues_log[i,j]
    }
  }
}

paletteLength <- 100
myBreaks <- c(seq(min(ctel_pvalues_log), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(ctel_pvalues_log)/paletteLength, max(ctel_pvalues_log), length.out=floor(paletteLength/2)))

myBreaks <- c(seq(-10, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(50/paletteLength, 10, length.out=floor(paletteLength/2)))

heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color[50] <- rgb(1,1,1)

pheatmap(t(ctel_pvalues_log),
         cluster_rows = F,
         cluster_cols = F,
         border_color = NA,
         color = heatmap_color,
         breaks = myBreaks,
         display_numbers = T)


##############################
#### 2. Owenia fusiformis ####
##############################


ofus_clusters <- read.table("01-Owenia_fusiformis_peak_clusters.txt", sep ='\t', header = T)
rownames(ofus_clusters) <- ofus_clusters$Peak_ID

ofus_motif_count <- read.table("02-Ofus_peakID_motifCount.txt", sep ='\t', header = T)
ofus_motif_count_labelled <- cbind(ofus_motif_count, ofus_clusters[rownames(ofus_motif_count),"Cluster"])
colnames(ofus_motif_count_labelled)[97] <- "Cluster"
ofus_motif_count_final <- aggregate(ofus_motif_count_labelled[-c(97)], list(ofus_motif_count_labelled$Cluster), FUN=sum) 
colnames(ofus_motif_count_final)[1] <- "Cluster"

ofus_rowsummed <- data.frame("cluster_motif_total" = rowSums(ofus_motif_count_final[-c(1)], na.rm = T))

ofus_pvalues <- data.frame(matrix(ncol = length(motif_list_ofus$motif_name), nrow = length(ofus_motif_count_final$Cluster)))
ofus_oddsratio <- data.frame(matrix(ncol = length(motif_list_ofus$motif_name), nrow = length(ofus_motif_count_final$Cluster)))

colnames(ofus_pvalues) <- motif_list_ofus$motif_name
colnames(ofus_oddsratio) <- motif_list_ofus$motif_name
rownames(ofus_pvalues) <- ofus_motif_count_final$Cluster
rownames(ofus_pvalues) <- ofus_motif_count_final$Cluster

for (i in motif_list_ofus$motif_name){
  for (j in ofus_motif_count_final$Cluster){
    test <- data.frame("motif_occurences" = c(ofus_motif_count_final[j,i],
                                              ofus_rowsummed[j,] - ofus_motif_count_final[j,i]),
                       "other_motif_occurences" = c(sum(ofus_motif_count_final[,i]) - ofus_motif_count_final[j,i],
                                                    sum(ofus_rowsummed) - (sum(ofus_motif_count_final[,i]) - ofus_motif_count_final[j,i]) - (ofus_rowsummed[j,] - ofus_motif_count_final[j,i])),
                       row.names = c("in_cluster", "not_in_cluster"),
                       stringsAsFactors = FALSE)
    ofus_oddsratio[j,i] <- fisher.test(test)$estimate[1]
    if (fisher.test(test)$estimate[1] > 1){
      ofus_pvalues[j,i] <- fisher.test(test)$p.value
    } else {
      ofus_pvalues[j,i] <- -fisher.test(test)$p.value
    }
  }
}

write.table(ofus_oddsratio, sep="\t", file="03-Odds_ratio_Owenia.txt", quote = F)


# Subet motifs of interest
ofus_pvalues <- ofus_pvalues[colnames(ofus_pvalues) %in% motif_list_ofus_reduced$motif_name]


for (i in motif_list_ofus_reduced$motif_name){
  ofus_pvalues[,i] <- p.adjust(ofus_pvalues[,i], method = "bonferroni")
}

ofus_pvalues[ofus_pvalues>0.05] <- 1
ofus_pvalues[ofus_pvalues<(-0.05)] <- 1
ofus_pvalues_log <- data.matrix(-log10(abs(ofus_pvalues)))

for (i in c(1:length(unique(ofus_motif_count_final$Cluster)))){
  for (j in c(1:length(unique(motif_list_ofus_reduced$motif_name)))){
    if (ofus_pvalues[i,j] < 0) {
      ofus_pvalues_log[i,j] <- -ofus_pvalues_log[i,j]
    }
  }
}

paletteLength <- 100
myBreaks <- c(seq(min(ofus_pvalues_log), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(ofus_pvalues_log)/paletteLength, max(ofus_pvalues_log), length.out=floor(paletteLength/2)))

myBreaks <- c(seq(-10, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(50/paletteLength, 10, length.out=floor(paletteLength/2)))

heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color[50] <- rgb(1,1,1)

pheatmap(t(ofus_pvalues_log),
         cluster_rows = F,
         cluster_cols = F,
         border_color = NA,
         color = heatmap_color,
         breaks = myBreaks,
         display_numbers = T)
