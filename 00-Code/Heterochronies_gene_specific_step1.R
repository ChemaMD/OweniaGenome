library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(topGO)
library(cowplot)
library(VennDiagram)
library(readxl)
library(writexl)

install.packages("writexl")

# Import clustered genomes
ofus_mfuzz <- read.table("/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/01-Owenia_fusiformis_mfuzz_clusters.txt", header=F, sep ='\t')
ctel_mfuzz <- read.table("/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/01-Capitella_teleta_mfuzz_clusters.txt", header=F, sep ='\t')
dgyr_mfuzz <- read.table("/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/01-Dimorphilus_gyrociliatus_mfuzz_clusters.txt", header=F, sep ='\t')

# Import 1-to-1 orthologues between Owenia and Capitella/Dimorphilus
ofus2ctel <- read.table("/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/02-Ofus2Ctel_orthologues.txt", header=T, sep ='\t')
ofus2dgyr <- read.table("/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/02-Ofus2Dgyr_orthologues.txt", header=T, sep ='\t')
ctel2dgyr <- read.table("/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/02-Ctel2Dgyr_orthologues.txt", header=T, sep ='\t')


# Subset genes from clusters 7-12 (Capitella) and 1-4 (Dimorphilus)
ctel_7to12 <- ctel_mfuzz[ctel_mfuzz$V2 == '7' | ctel_mfuzz$V2 == '8' | ctel_mfuzz$V2 == '9' | ctel_mfuzz$V2 == '10'| ctel_mfuzz$V2 == '11' | ctel_mfuzz$V2 == '12',] # 14,251 transcripts in clusters 1, 2, 4 and 5 in Capitella 
dgyr_23 <- dgyr_mfuzz[ dgyr_mfuzz$V2 == '2' | dgyr_mfuzz$V2 == '3',] # 4645 transcripts in clusters 2 and 3 in Dimorphilus

# Subset genes from clusters 8 and 12 (Owenia)
ofus_8to12 <- ofus_mfuzz[ofus_mfuzz$V2 == '8'| ofus_mfuzz$V2 == '9'| ofus_mfuzz$V2 == '10'| ofus_mfuzz$V2 == '11'|ofus_mfuzz$V2 == '12',] # 2,799 transcripts in cluster 7 in Owenia
#ofus_12 <- ofus_mfuzz[ofus_mfuzz$V2 == '12',] # 3,020 transcripts in cluster 12 in Owenia 

# Get the Owenia orthologs (where they exist) from the above subset genes 
ofus_ctel_7to12_orthologues <- data.frame("Gene_ID" = ofus2ctel[ofus2ctel$Ctel %in% ctel_7to12$V1,"Ofus"]) # 2999 1-to-1 orthologs out of the 20511 transcripts (14.62%)
ofus_dgyr_23_orthologues <- data.frame("Gene_ID" = ofus2dgyr[ofus2dgyr$Dgyr %in% dgyr_23$V1,"Ofus"]) # 2320 1-to-1 orthologs out of the 6,492 transcripts (35.73%)
ctel_dgyr_23_orthologues <- data.frame("Gene_ID" = ctel2dgyr[ctel2dgyr$Dgyr %in% dgyr_23$V1,"Ctel"]) # 2239 1-to-1 orthologs out of the 4645 transcripts (48.20%)

# Import transcription factors
# Owenia transcription Factors
ofus_tf_list <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/16-Transcription_factors_Owenia_Capitella_Dimorphilus/06-Owenia_fusiformis_TF_list_classified.txt", 
                           header = F, sep ='\t')
colnames(ofus_tf_list) <- c("Gene_ID", "PFAM_ID", "PFAM_annotation")
ofus_tf_stats <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/16-Transcription_factors_Owenia_Capitella_Dimorphilus/07-Owenia_fusiformis_TF_stats.txt", 
                            header = F, sep ='\t')
colnames(ofus_tf_stats) <- c("PFAM_ID", "PFAM_annotation","count","Gene_ID_list")

# Capitella Transcription Factors
ctel_tf_list <- read.table("~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/16-Transcription_factors_Owenia_Capitella_Dimorphilus/06-Capitella_teleta_TF_list_classified.txt", 
                           header = F, sep ='\t')
colnames(ctel_tf_list) <- c("Gene_ID", "PFAM_ID", "PFAM_annotation")

##################################################################################
###### Get the displaced genes for Owenia fusiformis and Capitella teleta ########
##################################################################################

# gene set 1. genes expressed in 7-12 in Capitella, expressed in 8-12 in Owenia
ofus_ctel_displaced_all <- data.frame("Gene_ID" = ofus_ctel_7to12_orthologues[ofus_ctel_7to12_orthologues$Gene_ID %in% ofus_8to12$V1,"Gene_ID"]) 
# 1602 displaced transcripts from Capitella in cluster 7to12 (1602/15941 = 10.05% of Owenia's cluster 8-12; 1602/20511 = 7.81% of Capitella's clusters 7-12 )

# gene set 2. genes expressed in 7-12 in Capitella, expressed in 2&3 D
ctel_dgyr_displaced_all <- data.frame("Gene_ID" = ctel_dgyr_23_orthologues[ctel_dgyr_23_orthologues$Gene_ID %in% ctel_7to12$V1,"Gene_ID"])
# 552 displaced transcripts from Capitella in cluster 7to12 (552/4645 = 11.88% of Dimorph's cluster 2&3; 552/20511 = 2.69% of Capitella's clusters 7-12 )

# gene set 3. genes expressed in 8-12 in Owenia, expressed in  2&3 D
ofus_dgyr_displaced_all <- data.frame("Gene_ID" = ofus_dgyr_23_orthologues[ofus_dgyr_23_orthologues$Gene_ID %in% ofus_8to12$V1,"Gene_ID"])
# 469 displaced transcripts from Capitella in cluster 7to12 (496/15941 = 3.11% of Owenia's cluster 8-12; 496/4645 = 10.68% of Dimorph's cluster 2&3 )

# gene set 4. genes expressed in Dgyr, Ctel and Ofus
ofus_ctel_dgyr_all_asOw <- data.frame("Gene_ID" = ofus2ctel[ofus2ctel$Ctel %in% ctel_dgyr_displaced_all$Gene_ID, "Ofus"])
ofus_ctel_dgyr_displaced_all <-data.frame(ofus_ctel_dgyr_all_asOw[ofus_ctel_dgyr_all_asOw$Gene_ID %in% ofus_dgyr_displaced_all$Gene_ID,])
# 163 common genes between Dgyr, Ofus and Ctel at given stages

ofus_dgyr_ctel_displaced_all <- data.frame

# extract all the TF from these data sets

tf_ofus_dgyr <- data.frame(ofus_dgyr_displaced_all[ofus_dgyr_displaced_all$Gene_ID %in% ofus_tf_list$Gene_ID,])
# 53/469 common orthologues for Dimorph and Owenia are TF

tf_ctel_dgyr <- data.frame(ctel_dgyr_displaced_all[ctel_dgyr_displaced_all$Gene_ID %in% ctel_tf_list$Gene_ID,])
# 70/552 common orthologues for Dimorph and Capitella are TF

tf_ctel_dgyr_asOw <- data.frame("Gene_ID" = ofus2ctel[ofus2ctel$Ctel %in% tf_ctel_dgyr$ctel_dgyr_displaced_all.ctel_dgyr_displaced_all.Gene_ID..in.., "Ofus"])
tf_dgyr_ctel_ofus <- data.frame(tf_ctel_dgyr_asOw[tf_ctel_dgyr_asOw$Gene_ID %in% tf_ofus_dgyr$ofus_dgyr_displaced_all.ofus_dgyr_displaced_all.Gene_ID..in..,])
# 28 Common TF between all 3


Ow_annotation <- read_excel('/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/Owenia_annotation_v060921.1_TrinoPanther.xlsx')
Ow_annotation <- data.frame(Ow_annotation)
CT_annotation <- read_excel('/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/Capitella_annotation_v300321_TrinoPanther_with_histones.xlsx')



################################# Extraction and Writing excel files##############################


# Commmon TF between Ofus and Dgyr
TF_DO_index <- c()

for (i in tf_ofus_dgyr$ofus_dgyr_displaced_all.ofus_dgyr_displaced_all.Gene_ID..in..) {
  output <- grep(i, Ow_annotation$transcript_id)
  TF_DO_index[i] <- output
}

  
TF_DO_subset <- Ow_annotation[TF_DO_index, ] 
TF_DO_subset <- data.frame(TF_DO_subset)

write_xlsx(TF_DO_subset, '/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/TF_Dgyr_Ofus.xlsx')


# Commmon TF between Ctel and Dgyr
TF_DC_index <- c()

for (i in tf_ctel_dgyr$ctel_dgyr_displaced_all.ctel_dgyr_displaced_all.Gene_ID..in..) {
  output <- grep(i, CT_annotation$Transcript_id)
  TF_DC_index[i] <- output
}

TF_DC_subset <- CT_annotation[TF_DC_index, ] 
TF_DC_subset <- data.frame(TF_DC_subset)

write_xlsx(TF_DC_subset, '/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/TF_Dgyr_Ctel.xlsx')

# Commmon TF between Ofus and Dgyr and Ctel
TF_DOC_index <- c()

for (i in tf_dgyr_ctel_ofus$tf_ctel_dgyr_asOw.tf_ctel_dgyr_asOw.Gene_ID..in..tf_ofus_dgyr.ofus_dgyr_displaced_all.ofus_dgyr_displaced_all.Gene_ID..in....) {
  output <- grep(i, Ow_annotation$transcript_id)
  TF_DOC_index[i] <- output
}

TF_DOC_subset <- Ow_annotation[TF_DOC_index, ] 
TF_DOC_subset <- data.frame(TF_DOC_subset)

write_xlsx(TF_DOC_subset, '/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/TF_Dgyr_Ofus_Ctel.xlsx')


# Commmon Genes between Ofus and Dgyr
G_DO_index <- c()

for (i in ofus_dgyr_displaced_all$Gene_ID){
  output <- grep(i, Ow_annotation$transcript_id)
  G_DO_index[i] <- output
}

G_DO_subset <- Ow_annotation[G_DO_index, ] 
G_DO_subset <- data.frame(G_DO_subset)

write_xlsx(G_DO_subset, '/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/G_Dgyr_Ofus.xlsx')


# Commmon Genes between Ctel and Dgyr
G_DC_index <- c()

for (i in ctel_dgyr_displaced_all$Gene_ID){
  output <- grep(i, CT_annotation$Transcript_id)
  G_DC_index[i] <- output
}

G_DC_subset <- Ow_annotation[G_DC_index, ] 
G_DC_subset <- data.frame(G_DC_subset)

write_xlsx(G_DC_subset, '/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/G_Dgyr_Ctel.xlsx')

# Common Genes between all 3
G_DCO_index <- c()

for (i in ofus_ctel_dgyr_displaced_all$ofus_ctel_dgyr_all_asOw.ofus_ctel_dgyr_all_asOw.Gene_ID..in..){
  output <- grep(i, Ow_annotation$transcript_id)
  G_DCO_index[i] <- output
}

G_DCO_subset <- Ow_annotation[G_DCO_index, ] 
G_DCO_subset <- data.frame(G_DCO_subset)

write_xlsx(G_DCO_subset, '/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/G_Ofus_Dgyr_Ctel.xlsx')


##################################### Subsetting Capitella and Owenia ############################


ctel_7to12 # Capitella Late
ofus_8to12 # Owenia Late
ctel_1to6 # Capitella Early
ofus_1to7 # Owenia Early


# Subset genes from clusters 1-6 (Capitella)
ctel_1to6 <- ctel_mfuzz[ctel_mfuzz$V2 == '1' | ctel_mfuzz$V2 == '2' | ctel_mfuzz$V2 == '3' | ctel_mfuzz$V2 == '4'| ctel_mfuzz$V2 == '5' | ctel_mfuzz$V2 == '6',] 

# Subset genes from clusters 1 and 7 (Owenia)
ofus_1to7 <- ofus_mfuzz[ofus_mfuzz$V2 == '1'| ofus_mfuzz$V2 == '2'| ofus_mfuzz$V2 == '3'| ofus_mfuzz$V2 == '4'|ofus_mfuzz$V2 == '5'|ofus_mfuzz$V2 == '6'|ofus_mfuzz$V2 == '7',]

# Get Orthologues from Early and Late Genes for Capitella
ofus_ctel_E_orthologues <- data.frame("Gene_ID" = ofus2ctel[ofus2ctel$Ctel %in% ctel_1to6$V1,"Ofus"])
ofus_ctel_L_orthologues <- data.frame("Gene_ID" = ofus2ctel[ofus2ctel$Ctel %in% ctel_7to12$V1,"Ofus"])


# Gene Comparisons for EE LE EL LL
ofus_ctel_EE <- data.frame("Gene_ID" = ofus_ctel_E_orthologues[ofus_ctel_E_orthologues$Gene_ID %in% ofus_1to7$V1,"Gene_ID"]) 
ofus_ctel_EL <- data.frame("Gene_ID" = ofus_ctel_L_orthologues[ofus_ctel_L_orthologues$Gene_ID %in% ofus_1to7$V1,"Gene_ID"]) 
ofus_ctel_LE <- data.frame("Gene_ID" = ofus_ctel_E_orthologues[ofus_ctel_E_orthologues$Gene_ID %in% ofus_8to12$V1,"Gene_ID"]) 
ofus_ctel_LL <- data.frame("Gene_ID" = ofus_ctel_L_orthologues[ofus_ctel_L_orthologues$Gene_ID %in% ofus_8to12$V1,"Gene_ID"]) 

# Extract Transcription Factors
tf_ofus_ctel_EE <- data.frame(ofus_ctel_EE[ofus_ctel_EE$Gene_ID %in% ofus_tf_list$Gene_ID,])
tf_ofus_ctel_LE <- data.frame(ofus_ctel_LE[ofus_ctel_LE$Gene_ID %in% ofus_tf_list$Gene_ID,])
tf_ofus_ctel_EL <- data.frame(ofus_ctel_EL[ofus_ctel_EL$Gene_ID %in% ofus_tf_list$Gene_ID,])
tf_ofus_ctel_LL <- data.frame(ofus_ctel_LL[ofus_ctel_LL$Gene_ID %in% ofus_tf_list$Gene_ID,])





###################################  Extract List into excel  #######################################################################

#########################################  All Genes  ###############################################################################

# EE
G_OC_EE_index <- c()

for (i in ofus_ctel_EE$Gene_ID) {
  output <- grep(i, Ow_annotation$transcript_id)
  G_OC_EE_index[i] <- output
}

G_OC_EE_subset <- Ow_annotation[G_OC_EE_index, ] 
G_OC_EE_subset <- data.frame(G_OC_EE_subset)

write_xlsx(G_OC_EE_subset, '/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/G_Ofus_Early_Ctel_Early.xlsx')

# EL
G_OC_EL_index <- c()

for (i in ofus_ctel_EL$Gene_ID) {
  output <- grep(i, Ow_annotation$transcript_id)
  G_OC_EL_index[i] <- output
}

G_OC_EL_subset <- Ow_annotation[G_OC_EL_index, ] 
G_OC_EL_subset <- data.frame(G_OC_EL_subset)

write_xlsx(G_OC_EL_subset, '/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/G_Ofus_Early_Ctel_Late.xlsx')


# LE
G_OC_LE_index <- c()

for (i in ofus_ctel_LE$Gene_ID) {
  output <- grep(i, Ow_annotation$transcript_id)
  G_OC_LE_index[i] <- output
}

G_OC_LE_subset <- Ow_annotation[G_OC_LE_index, ] 
G_OC_LE_subset <- data.frame(G_OC_LE_subset)

write_xlsx(G_OC_LE_subset, '/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/G_Ofus_Late_Ctel_Early.xlsx')


# LL
G_OC_LL_index <- c()

for (i in ofus_ctel_LL$Gene_ID) {
  output <- grep(i, Ow_annotation$transcript_id)
  G_OC_LL_index[i] <- output
}

G_OC_LL_subset <- Ow_annotation[G_OC_LL_index, ] 
G_OC_LL_subset <- data.frame(G_OC_LL_subset)

write_xlsx(G_OC_LL_subset, '/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/G_Ofus_Late_Ctel_Late.xlsx')

################################################# TF ########################################


# EE
TF_OC_EE_index <- c()

for (i in tf_ofus_ctel_EE$ofus_ctel_EE.ofus_ctel_EE.Gene_ID..in..ofus_tf_list.Gene_ID..) {
  output <- grep(i, Ow_annotation$transcript_id)
  TF_OC_EE_index[i] <- output
}

TF_OC_EE_subset <- Ow_annotation[TF_OC_EE_index, ] 
TF_OC_EE_subset <- data.frame(TF_OC_EE_subset)

write_xlsx(TF_OC_EE_subset, '/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/TF_Ofus_Early_Ctel_Early.xlsx')

# EL

TF_OC_EL_index <- c()

for (i in tf_ofus_ctel_EL$ofus_ctel_EL.ofus_ctel_EL.Gene_ID..in..ofus_tf_list.Gene_ID..) {
  output <- grep(i, Ow_annotation$transcript_id)
  TF_OC_EL_index[i] <- output
}

TF_OC_EL_subset <- Ow_annotation[TF_OC_EL_index, ] 
TF_OC_EL_subset <- data.frame(TF_OC_EL_subset)

write_xlsx(TF_OC_EL_subset, '/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/TF_Ofus_Early_Ctel_Late.xlsx')

# LE

TF_OC_LE_index <- c()

for (i in tf_ofus_ctel_LE$ofus_ctel_LE.ofus_ctel_LE.Gene_ID..in..ofus_tf_list.Gene_ID..) {
  output <- grep(i, Ow_annotation$transcript_id)
  TF_OC_LE_index[i] <- output
}

TF_OC_LE_subset <- Ow_annotation[TF_OC_LE_index, ] 
TF_OC_LE_subset <- data.frame(TF_OC_LE_subset)

write_xlsx(TF_OC_LE_subset, '/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/TF_Ofus_Late_Ctel_Early.xlsx')

# LL

TF_OC_LL_index <- c()

for (i in tf_ofus_ctel_LL$ofus_ctel_LL.ofus_ctel_LL.Gene_ID..in..ofus_tf_list.Gene_ID..) {
  output <- grep(i, Ow_annotation$transcript_id)
  TF_OC_LL_index[i] <- output
}

TF_OC_LL_subset <- Ow_annotation[TF_OC_LL_index, ] 
TF_OC_LL_subset <- data.frame(TF_OC_LL_subset)

write_xlsx(TF_OC_LL_subset, '/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/TF_Ofus_Late_Ctel_Late.xlsx')





###################################### GO Enrichment Analysis ####################################

geneID2GO <- readMappings(file = "/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/22-Gene_specific_heterochrony/02-New_clustering/Owenia_fusiformis_geneID2GO_topGO_onlyGOgenes.txt")
geneUniverse <- names(geneID2GO)

gene_EE_names <- as.character(ofus_ctel_EE$Gene_ID)
gene_EE_list_for_GO <- factor(as.integer(geneUniverse %in% gene_EE_names))
names(gene_EE_list_for_GO) <- geneUniverse

gene_EL_names <- as.character(ofus_ctel_EL$Gene_ID)
gene_EL_list_for_GO <- factor(as.integer(geneUniverse %in% gene_EL_names))
names(gene_EL_list_for_GO) <- geneUniverse

gene_LE_names <- as.character(ofus_ctel_LE$Gene_ID)
gene_LE_list_for_GO <- factor(as.integer(geneUniverse %in% gene_LE_names))
names(gene_LE_list_for_GO) <- geneUniverse

gene_LL_names <- as.character(ofus_ctel_LL$Gene_ID)
gene_LL_list_for_GO <- factor(as.integer(geneUniverse %in% gene_LL_names))
names(gene_LL_list_for_GO) <- geneUniverse

#################################################################################################


gene_EE_GOdata_BP <- new("topGOdata", description="gene_EE_program_BP",
                           ontology="BP", allGenes=gene_EE_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) 
gene_EE_resultFisher_BP <- runTest(gene_EE_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_gene_EE_BP <- GenTable(gene_EE_GOdata_BP, classicFisher = gene_EE_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)


gene_EL_GOdata_BP <- new("topGOdata", description="gene_EL_program_BP",
                         ontology="BP", allGenes=gene_EL_list_for_GO,
                         annot = annFUN.gene2GO, gene2GO = geneID2GO) 
gene_EL_resultFisher_BP <- runTest(gene_EL_GOdata_BP,
                                   algorithm="classic", statistic="fisher") 
results_gene_EL_BP <- GenTable(gene_EL_GOdata_BP, classicFisher = gene_EL_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)


gene_LE_GOdata_BP <- new("topGOdata", description="gene_LE_program_BP",
                         ontology="BP", allGenes=gene_LE_list_for_GO,
                         annot = annFUN.gene2GO, gene2GO = geneID2GO) 
gene_LE_resultFisher_BP <- runTest(gene_LE_GOdata_BP,
                                   algorithm="classic", statistic="fisher") 
results_gene_LE_BP <- GenTable(gene_LE_GOdata_BP, classicFisher = gene_LE_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)


gene_LL_GOdata_BP <- new("topGOdata", description="gene_LL_program_BP",
                         ontology="BP", allGenes=gene_LL_list_for_GO,
                         annot = annFUN.gene2GO, gene2GO = geneID2GO) 
gene_LL_resultFisher_BP <- runTest(gene_LL_GOdata_BP,
                                   algorithm="classic", statistic="fisher") 
results_gene_LL_BP <- GenTable(gene_LL_GOdata_BP, classicFisher = gene_LL_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)

############################################################################ Export GOterms #########################################################

write.table(results_gene_EE_BP, file='GOterms_results_gene_EE_BP.txt',sep='\t', quote=FALSE, col.names = T, row.names = F)
write.table(results_gene_EL_BP, file='GOterms_results_gene_EL_BP.txt',sep='\t', quote=FALSE, col.names = T, row.names = F)
write.table(results_gene_LE_BP, file='GOterms_results_gene_LE_BP.txt',sep='\t', quote=FALSE, col.names = T, row.names = F)
write.table(results_gene_LL_BP, file='GOterms_results_gene_LL_BP.txt',sep='\t', quote=FALSE, col.names = T, row.names = F)

####################################################################################################

gene_EE_GOdata_MF <- new("topGOdata", description="gene_EE_program_MF",
                           ontology="MF", allGenes=gene_EE_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) 
gene_EE_resultFisher_MF <- runTest(gene_EE_GOdata_MF,
                                     algorithm="classic", statistic="fisher") 
results_gene_EE_MF <- GenTable(gene_EE_GOdata_MF, classicFisher = gene_EE_resultFisher_MF, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)


gene_EL_GOdata_MF <- new("topGOdata", description="gene_EL_program_MF",
                         ontology="MF", allGenes=gene_EL_list_for_GO,
                         annot = annFUN.gene2GO, gene2GO = geneID2GO) 
gene_EL_resultFisher_MF <- runTest(gene_EL_GOdata_MF,
                                   algorithm="classic", statistic="fisher") 
results_gene_EL_MF <- GenTable(gene_EL_GOdata_MF, classicFisher = gene_EL_resultFisher_MF, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)


gene_LE_GOdata_MF <- new("topGOdata", description="gene_LE_program_MF",
                         ontology="MF", allGenes=gene_LE_list_for_GO,
                         annot = annFUN.gene2GO, gene2GO = geneID2GO) 
gene_LE_resultFisher_MF <- runTest(gene_LE_GOdata_MF,
                                   algorithm="classic", statistic="fisher") 
results_gene_LE_MF <- GenTable(gene_LE_GOdata_MF, classicFisher = gene_LE_resultFisher_MF, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)


gene_LL_GOdata_MF <- new("topGOdata", description="gene_LL_program_MF",
                         ontology="MF", allGenes=gene_LL_list_for_GO,
                         annot = annFUN.gene2GO, gene2GO = geneID2GO) 
gene_LL_resultFisher_MF <- runTest(gene_LL_GOdata_MF,
                                   algorithm="classic", statistic="fisher") 
results_gene_LL_MF <- GenTable(gene_LL_GOdata_MF, classicFisher = gene_LL_resultFisher_MF, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)

#########################################################################################################

results_gene_EE_BP$classicFisher[results_gene_EE_BP$classicFisher == '<1e-30'] <- 1e-30
results_gene_EL_BP$classicFisher[results_gene_EL_BP$classicFisher == '<1e-30'] <- 1e-30
results_gene_LE_BP$classicFisher[results_gene_LE_BP$classicFisher == '<1e-30'] <- 1e-30
results_gene_LL_BP$classicFisher[results_gene_LL_BP$classicFisher == '<1e-30'] <- 1e-30

results_gene_EE_MF$classicFisher[results_gene_EE_MF$classicFisher == '<1e-30'] <- 1e-30
results_gene_EL_MF$classicFisher[results_gene_EL_MF$classicFisher == '<1e-30'] <- 1e-30
results_gene_LE_MF$classicFisher[results_gene_LE_MF$classicFisher == '<1e-30'] <- 1e-30
results_gene_LL_MF$classicFisher[results_gene_LL_MF$classicFisher == '<1e-30'] <- 1e-30

############################################### BP PLot #################################################


goEnrichment <- results_gene_EE_BP
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

gene_EE_plot_BP <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle(" Owenia Early / Capitella Early (Biological Process)") +
  scale_y_continuous(limits=c(0,30),breaks=round(seq(0,30, by = 2), 1)) +
  theme_classic() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, vjust=0.5),
    axis.title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()



goEnrichment <- results_gene_EL_BP
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

gene_EL_plot_BP <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle(" Owenia Early / Capitella Late (Biological Process)") +
  scale_y_continuous(limits=c(0,30),breaks=round(seq(0,30, by = 2), 1)) +
  theme_classic() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, vjust=0.5),
    axis.title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()



goEnrichment <- results_gene_LE_BP
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

gene_LE_plot_BP <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle(" Owenia Late / Capitella Early (Biological Process)") +
  scale_y_continuous(limits=c(0,30),breaks=round(seq(0,30, by = 2), 1)) +
  theme_classic() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, vjust=0.5),
    axis.title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()



goEnrichment <- results_gene_LL_BP
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

gene_LL_plot_BP <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle(" Owenia Late / Capitella Late (Biological Process)") +
  scale_y_continuous(limits=c(0,30),breaks=round(seq(0,30, by = 2), 1)) +
  theme_classic() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, vjust=0.5),
    axis.title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

plot_grid(gene_EE_plot_BP, gene_EL_plot_BP,
          gene_LE_plot_BP, gene_LL_plot_BP,
          ncol = 2, align = "v")

plot_grid(gene_EE_plot_BP,
          gene_LL_plot_BP,
          ncol = 1, align = "v")


####################################### MF PLot ####################################################



goEnrichment <- results_gene_EE_MF
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

gene_EE_plot_MF <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Molecular Function") +
  ylab("-log10(p-value)") +
  ggtitle("Owenia Early / Capitella Early (Molecular Function)") +
  scale_y_continuous(limits=c(0,30),breaks=round(seq(0,30, by = 2), 1)) +
  theme_classic() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, vjust=0.5),
    axis.title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()


goEnrichment <- results_gene_EL_MF
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

gene_EL_plot_MF <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Owenia Early / Capitella Late (Molecular Function)") +
  scale_y_continuous(limits=c(0,30),breaks=round(seq(0,30, by = 2), 1)) +
  theme_classic() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, vjust=0.5),
    axis.title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()


goEnrichment <- results_gene_LE_MF
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

gene_LE_plot_MF <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Owenia Late / Capitella Early (Molecular Function)") +
  scale_y_continuous(limits=c(0,30),breaks=round(seq(0,30, by = 2), 1)) +
  theme_classic() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, vjust=0.5),
    axis.title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()


goEnrichment <- results_gene_LL_MF
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

gene_LL_plot_MF <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Molecular Function") +
  ylab("-log10(p-value)") +
  ggtitle("Owenia Late / Capitella Late (Molecular Function)") +
  scale_y_continuous(limits=c(0,30),breaks=round(seq(0,30, by = 2), 1)) +
  theme_classic() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, vjust=0.5),
    axis.title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()


plot_grid(gene_EE_plot_MF, gene_EL_plot_MF,
          gene_LE_plot_MF, gene_LL_plot_MF,
          ncol = 2, align = "v")

plot_grid(gene_EE_plot_MF,
          gene_LL_plot_MF,
          ncol = 1, align = "v")

######################################## Go Enrichment all genes at time points #################################

G_DCO_subset # Take transcript ID

gene_DCO_names <- as.character(G_DCO_subset$transcript_id)
gene_DCO_list_for_GO <- factor(as.integer(geneUniverse %in% gene_DCO_names))
names(gene_DCO_list_for_GO) <- geneUniverse

gene_DCO_GOdata_BP <- new("topGOdata", description="gene_DCO_program_BP",
                         ontology="BP", allGenes=gene_DCO_list_for_GO,
                         annot = annFUN.gene2GO, gene2GO = geneID2GO) 
gene_DCO_resultFisher_BP <- runTest(gene_DCO_GOdata_BP,
                                   algorithm="classic", statistic="fisher") 
results_gene_DCO_BP <- GenTable(gene_DCO_GOdata_BP, classicFisher = gene_DCO_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)


gene_DCO_GOdata_MF <- new("topGOdata", description="gene_DCO_program_MF",
                         ontology="MF", allGenes=gene_DCO_list_for_GO,
                         annot = annFUN.gene2GO, gene2GO = geneID2GO) 
gene_DCO_resultFisher_MF <- runTest(gene_DCO_GOdata_MF,
                                   algorithm="classic", statistic="fisher") 
results_gene_DCO_MF <- GenTable(gene_DCO_GOdata_MF, classicFisher = gene_DCO_resultFisher_MF, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)



results_gene_DCO_BP$classicFisher[results_gene_DCO_BP$classicFisher == '< 1e-30'] <- 1e-30
results_gene_DCO_MF$classicFisher[results_gene_DCO_MF$classicFisher == '< 1e-30'] <- 1e-30


goEnrichment <- results_gene_DCO_BP
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

gene_DCO_plot_BP <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle(" Dimorphilus / Capitella / Owenia (Biological Process)") +
  scale_y_continuous(limits=c(0,20),breaks=round(seq(0,20, by = 2), 1)) +
  theme_classic() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, vjust=0.5),
    axis.title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

goEnrichment <- results_gene_DCO_MF
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

gene_DCO_plot_MF <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Molecular Function") +
  ylab("-log10(p-value)") +
  ggtitle("Dimorphilus / Capitella / Owenia (Molecular Function)") +
  scale_y_continuous(limits=c(0,20),breaks=round(seq(0,20, by = 2), 1)) +
  theme_classic() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, vjust=0.5),
    axis.title=element_text(size=12),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

plot_grid(gene_DCO_plot_MF, gene_DCO_plot_BP,
          ncol = 1, align = "v")
