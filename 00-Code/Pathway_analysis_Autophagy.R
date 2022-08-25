# 1. Obtain ENSEMBL Gene and Transcript IDs for the HGNC symbols of the 
# genes we have fetched from the KEGG autophagy pathway for H. sapiens

library(biomaRt)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

transcript_ids <- read.table("02-Human_gene_IDs.txt", sep ='\t', header = F)
colnames(transcript_ids) <- c("HGNC_gene_symbol","Full_name")

res <- getBM(attributes = c('hgnc_symbol', 
                            'ensembl_gene_id',
                            'ensembl_transcript_id_version'),
             filters = 'hgnc_symbol', 
             values = transcript_ids$HGNC_gene_symbol,
             mart = mart)

length(unique(transcript_ids$HGNC_gene_symbol)) == length(unique(res$hgnc_symbol)) # continue if TRUE

colnames(res) <- c("HGNC_gene_symbol", "ENSEMBL_Gene_ID", "ENSEMBL_Transcript_ID")
res$ENSEMBL_Transcript_ID_custom <- paste("Hsap", gsub("\\..*","",res$ENSEMBL_Transcript_ID), sep ='_')

write.table(res, "03-Human_ENSEMBL_ID.txt", sep = '\t', quote = F, row.names = F)


# 2. Get orthologs for Ofus, and from this, for Ctel

ofus2hsap <- read.table("04-Ofus2Hsap.txt", sep = '\t', header = T)
ofus2ctel <- read.table("04-Ofus2Ctel.txt", sep = '\t', header = T)
hsap2ofus <- ofus2hsap[c(3,2)]
rownames(hsap2ofus) <- hsap2ofus$Hsap
ofus2ctel <- ofus2ctel[c(3,2)]
rownames(ofus2ctel) <- ofus2ctel$Ofus

res2 <- res
rownames(res2) <- res2$ENSEMBL_Transcript_ID_custom

hsap2ofus_orthology <- hsap2ofus[hsap2ofus$Hsap %in% res2$ENSEMBL_Transcript_ID_custom,]
hsap2ofus_orthology$HGNC_gene_symbol <- res2[rownames(hsap2ofus_orthology),"HGNC_gene_symbol"]
hsap2ofus_orthology$Ctel <- ofus2ctel[hsap2ofus_orthology$Ofus,"Ctel"]

res3 <- getBM(attributes = c('hgnc_symbol', 
                             'ensembl_gene_id'),
              filters = 'hgnc_symbol', 
              values = transcript_ids$HGNC_gene_symbol,
              mart = mart)

final_dataset <- res3
colnames(final_dataset) <- c("HGNC_gene_symbol", "ENSEMBL_Gene_ID")
final_dataset$orthology <- ifelse(final_dataset$HGNC_gene_symbol %in% hsap2ofus_orthology$HGNC_gene_symbol, "Yes", "No")
final_dataset$Hsap <- ifelse(final_dataset$HGNC_gene_symbol %in% hsap2ofus_orthology$HGNC_gene_symbol, hsap2ofus_orthology$Hsap, NA)
final_dataset$Ofus <- ifelse(final_dataset$HGNC_gene_symbol %in% hsap2ofus_orthology$HGNC_gene_symbol, hsap2ofus_orthology$Ofus, NA)
final_dataset$Ctel <- ifelse(final_dataset$HGNC_gene_symbol %in% hsap2ofus_orthology$HGNC_gene_symbol, hsap2ofus_orthology$Ctel, NA)

write.table(final_dataset, "05-Annelid_autophagy_orthology.txt", sep ='\t', quote = F, row.names = F)


# 3. Get quadrant for genes

library(readxl)

oece <- read_xlsx("G_Ofus_Early_Ctel_Early.xlsx")
oece$quadrant <- "early"
oecl <- read_xlsx("G_Ofus_Early_Ctel_Late.xlsx")
oecl$quadrant <- "delayed"
olce <- read_xlsx("G_Ofus_Late_Ctel_Early.xlsx")
olce$quadrant <- "accelerated"
olcl <- read_xlsx("G_Ofus_Late_Ctel_Late.xlsx")
olcl$quadrant <- "late"

all <- rbind(oece, oecl, olce, olcl)
all <- data.frame(all)
rownames(all) <- all$transcript_id

quadrant_analysis <- final_dataset
quadrant_analysis$quadrant <- all[quadrant_analysis$Ofus, "quadrant"]

write.table(quadrant_analysis, "06-Autophagy_annelids_by_quadrants.txt", sep = '\t', quote = F, row.names = F)



# 4. Try to get ALL orthologs, not just the 1-to-1 ones
# Then remove the ones that cannot be resolved for Owenia (more than 1 gene ID per entry)
# 

ofus2hsap_full <- read.table("07-Ofus2Hsap_full.txt", sep ='\t', header = T)
hsap2ofus_full <- ofus2hsap_full[c(3,2)]
hsap2ofus_full[,3] <- rep(NA, nrow(hsap2ofus_full))
hsap2ofus_full[,4] <- rep(NA, nrow(hsap2ofus_full))
colnames(hsap2ofus_full)[c(3,4)] <- c("Autophagy_check", "HGNC_gene_symbol")

for (i in rownames(hsap2ofus_full)){
  summary <- unlist(strsplit(hsap2ofus_full[i,"Hsap"], ","))
  for (j in summary){
    if (j %in% res2$ENSEMBL_Transcript_ID_custom){
      hsap2ofus_full[i,"Autophagy_check"] <- "Autophagy"
      hsap2ofus_full[i,"HGNC_gene_symbol"] <- res2[j,"HGNC_gene_symbol"]
    }
  }
}


hsap2ofus_full_orthology <- hsap2ofus_full[!is.na(hsap2ofus_full$Autophagy_check),]
hsap2ofus_full_orthology$uniqueness_check <- rep(NA, nrow(hsap2ofus_full_orthology))

for (i in rownames(hsap2ofus_full_orthology)){
  hsap2ofus_full_orthology[i,"uniqueness_check"] <- length(unlist(strsplit(hsap2ofus_full_orthology[i,"Ofus"], ",")))
}

copy_file <- transcript_ids
rownames(copy_file) <- copy_file$HGNC_gene_symbol

hsap2ofus_full_orthology_single <- hsap2ofus_full_orthology[hsap2ofus_full_orthology$uniqueness_check == 1,c(1,2,4)]
hsap2ofus_full_orthology_single$Ctel <- ofus2ctel[hsap2ofus_full_orthology_single$Ofus,"Ctel"]
hsap2ofus_full_orthology_final <- hsap2ofus_full_orthology_single[!is.na(hsap2ofus_full_orthology_single$Ctel),c(1,2,4,3)]
hsap2ofus_full_orthology_final$Full_name <- copy_file[hsap2ofus_full_orthology_final$HGNC_gene_symbol,"Full_name"]
hsap2ofus_full_orthology_final$quadrant <- all[hsap2ofus_full_orthology_final$Ofus, "quadrant"]

table(hsap2ofus_full_orthology_final$quadrant)

write.table(hsap2ofus_full_orthology_final, "08-Autophagy_annelids_complete_by_quadrants.txt", sep = '\t', quote = F, row.names = F)
