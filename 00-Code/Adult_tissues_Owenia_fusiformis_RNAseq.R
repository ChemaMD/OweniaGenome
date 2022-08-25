# Necessary libraries
library(tximport)
library(DESeq2)
library(rhdf5)

# Get the working directory: make sure it's the folder where the folders from the kallisto outputs are located
dir <- getwd()

# Open the table we created
samples <- read.table(file.path(dir, "owenia_tissue_analysis_samples_table.txt"), header = TRUE) # change <table> accordingly to the name you've given to the table we created above

# Create a vector named files containing the paths of the abundance.h5 files of each experiment, which we get from the "run" column (i.e., samples$run).
files <- file.path(dir, samples$run, "abundance.h5")

# Make sure that all the files exist
all(file.exists(files))

names(files) <- samples$experiment

# Import kallisto files, that is, the abundance.h5 files
kallisto_output <- tximport(files, type = "kallisto", txOut = TRUE)

# Create an empty data frame for each of the samples and replicates
df <- data.frame(stage = factor(rep(c("tail","testis","ovary","retractor_muscle",
                                      "head_chaetae", "gut", "head", "blood_vessel", 
                                      "body_wall"), each=1)),
                                    replicate = factor(c("rep-1","rep-1","rep-1","rep-1",
                                                         "rep-1","rep-1","rep-1","rep-1","rep-1")))
rownames(df) <- colnames(kallisto_output$counts)
dds <- DESeqDataSetFromTximport(kallisto_output, colData = df, design = ~1) #it works but won't work when passed to the DESeq function
# dds stands for DESeq data set
dds_DEseq2 <- DESeq(dds)

# Get the table with normalised counts
counts_DESeq2 <- counts(dds_DEseq2, normalized=TRUE) # change the name if you want
counts_DESeq2 <- as.data.frame(counts_DESeq2)
counts_DESeq2 <- rename(counts_DESeq2, tail = owe_tiss_1, testis = owe_tiss_2, ovary = owe_tiss_3,
                         retractor_muscle = owe_tiss_4, head_chaetae = owe_tiss_5, gut = owe_tiss_6, 
                         head = owe_tiss_7, blood_vessel = owe_tiss_8, body_wall = owe_tiss_9)
write.table(counts_DESeq2, "DESeq2_normalised_counts_tissues.txt", sep="\t", quote = F) # change the name if you want

#creating tpm matrix
#getting this done with the other scripts conventions and names and then renaming it all
library(tximport)
library(rhdf5)
library(dplyr)

dir <- getwd()
samples <- read.table(file.path(dir, "owenia_tissue_analysis_samples_table.txt"), header = TRUE)
files <- file.path(dir, samples$run, "abundance.h5")
all(file.exists(files))
names(files) <- samples$experiment
tximport_kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
head(tximport_kallisto$counts)
head(tximport_kallisto[["abundance"]])
tpm_tissues <- as.data.frame(tximport_kallisto[["abundance"]])
tpm_tissues <- rename(tpm_tissues, tail = owe_tiss_1, testis = owe_tiss_2, ovary = owe_tiss_3,
                        retractor_muscle = owe_tiss_4, head_chaetae = owe_tiss_5, gut = owe_tiss_6, 
                        head = owe_tiss_7, blood_vessel = owe_tiss_8, body_wall = owe_tiss_9)
write.table(tpm_tissues, "owenia_tpm_tissues.txt", sep="\t",quote=F)

#correlation matrices

# Necessary libraries
library(ggplot2)
library(corrplot)
library(RColorBrewer)

# Import gene expression matrix in TPM. This is an example for Owenia fusiformis' stage-wise RNA-seq time course from the blastula stage to the juvenile adult
df_TPM <- read.table("owenia_tpm_tissues.txt", header = TRUE)
colnames(df_TPM) <- c("tail", "testis",
                      "ovary", "retractor_muscle", "head_chaetae", "gut",
                      "head", "blood_vessel", "body_wall")

# Create a colour ramp based on RdBu. Pearson correlation coefficient is automatically displayed in a [-1,+1] range. 
#However, we are only interested in the [0,+1] range. For that matter, we colour-code the [0,+1] range with the RdBu palette, and the [-1,0) 
#with the Greys palette. In the end, we'll crop the [-1,0) range of the legend in Inkscape or in Adobe Illustrator.

# Create palette for the [0,+1] range
palet1 = colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(100)

# Create palette for the [-1,0) range
palet2 = colorRampPalette(rev(brewer.pal(n = 7, name = "Greys")))(100)

# Join both palettes into a single one
palet = cbind(palet2,palet1)

# Calculate the correlation using the Pearson correlation coefficient. This can be changed to others (e.g. Pearson, Euclidean, etc., look the documentation).
cc = cor(df_TPM, method = "pearson")

# Plot the correlation plot
corrplot(cc, tl.col = "black", order = "hclust", method = "color", col=palet, cl.lim=c(0,1), hclust())

#now for the normalised data set

df_normalised <- read.table("DESeq2_normalised_counts_tissues.txt", header = TRUE)
colnames(df_normalised) <- c("tail", "testis",
                      "ovary", "retractor_muscle", "head_chaetae", "gut",
                      "head", "blood_vessel", "body_wall")

# Create palette for the [0,+1] range
palet1 = colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(100)

# Create palette for the [-1,0) range
palet2 = colorRampPalette(rev(brewer.pal(n = 7, name = "Greys")))(100)

# Join both palettes into a single one
palet = cbind(palet2,palet1)

# Calculate the correlation using the Pearson correlation coefficient. This can be changed to others (e.g. Pearson, Euclidean, etc., look the documentation).
cc_norm = cor(df_normalised, method = "pearson")

# Plot the correlation plot
corrplot(cc_norm, tl.col = "black", order = "hclust", method = "color", col=palet, cl.lim=c(0,1))
heatmaply_cor(cc_norm)
library(heatmaply)

#subsetting genes, >2 tpm, by tissue
tail_tpm <- select(df_TPM, tail)
testis_tpm <- select(df_TPM, testis)  
ovary_tpm <- select(df_TPM, ovary)  
retractor_tpm <- select(df_TPM, retractor_muscle)  
chaetae_tpm <- select(df_TPM, head_chaetae)  
gut_tpm <- select(df_TPM, gut)
head_tpm <- select(df_TPM, head)
blood_vessel_tpm <- select(df_TPM, blood_vessel)
body_wall_tpm <- select(df_TPM, body_wall)
#tpm subsetted
write.table(tail_tpm, "tail_tpm_table.txt", sep="\t", quote = F)
write.table(testis_tpm, "testis_tpm_table.txt", sep="\t", quote = F)
write.table(ovary_tpm, "ovary_tpm_table.txt", sep="\t", quote = F)
write.table(retractor_tpm, "retractor_tpm_table.txt", sep="\t", quote = F)
write.table(chaetae_tpm, "chaetae_tpm_table.txt", sep="\t", quote = F)
write.table(gut_tpm, "gut_tpm_table.txt", sep="\t", quote = F)
write.table(head_tpm, "head_tpm_table.txt", sep="\t", quote = F)
write.table(blood_vessel_tpm, "blood_vessel_tpm_table.txt", sep="\t", quote = F)
write.table(body_wall_tpm, "body_wall_tpm_table.txt", sep="\t", quote = F)

#all tpm per tissue >2
tail_tpm2 <- filter(tail_tpm, tail > 2)
testis_tpm2 <- filter(testis_tpm, testis > 2)
ovary_tpm2 <- filter(ovary_tpm, ovary > 2)
retractor_tpm2 <- filter(retractor_tpm, retractor_muscle > 2)
chaetae_tpm2 <- filter(chaetae_tpm, head_chaetae > 2)
gut_tpm2 <- filter(gut_tpm, gut > 2)
head_tpm2 <- filter(head_tpm, head > 2)
blood_vessel_tpm2 <- filter(blood_vessel_tpm, blood_vessel > 2)
body_tpm2 <- filter(body_wall_tpm, body_wall > 2)
#tpm>2 subsetted
write.table(tail_tpm2, "tail_tpm2_table.txt", sep="\t", quote = F)
write.table(testis_tpm2, "testis_tpm2_table.txt", sep="\t", quote = F)
write.table(ovary_tpm2, "ovary_tpm2_table.txt", sep="\t", quote = F)
write.table(retractor_tpm2, "retractor_tpm2_table.txt", sep="\t", quote = F)
write.table(chaetae_tpm2, "chaetae_tpm2_table.txt", sep="\t", quote = F)
write.table(gut_tpm2, "gut_tpm2_table.txt", sep="\t", quote = F)
write.table(head_tpm2, "head_tpm2_table.txt", sep="\t", quote = F)
write.table(blood_vessel_tpm2, "blood_vessel_tpm2_table.txt", sep="\t", quote = F)
write.table(body_tpm2, "body_wall_tpm2_table.txt", sep="\t", quote = F)

#truly unique genes for the bodywall, tail, head and prostomium
#filter out all genes with tpm lower than 2, then only selecting rows which contain genes that experience only 
#tail related expression of genes
#these aren't necessary now but keeping it here anyway
unique_tail <- filter(df_TPM, tail > 2)
unique_tail <- filter(unique_tail, unique_tail$tail >= 2 & unique_tail$testis < 2 & unique_tail$ovary < 2 
                      & unique_tail$retractor_muscle < 2  & unique_tail$head_chaetae < 2  & unique_tail$gut < 2 
                      & unique_tail$head < 2  & unique_tail$blood_vessel < 2  & unique_tail$body_wall < 2 )
unique_bodywall <- filter(df_TPM, body_wall >2)
unique_bodywall <- filter(unique_bodywall, unique_bodywall$body_wall >= 2 & unique_bodywall$testis < 2  & unique_bodywall$ovary < 2 
                          & unique_bodywall$retractor_muscle< 2  & unique_bodywall$head_chaetae < 2  & unique_bodywall$gut < 2 
                          & unique_bodywall$head < 2  & unique_bodywall$blood_vessel < 2  & unique_bodywall$tail < 2 )
unique_head <- filter(df_TPM, head > 2)
unique_head <- filter(unique_head, unique_head$head >= 2 & unique_head$testis < 2  & unique_head$ovary < 2 
                      & unique_head$retractor_muscle < 2  & unique_head$head_chaetae < 2  & unique_head$gut < 2 
                      & unique_head$body_wall < 2  & unique_head$blood_vessel < 2  & unique_head$tail < 2 )
unique_prostomium <- filter(df_TPM, head_chaetae > 2)
unique_prostomium <- filter(unique_prostomium, unique_prostomium$head_chaetae >= 2 & unique_prostomium$testis < 2  & unique_prostomium$ovary < 2
                            & unique_prostomium$retractor_muscle < 2 & unique_prostomium$head < 2 & unique_prostomium$gut < 2 
                            & unique_prostomium$body_wall < 2 & unique_prostomium$blood_vessel < 2 & unique_prostomium$tail < 2 )
unique_anterior <- filter(df_TPM, head >2 & head_chaetae > 2)
unique_anterior <- filter(unique_anterior, unique_anterior$head_chaetae >= 2 & unique_anterior$head >= 2
                          & unique_anterior$testis < 2 & unique_anterior$ovary < 2
                          & unique_anterior$retractor_muscle < 2 & unique_anterior$gut < 2
                          & unique_anterior$body_wall < 2 & unique_anterior$blood_vessel < 2 & unique_anterior$tail < 2)


#venn diagrams
#these will compare the four interesting data bits: the anterior tissues, head and the Head_chaetae (the prostomium) 
#and the posterior tissues of interest, the tail and the body wall.

library(ggvenn)
library(tibble)
library(gplots)
library(VennDiagram)
head_tpm2 <- rownames_to_column(head_tpm2, "gene_id")
prostomium_tpm2 <- rownames_to_column(chaetae_tpm2, "gene_id")
tail_tpm2 <- rownames_to_column(tail_tpm2, "gene_id")
body_wall_tpm2 <- rownames_to_column(body_tpm2, "gene_id")

x <- list(Head = head_tpm2$gene_id, Prostomium = prostomium_tpm2$gene_id, 
          Tail = tail_tpm2$gene_id, Body_wall = body_wall_tpm2$gene_id)

ggvenn(x, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4)
v.table <- venn(x)
overlap <- calculate.overlap(x)

#posterior should be the tail, the bodywall and overlap
unique_genes_posterior <- c(overlap$a3, overlap$a2, overlap$a1)
#anterior should be the head, prostomium and the overlap
unique_anterior_genes <- c(overlap$a9, overlap$a14, overlap$a15)
unique_prostomium_genes <- overlap$a14
unique_tail_genes <- overlap$a1
unique_head_genes <- overlap$a9
unique_bodywall_genes <- overlap$a3
unique_posterior_overlap <- overlap$a2
unique_anterior_overlap <- overlap$a15
#subsetting and saving the names of the unique genes from each interesting segment of the venn diagram
write.table(unique_genes_posterior, "unique_posterior_genes.txt",sep="\t", quote = F)
write.table(unique_anterior_genes, "unique_anterior_genes.txt",sep="\t", quote = F)
write.table(unique_prostomium_genes, "unique_prostomium_genes.txt",sep="\t", quote = F)
write.table(unique_tail_genes, "unique_tail_genes.txt",sep="\t", quote = F)
write.table(unique_head_genes, "unique_head_genes.txt",sep="\t", quote = F)
write.table(unique_bodywall_genes, "unique_bodywall_genes.txt",sep="\t", quote = F)
write.table(unique_posterior_overlap, "unique_posterior_overlap_genes.txt",sep="\t", quote = F)
write.table(unique_anterior_overlap, "unique_anterior_overlap_genes.txt",sep="\t", quote = F)

#go term enrichment of the venn diagram slices above

library(ggplot2)
library(ggpubr)
library(topGO)
library(cowplot)

geneID2GO <- readMappings(file = "Ofus_geneID2GO_topGO_onlyGOgenes.txt") # change <GO-universe> according to the name that you gave to the original file
geneUniverse <- names(geneID2GO)

#posterior cluster first
cluster_1 <- read.table("unique_posterior_genes.txt", header=TRUE)
cluster_1_names <- as.character(cluster_1$x)
cluster_1_list_for_GO <- factor(as.integer(geneUniverse %in% cluster_1_names))
names(cluster_1_list_for_GO) <- geneUniverse

cluster_1_GOdata_BP <- new("topGOdata", description="cluster_1_program",
                           ontology="BP", allGenes=cluster_1_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) # change BP for MF or CC if needed
cluster_1_resultFisher_BP <- runTest(cluster_1_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_cluster_1 <- GenTable(cluster_1_GOdata_BP, classicFisher = cluster_1_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30) 

write.table(results_cluster_1, "unique_posterior_clusterGO.txt", quote=FALSE, sep='\t') # change <cluster-1> according to the name you wanna give the table
#for if you have to restart
results_cluster_1 <- read.table("unique_posterior_clusterGO.txt", header=TRUE)

goEnrichment <- results_cluster_1
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

cluster_1_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Posterior genes") +
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
cluster_1_plot

#anterior cluster
cluster_2 <- read.table("unique_anterior_genes.txt", header=TRUE)
cluster_2_names <- as.character(cluster_2$x)
cluster_2_list_for_GO <- factor(as.integer(geneUniverse %in% cluster_2_names))
names(cluster_2_list_for_GO) <- geneUniverse

cluster_2_GOdata_BP <- new("topGOdata", description="cluster_2_program",
                           ontology="BP", allGenes=cluster_2_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) # change BP for MF or CC if needed
cluster_2_resultFisher_BP <- runTest(cluster_2_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_cluster_2 <- GenTable(cluster_2_GOdata_BP, classicFisher = cluster_2_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30) 

write.table(results_cluster_2, "unique_anterior_clusterGO.txt", quote=FALSE, sep='\t') # change <cluster-1> according to the name you wanna give the table
#for if you have to restart
results_cluster_2 <- read.table("unique_anterior_clusterGO.txt", header=TRUE)

goEnrichment <- results_cluster_2
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

cluster_2_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Anterior genes") +
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
cluster_2_plot

#prostomium cluster
cluster_3 <- read.table("unique_prostomium_genes.txt", header=TRUE)
cluster_3_names <- as.character(cluster_3$x)
cluster_3_list_for_GO <- factor(as.integer(geneUniverse %in% cluster_3_names))
names(cluster_3_list_for_GO) <- geneUniverse

cluster_3_GOdata_BP <- new("topGOdata", description="cluster_3_program",
                           ontology="BP", allGenes=cluster_3_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) # change BP for MF or CC if needed
cluster_3_resultFisher_BP <- runTest(cluster_3_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_cluster_3 <- GenTable(cluster_3_GOdata_BP, classicFisher = cluster_3_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30) 

write.table(results_cluster_3, "unique_prostomium_clusterGO.txt", quote=FALSE, sep='\t') # change <cluster-1> according to the name you wanna give the table
#for if you have to restart
results_cluster_3 <- read.table("unique_prostomium_clusterGO.txt", header=TRUE)

goEnrichment <- results_cluster_3
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

cluster_3_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Prostomium genes") +
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
cluster_3_plot

#tail cluster
cluster_4 <- read.table("unique_tail_genes.txt", header=TRUE)
cluster_4_names <- as.character(cluster_4$x)
cluster_4_list_for_GO <- factor(as.integer(geneUniverse %in% cluster_4_names))
names(cluster_4_list_for_GO) <- geneUniverse

cluster_4_GOdata_BP <- new("topGOdata", description="cluster_4_program",
                           ontology="BP", allGenes=cluster_4_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) # change BP for MF or CC if needed
cluster_4_resultFisher_BP <- runTest(cluster_4_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_cluster_4 <- GenTable(cluster_4_GOdata_BP, classicFisher = cluster_4_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30) 

write.table(results_cluster_4, "unique_tail_clusterGO.txt", quote=FALSE, sep='\t') # change <cluster-1> according to the name you wanna give the table
#for if you have to restart
results_cluster_2 <- read.table("unique_tail_clusterGO.txt", header=TRUE)

goEnrichment <- results_cluster_4
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

cluster_4_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Tail genes") +
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
cluster_4_plot

#head cluster
cluster_5 <- read.table("unique_head_genes.txt", header=TRUE)
cluster_5_names <- as.character(cluster_5$x)
cluster_5_list_for_GO <- factor(as.integer(geneUniverse %in% cluster_5_names))
names(cluster_5_list_for_GO) <- geneUniverse

cluster_5_GOdata_BP <- new("topGOdata", description="cluster_5_program",
                           ontology="BP", allGenes=cluster_5_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) # change BP for MF or CC if needed
cluster_5_resultFisher_BP <- runTest(cluster_5_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_cluster_5 <- GenTable(cluster_5_GOdata_BP, classicFisher = cluster_5_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30) 

write.table(results_cluster_5, "unique_head_clusterGO.txt", quote=FALSE, sep='\t') # change <cluster-1> according to the name you wanna give the table
#for if you have to restart
results_cluster_5 <- read.table("unique_head_clusterGO.txt", header=TRUE)

goEnrichment <- results_cluster_5
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

cluster_5_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Head genes") +
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
cluster_5_plot

#bodywall cluster
cluster_6 <- read.table("unique_bodywall_genes.txt", header=TRUE)
cluster_6_names <- as.character(cluster_6$x)
cluster_6_list_for_GO <- factor(as.integer(geneUniverse %in% cluster_6_names))
names(cluster_6_list_for_GO) <- geneUniverse

cluster_6_GOdata_BP <- new("topGOdata", description="cluster_6_program",
                           ontology="BP", allGenes=cluster_6_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) # change BP for MF or CC if needed
cluster_6_resultFisher_BP <- runTest(cluster_6_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_cluster_6 <- GenTable(cluster_6_GOdata_BP, classicFisher = cluster_6_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30) 

write.table(results_cluster_6, "unique_bodywall_clusterGO.txt", quote=FALSE, sep='\t') # change <cluster-1> according to the name you wanna give the table
#for if you have to restart
results_cluster_6 <- read.table("unique_bodywall_clusterGO.txt", header=TRUE)

goEnrichment <- results_cluster_6
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

cluster_6_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Bodywall genes") +
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
cluster_6_plot

#plot of the enriched clusters all at once
plot_grid(cluster_1_plot, cluster_2_plot, cluster_3_plot, cluster_4_plot,cluster_5_plot,cluster_6_plot,
          ncol = 1, align = "v") +
  ggsave("anterior_posterior_genes.pdf", height = 50, width = 20, units = "cm")

#the thing to do now is
#tracking the expression values of the unique genes over time using the stages of development
development_stages <- read.table("new_7_stages.txt", header = TRUE)
development_stages <- rename(development_stages, gene_id = ID) 
#unique_anterior, unique_head, unique_tail, unique_prostomium, unique_bodywall
unique_prostomium <- rownames_to_column(unique_prostomium, "gene_id")
unique_head <- rownames_to_column(unique_head, "gene_id")
unique_bodywall <- rownames_to_column(unique_bodywall, "gene_id")
unique_tail <- rownames_to_column(unique_tail, "gene_id")
unique_anterior <- rownames_to_column(unique_anterior, "gene_id")

stages_unique_anterior <- left_join(unique_anterior, development_stages, by = "gene_id")
stages_unique_tail <- inner_join(unique_tail, development_stages, by = "gene_id")
stages_unique_bodywall <- inner_join(unique_bodywall, development_stages, by = "gene_id")
stages_unique_prostomium <- inner_join(unique_prostomium, development_stages, by = "gene_id")
stages_unique_head <- inner_join(unique_head, development_stages, by = "gene_id")

write.table(stages_unique_anterior, "anterior_unique_genes.txt", quote = FALSE, sep = '\t')
write.table(stages_unique_tail, "tail_unique_genes.txt", quote = FALSE, sep = '\t')
write.table(stages_unique_bodywall, "bodywall_unique_genes.txt", quote = FALSE, sep = '\t')
write.table(stages_unique_prostomium, "prostomium_unique_genes.txt", quote = FALSE, sep = '\t')
write.table(stages_unique_head, "head_unique_genes.txt", quote = FALSE, sep = '\t')

#BELOW IS USELESS
#go term enrichment of the unique genes
library(ggplot2)
library(ggpubr)
library(topGO)
library(cowplot)

geneID2GO <- readMappings(file = "Ofus_geneID2GO_topGO_onlyGOgenes.txt") # change <GO-universe> according to the name that you gave to the original file
geneUniverse <- names(geneID2GO)

#anterior cluster first
cluster_1 <- read.table("anterior_unique_genes.txt", header=TRUE)
cluster_1_names <- as.character(cluster_1$gene)
cluster_1_list_for_GO <- factor(as.integer(geneUniverse %in% cluster_1_names))
names(cluster_1_list_for_GO) <- geneUniverse

cluster_1_GOdata_BP <- new("topGOdata", description="cluster_1_program",
                           ontology="BP", allGenes=cluster_1_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) # change BP for MF or CC if needed
cluster_1_resultFisher_BP <- runTest(cluster_1_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_cluster_1 <- GenTable(cluster_1_GOdata_BP, classicFisher = cluster_1_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15) 

write.table(results_cluster_1, "anterior_unique_genes_clusterGO.txt", quote=FALSE, sep='\t') # change <cluster-1> according to the name you wanna give the table
#for if you have to restart
results_cluster_1 <- read.table("anterior_unique_genes_clusterGO.txt", header=TRUE)

goEnrichment <- results_cluster_1
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

cluster_1_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Anterior genes") +
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
cluster_1_plot

#tail cluster next
cluster_2 <- read.table("tail_unique_genes.txt", header=TRUE)
cluster_2_names <- as.character(cluster_2$gene)
cluster_2_list_for_GO <- factor(as.integer(geneUniverse %in% cluster_2_names))
names(cluster_2_list_for_GO) <- geneUniverse

cluster_2_GOdata_BP <- new("topGOdata", description="cluster_2_program",
                           ontology="BP", allGenes=cluster_2_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) # change BP for MF or CC if needed
cluster_2_resultFisher_BP <- runTest(cluster_2_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_cluster_2 <- GenTable(cluster_2_GOdata_BP, classicFisher = cluster_2_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15) 

write.table(results_cluster_2, "tail_unique_genes_clusterGO.txt", quote=FALSE, sep='\t') # change <cluster-1> according to the name you wanna give the table
#for if you have to restart
results_cluster_2 <- read.table("tail_unique_genes_clusterGO.txt", header=TRUE)

goEnrichment <- results_cluster_2
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

cluster_2_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Tail genes") +
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
cluster_2_plot

#bodywall genes cluster
cluster_3 <- read.table("bodywall_unique_genes.txt", header=TRUE)
cluster_3_names <- as.character(cluster_3$gene)
cluster_3_list_for_GO <- factor(as.integer(geneUniverse %in% cluster_3_names))
names(cluster_3_list_for_GO) <- geneUniverse

cluster_3_GOdata_BP <- new("topGOdata", description="cluster_3_program",
                           ontology="BP", allGenes=cluster_3_list_for_GO,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO) # change BP for MF or CC if needed
cluster_3_resultFisher_BP <- runTest(cluster_3_GOdata_BP,
                                     algorithm="classic", statistic="fisher") 
results_cluster_3 <- GenTable(cluster_3_GOdata_BP, classicFisher = cluster_3_resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15) 

write.table(results_cluster_3, "bodywall_unique_genes_clusterGO.txt", quote=FALSE, sep='\t') # change <cluster-1> according to the name you wanna give the table
#for if you have to restart
results_cluster_3 <- read.table("bodywall_unique_genes_clusterGO.txt", header=TRUE)

goEnrichment <- results_cluster_3
goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
goEnrichment <- goEnrichment[,c("GO.ID","Term","classicFisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

cluster_3_plot <- ggplot(goEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("Bodywall genes") +
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
cluster_3_plot

plot_grid(cluster_1_plot, cluster_2_plot, cluster_3_plot,
          ncol = 1, align = "v") +
  ggsave("tail_bodywall_anterior_genesGO.pdf", height = 50, width = 20, units = "cm")
#scratch all the stuf about the truly unique genes, sticking with just the venn diagram ones instead
#ABOVE IS USELESS

#loess regression analysis for the isolated genes 
library(ggplot2)
library(ggpubr)
library(gplots)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(matrixStats)
library(data.table)
library(mgcv)
library(scales)
# Create colour palette, in this case we use the viridis one
palette <- viridis_pal(option = "C",direction = 1)(4) # change 4 to the actual number of clusters of your gene clustering 

# Set an x-axis here named time with the number of time points of your time course. In this case it is 7 because we sampled from the blastula stage to the juvenile stage of Owenia fusiformis' development.
time <- c(1:7)

#creating the tables with the data from the averaged replicates and the names of the genes in each cluster
#load in data and conjoin genes with their relevant data
unique_posterior_genes <- read.table("unique_posterior_genes.txt", header = TRUE)
unique_anterior_genes <- read.table("unique_anterior_genes.txt", header = TRUE)
unique_prostomium_genes <- read.table("unique_prostomium_genes.txt", header = TRUE)
unique_tail_genes <- read.table("unique_tail_genes.txt", header = TRUE)
unique_head_genes <- read.table("unique_head_genes.txt", header = TRUE)
unique_bodywall_genes <- read.table("unique_bodywall_genes.txt", header = TRUE)
seven_stages <- read.table("new_7_stages.txt", header = TRUE)
anterior_overlap_genes <- read.table("unique_posterior_overlap_genes.txt", header = TRUE)
posterior_overlap_genes <- read.table("unique_anterior_overlap_genes.txt", header = TRUE)

unique_posterior_genes <- rename(unique_posterior_genes, ID = x)
unique_anterior_genes <- rename(unique_anterior_genes, ID = x)
unique_prostomium_genes <- rename(unique_prostomium_genes, ID = x)
unique_tail_genes <- rename(unique_tail_genes, ID = x)
unique_head_genes <- rename(unique_head_genes, ID = x)
unique_bodywall_genes <- rename(unique_bodywall_genes, ID = x)
anterior_overlap_genes <- rename(anterior_overlap_genes, ID = x)
posterior_overlap_genes <- rename(posterior_overlap_genes, ID = x)

unique_posterior_genes <- left_join(unique_posterior_genes, seven_stages, by = "ID")
unique_anterior_genes <- left_join(unique_anterior_genes, seven_stages, by = "ID")
unique_anterior_genes <- select(unique_anterior_genes, -x)
unique_prostomium_genes<- left_join(unique_prostomium_genes, seven_stages, by = "ID")
unique_prostomium_genes <- select(unique_prostomium_genes, -x)
unique_tail_genes <- left_join(unique_tail_genes, seven_stages, by = "ID")
unique_tail_genes <- select(unique_tail_genes, -x)
unique_head_genes <- left_join(unique_head_genes, seven_stages, by = "ID")
unique_head_genes <- select(unique_head_genes, -x)
unique_bodywall_genes <- left_join(unique_bodywall_genes, seven_stages, by = "ID")
unique_bodywall_genes <- select(unique_bodywall_genes, -x)
posterior_overlap_genes <- left_join(posterior_overlap_genes, seven_stages, by = "ID")
anterior_overlap_genes <- left_join(anterior_overlap_genes, seven_stages, by = "ID")

#exporting tables of the unique genes at each developmental stage

write.table(unique_posterior_genes, "posterior_genes_7stages.txt",sep="\t", quote = F)
write.table(unique_anterior_genes, "anterior_genes_7stages.txt",sep="\t", quote = F)
write.table(unique_prostomium_genes, "prostomium_genes_7stages.txt",sep="\t", quote = F)
write.table(unique_tail_genes, "tail_genes_7stages.txt",sep="\t", quote = F)
write.table(unique_head_genes, "head_genes_7stages.txt",sep="\t", quote = F)
write.table(unique_bodywall_genes, "bodywall_genes_7stages.txt",sep="\t", quote = F)
write.table(posterior_overlap_genes, "posterior_overlap_genes_7stages.txt", sep="\t", quote = F)
write.table(anterior_overlap_genes, "anterior_overlap_genes_7stages.txt", sep="\t", quote = F)

# For cluster 1
df_cluster1_raw <- unique_posterior_genes
df_cluster1_raw <- column_to_rownames(df_cluster1_raw, var = "ID")
mat_cluster1_raw <- as.matrix(df_cluster1_raw) # convert to matrix
mat_cluster1_norm <- t(scale(t(mat_cluster1_raw))) # normalise to z-score
df_cluster1_norm <- data.frame(mat_cluster1_norm) # convert back to data frame
df_cluster1_clean <- data.frame(cbind(time,t(df_cluster1_norm))) # merge with time

# Convert from short to long format
df_cluster1_clean_long <- gather(df_cluster1_clean, ID, expression, -time)
df_cluster1_clean_long$time <- as.numeric(df_cluster1_clean_long$time)
df_cluster1_clean_long$expression <- as.numeric(df_cluster1_clean_long$expression)

# Create a statistical summary dataframe containing the number of genes in each cluster (length, n), the mean 
summary_cluster1 <- data.frame(time=df_cluster1_clean$time, 
                               n=tapply(df_cluster1_clean_long$expression, 
                                        df_cluster1_clean_long$time, length), 
                               mean=tapply(df_cluster1_clean_long$expression, 
                                           df_cluster1_clean_long$time, mean, na.rm=TRUE))

# Plot data
cluster1_plot <- ggplot(summary_cluster1, aes(x=time, y=mean)) +
  geom_line(data = df_cluster1_clean_long, aes(x=time, y=expression, group=ID), 
            color="gray") +
  stat_smooth(colour = palette[1], # change 1 to match the number of the cluster, so that each cluster has its own colour of the continuous palette
              fill = palette[1], # change 1 to match the number of the cluster, so that each cluster has its own colour of the continuous palette
              size = 2, alpha = 0.5) +  # stat_smooth performs the loess regression by default unless a different method is specified
  geom_hline(yintercept = 0, linetype="dashed", color = "black") + # y-intercept at y=0 so that it's easy to see whether they are expressed less or more than average. Remove or comment this line if you don't want the line to show.
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(title = "Posterior genes", y = "Normalised gene expression (z-score)", x = "Development stage") + # Change "Cluster 1" accordingly to the cluster number
  scale_y_continuous(breaks = c(-2:2)) + # change the limits according to the range of your data, this worked for us, but there will be datasets in which you probably need to make them c(-5,5) for example, if the standard deviation of the expression values is higher
  scale_x_continuous(breaks = c(1:7), expand = c(0,0), limits = c(0.8,7.2),
                     labels = c("blastula", "gastrula","elongation", "early_larva", 
                                "mitraria_larva", "competent_larva", "juvenille")) + # change labels accordingly to the names of the stages you have, change 7 in breaks according to the number of stages you have
  theme(plot.title = element_text(color="black",size=14,hjust=0.5),
        axis.title.x = element_text(color="black",size=14),
        axis.title.y = element_text(color="black",size=14),
        axis.text.x = element_text(color="black",size=12,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=12,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm")) +
  ggsave("cluster_posterior_loess.pdf", width = 25, height = 35, units = "cm") # change <name-of-file> and size (height and width) according to your needs
cluster1_plot
#anterior cluster next
df_cluster2_raw <- unique_anterior_genes
df_cluster2_raw <- column_to_rownames(df_cluster2_raw, var = "ID")
mat_cluster2_raw <- as.matrix(df_cluster2_raw) # convert to matrix
mat_cluster2_norm <- t(scale(t(mat_cluster2_raw))) # normalise to z-score
df_cluster2_norm <- data.frame(mat_cluster2_norm) # convert back to data frame
df_cluster2_clean <- data.frame(cbind(time,t(df_cluster2_norm))) # merge with time

# Convert from short to long format
df_cluster2_clean_long <- gather(df_cluster2_clean, ID, expression, -time)

# Create a statistical summary dataframe containing the number of genes in each cluster (length, n), the mean 
summary_cluster2 <- data.frame(time=df_cluster2_clean$time, 
                               n=tapply(df_cluster2_clean_long$expression, 
                                        df_cluster2_clean_long$time, length), 
                               mean=tapply(df_cluster2_clean_long$expression, 
                                           df_cluster2_clean_long$time, mean, na.rm=TRUE))

# Plot data
cluster2_plot <- ggplot(summary_cluster2, aes(x=time, y=mean)) +
  geom_line(data = df_cluster2_clean_long, aes(x=time, y=expression, group=ID), 
            color="gray") +
  stat_smooth(colour = palette[1], # change 1 to match the number of the cluster, so that each cluster has its own colour of the continuous palette
              fill = palette[1], # change 1 to match the number of the cluster, so that each cluster has its own colour of the continuous palette
              size = 2, alpha = 0.5) +  # stat_smooth performs the loess regression by default unless a different method is specified
  geom_hline(yintercept = 0, linetype="dashed", color = "black") + # y-intercept at y=0 so that it's easy to see whether they are expressed less or more than average. Remove or comment this line if you don't want the line to show.
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(title = "Anterior genes", y = "Normalised gene expression (z-score)", x = "Developmental stage") + # Change "Cluster 1" accordingly to the cluster number
  scale_y_continuous(breaks = c(-2:2), limits = c(-3,4)) + # change the limits according to the range of your data, this worked for us, but there will be datasets in which you probably need to make them c(-5,5) for example, if the standard deviation of the expression values is higher
  scale_x_continuous(breaks = c(1:7), expand = c(0,0), limits = c(0.8,7.2),
                     labels = c("blastula", "gastrula","elongation", "early_larva", 
                                "mitraria_larva", "competent_larva", "juvenille")) + # change labels accordingly to the names of the stages you have, change 7 in breaks according to the number of stages you have
  theme(plot.title = element_text(color="black",size=14,hjust=0.5),
        axis.title.x = element_text(color="black",size=14),
        axis.title.y = element_text(color="black",size=14),
        axis.text.x = element_text(color="black",size=12,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=12,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm")) +
  ggsave("anterior_cluster_loess.pdf", width = 25, height = 35, units = "cm") # change <name-of-file> and size (height and width) according to your needs
cluster2_plot
#tail cluster
df_cluster3_raw <- unique_tail_genes
df_cluster3_raw <- column_to_rownames(df_cluster3_raw, var = "ID")
mat_cluster3_raw <- as.matrix(df_cluster3_raw) # convert to matrix
mat_cluster3_norm <- t(scale(t(mat_cluster3_raw))) # normalise to z-score
df_cluster3_norm <- data.frame(mat_cluster3_norm) # convert back to data frame
df_cluster3_clean <- data.frame(cbind(time,t(df_cluster3_norm))) # merge with time

# Convert from short to long format
df_cluster3_clean_long <- gather(df_cluster3_clean, ID, expression, -time)

# Create a statistical summary dataframe containing the number of genes in each cluster (length, n), the mean 
summary_cluster3 <- data.frame(time=df_cluster3_clean$time, 
                               n=tapply(df_cluster3_clean_long$expression, 
                                        df_cluster3_clean_long$time, length), 
                               mean=tapply(df_cluster3_clean_long$expression, 
                                           df_cluster3_clean_long$time, mean, na.rm = TRUE))

# Plot data
cluster3_plot <- ggplot(summary_cluster3, aes(x=time, y=mean)) +
  geom_line(data = df_cluster3_clean_long, aes(x=time, y=expression, group=ID), 
            color="gray") +
  stat_smooth(colour = palette[1], # change 1 to match the number of the cluster, so that each cluster has its own colour of the continuous palette
              fill = palette[1], # change 1 to match the number of the cluster, so that each cluster has its own colour of the continuous palette
              size = 2, alpha = 0.5) +  # stat_smooth performs the loess regression by default unless a different method is specified
  geom_hline(yintercept = 0, linetype="dashed", color = "black") + # y-intercept at y=0 so that it's easy to see whether they are expressed less or more than average. Remove or comment this line if you don't want the line to show.
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(title = "Tail genes", y = "Normalised gene expression (z-score)", x = "Developmental stage") + # Change "Cluster 1" accordingly to the cluster number
  scale_y_continuous(breaks = c(-2:2), limits = c(-3,4)) + # change the limits according to the range of your data, this worked for us, but there will be datasets in which you probably need to make them c(-5,5) for example, if the standard deviation of the expression values is higher
  scale_x_continuous(breaks = c(1:7), expand = c(0,0), limits = c(0.8,7.2),
                     labels = c("blastula", "gastrula","elongation", "early_larva", 
                                "mitraria_larva", "competent_larva", "juvenille")) + # change labels accordingly to the names of the stages you have, change 7 in breaks according to the number of stages you have
  theme(plot.title = element_text(color="black",size=14,hjust=0.5),
        axis.title.x = element_text(color="black",size=14),
        axis.title.y = element_text(color="black",size=14),
        axis.text.x = element_text(color="black",size=12,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=12,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm")) +
  ggsave("tail_cluster_loess.pdf", width = 25, height = 35, units = "cm") # change <name-of-file> and size (height and width) according to your needs

#unique bodywall genes
df_cluster4_raw <- unique_bodywall_genes
df_cluster4_raw <- column_to_rownames(df_cluster4_raw, var = "ID")
mat_cluster4_raw <- as.matrix(df_cluster4_raw) # convert to matrix
mat_cluster4_norm <- t(scale(t(mat_cluster4_raw))) # normalise to z-score
df_cluster4_norm <- data.frame(mat_cluster4_norm) # convert back to data frame
df_cluster4_clean <- data.frame(cbind(time,t(df_cluster4_norm))) # merge with time

# Convert from short to long format
df_cluster4_clean_long <- gather(df_cluster4_clean, ID, expression, -time)

# Create a statistical summary dataframe containing the number of genes in each cluster (length, n), the mean 
summary_cluster4 <- data.frame(time=df_cluster4_clean$time, 
                               n=tapply(df_cluster4_clean_long$expression, 
                                        df_cluster4_clean_long$time, length), 
                               mean=tapply(df_cluster4_clean_long$expression, 
                                           df_cluster4_clean_long$time, mean, na.rm = TRUE))

# Plot data
cluster4_plot <- ggplot(summary_cluster4, aes(x=time, y=mean)) +
  geom_line(data = df_cluster4_clean_long, aes(x=time, y=expression, group=ID), 
            color="gray") +
  stat_smooth(colour = palette[1], # change 1 to match the number of the cluster, so that each cluster has its own colour of the continuous palette
              fill = palette[1], # change 1 to match the number of the cluster, so that each cluster has its own colour of the continuous palette
              size = 2, alpha = 0.5) +  # stat_smooth performs the loess regression by default unless a different method is specified
  geom_hline(yintercept = 0, linetype="dashed", color = "black") + # y-intercept at y=0 so that it's easy to see whether they are expressed less or more than average. Remove or comment this line if you don't want the line to show.
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(title = "Bodywall genes", y = "Normalised gene expression (z-score)", x = "Developmental stage") + # Change "Cluster 1" accordingly to the cluster number
  scale_y_continuous(breaks = c(-2:2), limits = c(-3,4)) + # change the limits according to the range of your data, this worked for us, but there will be datasets in which you probably need to make them c(-5,5) for example, if the standard deviation of the expression values is higher
  scale_x_continuous(breaks = c(1:7), expand = c(0,0), limits = c(0.8,7.2),
                     labels = c("blastula", "gastrula","elongation", "early_larva", 
                                "mitraria_larva", "competent_larva", "juvenille")) + # change labels accordingly to the names of the stages you have, change 7 in breaks according to the number of stages you have
  theme(plot.title = element_text(color="black",size=14,hjust=0.5),
        axis.title.x = element_text(color="black",size=14),
        axis.title.y = element_text(color="black",size=14),
        axis.text.x = element_text(color="black",size=12,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=12,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm")) +
  ggsave("bodywall_cluster_loess.pdf", width = 25, height = 35, units = "cm") # change <name-of-file> and size (height and width) according to your needs

#unique head genes
df_cluster5_raw <- unique_head_genes
df_cluster5_raw <- column_to_rownames(df_cluster5_raw, var = "ID")
mat_cluster5_raw <- as.matrix(df_cluster5_raw) # convert to matrix
mat_cluster5_norm <- t(scale(t(mat_cluster5_raw))) # normalise to z-score
df_cluster5_norm <- data.frame(mat_cluster5_norm) # convert back to data frame
df_cluster5_clean <- data.frame(cbind(time,t(df_cluster5_norm))) # merge with time

# Convert from short to long format
df_cluster5_clean_long <- gather(df_cluster5_clean, ID, expression, -time)

# Create a statistical summary dataframe containing the number of genes in each cluster (length, n), the mean 
summary_cluster5 <- data.frame(time=df_cluster5_clean$time, 
                               n=tapply(df_cluster5_clean_long$expression, 
                                        df_cluster5_clean_long$time, length), 
                               mean=tapply(df_cluster5_clean_long$expression, 
                                           df_cluster5_clean_long$time, mean, na.rm = TRUE))

# Plot data
cluster5_plot <- ggplot(summary_cluster5, aes(x=time, y=mean)) +
  geom_line(data = df_cluster5_clean_long, aes(x=time, y=expression, group=ID), 
            color="gray") +
  stat_smooth(colour = palette[1], # change 1 to match the number of the cluster, so that each cluster has its own colour of the continuous palette
              fill = palette[1], # change 1 to match the number of the cluster, so that each cluster has its own colour of the continuous palette
              size = 2, alpha = 0.5) +  # stat_smooth performs the loess regression by default unless a different method is specified
  geom_hline(yintercept = 0, linetype="dashed", color = "black") + # y-intercept at y=0 so that it's easy to see whether they are expressed less or more than average. Remove or comment this line if you don't want the line to show.
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(title = "head genes", y = "Normalised gene expression (z-score)", x = "Developmental stage") + # Change "Cluster 1" accordingly to the cluster number
  scale_y_continuous(breaks = c(-2:2), limits = c(-3,4)) + # change the limits according to the range of your data, this worked for us, but there will be datasets in which you probably need to make them c(-5,5) for example, if the standard deviation of the expression values is higher
  scale_x_continuous(breaks = c(1:7), expand = c(0,0), limits = c(0.8,7.2),
                     labels = c("blastula", "gastrula","elongation", "early_larva", 
                                "mitraria_larva", "competent_larva", "juvenille")) + # change labels accordingly to the names of the stages you have, change 7 in breaks according to the number of stages you have
  theme(plot.title = element_text(color="black",size=14,hjust=0.5),
        axis.title.x = element_text(color="black",size=14),
        axis.title.y = element_text(color="black",size=14),
        axis.text.x = element_text(color="black",size=12,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=12,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm")) +
  ggsave("head_cluster_loess.pdf", width = 25, height = 35, units = "cm") # change <name-of-file> and size (height and width) according to your needs

#prostomium cluster
df_cluster6_raw <- unique_prostomium_genes
df_cluster6_raw <- column_to_rownames(df_cluster6_raw, var = "ID")
mat_cluster6_raw <- as.matrix(df_cluster6_raw) # convert to matrix
mat_cluster6_norm <- t(scale(t(mat_cluster6_raw))) # normalise to z-score
df_cluster6_norm <- data.frame(mat_cluster6_norm) # convert back to data frame
df_cluster6_clean <- data.frame(cbind(time,t(df_cluster6_norm))) # merge with time

# Convert from short to long format
df_cluster6_clean_long <- gather(df_cluster6_clean, ID, expression, -time)

# Create a statistical summary dataframe containing the number of genes in each cluster (length, n), the mean 
summary_cluster6 <- data.frame(time=df_cluster6_clean$time, 
                               n=tapply(df_cluster6_clean_long$expression, 
                                        df_cluster6_clean_long$time, length), 
                               mean=tapply(df_cluster6_clean_long$expression, 
                                           df_cluster6_clean_long$time, mean, na.rm=TRUE))

# Plot data
cluster6_plot <- ggplot(summary_cluster6, aes(x=time, y=mean)) +
  geom_line(data = df_cluster6_clean_long, aes(x=time, y=expression, group=ID), 
            color="gray") +
  stat_smooth(colour = palette[1], # change 1 to match the number of the cluster, so that each cluster has its own colour of the continuous palette
              fill = palette[1], # change 1 to match the number of the cluster, so that each cluster has its own colour of the continuous palette
              size = 2, alpha = 0.5) +  # stat_smooth performs the loess regression by default unless a different method is specified
  geom_hline(yintercept = 0, linetype="dashed", color = "black") + # y-intercept at y=0 so that it's easy to see whether they are expressed less or more than average. Remove or comment this line if you don't want the line to show.
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(title = "prostomium genes", y = "Normalised gene expression (z-score)", x = "Developmental stage") + # Change "Cluster 1" accordingly to the cluster number
  scale_y_continuous(breaks = c(-2:2), limits = c(-3,4)) + # change the limits according to the range of your data, this worked for us, but there will be datasets in which you probably need to make them c(-5,5) for example, if the standard deviation of the expression values is higher
  scale_x_continuous(breaks = c(1:7), expand = c(0,0), limits = c(0.8,7.2),
                     labels = c("blastula", "gastrula","elongation", "early_larva", 
                                "mitraria_larva", "competent_larva", "juvenille")) + # change labels accordingly to the names of the stages you have, change 7 in breaks according to the number of stages you have
  theme(plot.title = element_text(color="black",size=14,hjust=0.5),
        axis.title.x = element_text(color="black",size=14),
        axis.title.y = element_text(color="black",size=14),
        axis.text.x = element_text(color="black",size=12,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=12,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm")) +
  ggsave("prostomium_cluster_loess.pdf", width = 25, height = 35, units = "cm") # change <name-of-file> and size (height and width) according to your needs

allfigures <- ggarrange(cluster1_plot, cluster2_plot, cluster3_plot, cluster4_plot,
                        cluster5_plot, cluster6_plot,
                        ncol = 2, nrow = 3) +
  ggsave("tissue_analysis_6_loess_graphs.pdf", width = 35, height = 50, units = "cm") # change <name-of-file> and size (height and width) according to your needs 

ant_post <- ggarrange(cluster1_plot, cluster2_plot, ncol = 2, nrow = 1)+
  ggsave("anterior_posterior_loess_graphs.pdf", width = 50, height = 50, units = "cm") 

#boxplots and heatmaps


rescale_custom <- function(x) (x/(max(x)))

# Owenia fusiformis
raw <- read.table("anterior_overlap_genes_7stages.txt", header = TRUE)
rownames(raw) <- raw$ID

df_raw <- raw[,-c(1)]
df_norm <- t(apply(df_raw, 1, rescale_custom))
df_norm <- na.omit(df_norm)

heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color <- heatmap_color[50:100]

ComplexHeatmap::Heatmap(df_norm,
                        cluster_columns = FALSE,
                        cluster_rows = TRUE,
                        col = heatmap_color,
                        heatmap_legend_param = list(color_bar = "continuous"))

raw2 <- read.table("posterior_overlap_genes_7stages.txt", header = TRUE)
rownames(raw2) <- raw2$ID

df_raw2 <- raw2[,-c(1)]
df_norm2 <- t(apply(df_raw2, 1, rescale_custom))
df_norm2 <- na.omit(df_norm2)

heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color <- heatmap_color[50:100]

ComplexHeatmap::Heatmap(df_norm2,
                        cluster_columns = FALSE,
                        cluster_rows = TRUE,
                        col = heatmap_color,
                        heatmap_legend_param = list(color_bar = "continuous"))

#boxplots-violin plots

# Owenia fusiformis boxplot deseq2

stages <-  c("blastula", "gastrula","elongation",
             "early larva", "mitraria larva",
             "competent larva","juvenile")
df_DESeq2 <- read.table("combined_overlap_genes_7stages.txt", header = TRUE)
rownames(df_DESeq2) <- df_DESeq2$ID

df_copy <- df_DESeq2[,-c(1:2)]
df_clean <- data.frame(cbind(stages,t(df_copy)))
df_clean_long <- gather(df_clean, gene_ID, expression, -stages)
df_clean_long$family <- rep(df_DESeq2$Family, each = 7)
df_clean_long$expression <- as.numeric(df_clean_long$expression)
df_clean_long$family <- factor(df_clean_long$family, levels = c("Anterior", "Posterior"))
df_clean_long$stages <- factor(df_clean_long$stages, levels = c("blastula", "gastrula","elongation",
                                                                "early larva", "mitraria larva",
                                                                "competent larva","juvenile"))

ggplot(df_clean_long, aes(x=stages, y=expression, fill=family)) +
  geom_boxplot(outlier.shape = NA) +  
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "DESeq2 normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(0,1200)) +
  scale_y_continuous(breaks = c(seq(0,5000,500))) +
  theme(plot.title = element_text(color="black",size=14,hjust=0.5),
        axis.title.x = element_text(color="black",size=14),
        axis.title.y = element_text(color="black",size=14),
        axis.text.x = element_text(color="black",size=12,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=12,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"))

#owenia boxplots z score
stages <-  c("blastula", "gastrula","elongation",
             "early larva", "mitraria larva",
             "competent larva","juvenile")
df_DESeq2 <- read.table("combined_overlap_genes_7stages.txt", header = TRUE)
rownames(df_DESeq2) <- df_DESeq2$ID


df_copy <- df_DESeq2[,-c(1:2)]
df_norm <- data.frame(t(scale(t(data.matrix(df_copy)))))
df_clean <- data.frame(cbind(stages,t(df_norm)))
df_clean_long <- gather(df_clean, gene_ID, expression, -stages)
df_clean_long$family <- rep(df_DESeq2$Family, each = 7)
df_clean_long$expression <- as.numeric(df_clean_long$expression)
df_clean_long$family <- factor(df_clean_long$family, levels = c("Anterior", "Posterior", "Trunk"))
df_clean_long$stages <- factor(df_clean_long$stages, levels = c("blastula", "gastrula","elongation",
                                                                "early larva", "mitraria larva",
                                                                "competent larva","juvenile"))

ggplot(df_clean_long, aes(x=stages, y=expression, fill=family)) +
  geom_boxplot(outlier.shape = NA) +  
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "DESeq2 normalised expression", x = "developmental stage") +
  coord_cartesian(ylim = c(-3,3)) +
  scale_y_continuous(breaks = c(seq(-3,3,1))) +
  theme(plot.title = element_text(color="black",size=14,hjust=0.5),
        axis.title.x = element_text(color="black",size=14),
        axis.title.y = element_text(color="black",size=14),
        axis.text.x = element_text(color="black",size=12,angle=90,
                                   hjust=0.95,vjust=0.2,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.text.y = element_text(color="black",size=12,
                                   margin=margin(0.3,0.3,0.3,0.3,"cm")),
        axis.line = element_line(colour="black", size=0.5, linetype = "solid"),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"))

#paired t-tests of the z scores, bonferroni corrected i.e.multiply the p values produced by 7
blastula <- pairwise.t.test(df_clean_long$expression[df_clean_long$stages == 'blastula'], df_clean_long$family[df_clean_long$stages == 'blastula'])
gastrula <- pairwise.t.test(df_clean_long$expression[df_clean_long$stages == 'gastrula'], df_clean_long$family[df_clean_long$stages == 'gastrula'])
elongation <- pairwise.t.test(df_clean_long$expression[df_clean_long$stages == 'elongation'], df_clean_long$family[df_clean_long$stages == 'elongation'])
early <- pairwise.t.test(df_clean_long$expression[df_clean_long$stages == 'early larva'], df_clean_long$family[df_clean_long$stages == 'early larva'])
mitraria <- pairwise.t.test(df_clean_long$expression[df_clean_long$stages == 'mitraria larva'], df_clean_long$family[df_clean_long$stages == 'mitraria larva'])
competent <- pairwise.t.test(df_clean_long$expression[df_clean_long$stages == 'competent larva'], df_clean_long$family[df_clean_long$stages == 'competent larva'])
juvenile <- pairwise.t.test(df_clean_long$expression[df_clean_long$stages == 'juvenile'], df_clean_long$family[df_clean_long$stages == 'juvenile'])
pvalues <- c(blastula$p.value, gastrula$p.value, elongation$p.value, early$p.value, mitraria$p.value, 
             competent$p.value, juvenile$p.value)

#to bonferroni correct, the pvalues need to be multiplied by 7
pvalues
final_pvalues <- pvalues*7
final_pvalues
dev_stage <- c("blastula", "gastrula", "elongation", "early larva", "mitraria larva",
               "competent larva", "juvenile")
pvalues_stages <- data.frame(dev_stage, pvalues)
pvalues_stages <- rename(pvalues_stages, bonferroni_pvalue = final_pvalues)
write.table(pvalues_stages, "pvalues_stages.txt", sep="\t", quote = F)

