library(topGO)
library(ggplot2)
library(ggpubr)
library(cowplot)

# Import gene universe: whole (GO-annotated) genome
geneID2GO <- readMappings(file = "Owenia_GO_universe.txt") ### 21108 transcripts have GO annotation
geneUniverse <- names(geneID2GO)

# Import and transform gene lists from modules from WCNA

modules <- read.csv("moduleLabels.csv")

module1 <- modules[modules$x == "1",]
module1 <- as.character(module1$X) 
module1genelist <- factor(as.integer(geneUniverse %in% module1))
names(module1genelist) <- geneUniverse

module2 <- modules[modules$x == "2",]
module2 <- as.character(module2$X) 
module2genelist <- factor(as.integer(geneUniverse %in% module2))
names(module2genelist) <- geneUniverse

module3 <- modules[modules$x == "3",]
module3 <- as.character(module3$X) 
module3genelist <- factor(as.integer(geneUniverse %in% module3))
names(module3genelist) <- geneUniverse

module4 <- modules[modules$x == "4",]
module4 <- as.character(module4$X) 
module4genelist <- factor(as.integer(geneUniverse %in% module4))
names(module4genelist) <- geneUniverse

module5 <- modules[modules$x == "5",]
module5 <- as.character(module5$X) 
module5genelist <- factor(as.integer(geneUniverse %in% module5))
names(module5genelist) <- geneUniverse

module6 <- modules[modules$x == "6",]
module6 <- as.character(module6$X) 
module6genelist <- factor(as.integer(geneUniverse %in% module6))
names(module6genelist) <- geneUniverse

module7 <- modules[modules$x == "7",]
module7 <- as.character(module7$X) 
module7genelist <- factor(as.integer(geneUniverse %in% module7))
names(module7genelist) <- geneUniverse

module8 <- modules[modules$x == "8",]
module8 <- as.character(module8$X) 
module8genelist <- factor(as.integer(geneUniverse %in% module8))
names(module8genelist) <- geneUniverse

module9 <- modules[modules$x == "9",]
module9 <- as.character(module9$X) 
module9genelist <- factor(as.integer(geneUniverse %in% module9))
names(module9genelist) <- geneUniverse

module10 <- modules[modules$x == "10",]
module10 <- as.character(module10$X) 
module10genelist <- factor(as.integer(geneUniverse %in% module10))
names(module10genelist) <- geneUniverse

module11 <- modules[modules$x == "11",]
module11 <- as.character(module11$X) 
module11genelist <- factor(as.integer(geneUniverse %in% module11))
names(module11genelist) <- geneUniverse

module12 <- modules[modules$x == "12",]
module12 <- as.character(module12$X) 
module12genelist <- factor(as.integer(geneUniverse %in% module12))
names(module12genelist) <- geneUniverse

module13 <- modules[modules$x == "13",]
module13 <- as.character(module13$X) 
module13genelist <- factor(as.integer(geneUniverse %in% module13))
names(module13genelist) <- geneUniverse

module14 <- modules[modules$x == "14",]
module14 <- as.character(module14$X) 
module14genelist <- factor(as.integer(geneUniverse %in% module14))
names(module14genelist) <- geneUniverse

module15 <- modules[modules$x == "15",]
module15 <- as.character(module15$X) 
module15genelist <- factor(as.integer(geneUniverse %in% module15))
names(module15genelist) <- geneUniverse

module16 <- modules[modules$x == "0",]
module16 <- as.character(module16$X) 
module16genelist <- factor(as.integer(geneUniverse %in% module16))
names(module16genelist) <- geneUniverse



# fisher testing of GO term enrichment for Biological Process (BP)
#module 1
module1_GOdata_BP <- new("topGOdata", description="module1", 
                          ontology="BP", allGenes=module1genelist,  
                          annot = annFUN.gene2GO, gene2GO = geneID2GO)
module1_resultFisher_BP <- runTest(module1_GOdata_BP, 
                                    alGOrithm="classic", statistic="fisher")
module1_BP <- GenTable(module1_GOdata_BP, classicFisher = module1_resultFisher_BP, 
                        orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)


#module 2
module2_GOdata_BP <- new("topGOdata", description="module2", 
                          ontology="BP", allGenes=module2genelist,  
                          annot = annFUN.gene2GO, gene2GO = geneID2GO)
module2_resultFisher_BP <- runTest(module2_GOdata_BP, 
                                    alGOrithm="classic", statistic="fisher")
module2_BP <- GenTable(module2_GOdata_BP, classicFisher = module2_resultFisher_BP, 
                        orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)

#module 3
module3_GOdata_BP <- new("topGOdata", description="module3", 
                          ontology="BP", allGenes=module3genelist,  
                          annot = annFUN.gene2GO, gene2GO = geneID2GO)
module3_resultFisher_BP <- runTest(module3_GOdata_BP, 
                                    alGOrithm="classic", statistic="fisher")
module3_BP <- GenTable(module3_GOdata_BP, classicFisher = module3_resultFisher_BP, 
                        orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)

#module 4
module4_GOdata_BP <- new("topGOdata", description="module4", 
                          ontology="BP", allGenes=module4genelist,  
                          annot = annFUN.gene2GO, gene2GO = geneID2GO)
module4_resultFisher_BP <- runTest(module4_GOdata_BP, 
                                    alGOrithm="classic", statistic="fisher")
module4_BP <- GenTable(module4_GOdata_BP, classicFisher = module4_resultFisher_BP, 
                        orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)

#module 5
module5_GOdata_BP <- new("topGOdata", description="module5", 
                          ontology="BP", allGenes=module5genelist,  
                          annot = annFUN.gene2GO, gene2GO = geneID2GO)
module5_resultFisher_BP <- runTest(module5_GOdata_BP, 
                                    alGOrithm="classic", statistic="fisher")
module5_BP <- GenTable(module5_GOdata_BP, classicFisher = module5_resultFisher_BP, 
                        orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)

#module 6
module6_GOdata_BP <- new("topGOdata", description="module6", 
                          ontology="BP", allGenes=module6genelist,  
                          annot = annFUN.gene2GO, gene2GO = geneID2GO)
module6_resultFisher_BP <- runTest(module6_GOdata_BP, 
                                    alGOrithm="classic", statistic="fisher")
module6_BP <- GenTable(module6_GOdata_BP, classicFisher = module6_resultFisher_BP, 
                        orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)

#module 7
module7_GOdata_BP <- new("topGOdata", description="module7", 
                          ontology="BP", allGenes=module7genelist,  
                          annot = annFUN.gene2GO, gene2GO = geneID2GO)
module7_resultFisher_BP <- runTest(module7_GOdata_BP, 
                                    alGOrithm="classic", statistic="fisher")
module7_BP <- GenTable(module7_GOdata_BP, classicFisher = module7_resultFisher_BP, 
                        orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)

#module 8
module8_GOdata_BP <- new("topGOdata", description="module8", 
                         ontology="BP", allGenes=module8genelist,  
                         annot = annFUN.gene2GO, gene2GO = geneID2GO)
module8_resultFisher_BP <- runTest(module8_GOdata_BP, 
                                   alGOrithm="classic", statistic="fisher")
module8_BP <- GenTable(module8_GOdata_BP, classicFisher = module8_resultFisher_BP, 
                       orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)

#module 9
module9_GOdata_BP <- new("topGOdata", description="module9", 
                         ontology="BP", allGenes=module9genelist,  
                         annot = annFUN.gene2GO, gene2GO = geneID2GO)
module9_resultFisher_BP <- runTest(module9_GOdata_BP, 
                                   alGOrithm="classic", statistic="fisher")
module9_BP <- GenTable(module9_GOdata_BP, classicFisher = module9_resultFisher_BP, 
                       orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)

#module 10
module10_GOdata_BP <- new("topGOdata", description="module10", 
                          ontology="BP", allGenes=module10genelist,  
                          annot = annFUN.gene2GO, gene2GO = geneID2GO)
module10_resultFisher_BP <- runTest(module10_GOdata_BP, 
                                    alGOrithm="classic", statistic="fisher")
module10_BP <- GenTable(module10_GOdata_BP, classicFisher = module10_resultFisher_BP, 
                        orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)

#module 11
module11_GOdata_BP <- new("topGOdata", description="module11", 
                          ontology="BP", allGenes=module11genelist,  
                          annot = annFUN.gene2GO, gene2GO = geneID2GO)
module11_resultFisher_BP <- runTest(module11_GOdata_BP, 
                                    alGOrithm="classic", statistic="fisher")
module11_BP <- GenTable(module11_GOdata_BP, classicFisher = module11_resultFisher_BP, 
                        orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)

#module 12
module12_GOdata_BP <- new("topGOdata", description="module12", 
                          ontology="BP", allGenes=module12genelist,  
                          annot = annFUN.gene2GO, gene2GO = geneID2GO)
module12_resultFisher_BP <- runTest(module12_GOdata_BP, 
                                    alGOrithm="classic", statistic="fisher")
module12_BP <- GenTable(module12_GOdata_BP, classicFisher = module12_resultFisher_BP, 
                        orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)

#module 13
module13_GOdata_BP <- new("topGOdata", description="module13", 
                          ontology="BP", allGenes=module13genelist,  
                          annot = annFUN.gene2GO, gene2GO = geneID2GO)
module13_resultFisher_BP <- runTest(module13_GOdata_BP, 
                                    alGOrithm="classic", statistic="fisher")
module13_BP <- GenTable(module13_GOdata_BP, classicFisher = module13_resultFisher_BP, 
                        orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)

#module 14
module14_GOdata_BP <- new("topGOdata", description="module14", 
                          ontology="BP", allGenes=module14genelist,  
                          annot = annFUN.gene2GO, gene2GO = geneID2GO)
module14_resultFisher_BP <- runTest(module14_GOdata_BP, 
                                    alGOrithm="classic", statistic="fisher")
module14_BP <- GenTable(module14_GOdata_BP, classicFisher = module14_resultFisher_BP, 
                        orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)

#module 15
module15_GOdata_BP <- new("topGOdata", description="module15", 
                          ontology="BP", allGenes=module15genelist,  
                          annot = annFUN.gene2GO, gene2GO = geneID2GO)
module15_resultFisher_BP <- runTest(module15_GOdata_BP, 
                                    alGOrithm="classic", statistic="fisher")
module15_BP <- GenTable(module15_GOdata_BP, classicFisher = module15_resultFisher_BP, 
                        orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 15)

## export results
write.csv(module1_BP,"module1_GO.txt")
write.csv(module2_BP,"module2_GO.txt")
write.csv(module3_BP,"module3_GO.txt")
write.csv(module4_BP,"module4_GO.txt")
write.csv(module5_BP,"module5_GO.txt")
write.csv(module6_BP,"module6_GO.txt")
write.csv(module7_BP,"module7_GO.txt")
write.csv(module8_BP,"module8_GO.txt")
write.csv(module9_BP,"module9_GO.txt")
write.csv(module10_BP,"module10_GO.txt")
write.csv(module11_BP,"module11_GO.txt")
write.csv(module12_BP,"module12_GO.txt")
write.csv(module13_BP,"module13_GO.txt")
write.csv(module14_BP,"module14_GO.txt")
write.csv(module15_BP,"module15_GO.txt")

###plot results
#just change sequentially the number in module1_BP and use ggarrange at the end
#feeding the names given to the ggplots 

GOEnrichment <- module7_BP
GOEnrichment$classicFisher <- as.numeric(GOEnrichment$classicFisher)
GOEnrichment <- GOEnrichment[,c("GO.ID","Term","classicFisher")]
GOEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", GOEnrichment$Term)
GOEnrichment$Term <- gsub("\\.\\.\\.$", "", GOEnrichment$Term)
GOEnrichment$Term <- paste(GOEnrichment$GO.ID, GOEnrichment$Term, sep=", ")
GOEnrichment$Term <- factor(GOEnrichment$Term, levels=rev(GOEnrichment$Term))

module_7_plot <- ggplot(GOEnrichment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p-value)") +
  ggtitle("module 7") +
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

module_7_plot

plot_grid(module_1_plot + rremove("y.title") + rremove("x.title"), 
          module_2_plot + rremove("y.title") + rremove("x.title"),
          module_3_plot + rremove("y.title") + rremove("x.title"),
          module_4_plot + rremove("y.title") + rremove("x.title"),
          ncol = 1, align="v")

plot_grid(module_5_plot + rremove("y.title") + rremove("x.title"), 
          module_6_plot + rremove("y.title") + rremove("x.title"),
          module_7_plot + rremove("y.title") + rremove("x.title"),
          module_7_plot + rremove("y.title") + rremove("x.title"),
          ncol = 1, align="v")