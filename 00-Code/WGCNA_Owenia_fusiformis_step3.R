#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "~/Dropbox/02-OweniaGenome/02-Figures/01-DataForFigures/00-RNAseq/WCNA";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression data saved in the first part
lnames = load(file = "step1.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "step2.RData");
lnames

names <- read.csv("moduleLabels.csv")


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


enableWGCNAThreads()
# Recalculate topological overlap if needed
TOM <- TOMsimilarityFromExpr(datExpr, power = 16);
# Read in the annotation file
annot <- read.csv(file = "Owenia_annotation_v250920.csv");
# Select modules (in this case, the whole genome)
colors <- read.csv(file = "module_colors.csv");
colors <- colors[2]
modules <- t(colors)
# Select module IDs
Gene_IDs <- names(datExpr)
inModule <- is.finite(match(moduleColors, modules));
modProbes <- Gene_IDs[inModule];
modGenes <- annot$gene_symbol[match(modProbes, annot$transcript_id)];
# Select the corresponding Topological Overlap
modTOM <- TOM[inModule, inModule];
dimnames(modTOM) <- list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt <- exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.25,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


# Export MEs and moduleLabels for further analysis
write.csv(MEs,"modules_by_stage.csv")
write.csv(moduleLabels,"moduleLabels.csv")
