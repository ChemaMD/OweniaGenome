#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
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


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Recalculate topological overlap if needed
TOM <- TOMsimilarityFromExpr(datExpr, power = 8);
# Read in the annotation file
annot <- read.csv(file = "00-Capitella_annotation_v050321_TrinoPanther.csv");
# Select modules (in this case, the whole genome)
colors <- read.csv(file = "moduleColors.csv");
colors <- unique(colors[2])
modules <- t(colors)
# Select module IDs
Gene_IDs <- names(datExpr)
inModule <- is.finite(match(moduleColors, modules));
modProbes <- Gene_IDs[inModule];
modGenes <- annot$Transcript_id[match(modProbes, annot$Transcript_id)];
# Select the corresponding Topological Overlap
modTOM <- TOM[inModule, inModule];
dimnames(modTOM) <- list(modProbes, modProbes)
write.table(modules, paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""))
# Export the network into edge and node list files Cytoscape can read
cyt <- exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.25,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])
