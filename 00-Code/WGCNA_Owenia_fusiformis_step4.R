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
# Import the file
names <- read.csv("04-module_by_stage.csv")
rownames(names) <- names$X
names <- names[-c(1)]

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Load the heatmap package and the colour package
library(ComplexHeatmap)
library(RColorBrewer)
# Create vector with colours
colours <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
# Plot heatmap
ComplexHeatmap::Heatmap(t(names),
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        col = colours)