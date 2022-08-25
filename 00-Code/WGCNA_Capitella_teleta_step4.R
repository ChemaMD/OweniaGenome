#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "~/Dropbox/02-OweniaGenome/06-Revision/00-DATA/06-Capitella_WGCNA";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# Import the file
names <- read.csv("modules_by_stage.csv")
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
