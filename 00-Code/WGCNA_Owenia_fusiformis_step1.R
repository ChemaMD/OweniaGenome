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
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the dataset
raw_data <- read.csv("normalized_counts_average.csv");
rownames(raw_data) <- t(raw_data[2])
raw_data <- raw_data[-c(1:9)]
colnames(raw_data) <- c("blastulae","gastrulae","elongation","early_larvae",
                        "mature_larvae","competent_larvae","juvenile")

# Take a quick look at what is in the data set:
dim(raw_data);
names(raw_data);


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


datExpr0 = as.data.frame(t(raw_data));


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


gsg <- goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


sampleTree <- hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================

# NOTE: Because there are no outlier samples and they agree in biological terms, 
# in the sense that they cluster according to developmental stage, we will keep
# all the samples (quite logical, but it's a step of the pipeline)

# Plot a line to show the cut
abline(h = 900000, col = "red");
# Determine cluster under the line
clust <- cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
# clust 0 contains the samples we want to keep.
keepSamples <- (clust==0)
datExpr <- datExpr0[keepSamples, ]
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


save(datExpr, file = "step1.RData")
