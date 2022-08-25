library(tidyverse)
library(dplyr)
library(ggplot2)
library(data.table)
library(gplots)
library(pheatmap)
library(dendextend)
library(factoextra)
library(ComplexHeatmap)
library(RColorBrewer)
library(NbClust)
library(scales)


raw_data <- read.delim('/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/12-Hox_expression_JSD_species/02-Blan_longest_isoform_Hox_DESeq2_average.txt')
raw_data
#.........................................Subsetting...............................................

#extract all geneIDs
names_row <- raw_data[,1]
names_row
#extract wanted gene ids (Chema excel) and create list
CT_IDs <- c('BL12289_cuf1',
            'BL01409_evm0',
            'EU921831.1',
            'EU921832.1',
            'BL11265_evm0',
            'BL01497_evm0',
            'BL02690_evm0',
            'BL01142_evm0',
            'EU921834.1',
            'BL02721_evm0',
            'BL22794_evm0',
            'BL01764_evm0',
            'BL11259_evm0',
            'BL11262_cuf5',
            'BL06042_cuf1')



CT_Names <- c('Hox1',
              'Hox2',
              'Hox3',
              'Hox4',
              'Hox5',
              'Hox6',
              'Hox7',
              'Hox8',
              'Hox9',
              'Hox10',
              'Hox11',
              'Hox12',
              'Hox13',
              'Hox14',
              'Hox15')


#Create emtpy list to add index of wanted rows
CT_index <- c()

for (i in CT_IDs) {
  output <- grep(i, names_row)
  CT_index[i] <- output
}

CT_index

#extract any colums not wanted

CTsubset <- raw_data[CT_index, ] 
CTsubset

#.....................................................................................................

df_curated <- data.matrix(CTsubset[,2:ncol(CTsubset)])
df_curated

rownames(df_curated) <- CT_Names


# 0 to 1 relative expression
rescale_custom <- function(x) (x/(max(x)))

df_normalised <- t(apply(df_curated, 1, rescale_custom))
df_normalised <- na.omit(df_normalised)
df_normalised

df_sorted <- df_normalised[match(CT_Names, rownames(df_normalised)),]
df_sorted

#1/2 rdbu (substitute of reds for 0 to 1)
heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color <- heatmap_color[50:100]

ComplexHeatmap::Heatmap(df_sorted,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        column_title = "Blan_Heatmap",
                        col = heatmap_color,
                        heatmap_legend_param = list(color_bar = "continuous"))

                       # column_labels = c('oocyte', 'fertilized embryo', 'polar body',	
                       #'2 cells', '4 cells', '8 cells', '16 cells',	'32 cells',	
                       #'blastula',	'emerged cilia', 'early trochophore',	'middle trochophore',
                       #'late trochophore',	'early segmentation', 'SL_mean',	'WL_mean'))library(tidyverse)
library(dplyr)
library(ggplot2)
library(data.table)
library(gplots)
library(pheatmap)
library(dendextend)
library(factoextra)
library(ComplexHeatmap)
library(RColorBrewer)
library(NbClust)
library(scales)


raw_data <- read.delim('/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/12-Hox_expression_JSD_species/02-Cele_longest_isoform_Hox_DESeq2_average.txt')
raw_data
#.........................................Subsetting...............................................

#extract all geneIDs
names_row <- raw_data[,1]
names_row
#extract wanted gene ids (Chema excel) and create list
CE_IDs <- c('R13A5.5.2', 'C07H6.7.2', 'C08C3.3.1', 'C08C3.1b.1', 'Y75B8A.1.1', 'Y75B8A.2a.1')
CE_Names <- c('ceh-13', 'lin-39', 'mab-5', 'egl-5', 'php-3', 'nob-1')

#Create emtpy list to add index of wanted rows
CE_index <- c()

for (i in CE_IDs) {
  output <- grep(i, names_row)
  CE_index[i] <- output
}

CE_index

#extract genes in hox_genesIDs

CEsubset <- raw_data[CE_index, ] 
CEsubset

#.....................................................................................................

df_curated <- data.matrix(CEsubset[,2:ncol(CEsubset)])
df_curated

rownames(df_curated) <- CE_Names


# 0 to 1 relative expression
rescale_custom <- function(x) (x/(max(x)))

df_normalised <- t(apply(df_curated, 1, rescale_custom))
df_normalised <- na.omit(df_normalised)
df_normalised

df_sorted <- df_normalised[match(CE_Names, rownames(df_normalised)),]
df_sorted

#1/2 rdbu (substitute of reds for 0 to 1)
heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color <- heatmap_color[50:100]

ComplexHeatmap::Heatmap(df_sorted,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        column_title = "Cele_Heatmap",
                        col = heatmap_color,
                        heatmap_legend_param = list(color_bar = "continuous"))

                       # column_labels = c('oocyte', 'fertilized embryo', 'polar body',	
                       #'2 cells', '4 cells', '8 cells', '16 cells',	'32 cells',	
                       #'blastula',	'emerged cilia', 'early trochophore',	'middle trochophore',
                       #'late trochophore',	'early segmentation', 'SL_mean',	'WL_mean'))library(tidyverse)
library(dplyr)
library(ggplot2)
library(data.table)
library(gplots)
library(pheatmap)
library(dendextend)
library(factoextra)
library(ComplexHeatmap)
library(RColorBrewer)
library(NbClust)
library(scales)


raw_data <- read.delim('/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/12-Hox_expression_JSD_species/02-Cgig_longest_isoform_Hox_DESeq2_average.txt')
raw_data
#.........................................Subsetting...............................................

#extract all geneIDs
names_row <- raw_data[,1]
names_row

#extract wanted gene ids (Chema excel) and create list
CT_IDs <-c('EKC32705',
'EKC32708',
'EKC32709',
'EKC32713',
'EKC32714',
'EKC41105',
'EKC41102',
'EKC39601',
'EKC29602',
'EKC29599')


CT_Names <- c('Hox1',
'Hox2',
'Hox3',
'Hox4',
'Hox5',
'Lox5',
'Lox4',
'Lox2',
'Post2',
'Post1')


#Create emtpy list to add index of wanted rows
CT_index <- c()

for (i in CT_IDs) {
  output <- grep(i, names_row)
  CT_index[i] <- output
}

CT_index

#extract genes in hox_genesIDs

CTsubset <- raw_data[CT_index, ] 
CTsubset

#.....................................................................................................

df_curated <- data.matrix(CTsubset[,2:ncol(CTsubset)])
df_curated

rownames(df_curated) <- CT_Names


# 0 to 1 relative expression
rescale_custom <- function(x) (x/(max(x)))

df_normalised <- t(apply(df_curated, 1, rescale_custom))
df_normalised <- na.omit(df_normalised)
df_normalised

df_sorted <- df_normalised[match(CT_Names, rownames(df_normalised)),]
df_sorted

#1/2 rdbu (substitute of reds for 0 to 1)
heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color <- heatmap_color[50:100]

ComplexHeatmap::Heatmap(df_sorted,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        column_title = "Cgig_Heatmap",
                        col = heatmap_color,
                        heatmap_legend_param = list(color_bar = "continuous"))

                       # column_labels = c('oocyte', 'fertilized embryo', 'polar body',	
                       #'2 cells', '4 cells', '8 cells', '16 cells',	'32 cells',	
                       #'blastula',	'emerged cilia', 'early trochophore',	'middle trochophore',
                       #'late trochophore',	'early segmentation', 'SL_mean',	'WL_mean'))library(tidyverse)
library(dplyr)
library(ggplot2)
library(data.table)
library(gplots)
library(pheatmap)
library(dendextend)
library(factoextra)
library(ComplexHeatmap)
library(RColorBrewer)
library(NbClust)
library(scales)


raw_data <- read.delim('/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/12-Hox_expression_JSD_species/02-Chem_longest_isoform_Hox_DESeq2_average.txt')
raw_data
#.........................................Subsetting...............................................

#extract all geneIDs
names_row <- raw_data[,1]
names_row
#extract wanted gene ids (Chema excel) and create list
CT_IDs <- c('CLYHEMT022475.1',
'CLYHEMT014275.1',
'CLYHEMT022840.1',
'CLYHEMT023651.1')



CT_Names <- c('Hox1 CLHE',
'Hox9-14A CLHE',
'Hox9-14B CLHE',
'Hox9-14C CLHE')

#Create emtpy list to add index of wanted rows
CT_index <- c()

for (i in CT_IDs) {
  output <- grep(i, names_row)
  CT_index[i] <- output
}

CT_index

#extract genes in hox_genesIDs

CTsubset <- raw_data[CT_index, ] 
CTsubset

#.....................................................................................................

df_curated <- data.matrix(CTsubset[,2:ncol(CTsubset)])
df_curated

rownames(df_curated) <- CT_Names


# 0 to 1 relative expression
rescale_custom <- function(x) (x/(max(x)))

df_normalised <- t(apply(df_curated, 1, rescale_custom))
df_normalised <- na.omit(df_normalised)
df_normalised

df_sorted <- df_normalised[match(CT_Names, rownames(df_normalised)),]
df_sorted

#1/2 rdbu (substitute of reds for 0 to 1)
heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color <- heatmap_color[50:100]

ComplexHeatmap::Heatmap(df_sorted,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        column_title = "Chem_Heatmap",
                        col = heatmap_color,
                        heatmap_legend_param = list(color_bar = "continuous"))

                       # column_labels = c('oocyte', 'fertilized embryo', 'polar body',	
                       #'2 cells', '4 cells', '8 cells', '16 cells',	'32 cells',	
                       #'blastula',	'emerged cilia', 'early trochophore',	'middle trochophore',
                       #'late trochophore',	'early segmentation', 'SL_mean',	'WL_mean'))library(tidyverse)
library(dplyr)
library(ggplot2)
library(data.table)
library(gplots)
library(pheatmap)
library(dendextend)
library(factoextra)
library(ComplexHeatmap)
library(RColorBrewer)
library(NbClust)
library(scales)


raw_data <- read.delim('/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/12-Hox_expression_JSD_species/02-Ctel_longest_isoform_Hox_DESeq2_average.txt')
raw_data
#.........................................Subsetting...............................................

#extract all geneIDs
names_row <- raw_data[,1]
names_row
#extract wanted gene ids (Chema excel) and create list
CT_IDs <- c('CTELG00000005701.1',
'CTELG00000005702.1',
'CTELG00000005703.1',
'CTELG00000005705.1',
'CTELG00000005706.1',
'CTELG00000005707.1',
'CTELG00000005708.1',
'CTELG00000005709.1',
'CTELG00000014784.1',
'CTELG00000014785.1',
'CTELG00000037722.1_post1')
#post1 not found in grep search

grep('ABY67961', names_row)
CT_Names <- c('lab',
'pb',
'Hox3',
'Dfd',
'Scr',
'Lox5',
'Antp',
'Lox4',
'Lox2',
'Post2',
'Post1')

#Create emtpy list to add index of wanted rows
CT_index <- c()

for (i in CT_IDs) {
  output <- grep(i, names_row)
  CT_index[i] <- output
}

CT_index

#extract genes in hox_genesIDs

CTsubset <- raw_data[CT_index,] 
CTsubset

#.....................................................................................................

df_curated <- data.matrix(CTsubset[,2:ncol(CTsubset)])
df_curated

rownames(df_curated) <- CT_Names


# 0 to 1 relative expression
rescale_custom <- function(x) (x/(max(x)))

df_normalised <- t(apply(df_curated, 1, rescale_custom))
df_normalised <- na.omit(df_normalised)
df_normalised

df_sorted <- df_normalised[match(CT_Names, rownames(df_normalised)),]
df_sorted

#1/2 rdbu (substitute of reds for 0 to 1)
heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color <- heatmap_color[50:100]

ComplexHeatmap::Heatmap(df_sorted,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        column_title = "Ctel_Heatmap",
                        col = heatmap_color,
                        heatmap_legend_param = list(color_bar = "continuous"))

                       # column_labels = c('oocyte', 'fertilized embryo', 'polar body',	
                       #'2 cells', '4 cells', '8 cells', '16 cells',	'32 cells',	
                       #'blastula',	'emerged cilia', 'early trochophore',	'middle trochophore',
                       #'late trochophore',	'early segmentation', 'SL_mean',	'WL_mean'))library(tidyverse)
library(dplyr)
library(ggplot2)
library(data.table)
library(gplots)
library(pheatmap)
library(dendextend)
library(factoextra)
library(ComplexHeatmap)
library(RColorBrewer)
library(NbClust)
library(scales)


raw_data <- read.delim('/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/12-Hox_expression_JSD_species/02-Dmel_longest_isoform_Hox_DESeq2_average.txt')
raw_data
#.........................................Subsetting...............................................

#extract all geneIDs
names_row <- raw_data[,1]
names_row
#extract wanted gene ids (Chema excel) and create list
CT_IDs <- c('FBtr0081696',
            'FBtr0089276',
            'FBtr0081670',
            'FBtr0081671',
            'FBtr0081621',
            'U10506.1',
            'FBtr0081625',
            'FBtr0081652',
            'FBtr0083347',
            'FBtr0083388',
            'FBtr0415463')

grep ('U10506.1', names_row)

CT_Names <- c('lab',
              'pb',
              'zen',
              'zen2',
              'Dfd',
              'scr',
              'ftz',
              'antp',
              'ubx',
              'abd-a',
              'abd-b')


#Create emtpy list to add index of wanted rows
CT_index <- c()

for (i in CT_IDs) {
  output <- grep(i, names_row)
  CT_index[i] <- output
}

CT_index

#extract genes in hox_genesIDs

CTsubset <- raw_data[CT_index, ] 
CTsubset

#.....................................................................................................

df_curated <- data.matrix(CTsubset[,2:ncol(CTsubset)])
df_curated

rownames(df_curated) <- CT_Names


# 0 to 1 relative expression
rescale_custom <- function(x) (x/(max(x)))

df_normalised <- t(apply(df_curated, 1, rescale_custom))
df_normalised <- na.omit(df_normalised)
df_normalised

df_sorted <- df_normalised[match(CT_Names, rownames(df_normalised)),]
df_sorted

#1/2 rdbu (substitute of reds for 0 to 1)
heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color <- heatmap_color[50:100]

ComplexHeatmap::Heatmap(df_sorted,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        column_title = "Dmel_Heatmap",
                        col = heatmap_color,
                        heatmap_legend_param = list(color_bar = "continuous"))

                       # column_labels = c('oocyte', 'fertilized embryo', 'polar body',	
                       #'2 cells', '4 cells', '8 cells', '16 cells',	'32 cells',	
                       #'blastula',	'emerged cilia', 'early trochophore',	'middle trochophore',
                       #'late trochophore',	'early segmentation', 'SL_mean',	'WL_mean'))library(tidyverse)
library(dplyr)
library(ggplot2)
library(data.table)
library(gplots)
library(pheatmap)
library(dendextend)
library(factoextra)
library(ComplexHeatmap)
library(RColorBrewer)
library(NbClust)
library(scales)


raw_data <- read.delim('/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/12-Hox_expression_JSD_species/02-Drer_longest_isoform_Hox_DESeq2_average.txt')
raw_data
#.........................................Subsetting...............................................

#extract all geneIDs
names_row <- raw_data[,1]
names_row
#extract wanted gene ids (Chema excel) and create list
CT_IDs <- c('ENSDART00000167757',
            'ENSDART00000163611',
            'NM_131535.1',
            'ENSDART00000169480',
            'ENSDART00000166129',
            'ENSDART00000169017',
            'ENSDART00000171610',
            'ENSDART00000009827',
            'ENSDART00000079383',
            'ENSDART00000046766',
            'ENSDART00000006043',
            'ENSDART00000052662',
            'ENSDART00000110682',
            'ENSDART00000146636',
            'ENSDART00000184299',
            'ENSDART00000012470',
            'ENSDART00000025449',
            'ENSDART00000024256',
            'ENSDART00000078453',
            'ENSDART00000046638',
            'ENSDART00000023674',
            'ENSDART00000111227',
            'ENSDART00000078428',
            'ENSDART00000136415',
            'ENSDART00000076161',
            'ENSDART00000025966',
            'ENSDART00000076154',
            'ENSDART00000103131',
            'ENSDART00000154825',
            'ENSDART00000103132',
            'ENSDART00000139319',
            'ENSDART00000130090',
            'ENSDART00000103139',
            'ENSDART00000127384',
            'ENSDART00000179992',
            'ENSDART00000103146',
            'ENSDART00000103147',
            'ENSDART00000103149',
            'ENSDART00000190008',
            'ENSDART00000164839',
            'ENSDART00000170593',
            'ENSDART00000182896',
            'ENSDART00000146131',
            'ENSDART00000082355',
            'ENSDART00000082354',
            'ENSDART00000080608',
            'ENSDART00000082344',
            'ENSDART00000082339',
            'ENSDART00000082332')

CT_Names <- c('Hox1Aa',
              'Hox3Aa',
              'Hox4Aa',
              'Hox5Aa',
              'Hox9Aa',
              'Hox11Aa',
              'Hox13Aa',
              'Hox2Ab',
              'Hox9Ab',
              'Hox10Ab',
              'Hox11Ab',
              'Hox13Ab',
              'Hox1Ba',
              'Hox2Ba',
              'Hox3Ba',
              'Hox4Ba',
              'Hox5Ba',
              'Hox6Ba',
              'Hox7Ba',
              'Hox8Ba',
              'Hox9Ba',
              'Hox10Ba',
              'Hox13Ba',
              'Hox1Bb',
              'Hox5Bb',
              'Hox6Bb',
              'Hox8Bb',
              'Hox1Ca',
              'Hox3Ca',
              'Hox4Ca',
              'Hox5Ca',
              'Hox6Ca',
              'Hox8Ca',
              'Hox9Ca',
              'Hox10Ca',
              'Hox11Ca',
              'Hox12Ca',
              'Hox13Ca',
              'Hox6Cb',
              'Hox11Cb',
              'Hox12Cb',
              'Hox13Cb',
              'Hox3Da',
              'Hox4Da',
              'Hox9Da',
              'Hox10Da',
              'Hox11Da',
              'Hox12Da',
              'Hox13Da')


#Create emtpy list to add index of wanted rows
CT_index <- c()

for (i in CT_IDs) {
  output <- grep(i, names_row)
  CT_index[i] <- output
}

CT_index

#extract genes in hox_genesIDs

CTsubset <- raw_data[CT_index, ] 
CTsubset

#.....................................................................................................

df_curated <- data.matrix(CTsubset[,2:ncol(CTsubset)])
df_curated

rownames(df_curated) <- CT_Names


# 0 to 1 relative expression
rescale_custom <- function(x) (x/(max(x)))

df_normalised <- t(apply(df_curated, 1, rescale_custom))
df_normalised <- na.omit(df_normalised)
df_normalised

df_sorted <- df_normalised[match(CT_Names, rownames(df_normalised)),]
df_sorted

#1/2 rdbu (substitute of reds for 0 to 1)
heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color <- heatmap_color[50:100]

ComplexHeatmap::Heatmap(df_sorted,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        column_title = "Drer_Heatmap",
                        col = heatmap_color,
                        heatmap_legend_param = list(color_bar = "continuous"))

                       # column_labels = c('oocyte', 'fertilized embryo', 'polar body',	
                       #'2 cells', '4 cells', '8 cells', '16 cells',	'32 cells',	
                       #'blastula',	'emerged cilia', 'early trochophore',	'middle trochophore',
                       #'late trochophore',	'early segmentation', 'SL_mean',	'WL_mean'))library(tidyverse)
library(dplyr)
library(ggplot2)
library(data.table)
library(gplots)
library(pheatmap)
library(dendextend)
library(factoextra)
library(ComplexHeatmap)
library(RColorBrewer)
library(NbClust)
library(scales)


raw_data <- read.delim('/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/12-Hox_expression_JSD_species/02-Nvec_longest_isoform_Hox_DESeq2_average.txt')
raw_data
#.........................................Subsetting...............................................

#extract all geneIDs
names_row <- raw_data[,1]
names_row
#extract wanted gene ids (Chema excel) and create list
CT_IDs <- c('EDO42189',
'EDO45896',
'EDO42185',
'EDO42186',
'EDO49562',
'EDO49445')


CT_Names <- c('Anthox6',
'Anthox6a',
'Anthox7',
'Anthox8',
'Anthox1a',
'Anthox9')


#Create emtpy list to add index of wanted rows
CT_index <- c()

for (i in CT_IDs) {
  output <- grep(i, names_row)
  CT_index[i] <- output
}

CT_index

#extract genes in hox_genesIDs

CTsubset <- raw_data[CT_index, ] 
CTsubset

#.....................................................................................................

df_curated <- data.matrix(CTsubset[,2:ncol(CTsubset)])
df_curated

rownames(df_curated) <- CT_Names


# 0 to 1 relative expression
rescale_custom <- function(x) (x/(max(x)))

df_normalised <- t(apply(df_curated, 1, rescale_custom))
df_normalised <- na.omit(df_normalised)
df_normalised

df_sorted <- df_normalised[match(CT_Names, rownames(df_normalised)),]
df_sorted

#1/2 rdbu (substitute of reds for 0 to 1)
heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color <- heatmap_color[50:100]

ComplexHeatmap::Heatmap(df_sorted,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        column_title = "Nvec_Heatmap",
                        col = heatmap_color,
                        heatmap_legend_param = list(color_bar = "continuous"))

                       # column_labels = c('oocyte', 'fertilized embryo', 'polar body',	
                       #'2 cells', '4 cells', '8 cells', '16 cells',	'32 cells',	
                       #'blastula',	'emerged cilia', 'early trochophore',	'middle trochophore',
                       #'late trochophore',	'early segmentation', 'SL_mean',	'WL_mean'))library(tidyverse)
library(dplyr)
library(ggplot2)
library(data.table)
library(gplots)
library(pheatmap)
library(dendextend)
library(factoextra)
library(ComplexHeatmap)
library(RColorBrewer)
library(NbClust)
library(scales)


raw_data <- read.delim('/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/12-Hox_expression_JSD_species/02-Ofus_longest_isoform_Hox_DESeq2_average.txt')
raw_data
#.........................................Subsetting...............................................

#extract all geneIDs
names_row <- raw_data[,1]
names_row
#extract wanted gene ids (Chema excel) and create list
CT_IDs <- c('OFUSG13332.1',
            'OFUSG13333.1',
            'OFUSG13334.1',
            'OFUSG13335.1',
            'OFUSG13336.1',
            'OFUSG13337.1',
            'OFUSG13338.1',
            'OFUSG13339.1',
            'OFUSG13340.1',
            'OFUSG13341.1',
            'OFUSG23416.1')


CT_Names <- c('Hox1',
              'Hox2',
              'Hox3',
              'Hox4',
              'Hox5',
              'Lox5',
              'Antp',
              'Lox4',
              'Lox2',
              'Post2',
              'Post1')


#Create emtpy list to add index of wanted rows
CT_index <- c()

for (i in CT_IDs) {
  output <- grep(i, names_row)
  CT_index[i] <- output
}

CT_index

#extract genes in hox_genesIDs

CTsubset <- raw_data[CT_index, ] 
CTsubset

#.....................................................................................................

df_curated <- data.matrix(CTsubset[,2:ncol(CTsubset)])
df_curated

rownames(df_curated) <- CT_Names


# 0 to 1 relative expression
rescale_custom <- function(x) (x/(max(x)))

df_normalised <- t(apply(df_curated, 1, rescale_custom))
df_normalised <- na.omit(df_normalised)
df_normalised

df_sorted <- df_normalised[match(CT_Names, rownames(df_normalised)),]
df_sorted

#1/2 rdbu (substitute of reds for 0 to 1)
heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color <- heatmap_color[50:100]

ComplexHeatmap::Heatmap(df_sorted,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        column_title = "Ofus_Heatmap",
                        col = heatmap_color,
                        heatmap_legend_param = list(color_bar = "continuous"))

# column_labels = c('oocyte', 'fertilized embryo', 'polar body',	
#'2 cells', '4 cells', '8 cells', '16 cells',	'32 cells',	
#'blastula',	'emerged cilia', 'early trochophore',	'middle trochophore',
#'late trochophore',	'early segmentation', 'SL_mean',	'WL_mean'))
library(tidyverse)
library(dplyr)
library(ggplot2)
library(data.table)
library(gplots)
library(pheatmap)
library(dendextend)
library(factoextra)
library(ComplexHeatmap)
library(RColorBrewer)
library(NbClust)
library(scales)


raw_data <- read.delim('/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/12-Hox_expression_JSD_species/02-Spur_longest_isoform_Hox_DESeq2_average.txt')
raw_data
#.........................................Subsetting...............................................

#extract all geneIDs
names_row <- raw_data[,1]
names_row
#extract wanted gene ids (Chema excel) and create list
CT_IDs <- c('SPU_017352-tr',
            'SPU_012252-tr',
            'SPU_027568-tr',
            'SPU_005169-tr',
            'SPU_005171-tr',
            'SPU_002634-tr',
            'SPU_021309-tr',
            'SPU_002630-tr',
            'SPU_002633-tr',
            'SPU_002632-tr',
            'SPU_002631-tr',
            'SPU_000388-tr')


CT_Names <- c('Hox1',
              'Hox2',
              'Hox3',
              'Hox5',
              'Hox6',
              'Hox7',
              'Hox8',
              'Hox8',
              'Hox9/10',
              'Hox11/13a',
              'Hox11/13b',
              'Hox11/13c')


#Create emtpy list to add index of wanted rows
CT_index <- c()

for (i in CT_IDs) {
  output <- grep(i, names_row)
  CT_index[i] <- output
}

CT_index

#extract genes in hox_genesIDs

CTsubset <- raw_data[CT_index, ] 
CTsubset

#.....................................................................................................

df_curated <- data.matrix(CTsubset[,2:ncol(CTsubset)])
df_curated

rownames(df_curated) <- CT_Names


# 0 to 1 relative expression
rescale_custom <- function(x) (x/(max(x)))

df_normalised <- t(apply(df_curated, 1, rescale_custom))
df_normalised <- na.omit(df_normalised)
df_normalised

df_sorted <- df_normalised[match(CT_Names, rownames(df_normalised)),]
df_sorted

#1/2 rdbu (substitute of reds for 0 to 1)
heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color <- heatmap_color[50:100]

ComplexHeatmap::Heatmap(df_sorted,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        column_title = "Spur_Heatmap",
                        col = heatmap_color,
                        heatmap_legend_param = list(color_bar = "continuous"))

                       # column_labels = c('oocyte', 'fertilized embryo', 'polar body',	
                       #'2 cells', '4 cells', '8 cells', '16 cells',	'32 cells',	
                       #'blastula',	'emerged cilia', 'early trochophore',	'middle trochophore',
                       #'late trochophore',	'early segmentation', 'SL_mean',	'WL_mean'))library(tidyverse)
library(dplyr)
library(ggplot2)
library(data.table)
library(gplots)
library(pheatmap)
library(dendextend)
library(factoextra)
library(ComplexHeatmap)
library(RColorBrewer)
library(NbClust)
library(scales)

# go to extended data figure 10 to change names

raw_data <- read.delim('/Users/billie/Dropbox/02-OweniaGenome/06-Revision/00-DATA/12-Hox_expression_JSD_species/02-Uuni_longest_isoform_Hox_DESeq2_average_Hox_genes_ONLY.txt')
raw_data
#.........................................Subsetting...............................................

#extract all geneIDs
names_row <- raw_data[,1]
names_row
#extract wanted gene ids (Chema excel) and create list
CT_IDs <- c('lab',
              'pb',
              'Hox3',
              'Dfd',
              'Scr',
              'Lox5',
              'Antp',
              'Lox4',
              'Lox2',
              'Post2')


CT_Names <- c('lab',
              'pb',
              'Hox3',
              'Dfd',
              'Scr',
              'Lox5',
              'Antp',
              'Lox4',
              'Lox2',
              'Post2'
              )


#Create emtpy list to add index of wanted rows
CT_index <- c()

for (i in CT_IDs) {
  output <- grep(i, names_row)
  CT_index[i] <- output
}

CT_index
CT_cols <- c(1, )

#extract genes in hox_genesIDs and columns

CTsubset <- raw_data[CT_index,] 
CTsubset

#.....................................................................................................

df_curated <- data.matrix(CTsubset[,2:ncol(CTsubset)])
df_curated

rownames(df_curated) <- CT_Names


# 0 to 1 relative expression
rescale_custom <- function(x) (x/(max(x)))

df_normalised <- t(apply(df_curated, 1, rescale_custom))
df_normalised <- na.omit(df_normalised)
df_normalised

df_sorted <- df_normalised[match(CT_Names, rownames(df_normalised)),]
df_sorted

#1/2 rdbu (substitute of reds for 0 to 1)
heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
heatmap_color <- heatmap_color[50:100]

ComplexHeatmap::Heatmap(df_sorted,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        column_title = "Uuni_Heatmap",
                        col = heatmap_color,
                        heatmap_legend_param = list(color_bar = "continuous"))

# column_labels = c('oocyte', 'fertilized embryo', 'polar body',	
#'2 cells', '4 cells', '8 cells', '16 cells',	'32 cells',	
#'blastula',	'emerged cilia', 'early trochophore',	'middle trochophore',
#'late trochophore',	'early segmentation', 'SL_mean',	'WL_mean'))
