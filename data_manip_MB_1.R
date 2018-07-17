
source("http://bioconductor.org/biocLite.R") 
biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
install.packages("WGCNA")
library(WGCNA)


#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
#workingDir = "/Users/anni_djurhuus/Documents/Projects/R/WGCNA/FemaleLiver-Data";
#setwd(workingDir); 
# Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the female liver data set
MB_Data = read.csv("/Users/anni_djurhuus/Documents/Projects/MBON/m2w/google_drive/data/MEGAN6/post_decont/MB/MB_all_taxa.csv")
MB_Data <- MB_Data()
# Take a quick look at what is in the data set:
dim(MB_Data);
names(MB_Data);

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


datExpr0_MB = as.data.frame(t(MB_Data))

names(datExpr0_MB) = MB_Data$X
rownames(datExpr0_MB) = names(MB_Data)
datExpr0_MB <- datExpr0_MB[-1,]
datExpr0_MB[, c(1:659)] <- sapply(datExpr0_MB[, c(1:659)], as.numeric)

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


gsg = goodSamplesGenes(datExpr0_MB, verbose = 3);
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
    printFlush(paste("Removing genes:", paste(names(datExpr0_MB)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0_MB)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0_MB = datExpr0_MB[gsg$goodSamples, gsg$goodGenes]
}


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


sampleTree_MB = hclust(dist(datExpr0_MB), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree_MB, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


# Plot a line to show the cut
abline(h = 15, col = "red");
# Determine cluster under the line
clust_MB = cutreeStatic(sampleTree_MB, cutHeight = 15, minSize = 10)
table(clust_MB)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr_MB = datExpr0_MB[keepSamples, ]
nGenes_MB = ncol(datExpr_MB)
nSamples_MB = nrow(datExpr_MB)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


traitData_MB = read.csv("/Users/anni_djurhuus/Documents/Projects/MBON/m2w/google_drive/data/MEGAN6/metadata/m2w_metadata_MB.csv")
dim(traitData_MB)
names(traitData_MB)

# remove columns that hold information we do not need.
allTraits_MB = traitData_MB[, -c(1:2,4:14,17:43,47,50:61)];
#allTraits_MB = allTraits_MB[, c(2, 11:36) ];
dim(allTraits_MB)
names(allTraits_MB)

# Form a data frame analogous to expression data that will hold the clinical traits.
Samples_MB = rownames(datExpr_MB)
traitRows_MB = match(Samples_MB, allTraits_MB$correct_sample_name);
datTraits_MB = allTraits_MB[traitRows_MB, -1];
rownames(datTraits_MB) = allTraits_MB[traitRows_MB, 1];

collectGarbage();


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


# Re-cluster samples
sampleTree2_MB = hclust(dist(datExpr_MB), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors_MB = numbers2colors(datTraits_MB, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2_MB, traitColors_MB,
                    groupLabels = names(datTraits_MB), 
                    main = "Sample dendrogram and trait heatmap")


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


save(datExpr_MB, datTraits_MB, file = "/Users/anni_djurhuus/Documents/Projects/MBON/m2w/google_drive/data/MEGAN6/network_data/MB_data_network.RData")

