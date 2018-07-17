#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames_MB = load(file = "/Users/anni_djurhuus/Documents/Projects/MBON/m2w/google_drive/data/MEGAN6/network_data/MB_data_network.RData");
#The variable lnames contains the names of loaded variables.
lnames_MB


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft_MB = pickSoftThreshold(datExpr_MB, powerVector = powers, verbose = 5)
# Plot the results:

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft_MB$fitIndices[,1], -sign(sft_MB$fitIndices[,3])*sft_MB$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft_MB$fitIndices[,1], -sign(sft_MB$fitIndices[,3])*sft_MB$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft_MB$fitIndices[,1], sft_MB$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft_MB$fitIndices[,1], sft_MB$fitIndices[,5], labels=powers, cex=cex1,col="red")


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


net_MB = blockwiseModules(datExpr_MB, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "MB_TOM", 
                       verbose = 3)


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors_MB = labels2colors(net_MB$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net_MB$dendrograms[[1]], mergedColors_MB[net_MB$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


moduleLabels_MB = net_MB$colors
moduleColors_MB = labels2colors(net_MB$colors)
MEs_MB = net_MB$MEs;
geneTree_MB = net_MB$dendrograms[[1]];
save(MEs_MB, moduleLabels_MB, moduleColors_MB, geneTree_MB, 
     file = "/Users/anni_djurhuus/Documents/Projects/MBON/m2w/google_drive/data/MEGAN6/network_data/MB_network_construction.RData")

