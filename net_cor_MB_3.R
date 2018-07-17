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
# Load the expression and trait data saved in the first part
lnames_MB = load(file = "/Users/anni_djurhuus/Documents/Projects/MBON/m2w/google_drive/data/MEGAN6/network_data/MB_data_network.RData");
#The variable lnames contains the names of loaded variables.
lnames_MB
# Load network data saved in the second part.
lnames_MB = load(file = "/Users/anni_djurhuus/Documents/Projects/MBON/m2w/google_drive/data/MEGAN6/network_data/MB_network_construction.RData");
lnames_MB


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Define numbers of genes and samples
nGenes_MB = ncol(datExpr_MB);
nSamples_MB = nrow(datExpr_MB);
# Recalculate MEs with color labels
MEs0_MB = moduleEigengenes(datExpr_MB, moduleColors_MB)$eigengenes
MEs_MB = orderMEs(MEs0_MB)
moduleTraitCor_MB = cor(MEs_MB, datTraits_MB, use = "p");
moduleTraitPvalue_MB = corPvalueStudent(moduleTraitCor_MB, nSamples_MB);


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix_MB =  paste(signif(moduleTraitCor_MB, 2), "\n(",
                    signif(moduleTraitPvalue_MB, 1), ")", sep = "");
dim(textMatrix_MB) = dim(moduleTraitCor_MB)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor_MB,
               xLabels = names(datTraits_MB),
               yLabels = names(MEs_MB),
               ySymbols = names(MEs_MB),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix_MB,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Define variable weight containing the weight column of datTrait
temp_MB = as.data.frame(datTraits_MB$temp);
names(temp_MB) = "temperature"
# names (colors) of the modules
modNames_MB = substring(names(MEs_MB), 3)

geneModuleMembership_MB = as.data.frame(cor(datExpr_MB, MEs_MB, use = "p"));
MMPvalue_MB = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_MB), nSamples_MB));

names(geneModuleMembership_MB) = paste("MM", modNames_MB, sep="");
names(MMPvalue_MB) = paste("p.MM", modNames_MB, sep="");

geneTraitSignificance_MB = as.data.frame(cor(datExpr_MB, temp_MB, use = "p"));
GSPvalue_MB = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_MB), nSamples_MB));

names(geneTraitSignificance_MB) = paste("GS.", names(temp_MB), sep="");
names(GSPvalue_MB) = paste("p.GS.", names(temp_MB), sep="");


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


module = "blue"
column = match(module, modNames_MB);
moduleGenes_MB = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership_MB[moduleGenes_MB, column]),
                   abs(geneTraitSignificance_MB[moduleGenes_MB, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for temperature",
                   main = paste("Module membership vs. taxa significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, pch=19,col = module)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


names(datExpr_MB)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


names(datExpr_MB)[moduleColors_MB=="brown"]


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


annot = read.csv(file = "GeneAnnotation.csv");
dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = probes,
                       geneSymbol = annot$gene_symbol[probes2annot],
                       LocusLinkID = annot$LocusLinkID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
geneInfo = geneInfo0[geneOrder, ]


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


write.csv(geneInfo, file = "geneInfo.csv")

