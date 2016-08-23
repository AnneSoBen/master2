###################################
##### DATA INPUT AND CLEANING #####
###################################



# from Langfelder and Horvath 2008

###########
### HOW TO:
# R --vanilla --args abundancesfile.tsv traitsfile.tsv < wgcna.R

# both data files are tab separated files (.tsv)

# structure of file abundancesfile.tsv:
#           gene1   gene2   gene3
# sample1   ab11    ab12    ab13
# sample2   ab21    ab22    ab23
# sample2   ab31    ab32    ab33

# genes can also be species or any other type of individuals we want to analyze

# structure of file traitsfile.tsv
#           trait1    trait2    trait3
# sample1   val11     val12     val13
# sample2   val21     val22     val23
# sample3   val31     val32     val33

# traits, in our study, are environmental parameters measured in each sample
#########


# the function commandArgs scans the arguments which have been supplied when the current R session was invoked
#args <- commandArgs(trailingOnly=TRUE)
#abundancesfile <- args[1]
#traitsfile <- args[2]

# create a "basename" for plots file names from abundancesfile

# remove extension
library(tools)
basename<-file_path_sans_ext(abundancesfile)

# create a directory for the plots
plotsDirectory = paste(basename,"Plots", sep="")
if (!file.exists(plotsDirectory)){
  dir.create(plotsDirectory)
}

## LOAD PACKAGES
library(WGCNA)
library(vegan)
library(gtools)

abundancesfile = "eukNormAbArcticFP.tsv"
traitsfile = "envdataok.tsv"

#=====================================================================================
#
#  Part 1: read files, clean them and Hellinger transformation
#
#=====================================================================================


# The following setting is important, do not omit
options(stringsAsFactors = FALSE);
# Read in the abundances data set
abundances = read.delim(abundancesfile);

# for eukaryotes : remove samples with missing size fractions
#abundances = abundances[-c(1,2,5,6,13,25,31,36,38,41,64),]

# assign the rownames to the data frame taxd
taxd = as.data.frame(abundances[,2:ncol(abundances)])
rownames(taxd) = abundances[,1]

## taxd -> species or genes relative abundances (norm. to 1) - samples as rows and species/genes as columns
## md   -> environmental data (no transformation) - samples as rows and parameters as columns 

env = read.delim(traitsfile)
md = as.data.frame(env[,2:ncol(env)])
rownames(md) = env[,1]

## intersecting samples on species and env. matrices
stations = intersect(rownames(taxd), rownames(md))
taxdma = taxd[which(rownames(taxd) %in% stations),]
mdma = md[which(rownames(md) %in% stations),]

taxdma = taxdma[mixedsort(rownames(taxdma)),]
mdma = mdma[mixedsort(rownames(mdma)),]

# remove null abundance species or genes
taxdma = taxdma[,-which(colSums(taxdma)==0)]

# Hellinger transformation
taxdma = decostand(taxdma, method="hellinger")


#=====================================================================================
#
#  Part 2: sample clustering and trait heatmap
#
#=====================================================================================


# Re-cluster samples
sampleTree2 = hclust(dist(taxdma), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(mdma, signed = TRUE);
# Plot the sample dendrogram and the colors underneath.
pdf(file = paste(plotsDirectory, "/", basename, "SampleClusteringAndTraitHeatmap.pdf", sep=""), width = 12, height = 9);
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(mdma), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

tiff(file = paste(plotsDirectory, "/", basename, "SampleClusteringAndTraitHeatmap.tiff", sep=""), bg="white", res=300, width = 12, height = 9, units='in');
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(mdma), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

#####################################################
##### NETWORK CONSTRUCTION AND MODULE DETECTION #####
#####################################################

#=====================================================================================
#
#  Part 3: select the soft-thesholding
#
#=====================================================================================

threshold = 0.9

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
# RsquaredCut corresponds to the threshold above which the power is selected
sft = pickSoftThreshold(taxdma, powerVector = powers, verbose = 5, RsquaredCut = threshold)
# Plot the results:
# pdf
pdf(file = paste(plotsDirectory, "/", basename, "NetworkTopologyAnalysis.pdf", sep = ""), width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=threshold,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# tiff
tiff(file = paste(plotsDirectory, "/", basename, "NetworkTopologyAnalysis.tiff", sep = ""), bg="white", res=300, width = 9, height = 5, units='in');
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=threshold,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# sft$powerEstimate is the first power that reaches the threshold
sft$powerEstimate


max(sft$fitIndices[,2])

if (is.na(sft$powerEstimate)){
    print("THE POWER THRESHOLD IS TOO HIGH! The highest scale-free topology fit index is:")
    print(max(sft$fitIndices[,2]))
    q()
}

#=====================================================================================
#
#  Part 4: Modules detection
#
#=====================================================================================


# calculate the adjacencies, using the right soft-thresholding
softPower = sft$powerEstimate;
adjacency = adjacency(log(taxdma+1), power = softPower, type = "signed", corOptions= "method = 'pearson'");

# Turn adjacency into topological overlap matrix, and calculate the corresponding dissimilarity
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
# pdf
pdf(file = paste(plotsDirectory, "/", basename, "GeneClusteringOnTOMbasedDissimilarity.pdf", sep = ""), width = 12, height = 9);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()

# tiff
tiff(file = paste(plotsDirectory, "/", basename, "GeneClusteringOnTOMbasedDissimilarity.tiff", sep = ""), width = 12, height = 9, bg="white", res=300, units='in');
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 20;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                deepSplit = 2, pamRespectsDendro = FALSE,
                minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
# pdf
pdf(file = paste(plotsDirectory, "/", basename, "GeneDendrogramAndModuleColors.pdf", sep = ""), width = 12, height = 9);
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

tiff(file = paste(plotsDirectory, "/", basename, "GeneDendrogramAndModuleColors.tiff", sep = ""), width = 12, height = 9, bg="white", res=300, units='in');
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(taxdma, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
# pdf
pdf(file = paste(plotsDirectory, "/", basename, "ClusteringOfModuleEigengenes.pdf", sep = ""), width = 7, height = 6);
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
     
# height cut of 0.25 corresponding to correlation of 0.75 to merge modules whose expression profiles are very similar
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

# tiff
tiff(file = paste(plotsDirectory, "/", basename, "ClusteringOfModuleEigengenes.tiff", sep = ""), width = 7, height = 6, bg="white", res=300, units='in');
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# height cut of 0.25 corresponding to correlation of 0.75 to merge modules whose expression profiles are very similar
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

# Call an automatic merging function
merge = mergeCloseModules(taxdma, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

# pdf
pdf(file = paste(plotsDirectory, "/", basename, "MergedClusteringOfModuleEigengenes.pdf", sep = ""), wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# tiff
tiff(file = paste(plotsDirectory, "/", basename, "MergedClusteringOfModuleEigengenes.tiff", sep = ""), wi = 9, he = 6, bg="white", res=300, units='in')
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

###################################################
##### CORRELATION TO ENVIRONMENTAL PARAMETERS #####
###################################################

# Define numbers of genes and samples
nGenes = ncol(taxdma);
nSamples = nrow(taxdma);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(taxdma, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, mdma, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# take a look at p-values before correction
typeof(moduleTraitPvalue)
moduleTraitPvalue

# function that allows to use the BH method for correction of FDR
p.adjust.BH <- function(x){
  p.adjust(x, method = "BH")
}

y = apply(moduleTraitPvalue, 2, p.adjust.BH)
moduleTraitPvalue = y

# take a look at p-values after correction
typeof(moduleTraitPvalue)
moduleTraitPvalue

mtpv = as.data.frame(moduleTraitPvalue)
modules = rownames(mtpv[mtpv[,interestVariable]<0.05,])
modules
typeof(modules)

# pdf
pdf(file = paste(plotsDirectory, "/", basename, "ModuleTraitRelationships.pdf", sep = ""), width = 10, height = 9);
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(mdma),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix, #comment to remove Pearson correlation and p-value
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# tiff
tiff(file = paste(plotsDirectory, "/", basename, "ModuleTraitRelationships.tiff", sep = ""), width = 10, height = 9, bg="white", res=300, units='in');
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(mdma),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix, #comment to remove Pearson correlation and p-value
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


