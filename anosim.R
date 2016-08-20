# Anne-Sophie Benoiston --- may/june/july 2016
#     /\ /\
#    (=^_^=)

###

# load the "vegan" library
library(vegan)

# load abundance files (prepared with "datasetFiltering.R")
data = readRDS("eukNormAbOceanArcticFP.rds")
dataRegions = readRDS("eukNormAbOceanArcticFPR.rds") # contains information about the regions (e.g. Arctic Ocean, South Pacific Ocean)

# create a new variable (groups arctic and antarctic stations together)
dataRegions$anosim.group = ifelse(dataRegions$Region %in% c("AO","SO"), "PO", "OO")
dataRegions$anosim.group = factor(dataRegions$anosim.group)

dissiMatrix = vegdist(data, method = "bray")
data.ano = anosim(dissiMatrix, dataRegions$anosim.group)

fileName = "eukmetaB"

# plot de results
tiff(file = paste(fileName, "anosim300.tiff", sep = ""), bg="white", res=300, width = 9, height = 9, units='in');
plot(data.ano, ylab = "Dissimilarity ranks", main = "ANOSIM / all fractions (-S67SUR)", col = "grey")
dev.off()
