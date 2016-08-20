# Anne-Sophie Benoiston --- may/june/july 2016
#     /\ /\
#    (=^_^=)

###

library(vegan)

data = readRDS("eukNormAbOceanArcticFP.rds")
dataRegions = readRDS("eukNormAbOceanArcticFPR.rds")

dataRegions$anosim.group = ifelse(dataRegions$Region %in% c("AO","SO"), "PO", "OO")

dataRegions$anosim.group = factor(dataRegions$anosim.group)

dissiMatrix = vegdist(data, method = "bray")
data.ano = anosim(dissiMatrix, dataRegions$anosim.group)

fileName = "eukmetaB"

tiff(file = paste(fileName, "anosim300.tiff", sep = ""), bg="white", res=300, width = 9, height = 9, units='in');
plot(data.ano, ylab = "Dissimilarity ranks", main = "ANOSIM / all fractions (-S67SUR)", col = "grey")
dev.off()
