# Anne-Sophie Benoiston --- may/june 2016
#     /\ /\
#    (=^_^=)

###

# load vegan library
library(vegan)

# read the abundances files
data = readRDS("eukNormAbOceanArcticFP.rds")
dataRegions = readRDS("eukNormAbOceanArcticFPR.rds")

dataRegions$Depth = factor(dataRegions$Depth)

dataRegions$anosim.group = ifelse(dataRegions$Region %in% c("AO", "SO"), "PO", "OO")
dataRegions$anosim.group = factor(dataRegions$anosim.group)

dissiMatrix = vegdist(data, method = "bray")
data.nmds = metaMDS(data)

fileName = "AOvsOO"

# shape of the points to distinguish DCM and SUR
pts = c(16,17)

# plot the results
tiff(file = paste(fileName, "nmds.tiff", sep = ""), bg="white", res=300, width = 9, height = 9, units='in');
cols = c("cadetblue","darkblue")
plot(data.nmds, type = "n", display = "sites", main = paste("NMDS/Bray - Stress =", round(data.nmds$stress,3))) # type = "t" to print samples names on the plot
points(data.nmds, col = cols[dataRegions$anosim.group], pch = pts[dataRegions$Depth])
legend("topleft", legend = levels(dataRegions$anosim.group), col = cols, pch = 15)
legend("topright", legend = levels(dataRegions$Depth), col = "black", pch = pts)
dev.off()
