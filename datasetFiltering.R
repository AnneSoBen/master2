#########################
### DATASET FILTERING ###
#########################

# Anne-Sophie Benoiston --- may/june/july 2016
#     /\ /\
#    (=^_^=)

gdBarEukRaw = read.delim("globaldataset.barcode")
# 6 294 617 barcodes, 1085 samples

saveRDS(gdBarEukRaw, "globaldataset.barcode.v20151126.rds")

# remove columns we don't need :
# column 1 = md5sum (barcode id)
# column 2 = totab (total abundance of the barcode) - as we remove samples, the total is not the same
# column 3 = cid (identity of the swarm to which this sequence belongs to)
# column 5 = reference sequences
# column 8 = sequence
# columns 9 to 48 = BV9_X (samples from Biomarks)

gdBarEuk = gdBarEukRaw[,-c(1:3,5,8:48)]
# 6 294 617 barcodes, 1046 samples

# remove bacteria and archaea
gdBarEuk = gdBarEuk[!gdBarEuk$taxogroup %in% c("Bacteria","Archaea"),]
# 5 625 064 barcodes, 1046 samples

# barcodes with pid < 97% are considered as unknown
levels(gdBarEuk$lineage) <- c(levels(gdBarEuk$lineage), "unknown")
gdBarEuk$lineage[gdBarEuk$pid < 97] = "unknown"
gdBarEuk$lineage[gdBarEuk$lineage %in% c("Eukaryota", "Eukaryota|environmental samples", "Eukaryota|environmental samples|uncultured eukaryote", "Eukaryota|environmental samples|uncultured marine eukaryote", "Eukaryota|environmental samples|uncultured rumen protozoa", "Eukaryota|Eukaryota_X|Eukaryota_X+sp.", "Eukaryota|Orphans")] = "unknown"

# remove taxogroup and pid
gdBarEuk = gdBarEuk[,-c(1,3)]

saveRDS(gdBarEuk, "gdBarEuk.rds")

# sum the abundances of each lineage
gdBarEukPool = aggregate(.~lineage, data = gdBarEuk, sum)
# 3948 lineages, 1046 samples
saveRDS(gdBarEukPool, "gdBarEukPool.rds")

# transpose the data frame
gdBarEukPoolT = as.data.frame(t(gdBarEukPool))
colnames(gdBarEukPoolT) = gdBarEukPool$lineage
gdBarEukPoolT = gdBarEukPoolT[2:nrow(gdBarEukPoolT),]

saveRDS(gdBarEukPoolT, "gdBarEukPoolT.rds")

# load sampleIdToPangaeaId and sampleDepthStationRegion
sampleIdToPangaeaId = read.delim("sampleIdToPangaeaId.tsv")
sampleDepthStationRegion = read.csv("sampleDepthStationRegion.csv", sep = " ")

# merge the two data frames (both data frames has a column with the sample name, e.g. "S1SUR")

# we get a data frame with:
# Sample_seq_id Sample Fraction Template Pangaea_sample_id Depth Station Region 
addInfos = merge(sampleIdToPangaeaId, sampleDepthStationRegion, by.x = "Sample", by.y = "Sample")
saveRDS(addInfos, "addInfos.rds")

# add the corresponding columns with stations names, depth and fraction size
eukRawAb = merge(gdBarEukPoolT, addInfos, by.x = "row.names", by.y = "Sample_seq_id")

# move the columns at the begining of the data frame
eukRawAb = eukRawAb[,c(1,3948:3954,2:3947)]
saveRDS(eukRawAb, "eukRawAb.rds")

# keep the DNA samples (so remove the RNA and WGA/DNA samples)
eukRawAb = subset(eukRawAb, Template == "DNA")
# 3946 lineages, 917 samples

###### COUNT THE NUMBER OF FRACTIONS AVAILABLE FOR EACH SAMPLE-DEPTH

# make a table with fractions available for each station/sample
# make a list of stations
listSamples = as.character(unique(eukRawAb$Sample))

# create a matrix with the fractions as columns and the samples as rows

mat = matrix(0L, nrow = length(listSamples), ncol = 13) # 0L means it is filled with zeros
rownames(mat) = listSamples
colnames(mat) = c("0.22-3","0.8-inf","0.8-3","0.8-5","0.8-20","0.8-180","0.8->","3-inf","3-20","5-20","20-180","20-200","180-2000")

# fill the matrix with the number of fractions per sample
for (i in listSamples) {
	for (j in colnames(mat)) {
		occ = as.data.frame(table(eukRawAb[eukRawAb$Sample == i,]$Fraction))
		mat[i,j] = occ[occ$Var1 == j,]$Freq
	}
}

# 3/5-2000
mat2 = as.data.frame(mat)
mat2$sum1 = rowSums(mat2[,c(9,11,13)])
mat2$sum2 = rowSums(mat2[,c(10,11,13)])

list3_2000 = rownames(mat2[mat2$sum1 >= 3,])
list5_2000 = rownames(mat2[mat2$sum2 >= 3,])

# "S124SUR" and "S153MIX" are in the two lists, I choose to remove them from list3_2000
rm = c("S124SUR", "S153MIX")
list3_2000 = list3_2000[!list3_2000 %in% rm]

list3_5_2000 = c(list3_2000, list5_2000)

saveRDS(list3_5_2000, "list3_5_2000.rds")

# filter unknown depths
# remove FSW, ZZZ and INT samples
eukRawAb = eukRawAb[!eukRawAb$Depth %in% c("ZZZ","FSW","INT"),]

# remove S152MIX, S153MIX, S175MIX, S168MIX
eukRawAb = eukRawAb[!eukRawAb$Sample %in% c("S152MIX", "S153MIX", "S175MIX", "S168MIX"),]
# 3946 lineages, 881 samples
saveRDS(eukRawAb, "eukRawAb2.rds")

# prepare the abundance table with fractions

# keep the fractions [3-20],[5-20],[20-180],[180-2000]
eukRawAbOceanArctic = subset(eukRawAb, Fraction %in% c("3-20","5-20","20-180","180-2000"))

eukRawAbOceanArctic = subset(eukRawAbOceanArctic, Sample %in% list3_5_2000)

######

# remove the columns we don't need to aggregate abundances of the three fractions (Row.names, Fraction, Template, Pangaea_sample_id, depth, Station, Region)
eukRawAbOceanArcticFP1 = eukRawAbOceanArctic[,-c(1,3:8)]

# change the type of columns to sum the abundances
for (i in colnames(eukRawAbOceanArcticFP1)){
	if (i=="Sample"){
		eukRawAbOceanArcticFP1[,i] = as.character(eukRawAbOceanArcticFP1[,i])
	}
	else{
		eukRawAbOceanArcticFP1[,i] = as.numeric(as.character(eukRawAbOceanArcticFP1[,i]))
	}

} 

# sum the abundances of each sample (=station + depth)
eukRawAbOceanArcticFP = aggregate(.~Sample, data = eukRawAbOceanArcticFP1, sum)
rownames(eukRawAbOceanArcticFP) = eukRawAbOceanArcticFP$Sample
eukRawAbOceanArcticFP = eukRawAbOceanArcticFP[,-1]
# 3946 lineages, 105 samples (with fractions pooled)

# normalization
eukNormAbOceanArcticFP = eukRawAbOceanArcticFP/rowSums(eukRawAbOceanArcticFP)

# remove lineages with abundance = 0
eukNormAbOceanArcticFP = eukNormAbOceanArcticFP[,-which(colSums(eukNormAbOceanArcticFP)==0)]
# 2919 lineages, 105 samples

saveRDS(eukNormAbOceanArcticFP, "eukNormAbOceanArcticFP.rds")

eukNormAbOceanArcticFPR = merge(eukNormAbOceanArcticFP, sampleDepthStationRegion, by.x = "row.names", by.y = "Sample", all.x = FALSE, all.y = FALSE)

saveRDS(eukNormAbOceanArcticFPR, "eukNormAbOceanArcticFPR.rds")

listArctic = as.character(sampleDepthStationRegion[sampleDepthStationRegion$Region == "AO",]$Sample)

remove = c("S155DCM", "S145SUR", "S155SUR", "S81DCM", "S82SUR", "S90SUR", "S87SUR", "S86SUR")
listArcticPlus = c(listArctic, remove)

saveRDS(listArcticPlus, "listArcticPlus.rds")
eukNormAbArcticPlusFP = subset(eukNormAbOceanArcticFP, rownames(eukNormAbOceanArcticFP) %in% listArcticPlus)

saveRDS(eukNormAbArcticPlusFP, "eukNormAbArcticPlusFP.rds")

listOcean = as.character(sampleDepthStationRegion[sampleDepthStationRegion$Region != "AO",]$Sample)

listOceanLess = listOcean[!listOcean %in% remove]
saveRDS(listOceanLess, "listOceanLess.rds")
eukNormAbOceanLessFP = subset(eukNormAbOceanArcticFP, rownames(eukNormAbOceanArcticFP) %in% listOceanLess)

saveRDS(eukNormAbOceanLessFP, "eukNormAbOceanLessFP.rds")




