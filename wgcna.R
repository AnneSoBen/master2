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
