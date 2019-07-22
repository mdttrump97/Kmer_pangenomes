library("FindMyFriends")

# Read in .faa files (same input as kmerCounter.R)
setwd("/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/FAA_files/medioid_dataset")
genomeFiles <- list.files(getwd(), full.names=TRUE, pattern='*.faa')
pan <- pangenome(genomeFiles[1:length(genomeFiles)], translated=TRUE, geneLocation='prodigal', lowMem=FALSE)

# Set core gene group threshold
#defaults(pan)$coreThreshold <- 0.9

setwd("/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/medioid_3mers/")

# Read in cluster membership from canopyClustering.py
groups <- read.csv('medioid_3mer_top_9_0.675_grouping_list.csv', sep = ',')
groupingList <- groups[,"clust"]
pan <- manualGrouping(pan, groupingList)

# Plot pie chart of pan-genome composition
plotStat(pan, type='qual', palette = 6)
# Plot pan-genome composition as genomes are added
plotEvolution(pan)
pan