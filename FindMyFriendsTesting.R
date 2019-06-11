library(FindMyFriends)
library(igraph)

setwd("/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/FAA_files")
genomeFiles <- list.files(getwd(), full.names=TRUE, pattern='*.faa')
pan <- pangenome(genomeFiles[1:8], translated=TRUE, geneLocation='prodigal', lowMem=FALSE)
pan <- gpcGrouping(pan, lowerLimit=0.5)
pan

plotStat(pan, type='qual', palette=6)

pan <- pangenome(genomeFiles[1:8], translated=TRUE, geneLocation='prodigal', lowMem=FALSE)
panSim <- kmerSimilarity(pan, lowerLimit=0.8, rescale=FALSE)
pan <- graphGrouping(pan, panSim)
pan
plotStat(pan, type='qual', palette=6)

pan2 <- pangenome(genomeFiles[1:8], translated=TRUE, geneLocation='prodigal', lowMem=FALSE)
panSim2<- kmerSimilarity(pan2, lowerLimit=0.8, rescale=FALSE)
pan2 <- graphGrouping(pan2, panSim2)
pan2
plotStat(pan2, type='qual', palette=6)

pan3 <- pangenome(genomeFiles[1:8], translated=TRUE, geneLocation='prodigal', lowMem=FALSE)
panSim3<- kmerSimilarity(pan3, lowerLimit=0.8, rescale=FALSE)
pan3 <- graphGrouping(pan3, panSim3)
pan3
plotStat(pan3, type='qual', palette=6)