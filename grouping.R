library("FindMyFriends")
library("kebabs")
library("Matrix")

setwd("/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/FAA_files/20_genomes")
genomeFiles <- list.files(getwd(), full.names=TRUE, pattern='*.faa')
pan <- pangenome(genomeFiles[1:20], translated=TRUE, geneLocation='prodigal', lowMem=FALSE)
setwd("/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/20_genomes/")
groupingList <- as.integer(read.csv("20_genome_4mer_top_5_all_FYWH_CM_grouping_list.csv", header=FALSE, sep=',', stringsAsFactors = FALSE))
pan <- manualGrouping(pan, groupingList)
defaults(pan)$coreThreshold <- 0.9
plotStat(pan, type='qual', palette = 6)
plotEvolution(pan)
pan