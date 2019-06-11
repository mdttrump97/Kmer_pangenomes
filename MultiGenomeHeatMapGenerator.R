library("ape")
library("kmer")
library("MASS")
setwd("C:/Users/Matthew/Documents/UAMS SURF/K-mer testing/FAA files")

files <- scan(file = "listOfFAAFiles.txt", what = "character")

i <- length(files) - 1
for(filename in files) {
	sequences[i] <- read.FASTA(filename, type = "AA")
}

DistanceMatrix <- kdistance(sequences, k = 4)
Tree <- as.dendrogram(hclust(DistanceMatrix), "average")
plot(heatmap(as.matrix(DistanceMatrix), Rowv=Tree, Colv = "Rowv", scale = "none"))