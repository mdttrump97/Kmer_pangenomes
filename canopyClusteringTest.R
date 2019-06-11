library(FindMyFriends)
library(kebabs)
library(rlist)

setwd("/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/FAA_files")
genomeFiles <- list.files(getwd(), full.names=TRUE, pattern='*.faa')
pan <- pangenome(genomeFiles[1:3], translated=TRUE, geneLocation='prodigal', lowMem=FALSE)
distance_matrix <- kmerSimilarity(pan, kmerSize = 4, lowerLimit=0, rescale=TRUE)

cleaned_names <- c()
for(name in rownames(distance_matrix)){
  cleaned_names <- c(cleaned_names, strsplit(name, " ")[[1]][1])
}
rownames(distance_matrix) <- cleaned_names
colnames(distance_matrix) <- cleaned_names

points <- rownames(distance_matrix)
canopies <- list()
while(length(points) > 0){
  center_point = tail(points, 1)
  print(c("canopy center point: ", center_point))
  points <- points[!(points == center_point)]
  
  T1 <- 0.002
  T2 <- 0.001
  canopy_points <- c()
  
  for(entry in names(distance_matrix[center_point,])){
    if(distance_matrix[center_point, entry] < T1){
      canopy_points <- c(canopy_points, entry)
      if(distance_matrix[center_point, entry] < T2){
        points <- points[names(points) != entry] 
        #distance_matrix <- distance_matrix[,!colnames(distance_matrix) == entry, ]
      }
    }
  }
  print(c("points not in canopy at end", length(points)))
  canopies[[center_point]] <- canopy_points[]
}