library("ape")
library("kmer")
library("Matrix")

input_folder = "/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/FAA_files/medioid_dataset"
output_folder = "/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/medioid_3mers/"
kmer_length = 3
kmer_string = paste(toString(kmer_length), "mer", sep = '')
  
setwd(input_folder)
genome_files <- list.files(getwd(), full.names=TRUE, pattern='*.faa')

# Output gene ordering from FindMyFriends for canopyClustering.py output
pan <- pangenome(genomeFiles[1:length(genomeFiles)], translated=TRUE, geneLocation='prodigal', lowMem=FALSE)
grouping_list <- c()
for(x in pan@sequences@ranges@NAMES)
{
  grouping_list <- c(grouping_list, strsplit(x, "#")[[1]][1])
}
write.csv(grouping_list, file = paste(output_folder,"find_my_friends_gene_ordering_list.csv", sep=''))

for (file in genome_files)
{
  proteins <- read.FASTA(file, type = "AA")
  print(paste('Genome ',file,' proteins are loaded'))
  kmerCounts <- kcount(proteins, k = kmer_length, compress = FALSE)
  print(paste('Genome ',file, ' kmerCounts is finished'))
  sparseKmerCounts <- Matrix(kmerCounts, sparse = TRUE)
  
  file_name <- strsplit(file, ".faa")
  file_name <- strsplit(file_name[[1]][1], "FAA_files/medioid_dataset/")
  file_name <- paste(output_folder, file_name[[1]][2], sep = '')
  out_path <- paste(file_name,paste('_', kmer_string, sep = ''),sep = '')
  out_path <- paste(out_path, '_count_matrix_full_alphabet.csv',sep = '')
  
  protein_out_path <- paste(strsplit(file_name, paste(kmer_string, "_count", sep = ''))[[1]][1], "_protein_list.csv", sep = '')

  setwd(output_folder)
  writeMM(sparseKmerCounts, file=out_path)
  write.csv(rownames(kmerCounts), file = protein_out_path)
  write.csv(colnames(kmerCounts), file = paste(kmer_string, "_list.csv", sep = ''))
  print(paste('Genome ', file, ' matrix is output'))
  print(paste('Genome ',file,' is complete', sep = ''))
}