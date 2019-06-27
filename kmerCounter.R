library("ape")
library("kmer")
library("MASS")
library("Matrix")
setwd("/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/FAA_files/medioid_dataset")

genome_files <- list.files(getwd(), full.names=TRUE, pattern='*.faa')

for (file in genome_files)
{
  file_name <- strsplit(file, ".faa")[1]
  out_path <- paste(file,'_4mer_count_matrix_full_alphabet.csv',sep = '')
  proteins <- read.FASTA(file, type = "AA")
  print(paste('Genome ',file,' proteins are loaded'))
  kmerCounts <- kcount(proteins, k = 4, compress = FALSE)
  print(paste('Genome ',file, ' kmerCounts is finished'))
  sparseKmerCounts <- Matrix(kmerCounts, sparse = TRUE)
  setwd("/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/medioid_4mers")
  writeMM(sparseKmerCounts, file=out_path)
  write.csv(rownames(kmerCounts), file = paste(file, "_protein_list.csv"), sep = '')
  write.csv(colnames(kmerCounts), file = "4mer_list.csv")
  print(paste('Genome ', file, ' matrix is output'))
  print(paste('Genome ',file,' is complete', sep = ''))
}