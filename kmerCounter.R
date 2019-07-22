library("ape")
library("kmer")
library("Matrix")
library("FindMyFriends")

# Format input and output
input_folder = "/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/FAA_files/phylotypeA/"
folder_name = "FAA_files/phylotypeA/"
output_folder = "/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/phylotypeA/"

# Set k-mer length
kmer_length = 3
kmer_string = paste(toString(kmer_length), "mer", sep = '')

# Read in .faa files after running Prodigal on genomic FASTA (.fna) files
setwd(input_folder)
genome_files <- list.files(getwd(), full.names=TRUE, pattern='*.faa')

# Output ordering of genes from FindMyFriends to combine with canopyClustering.py output
pan <- pangenome(genome_files[1:length(genome_files)], translated=TRUE, geneLocation='prodigal', lowMem=FALSE)
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
  
  # Count k-mers in each protein sequence
  kmerCounts <- kcount(proteins, k = kmer_length, compress = FALSE)
  print(paste('Genome ',file, ' kmerCounts is finished'))
  # Convert to a sparse storage format to reduce file size when storing
  sparseKmerCounts <- Matrix(kmerCounts, sparse = TRUE)
  
  # Format output file name
  file_name <- strsplit(file, ".faa")
  file_name <- strsplit(file_name[[1]][1], folder_name)
  file_name <- paste(output_folder, file_name[[1]][2], sep = '')
  out_path <- paste(file_name,paste('_', kmer_string, sep = ''),sep = '')
  out_path <- paste(out_path, '_count_matrix_full_alphabet.mtx',sep = '')
  
  protein_out_path <- paste(strsplit(file_name, paste(kmer_string, "_count", sep = ''))[[1]][1], "_protein_list.csv", sep = '')

  setwd(output_folder)
  # Output to be used in kmerSelector.py
  # Write matrix out in matrix market format
  writeMM(sparseKmerCounts, file=out_path)
  write.csv(rownames(kmerCounts), file = protein_out_path)
  write.csv(colnames(kmerCounts), file = paste(kmer_string, "_list.csv", sep = ''))
  print(paste('Genome ', file, ' matrix is output'))
  print(paste('Genome ',file,' is complete', sep = ''))
}