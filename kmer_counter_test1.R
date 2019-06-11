library("ape")
library("kmer")
library("MASS")
library("Matrix")
setwd("/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/FAA_files")


file_list <- c("GCA_000005845.2_ASM584v2","GCA_000006665.1_ASM666v1","GCA_000006925.2_ASM692v2","GCA_000007445.1_ASM744v1","GCA_000008865.2_ASM886v2","GCA_000010245.1_ASM1024v1","GCA_000010485.1_ASM1048v1","GCA_000012005.1_ASM1200v1","GCA_000012025.1_ASM1202v1","GCA_000019645.1_ASM1964v1","GCA_000026325.2_ASM2632v2","GCA_000091005.1_ASM9100v1","GCA_000092525.1_ASM9252v1","GCA_000167875.1_ASM16787v1","GCA_000176575.2_ASM17657v2","GCA_000176655.2_ASM17665v2","GCA_000176695.2_ASM17669v2","GCA_000333215.1_ASM33321v1","GCA_000350785.1_Esch_coli_KTE21_V1","GCA_000613265.1_ASM61326v1","GCA_000690815.1_ASM69081v1","GCA_000714595.1_ASM71459v1","GCA_000734955.1_ASM73495v1","GCA_002007705.1_ASM200770v1","GCA_003697165.2","GCA_900092615.1_PRJEB14041")

for (file in file_list)
{
  setwd("/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/FAA_files")
  file_path <- paste(file,'.faa',sep = '')
  out_path <- paste(file,'_5mer_count_matrix_full_alphabet.csv',sep = '')
  proteins <- read.FASTA(file_path, type = "AA")
  print(paste('Genome ',file,' proteins are loaded'))
  kmerCounts <- kcount(proteins, k = 5, compress = FALSE)
  print(paste('Genome ',file, ' kmerCounts is finished'))
  sparseKmerCounts <- Matrix(kmerCounts, sparse = TRUE)
  setwd("/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/5mer_counts_sparse")
  writeMM(sparseKmerCounts, file=out_path)
  write.csv(rownames(kmerCounts), file = paste(file, "_protein_list.csv"), sep = '')
  write.csv(colnames(kmerCounts), file = "5mer_list.csv")
  print(paste('Genome ', file, ' matrix is output'))
  print(paste('Genome ',file,' is complete', sep = ''))
}