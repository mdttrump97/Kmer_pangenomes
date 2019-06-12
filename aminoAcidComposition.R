library("protr")
library("ape")
library("fmsb")


#file_list <- c("GCA_000005845.2_ASM584v2","GCA_000006665.1_ASM666v1","GCA_000006925.2_ASM692v2","GCA_000007445.1_ASM744v1","GCA_000008865.2_ASM886v2","GCA_000010245.1_ASM1024v1","GCA_000010485.1_ASM1048v1","GCA_000012005.1_ASM1200v1","GCA_000012025.1_ASM1202v1","GCA_000019645.1_ASM1964v1","GCA_000026325.2_ASM2632v2","GCA_000091005.1_ASM9100v1","GCA_000092525.1_ASM9252v1","GCA_000167875.1_ASM16787v1","GCA_000176575.2_ASM17657v2","GCA_000176655.2_ASM17665v2","GCA_000176695.2_ASM17669v2","GCA_000333215.1_ASM33321v1","GCA_000350785.1_Esch_coli_KTE21_V1","GCA_000613265.1_ASM61326v1","GCA_000690815.1_ASM69081v1","GCA_000714595.1_ASM71459v1","GCA_000734955.1_ASM73495v1","GCA_002007705.1_ASM200770v1","GCA_003697165.2","GCA_900092615.1_PRJEB14041")
file_list <- c("GCA_000005845.2_ASM584v2","GCA_000006665.1_ASM666v1","GCA_000006925.2_ASM692v2","GCA_000007445.1_ASM744v1","GCA_000008865.2_ASM886v2")
for (file in file_list)
{
  setwd("/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/FAA_files")
  file_path <- paste(file,'.faa',sep = '')
  
  x <- readFASTA(file_path)[]
  concatenated_proteins <- ""
  for(protein in x){
    concatenated_proteins <- paste(concatenated_proteins, protein, sep = "")
  }
  
  concatenated_proteins <- gsub("[^AVYWTSPFMKLIHGQECDNR]", "", concatenated_proteins)
  concatenated_proteins[0:10]
  aac <- extractAAC(concatenated_proteins)
  aac <- data.frame(t(as.matrix(aac)))
  aac <- rbind(rep(max(aac), 20), rep(min(aac),20), aac)
  radarchart(aac, title = file_path)
}