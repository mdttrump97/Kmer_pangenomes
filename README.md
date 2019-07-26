# **Kmer_pangenomes**

## **Alignment-free Construction of Pan-genomes**
This repository stores code created during my time at the University of Arkansas for Medical Sciences Summer Undergraduate Research Program (UAMS SURF). 

### **What is a pan-genome?**
A pan-genome is the total collection of genes present in a group of organisms. As bacteria can easily transfer and receive genetic material via mechanisms of horizontal gene transfer, 
a diverse set of genes may be present among organisms of the same species. Pan-genomes are composed of core genes (genes present in all genomes analyzed), accessory genes (genes present in some genomes analyzed), 
and singleton genes (genes present in only one genome). These divisions carry biological signficance, with core genes typically being associated with essential housekeeping functions, and accessory and singleton genes
being associated with dispensable, adaptable functions.

### **Alignment-free?**
Alignment-free methods are used to compare sequences without using sequence alignment. Alignment-free methods scale better than sequence alignment methods as input data sizes increase, so development using 
alignment-free methods has increased as more genomes have been sequenced and have become available. 

### **K-mer counting?**
A *k-mer* is a subsequence of length k of a larger sequence. K-mers are recorded by finding all overlapping subsequences by sliding a small window of length k across a sequence. By recording the occurrence of k-mers
in a sequence, information about the original sequence is stored. K-mer occurrences can then be converted to vectors, and multiple vectors can be compared using metrics from linear algebra, such as the cosine similarity or the euclidean distance.

## **Workflow**
GenBank -> Prodigal -> kmerCounter.R -> kmerSelector.py -> canopyClustering.py -> grouping.R

### **Step 1:** Acquiring genomes
Genomes are downloaded using [Batch Entrez](https://www.ncbi.nlm.nih.gov/sites/batchentrez) by uploading a text file containing the accession numbers of the desired assemblies. Choose database "Assembly", and then download genomes
as GenBank Genomic FASTA files (.fna files).

### **Step 2:** Translating genomes
Genomes are translated using [Prodigal](https://github.com/hyattpd/Prodigal). 
`prodigal -i ___.fna -a ___.faa -q`

### **Step 3:** Counting k-mers with kmerCounter.R
This program counts k-mers of a desired length from each protein in each genome and outputs the k-mers, the proteins, and the k-mer counts for each protein

**Input:** 
* .faa files from [Prodigal](https://github.com/hyattpd/Prodigal) of genomes that were downloaded from GenBank.

**Output:** 
* One .mtx file per genome (matrix market format) of k-mer counts for every protein in the genomes
* One .csv file per genome of protein IDs
* One .csv file of k-mers
* An ordered list of proteins to organize clustering output in canopyClustering.py

### **Step 4:** Selecting highly occurrent k-mers with kmerSelector.R
This program reads in the k-mer counts for each protein and chooses a small selection of the most occurrent k-mers to represent each protein. 

**Input:**
* Output from kmerCounter.R

**Output:**
* One large matrix .mtx file (matrix market format) of the k-mer counts for each protein in each genome
* One .csv file of the k-mers that correspond to the matrix rows
* One .csv file of the protein IDs that correspond to the matrix columns

### **Step 5:** Grouping similar genes with canopyClustering.py
This program reads in the large k-mer count matrix and calculates the pairwise distances between every pair of proteins in the matrix. The proteins are then clustered together if their cosine similarity distance 
is less than a given threshold

**Input:**
* Output from kmerSelector.R

**Output:**
* One .csv file recording the cluster membership for each protein 
* One .csv files recording the contents of each cluster

### **Step 6:** Viewing pan-genome graphs with grouping.R
This program reads in the cluster membership for each protein, applies it to each protein, and then divides each cluster into core, accessory, or singleton categories to describe the pan-genomes of the samples. To have updated graphs from FindMyFriends, download the [updated version](https://github.com/mdttrump97/FindMyFriends) and the build and install the package from source using R Studio. 

**Input:**
* Output from canopyClustering.R
* Ordered list of proteins from kmerCounter.R

**Output:**
* No direct file output, but pan-genome graphs can be created.

## Documentation:
Program creation and documentation across 10 weeks is documented in a [jupyter notebook](https://github.com/mdttrump97/UAMS-SURF-Jupyter-Notebook).
