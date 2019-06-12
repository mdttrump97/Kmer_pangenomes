from scipy import io, sparse
import pandas as pd
import glob
import csv

kmer_list = pd.read_csv("/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/10_genome_4mer_counts_sparse/4mer_list.csv")
del kmer_list["Unnamed: 0"]
kmer_list = kmer_list["x"]
kmer_dict = dict(zip(list(range(0, len(kmer_list))), list(kmer_list)))

data_folder = '/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/10_genome_4mer_counts_sparse/'
file_list = glob.glob(data_folder + '*4mer_count_matrix_full_alphabet*')
iteration = 0
dicty = dict()
for file in file_list: 
    iteration = iteration + 1
    print("Reading in " + str(file))
    sparse_matrix = io.mmread(str(file))
    df = pd.DataFrame(sparse_matrix.toarray())   
    df = df.T
    
    print("Reading in protein list")
    protein_list = pd.read_csv(file.split("4mer_count_")[0] + "protein_list.csv")
    del protein_list["Unnamed: 0"]
    protein_list = protein_list["x"]
    protein_list = [x.split(' ')[0] for x in protein_list]
    protein_dict = dict(zip(list(range(0, len(protein_list))), list(protein_list)))
    
    df.index = df.index.to_series().map(kmer_dict)
    df.columns = df.columns.to_series().map(protein_dict)
    
    print("Getting top 10 kmers")

    for col in list(df.columns):
        loop = df[df[col] != 0]
        loop = loop[col]
        loop = loop.sort_values(ascending = False)
        # how many top occuring kmers?
        loop = loop[0:5]
        dicty[col] = loop

top_x_df = pd.DataFrame(dicty)
top_x_df = top_x_df.fillna(0)
top_x_coo = sparse.coo_matrix(top_x_df)
print("Writing to file")
io.mmwrite('/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/10_genome_4mer_counts_sparse/4mer_top_10_table_full_alphabet_compressed.mtx', top_x_coo)
with open('/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/10_genome_4mer_counts_sparse/top_10_kmers.csv', 'w') as writeFile:
    writer = csv.writer(writeFile)
    writer.writerow(top_x_df.index)
with open('/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/10_genome_4mer_counts_sparse/top_10_protein_list.csv', 'w') as writeFile:
    writer = csv.writer(writeFile)
    writer.writerow(top_x_df.columns)