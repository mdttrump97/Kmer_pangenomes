from scipy import io, sparse
import pandas as pd
import glob
import csv

def get_kmer_subset(amino_acid_subset, kmer_list):
    kmer_subset = list()

    for kmer in kmer_list:
        if any(amino_acid in kmer for amino_acid in amino_acid_subset):
            kmer_subset.append(kmer)
    return kmer_subset

def sort_kmer_counts(kmer_counts):
    kmer_count_tuple = list(zip(kmer_counts.index,kmer_counts))
    sorted_kmer_counts = sorted(kmer_count_tuple, key=lambda item: (-item[1], item[0]))

    sorted_kmer_counts = pd.DataFrame(sorted_kmer_counts)
    sorted_kmer_counts.index = sorted_kmer_counts[0]
    del sorted_kmer_counts[0]
    return sorted_kmer_counts

def append_kmer_subset(amino_acid_subset, kmer_counts, protein_vector, top_x):
    selected_kmers = get_kmer_subset(amino_acid_subset, kmer_counts.index)
    #print("Selected kmers: " + str(len(selected_kmers)))
    if(len(selected_kmers) > top_x):
        selected_kmer_counts = kmer_counts.loc[selected_kmers,]
    
        for kmer in selected_kmer_counts.index:
            if kmer in protein_vector.index:
                selected_kmer_counts = selected_kmer_counts.drop(kmer)
        
        sorted_selected_kmer_counts = sort_kmer_counts(selected_kmer_counts)[1]
        
        #print("selected kmers: ")
        #print(str(sorted_selected_kmer_counts[0:3]))
        # how many top occuring selected kmers?
        protein_vector = protein_vector.append(sorted_selected_kmer_counts[0:top_x])
    return protein_vector

kmer_list = pd.read_csv("/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/10_genome_4mer_counts_sparse/4mer_list.csv")
del kmer_list["Unnamed: 0"]
kmer_list = kmer_list["x"]
kmer_dict = dict(zip(list(range(0, len(kmer_list))), list(kmer_list)))

data_folder = '/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/10_genome_4mer_counts_sparse/'
file_list = glob.glob(data_folder + '*4mer_count_matrix_full_alphabet*')
dicty = dict()
for file in file_list: 
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
    
    print("Generating protein vectors")
    
    top_x = 3
    for protein in list(df.columns):
        nonzero_kmer_df = df[df[protein] != 0]
        kmer_counts = nonzero_kmer_df[protein]
        
        kmer_count_tuple = list(zip(kmer_counts.index,kmer_counts))
        sorted_kmer_counts = sort_kmer_counts(kmer_counts)[1]
        
        # how many top occuring kmers?
        protein_vector = sorted_kmer_counts[0:top_x]
        #print("top 3 vector:")
        #print(protein_vector)
        
        #protein_vector = append_kmer_subset(["F","Y","W","H"], kmer_counts, protein_vector, top_x)
        #print("top 3 + aromatics:")
        #print(protein_vector)
        protein_vector = append_kmer_subset(["C","M"], kmer_counts, protein_vector, top_x)
        
        dicty[protein] = protein_vector
        #print("top 3 + aromatics + CM")
        #print(protein_vector)
        #print(protein)
        #print("------------------------------------------------")

print("Converting to dataframe")
top_x_df = pd.DataFrame(dicty)
top_x_df = top_x_df.fillna(0)
print("Converting to matrix market format")
top_x_coo = sparse.coo_matrix(top_x_df)
print("Writing to file")
io.mmwrite('/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/10_genome_4mer_counts_sparse/4mer_top_3_all_CM_table.mtx', top_x_coo)
with open('/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/10_genome_4mer_counts_sparse/top_3_all_CM_kmers.csv', 'w') as writeFile:
    writer = csv.writer(writeFile)
    writer.writerow(top_x_df.index)
with open('/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/10_genome_4mer_counts_sparse/top_3_all_CM_protein_list.csv', 'w') as writeFile:
    writer = csv.writer(writeFile)
    writer.writerow(top_x_df.columns)