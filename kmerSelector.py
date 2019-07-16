from scipy import io, sparse
import pandas as pd
import glob
import csv

# get kmers from "kmer_list" containing any amino acid listed in "amino_acid_subset"
def get_kmer_subset(amino_acid_subset, kmer_list):
    kmer_subset = list()
    for kmer in kmer_list:
        if any(amino_acid in kmer for amino_acid in amino_acid_subset):
            kmer_subset.append(kmer)
    return kmer_subset

# sort kmers by occurrence, then alphabetically once occurrences tie
def sort_kmer_counts(kmer_counts):
    kmer_count_tuple = list(zip(kmer_counts.index,kmer_counts))
    sorted_kmer_counts = sorted(kmer_count_tuple, key=lambda item: (-item[1], item[0]))

    sorted_kmer_counts = pd.DataFrame(sorted_kmer_counts)
    sorted_kmer_counts.index = sorted_kmer_counts[0]
    del sorted_kmer_counts[0]
    return sorted_kmer_counts

# append "top_x" kmers containing amino acids from "amino_acid_subset" to "protein_vector"
def append_kmer_subset(amino_acid_subset, kmer_counts, protein_vector, top_x):
    selected_kmers = get_kmer_subset(amino_acid_subset, kmer_counts.index)
    #print("Selected kmers: " + str(len(selected_kmers)))
    if(len(selected_kmers) > top_x):
        selected_kmer_counts = kmer_counts.loc[selected_kmers,]
    
        for kmer in selected_kmer_counts.index:
            if kmer in protein_vector.index:
                selected_kmer_counts = selected_kmer_counts.drop(kmer)
        
        sorted_selected_kmer_counts = sort_kmer_counts(selected_kmer_counts)[1]
        
        #print("Top 3 most occurrent 5-mers containing C or M: ")
        #print(str(sorted_selected_kmer_counts[0:3]))
        #print()
        # how many top occuring selected kmers?
        protein_vector = protein_vector.append(sorted_selected_kmer_counts[0:top_x])
    return protein_vector

# Format input and output
kmer_length = '3'
data_folder = '/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/medioid_3mers/'
output_description = "medioid_3mer_test"

# How many kmers to represent each gene?
top_x_most_occurrent = 6

# Get kmers from file
kmer_list = pd.read_csv(data_folder + kmer_length + "mer_list.csv")
del kmer_list["Unnamed: 0"]
kmer_list = kmer_list["x"]
kmer_dict = dict(zip(list(range(0, len(kmer_list))), list(kmer_list)))

kmer_count_dictionary = dict()

file_list = glob.glob(data_folder + '*' + kmer_length + 'mer_count_matrix_full_alphabet*')
for x in range(len(file_list)):
    print("Iteration: " + str(x))
    file = file_list[x]
    print("Reading in " + str(file))
    sparse_matrix = io.mmread(str(file))
    df = pd.DataFrame(sparse_matrix.toarray())   
    df = df.T
    
    # Get proteins from file
    print("Reading in protein list")
    protein_list = pd.read_csv(file.split("_" + kmer_length + "mer_count_")[0] + "_protein_list.csv")
    del protein_list["Unnamed: 0"]
    protein_list = protein_list["x"]
    protein_list = [x.split(' ')[0] for x in protein_list]
    protein_dict = dict(zip(list(range(0, len(protein_list))), list(protein_list)))
    
    # Add kmers and proteins to dataframe
    df.index = df.index.to_series().map(kmer_dict)
    df.columns = df.columns.to_series().map(protein_dict)
    
    print("Generating protein vectors")
    for protein in list(df.columns):
        nonzero_kmer_df = df[df[protein] != 0]
        kmer_counts = nonzero_kmer_df[protein]
        
        kmer_count_tuple = list(zip(kmer_counts.index,kmer_counts))
        sorted_kmer_counts = sort_kmer_counts(kmer_counts)[1]
        
        # How many top occuring kmers?
        protein_vector = sorted_kmer_counts[0:top_x_most_occurrent]
        #print(protein + str(":"))
        #print("Top 3 most occurent 4-mers:")
        #print(protein_vector)
        #print()
        
        #protein_vector = append_kmer_subset(["F","Y","W","H"], kmer_counts, protein_vector, top_x_most_occurrent)
        #protein_vector = append_kmer_subset(["C","M"], kmer_counts, protein_vector, top_x_most_occurrent)
        
        kmer_count_dictionary[protein] = protein_vector
        #print("Final protein vector:")
        #print(protein_vector)
        #print("------------------------------------------------")

print("Converting to dataframe")
top_x_most_occurrent_df = pd.DataFrame(kmer_count_dictionary)
top_x_most_occurrent_df = top_x_most_occurrent_df.fillna(0)

print("Converting to matrix market format")
top_x_most_occurrent_coo = sparse.coo_matrix(top_x_most_occurrent_df)

print("Writing to file")
io.mmwrite(data_folder + output_description + ".mtx", top_x_most_occurrent_coo)
with open(data_folder + output_description + '_kmers.csv', 'w') as writeFile:
    writer = csv.writer(writeFile)
    writer.writerow(top_x_most_occurrent_df.index)
with open(data_folder + output_description + '_protein_list.csv', 'w') as writeFile:
    writer = csv.writer(writeFile)
    writer.writerow(top_x_most_occurrent_df.columns)