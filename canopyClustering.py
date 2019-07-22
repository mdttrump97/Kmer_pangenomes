import pandas as pd
from sklearn.metrics.pairwise import pairwise_distances
from scipy import io, sparse
import csv

# Get organism ID from protein id 
def get_organism(s):
    head = s.split('.')[0].rstrip('0123456789')
    return head

# If cosine distance less than this threshold, cluster proteins together
clustering_threshold = 0.6

# Format input and output
data_folder = '/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/staph/'
input_description = 'staph_3mer_top_9'
output_description = input_description + '_' + str(clustering_threshold)

# Read in k-mer count matrix (output from kmerSelector.py)
sparse_matrix = io.mmread(data_folder + input_description + '.mtx')
print("Converting to csr_matrix")
df = sparse.csr_matrix(sparse_matrix)

# Proteins are index, kmers are columns
df = df.transpose()

# Read in proteins (output from kmerSelector.py)
proteins = pd.read_csv(data_folder + input_description + '_protein_list.csv')

# Calculate and format distance matrix
print("Pairwise distance calculation")
pairwise_distances = pd.DataFrame(pairwise_distances(df, metric='cosine'))
pairwise_distances.columns = proteins.columns
pairwise_distances.index = proteins.columns

elligible = list(pairwise_distances.columns)

# Make a list of precluster centers
# Make a sub dataframe with 100 proteins
subset_columns = elligible[0:101]
# Reduce the distance matrix to the 100 proteins
subset_distances = pairwise_distances[subset_columns]
# Isolate the 100% duplicate rows across 100 proteins
subset_duplicates = subset_distances[subset_distances.duplicated(keep = False)]
# List of 100% identical proteins
eligible_proteins = list(subset_duplicates.index)
# Reduce duplicate rows to unique points 
cluster_centers = subset_duplicates.drop_duplicates()

used_points = []
preclusters = dict()
cluster_centers = list(cluster_centers.index)
i = 1
print("Preclustering exact duplicates")
while len(cluster_centers) != 0:
    cluster_center = cluster_centers.pop()
    sorted_distances = pairwise_distances.loc[cluster_center].sort_values()
    available_points = sorted_distances[~sorted_distances.index.isin(used_points)]
    identical_points = list(available_points[available_points == 0].index)
    points_under_threshold = list(available_points[available_points < clustering_threshold].index)
    used_points = used_points + points_under_threshold
    used_points.append(cluster_center)
    used_points = list(set(used_points))
    cluster_centers = list(set(cluster_centers) - set(used_points))
    preclusters['clus_' + str(i)] = points_under_threshold
    i+=1

points_not_in_precluster = list(set(elligible) - set(used_points))

# Store these variables for reference
points_without_duplicates = points_not_in_precluster
points_with_duplicates = used_points

clusters = dict()
possible_cluster_centers = points_not_in_precluster
print("Clustering non-duplicate genes")
while len(possible_cluster_centers) != 0:
    cluster_center = possible_cluster_centers.pop()
    sorted_distances = pairwise_distances.loc[cluster_center].sort_values()
    available_points = sorted_distances[~sorted_distances.index.isin(used_points)]
    points_under_clustering_threshold = list(available_points[available_points < clustering_threshold].index)
    used_points = used_points + points_under_clustering_threshold
    #used_points.append(cluster_center)
    used_points = list(set(used_points))
    possible_cluster_centers = list(set(possible_cluster_centers) - set(used_points))
    clusters['clus_' + str(i)] = points_under_clustering_threshold
    i+=1

# Remove empty clusters
non_empty_clusters = dict()
for cluster_center in list(clusters.keys()):
    cluster = clusters[cluster_center]
    if len(cluster) > 0:
        non_empty_clusters[cluster_center] = cluster

total_clusters = {**preclusters, **non_empty_clusters}

'''
organism_count = dict()
for cluster_center in total_clusters.keys():
    cluster = total_clusters[cluster_center]
    org_list = [get_organism(protein) for protein in cluster]
    org_list = list(set(org_list))
    organism_count[cluster_center] = pd.Series(org_list).value_counts()
    
organism_count = pd.DataFrame(organism_count)

df_100 = organism_count.dropna(axis = 1, thresh = 14)
df_93 = organism_count.dropna(axis = 1, thresh = 13)
df_85 = organism_count.dropna(axis = 1, thresh = 12) 
'''

grouping_list = dict()
i = 1

# Reformat clusters for output
print("Formatting clusters")
for cluster_center in total_clusters.keys():
    formatted_clusters = dict()
    formatted_clusters['c'] = i
    formatted_clusters['points'] = total_clusters[cluster_center]
    clusters[i] = formatted_clusters
    for protein in formatted_clusters['points']:
        grouping_list[protein] = i
    i += 1

# Get ordering of genes (output from kmerCounter.R) to assist with data import into FindMyFriends (grouping.R)
ordered_gene_list = list(pd.read_csv(data_folder + 'find_my_friends_gene_ordering_list.csv')['x'])

# Create list of cluster membership for each protein for output
out_list = []
for entry in ordered_gene_list:
    out_list.append(str(grouping_list[str(entry).strip()]))
df_out = pd.DataFrame()
df_out['gene'] = ordered_gene_list
df_out['clust'] = out_list

# Write cluster membership for each protein (input for grouping.R)
df_out.to_csv(data_folder + output_description + '_grouping_list.csv', index = None)

# Write clusters out for comparison with UCLUST
with open(data_folder + output_description + '_clusters.txt', 'w') as csv_file:
    writer = csv.writer(csv_file)
    for key, value in total_clusters.items():
        writer.writerow([key, value])