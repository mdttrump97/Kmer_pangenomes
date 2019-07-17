import pandas as pd
from sklearn.metrics.pairwise import pairwise_distances
from scipy import io, sparse
import csv

def get_organism(s):
    head = s.split('.')[0].rstrip('0123456789')
    return head

threshold = 0.5
data_folder = '/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/medioid_3mers/'
input_description = 'medioid_3mer_top_9'
output_description = 'medioid_3mer_top_9_0.675'

# read in and format data from .csv file
#df = pd.read_csv('/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/10_genome_4mer_top_3_table_full_alphabet.csv')
sparse_matrix = io.mmread(data_folder + input_description + '.mtx')
#df = pd.DataFrame(sparse_matrix.toarray()) 
print("Converting to csr_matrix")
df = sparse.csr_matrix(sparse_matrix)

# proteineins are index, kmers are columns
df = df.transpose()

proteins = pd.read_csv(data_folder + input_description + '_protein_list.csv')
# calculate and format distance matrix
print("Pairwise distance calculation")
pairwise_distances = pd.DataFrame(pairwise_distances(df, metric='cosine'))
pairwise_distances.columns = proteins.columns
pairwise_distances.index = proteins.columns

elligible = list(pairwise_distances.columns)

# make a sub dataframe with 100 proteineins
subset_columns = elligible[0:101]
# reduce the dist mat to the 100 proteineins
subset_distances = pairwise_distances[subset_columns]
# isolate the 100% duplicate rows across 100 proteineins
subset_duplicates = subset_distances[subset_distances.duplicated(keep = False)]
# list of 100% identical proteineins
eligible_proteins = list(subset_duplicates.index)
                            
# reduce duplicate rows to uniques 
# basically list(set(list(meh)))
cluster_centers = subset_duplicates.drop_duplicates()

used_points = []
preclusters = dict()
cluster_centers = list(cluster_centers.index)
i = 1
print("Preclustering exact duplicates")
while len(cluster_centers) != 0:
    protein = cluster_centers.pop()
    sorted_distances = pairwise_distances.loc[protein].sort_values()
    available_points = sorted_distances[~sorted_distances.index.isin(used_points)]
    identical_points = list(available_points[available_points == 0].index)
    points_under_threshold = list(available_points[available_points < threshold].index)
    if(len(points_under_threshold) == 0):
        print("i: " + str(i))
        print("Center point: " + str(protein))
        print("Identical points: ")
        print(identical_points)
        print(protein in used_points)
    non_identical_points = list(set(points_under_threshold) - set(identical_points))
    used_points = used_points + points_under_threshold
    used_points.append(protein)
    used_points = list(set(used_points))
    cluster_centers = list(set(cluster_centers) - set(used_points))
    preclusters['clus_' + str(i)] = points_under_threshold
    i+=1

points_not_in_precluster = list(set(elligible) - set(used_points))
test2 = points_not_in_precluster
test3 = used_points
clusters = dict()
#i = 1
print("Clustering non-duplicate genes")
while len(points_not_in_precluster) != 0:
    protein = points_not_in_precluster.pop()
    sorted_distances = pairwise_distances.loc[protein].sort_values()
    available_points = sorted_distances[~sorted_distances.index.isin(used_points)]
    #identical = list(loop_series[loop_series == 0].index)
    points_under_threshold = list(available_points[available_points < threshold].index)
    #non_identical = list(set(merge) - set(identical))
    used_points = used_points + points_under_threshold
    #used_points.append(protein)
    used_points = list(set(used_points))
    points_not_in_precluster = list(set(points_not_in_precluster) - set(used_points))
    clusters['clus_' + str(i)] = points_under_threshold
    i+=1

non_empty_clusters = dict()
for cluster_center in list(clusters.keys()):
    cluster = clusters[cluster_center]
    if len(cluster) > 0:
        non_empty_clusters[cluster_center] = cluster

results = {**preclusters, **non_empty_clusters}

'''
organism_count = dict()
for cluster_center in results.keys():
    cluster = results[cluster_center]
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

print("Formatting clusters")
for cluster_center in results.keys():
    formatted_clusters = dict()
    formatted_clusters['c'] = i
    formatted_clusters['points'] = results[cluster_center]
    clusters[i] = formatted_clusters
    for protein in formatted_clusters['points']:
        grouping_list[protein] = i
    i += 1

ordered_gene_list = list(pd.read_csv(data_folder + 'find_my_friends_gene_ordering_list.csv')['x'])

out_list = []
for entry in ordered_gene_list:
    out_list.append(str(grouping_list[str(entry).strip()]))

df_out = pd.DataFrame()
df_out['gene'] = ordered_gene_list
df_out['clust'] = out_list
df_out.to_csv(data_folder + output_description + '_grouping_list2.csv', index = None)

with open(data_folder + output_description + '_clusters2.txt', 'w') as csv_file:
    writer = csv.writer(csv_file)
    for key, value in results.items():
        writer.writerow([key, value])