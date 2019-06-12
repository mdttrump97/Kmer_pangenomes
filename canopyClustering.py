import pandas as pd
import csv
from sklearn.metrics.pairwise import pairwise_distances
from scipy import io, sparse

# X shoudl be a numpy matrix, very likely sparse matrix: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.sparse.csr_matrix.html#scipy.sparse.csr_matrix
# T1 > T2 for overlapping clusters
# T1 = Distance to centroid point to not include in other clusters
# T2 = Distance to centroid point to include in cluster
# T1 > T2 for overlapping clusters
# T1 < T2 will have points which reside in no clusters
# T1 == T2 will cause all points to reside in mutually exclusive clusters
# Distance metric can be any from here: http://scikit-learn.org/stable/modules/generated/sklearn.metrics.pairwise.pairwise_distances.html
# filemap may be a list of point names in their order in X. If included, row numbers from X will be replaced with names from filemap. 
 
def main():
    # read in and format data from .csv file
    #df = pd.read_csv('/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/10_genome_4mer_top_3_table_full_alphabet.csv')
    sparse_matrix = io.mmread('/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/10_genome_4mer_counts_sparse/4mer_top_10_table_full_alphabet_compressed')
    #df = pd.DataFrame(sparse_matrix.toarray()) 
    print("Converting to csr_matrix")
    df = sparse.csr_matrix(sparse_matrix)
    
    '''
    kmers = pd.read_csv('/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/10_genome_4mer_counts_sparse/top_10_kmers.csv')
    
    df.index = kmers.columns
    df.columns = proteins.columns
    '''
    # proteins are index, kmers are columns
    df = df.transpose()
    
    proteins = pd.read_csv('/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/10_genome_4mer_counts_sparse/top_10_protein_list.csv')
    # calculate and format distance matrix
    print("Pairwise distance calculation")
    X1_dist = pd.DataFrame(pairwise_distances(df, metric='cosine'))
    X1_dist.columns = proteins.columns
    X1_dist.index = proteins.columns
    
    # distance threshold for clustering proteins together
    threshold = 0.5
    
    # set up and format data structures for algorithm
    canopies = dict()
    elligible_points = list(X1_dist.columns)
    points_in_threshold = []
    clustered_points = []
    
    # set up and format final grouping_list 
    grouping_list = pd.DataFrame(list([None] * X1_dist.shape[0]))
    grouping_list = grouping_list.T
    grouping_list.columns = X1_dist.columns

    iteration = 0
    while len(elligible_points) > 0:
        iteration = iteration + 1
        print("iteration: " + str(iteration))
        # record which organisms' proteins have been clustered into this cluster
        # only one protein per organism is allowed in this cluster
        clustered_organisms = []

        center_point = list(elligible_points)[-1]
        i = len(canopies)
        
        for point in elligible_points:
            if(X1_dist[center_point][point] < threshold):
                # check for organisms that have already been clustered
                if(point.split('_')[0] not in clustered_organisms):
                    points_in_threshold.append(point) 
                    clustered_organisms.append(point.split('_')[0])
                    
        # cluster points which have not already been clustered
        points_to_cluster = set(points_in_threshold).difference(set(clustered_points))
        clustered_points.extend(points_to_cluster)
        elligible_points = set(elligible_points).difference(set(points_to_cluster))
        
        canopies[i] = {"c":iteration, "points": points_to_cluster}
        
        for entry in canopies[i]["points"]:
            grouping_list[entry] = iteration
    
    with open('/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/10_genome_top_10_4mer_full_alphabet_cosine_clusters_0.5_ps.csv', 'w') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in canopies.items():
            writer.writerow([key, value])
    
    with open('/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/10_genome__4mer_full_alphabet_cosine_grouping_list_0.5_ps.csv', 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(list(grouping_list.iloc[0])) 
    
if __name__ == "__main__":
    main()