import pandas as pd
import csv
from sklearn.metrics.pairwise import pairwise_distances
from scipy import io, sparse
import re

def get_nonzero_minimum(distances, point, threshold):
    distances = pd.DataFrame(distances)
    distances = distances.loc[distances[point] < threshold]
    distances = distances[point].sort_values(ascending=True)
    
    for distance in distances:
        if distance > 0:
            return distance
    return 1
        
threshold = 0.675
data_folder = '/Users/matthewthompson/Documents/UAMS_SURF/K-mer_testing/CSV_files/medioid_3mers/'
input_description = 'medioid_3mer_top_9'
output_description = 'medioid_3mer_top_9_0.675'
    
# read in and format data from .csv file
sparse_matrix = io.mmread(data_folder + input_description + '.mtx')
#df = pd.DataFrame(sparse_matrix.toarray()) 
print("Converting to csr_matrix")
df = sparse.csr_matrix(sparse_matrix)

# proteins are index, kmers are columns
df = df.transpose()

proteins = pd.read_csv(data_folder + input_description + '_protein_list.csv')
# calculate and format distance matrix
print("Pairwise distance calculation")
pairwise_distance_matrix = pd.DataFrame(pairwise_distances(df, metric='cosine'))
pairwise_distance_matrix.columns = proteins.columns
pairwise_distance_matrix.index = proteins.columns

# set up and format data structures for algorithm
canopies = dict()

# set up and format final grouping_list 
grouping_list = pd.DataFrame(list([None] * pairwise_distance_matrix.shape[0]))
grouping_list = grouping_list.T
grouping_list.columns = pairwise_distance_matrix.columns

preclusters = dict()
elligible_points = list(proteins.columns)
preclustered_points = []
while len(elligible_points) > 0:
    points_in_precluster = []
    precluster_center = elligible_points.pop()
    print("Center: " + str(precluster_center))
    for point in elligible_points:
        if(pairwise_distance_matrix[precluster_center][point] == 0):
            points_in_precluster.append(point)
            preclustered_points.append(point)
            elligible_points.remove(point)
    points_in_precluster.append(precluster_center)
    preclusters[precluster_center] = points_in_precluster
    print("Points still available: " + str(len(elligible_points)))
    print("----------------------------------------------------------")

elligible_points = list(preclusters.keys())
canopies = dict()
reduced_distance_matrix = pairwise_distance_matrix.drop(list(preclustered_points))
reduced_distance_matrix = reduced_distance_matrix.drop(preclustered_points, axis = 1)

print("Finding nonzero mininum distances")
nonzero_min_distances = dict()
for point in elligible_points:
    nonzero_min_distances[point] = get_nonzero_minimum(reduced_distance_matrix, point, threshold)
            
iteration = 0
while len(elligible_points) > 0:
    #print(elligible_points[0:10])
    points_in_threshold = []
    # record which organisms' proteins have been clustered into this cluster
    # only one protein per organism is allowed in this cluster
    #clustered_organisms = []

    center_point = elligible_points.pop()
    print("Points still available: " + str(len(elligible_points)))
    print("----------------------------------------------------------")
    
    for point in elligible_points:
        distance_to_cluster_center = reduced_distance_matrix[center_point][point]
        
        nonzero_min_distance = nonzero_min_distances[point]
        if (distance_to_cluster_center == nonzero_min_distance) & (distance_to_cluster_center < threshold):
            points_in_threshold.append(point)
    #reduced_distance_matrix.drop(points_in_threshold, inplace = True)
    #reduced_distance_matrix.drop(points_in_threshold, axis = 1, inplace = True)
        
    '''
    if(pairwise_distance_matrix[center_point][point] < threshold):
        #organism = re.split("[0-9]", point)[0]
        # check for organisms that have already been clustered
        #if(organism not in clustered_organisms):
        points_in_threshold.append(point) 
        #clustered_organisms.append(organism)    
    '''
    points_in_threshold.append(center_point)
    print("points in threshold: " + str(points_in_threshold))

    points_in_canopy = []
    for point in points_in_threshold:
        points_in_canopy.extend(preclusters[point])
    
    canopies[iteration] = {"c":iteration, "points": points_in_canopy}
    
    points_in_canopy = list(set(points_in_canopy))
    
    for entry in canopies[iteration]["points"]:
        grouping_list[entry] = iteration
        
    elligible_points = [x for x in elligible_points if x not in points_in_threshold]
    iteration = iteration + 1

with open(data_folder + output_description + '_clusters.txt', 'w') as csv_file:
    writer = csv.writer(csv_file)
    for key, value in canopies.items():
        writer.writerow([key, value])

with open(data_folder + output_description + '_grouping_list.csv', 'w') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(list(grouping_list.iloc[0])) 