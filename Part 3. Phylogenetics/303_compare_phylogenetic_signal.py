print("\n\n\t**Running Script: 303_compare_phylogenetic_signal.py**\n\n")

#---------------------------------------------------------------------------
# 1. import modules
#---------------------------------------------------------------------------
from pathlib import Path
import itertools
import pandas as pd
import numpy as np
from skbio.stats.distance import DistanceMatrix
from skbio.stats.distance import mantel
import matplotlib.pyplot as plt
import seaborn as sns


#---------------------------------------------------------------------------
# 2. define file paths
#---------------------------------------------------------------------------
#input
p3_folder = Path(__file__).resolve().parent
PA_folder = p3_folder.parent
p3_outfolder = p3_folder / "output_folder"
distance_matrix_folder = p3_outfolder / "distance matrices used to calcualte trees"


#verify nput folders were found
folders_to_check = [p3_outfolder, distance_matrix_folder]
for folder_path in folders_to_check:
    if not folder_path.is_dir():
        print(f"\n\nfolder not found: {folder_path}\n\n")
        raise FileNotFoundError

#output
phylogenetic_signal_folder = p3_outfolder / "phylogenetic signal comparison"
phylogenetic_signal_folder.mkdir(parents=True, exist_ok=True)


#---------------------------------------------------------------------------
# 3. helper functiont o laod a distance matrix
#---------------------------------------------------------------------------
def load_distance_matrix(distance_matrix_file):
    names = []
    matrix = []
    with open(distance_matrix_file) as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                float(parts[1])
            except ValueError:
                continue
            names.append(parts[0])
            current_distances = [float(x) for x in parts[1: ]]
            matrix.append(current_distances)
    return names, matrix
####does this rely on the names being in the same order for all frequency matrices to do that stat tests accurately?

#---------------------------------------------------------------------------
# 4. find, validate, and load distance matrices
#---------------------------------------------------------------------------
#locate the files
distance_matrix_files_list = list(distance_matrix_folder.glob("*distance_matrix"))


#verify the files
EXPECTED_FILE_COUNT = 8 #update this if you change how many phylogenetic trees you have

for file in distance_matrix_files_list:
    if not file.is_file():
        raise FileNotFoundError(f"file not found for:\n  {file}")
if len(distance_matrix_files_list) != EXPECTED_FILE_COUNT:
    raise RuntimeError(
        f"Expected {EXPECTED_FILE_COUNT} input files "
        f"but found {len(distance_matrix_files_list)}")


#load them into a new dictionary, where the relevant part of the file name and the matrix are stored together
distance_matrix_dictionary = {}
for file in distance_matrix_files_list:
    matrix_name = file.stem.replace("_distance_matrix","")
    names, triangle_distance_matrix = load_distance_matrix(file)
    number_of_sequences = len(names)
    if len(triangle_distance_matrix) !=number_of_sequences:
        raise RuntimeError(f"\nUnexpected matrix shapre for {matrix_name}\n")

    #turn the triangle matrix into a square matrix
    square_distance_matrix = np.zeros(
        (number_of_sequences,
        number_of_sequences))
    for row_index, row in enumerate(triangle_distance_matrix):
        for column_index, distance in enumerate(row):
            square_distance_matrix[row_index, column_index] = distance
            square_distance_matrix[column_index, row_index] = distance
    distance_matrix_dictionary[matrix_name] = {
        "names": names,
        "triangle matrix": triangle_distance_matrix,
        "square matrix": square_distance_matrix,
        }
    print(f"{matrix_name}: {len(names)} plasmids")
            

#list which ifles are in the dictionary

print("\n\nfiles in the defined dictionary:")
for matrix_name in distance_matrix_dictionary:
    print(f"\t{matrix_name}")


#verify the names correspond to each other
reference_set = None
for matrix_name in distance_matrix_dictionary:
    current_set = set(distance_matrix_dictionary[matrix_name]["names"])
    if reference_set is None:
        reference_set = current_set
    elif current_set != reference_set:
        print(f"\n\nERROR: {matrix_name} has different set of plasmids")
        raise RuntimeError
reference_names = None
for matrix_name in distance_matrix_dictionary:
    current_names = (distance_matrix_dictionary[matrix_name]["names"])
    if reference_names is None:
        reference_names = current_names
    elif current_names != reference_names:
        print(f"\nERROR:\n{matrix_name}\nhas plasmids in different order\n")

        raise RuntimeError

#---------------------------------------------------------------------------
# 5. perform pairwisw Mantel tests and save results
#---------------------------------------------------------------------------

print("\n\nrunning Mantel tests. might take a few minutes")
#make a list of all pairwise combinations in the dictionary
pairwise_comparisons = itertools.combinations(distance_matrix_dictionary.keys(), 2)


#make a results list

mantel_results = []
for matrix_A_name, matrix_B_name in pairwise_comparisons: # go through every pairwise combination
    
    #fetch the corresponding matrices
    square_distance_matrix_A = DistanceMatrix(
        distance_matrix_dictionary[matrix_A_name]["square matrix"],
        ids= distance_matrix_dictionary[matrix_A_name]["names"])
    square_distance_matrix_B = DistanceMatrix(
        distance_matrix_dictionary[matrix_B_name]["square matrix"],
        ids= distance_matrix_dictionary[matrix_B_name]["names"])

    #run mantel
    r_value, p_value, n = mantel(
        square_distance_matrix_A, 
        square_distance_matrix_B, 
        method="pearson", 
        permutations = 999)

    mantel_results.append({"matrix_A" : matrix_A_name, "matrix_B" : matrix_B_name, "mantel_r" : r_value, "p_value" : p_value})


#write a reuslts table and save an outfile
results_table = pd.DataFrame(mantel_results)

mantel_results_file = (phylogenetic_signal_folder / "pairwise_mantel_results.tsv")
results_table.to_csv(mantel_results_file,sep="\t",index=False)

#---------------------------------------------------------------------------
# 6. heatmap
#---------------------------------------------------------------------------

#make a table with all the r's
matrix_names = sorted(distance_matrix_dictionary.keys())
heatmap_table = pd.DataFrame(
    1.0,
    index = matrix_names,
    columns = matrix_names)
for _, row in results_table.iterrows():
    a = row["matrix_A"]
    b = row["matrix_B"]
    r = row["mantel_r"]
    heatmap_table.loc[a, b] = r
    heatmap_table.loc[b, a] = r
heatmap_table.to_csv(phylogenetic_signal_folder / "mantel_heatmap_values.tsv", sep="\t")

# make it into a heatmap
plt.figure(figsize = (10,8))
sns.heatmap(heatmap_table, annot = True, cmap= "viridis", vmin=-1, vmax=1)
plt.tight_layout()
plt.savefig(phylogenetic_signal_folder / "mantel_heatmap.png", dpi=300)




#---------------------------------------------------------------------------
# 7. finished!
#---------------------------------------------------------------------------

print("\n\n\t***inished running script 303!!***\n\n")
