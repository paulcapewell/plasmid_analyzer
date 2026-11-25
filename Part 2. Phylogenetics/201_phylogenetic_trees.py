print("\n\n\n**Running Script: 201_phylogenetic_trees.py**\n\n")



#---------------------------------------------------------------------------
# 1. set up
#---------------------------------------------------------------------------

#import modules
from pathlib import Path
import subprocess
import os
from Bio import SeqIO
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import (DistanceCalculator, DistanceTreeConstructor)
from Bio import Phylo

# define folders
p2_folder = Path(__file__).resolve().parent
PA_folder = p2_folder.parent
p1_folder = PA_folder / "Part 1. emboss Scan"
p1_outfolder = p1_folder / "output_folder"
comp_seq_folder = p1_outfolder / "compiled_sequences_near_xerCD"
input_sequences_folder = comp_seq_folder / "218+ bp"

p2_outfolder = p2_folder / "output_folder"
p2_outfolder.mkdir(parents=True, exist_ok=True)

aligned_sequences_folder = p2_outfolder / "MAFFT aligned sequences"
aligned_sequences_folder.mkdir(parents=True, exist_ok=True)

distance_matrix_folder = p2_outfolder / "distance matrices (used to calcualte trees)"
distance_matrix_folder.mkdir(parents=True, exist_ok=True)

tree_folder = p2_outfolder / "trees"
tree_folder.mkdir(parents=True, exist_ok=True)


#define input files
input_sequences_files = [
    input_sequences_folder / "central_218bp_sequences_collapsed.txt",
    input_sequences_folder / "XerC_and_up_1000bp_sequences.txt",
    input_sequences_folder / "XerD_and_down_1000bp_sequences.txt"]

EXPECTED_FILE_COUNT = 3 #i.e. the number of input_sequnece_Files


#validate folders and files exist
folders_to_check = [input_sequences_folder, p1_outfolder, p1_folder, comp_seq_folder]
for folder_path in folders_to_check:
    if not folder_path.is_dir():
        print(f"folder not found: {folder_path}")
        raise FileNotFoundError
for file_path in input_sequences_files:
    if not file_path.is_file():
        raise FileNotFoundError(f"file not fount for:\n  {file_path}")

print("\nInput files located successfully:\n")
for file_path in input_sequences_files:
    print(f"  ✓ {file_path.name}")
print(f"\nTotal files found: {len(input_sequences_files)}\n")
if len(input_sequences_files) != EXPECTED_FILE_COUNT:
    raise RuntimeError(
        f"Expected {EXPECTED_FILE_COUNT} input files "
        f"but found {len(input_sequences_files)}")

#---------------------------------------------------------------------------
# 2. loop through input sequence files
#---------------------------------------------------------------------------

for input_sequences_file in input_sequences_files:
    file_name_root = input_sequences_file.stem
    print(f"processing {file_name_root}. \n\n\tRunning MAFFT. probably less than 1 minute")

    #---------------------------------------------------------------------------
    # 3. run multiple sequence alignment
    #---------------------------------------------------------------------------
    aligned_sequence_file = (aligned_sequences_folder / f"{file_name_root}_aligned.fasta")

    command = ["mafft", "--auto", str(input_sequences_file)]

    with open(aligned_sequence_file, "w") as outfile:
        result = subprocess.run(command, stdout=outfile, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print(result.stderr)
        raise RuntimeError("MAFFT alignment failed")

    alignment = AlignIO.read(aligned_sequence_file, "fasta")
    print("\t\tfinished\n\n"
        f"\t\t{len(alignment)} sequences | {alignment.get_alignment_length()} columns\n\n")

    #---------------------------------------------------------------------------
    # 4. create pairwise distance matrix
    #---------------------------------------------------------------------------

    print("\tcalculating distance matrix. Might take 2-20 minutes")
    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)
    distance_matrix_file = (distance_matrix_folder / f"{file_name_root}_distance_matrix")

    with open(distance_matrix_file, "w") as outfile:
        outfile.write(str(distance_matrix))

    print("\t\tfinished\n\n")


    #---------------------------------------------------------------------------
    # 5. create phylogenetic trees
    #---------------------------------------------------------------------------
    constructor = DistanceTreeConstructor()

    print("\tCreating NJ tree. Might take 15-45 minutes")
    tree = constructor.nj(distance_matrix)

    tree_file = (tree_folder / f"{file_name_root}_tree.nwk")
    Phylo.write(tree, tree_file, "newick")
    print("\t\tfinished\n\n")

print("\n\n\t***Finished running script 201!!***\n\n")