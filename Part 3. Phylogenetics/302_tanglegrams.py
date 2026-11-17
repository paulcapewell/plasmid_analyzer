print("\n\n\n**Running Script: 302_tanglegrams.py**\n\n")





from pathlib import Path
import subprocess
import os
from Bio import Phylo
import matplotlib.pyplot as plt

# define folders
p3_folder = Path(__file__).resolve().parent
PA_folder = p3_folder.parent
p2_folder = PA_folder / "Part 2. XerCD Sequences Metadata"
p1_folder = PA_folder / "Part 1. emboss Scan"
p2_outfolder = p2_folder / "output_folder"
p3_outfolder = p3_folder / "output_folder"
phylogenetic_folder = p3_outfolder / "phylogenetic tree files"

#validate folders exist
folders_to_check = {
    "Part 2 folder": p2_folder,                
    "Part 1 folder": p1_folder,
    "Part 2 output folder": p2_outfolder,
    "Part 3 output fiolder": p3_outfolder,
    "phylogenetic tree files": phylogenetic_folder}
for name, folder_path in folders_to_check.items():
    if not folder_path.is_dir():
        print(f"{name} not found for {folder_path}")
        raise FileNotFoundError

#define files
xerC_tree_file = phylogenetic_folder / "xerC and upstream sequences tree"
xerD_tree_file = phylogenetic_folder / "xerD and downstream sequences tree"

#validate files exist
files_to_check = {
    "xerC_tree_file": xerC_tree_file,
    "xerD_tree_file": xerD_tree_file}
for name, file_path in files_to_check.items():
    if not file_path.is_file():
        print(f"{name} not found for {file_path}")
        raise(FileNotFoundError)
    if not file_path.exists():
        print(f"{name}, i.e. {file_path}, does not exist")
        raise(FileNotFoundError)


#read the two Newick trees
tree_xerC = Phylo.read(xerC_tree_file, "newick")
tree_xerD = Phylo.read(xerD_tree_file, "newick")

#get coordinates for each tree
def get_tree_coordinates(tree):
    terminals = tree.get_terminals() 
    y_positions = {tip: i for i, tip in enumerate(terminals)}
    coordinates = {}
    def calculate(clade, parent_x=0):
        branch_length = clade.branch_length or 0
        x = parent_x + branch_length
        if clade.is_terminal():
            y = y_positions[clade]
        else:  
            child_ys = [
                calculate(child, x) [1]
                for child in clade.clades]
            y = sum(child_ys) / len(child_ys)
        coordinates[clade] = (x, y)
        return x, y
    calculate(tree.root) 
    max_x = max (x for x, y in coordinates.values())

    if max_x > 0:
        coordinates = {
            clade: (x / max_x, y)
            for clade, (x, y) in coordinates.items()}
    return coordinates
coords_xerC_tree = get_tree_coordinates(tree_xerC)
coords_xerD_tree = get_tree_coordinates(tree_xerD)

#get plasmidnames and confirm they match
names_xerC_plasmids = [tip.name for tip in tree_xerC.get_terminals()]
set_names_xerC_plasmids = set(names_xerC_plasmids)

names_xerD_plasmids = [tip.name for tip in tree_xerD.get_terminals()]
set_names_xerD_plasmids = set(names_xerD_plasmids)

if set_names_xerC_plasmids != set_names_xerD_plasmids:
    print(f"WARNING: Trees contain different plasmid names!\n  Only in xerC tree: {set_names_xerC_plasmids - set_names_xerD_plasmids} \n  Only in xerD tree: {set_names_xerD_plasmids - set_names_xerC_plasmids}")
    raise ValueError("Plasmid names need to be matched first.")
if len(names_xerC_plasmids) != len(set_names_xerC_plasmids) or len(names_xerD_plasmids) != len(set_names_xerD_plasmids):
    raise ValueError("Duplicate Plasmid Names Detected")


#set up figure
def draw_tree(ax, tree, coordinates, offset = 0, mirror = False):
    for clade in tree.find_clades(order = "level"):
        x, y = coordinates[clade]
        if mirror:
            x = offset - x
        else:
            x = offset + x
        child_coords = [] #what does clade and child mean in these commands and in phylogenetic trees?
        for child in clade.clades:
            cx, cy = coordinates[child]
            if mirror:
                cx = offset - cx

            else:
                cx = offset + cx
            child_coords.append((cx,cy))

            ax.plot([x,cx], [cy, cy], color = "black", lw=1) #wat's the purpose of this line? why is it xcx cycy?
        if child_coords:
            ys = [cy for cx, cy in child_coords]
            ax.plot(
                [x,x],
                [min(ys), max(ys)],
                color = "black",
                lw=1)



#draw the figure
fig, ax = plt.subplots(figsize=(14, 9))
draw_tree(ax, tree_xerC, coords_xerC_tree, offset = 0)
draw_tree(ax, tree_xerD, coords_xerD_tree, offset=3, mirror = True)

for tipC in tree_xerC.get_terminals(): #explain this code block to me
    name = tipC.name
    tipD=next(
        candidate for candidate in tree_xerD.get_terminals()
        if candidate.name == name)
    xC, yC = coords_xerC_tree[tipC]
    xD, yD = coords_xerD_tree[tipD]
    #explain the part above, I understand the part below

    ax.plot( #this block adds the line connecting the two trees
        [1, 2],
        [yC, yD],
        color = "red",
        alpha=0.5,
        lw=0.8)


#formatting
ax.set_title("Phylogenetic comparison of XerD and upstream vs XerD and downstream")
ax.set_xlim(-0.1, 3.1)
ax.axis("off")
plt.tight_layout()
plt.savefig(p3_outfolder / "plasmid_tangelgram.png", dpi=300, bbox_inches="tight")
plt.show()
    