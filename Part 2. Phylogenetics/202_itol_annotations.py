print("\n\n\n\t**Running Script: 202_annotations.py**\n\n")



#---------------------------------------------------------------------------
# 1. set up
#---------------------------------------------------------------------------
#import modules
from pathlib import Path
import subprocess
import os
import glob
from Bio import SeqIO
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import (DistanceCalculator, DistanceTreeConstructor)
from Bio import Phylo
import statistics


#define folder structure
p2_folder = Path(__file__).resolve().parent
PA_folder = p2_folder.parent
p1_folder = PA_folder / "Part 1. emboss Scan"
p1_outfolder = p1_folder / "output_folder"
p2_outfolder = p2_folder / "output_folder"
aligned_sequences_folder = p2_outfolder / "MAFFT aligned sequences"
distance_matrix_folder = p2_outfolder / "distance matrices (used to calcualte trees)"
tree_folder = p2_outfolder / "trees"


#verify input folder exist
folders_to_check = [p1_folder, p1_outfolder, p2_outfolder, aligned_sequences_folder, distance_matrix_folder, tree_folder]
for folder_path in folders_to_check:
    if not folder_path.is_dir():
        print(f"folder not found: {folder_path}")
        raise FileNotFoundError


#define input files
tree_files_list = list(tree_folder.glob("*_tree.nwk"))
distance_matrix_files_list = list(distance_matrix_folder.glob("*distance_matrix*"))


#verify input files
EXPECTED_FILE_COUNT = 3 #i.e. the number of input_sequnece_Files
if len(tree_files_list) != EXPECTED_FILE_COUNT:
    raise RuntimeError(
        f"Expected {EXPECTED_FILE_COUNT} input files, but found {len(tree_files_list)}")
if len(distance_matrix_files_list) != EXPECTED_FILE_COUNT:
    raise RuntimeError(
        f"Expected {EXPECTED_FILE_COUNT} input files, but found {len(tree_files_list)}")


#define new outfolders
annotation_folder = p2_outfolder / "annotations"
annotation_folder.mkdir(parents=True, exist_ok=True)
annotation_218bp_folder = (annotation_folder / "for_collapsed_218bp_tree")
annotation_218bp_folder.mkdir(parents=True, exist_ok=True)
annotation_1kb_folder = (annotation_folder /"for_1kb_trees")
annotation_1kb_folder.mkdir(parents=True, exist_ok=True)

#---------------------------------------------------------------------------
# 2. host lookup table
#---------------------------------------------------------------------------

print("finding host species of each plasmid... shouldn't take more than ~1 minutes\n")
host_lookup = {}
gbk_folder = (
    p1_outfolder /
    "individual_GBKs" /
    "full sequence")
for gbk_file in gbk_folder.glob("*.gbk"):
    record = SeqIO.read(gbk_file, "genbank")
    accession = record.id
    organism = record.annotations.get(
        "organism",
        "Unknown")
    host_lookup[accession] = organism
genus_lookup = {}
for accession, organism in host_lookup.items():
    genus = organism.split()[0]
    genus_lookup[accession] = genus



#---------------------------------------------------------------------------
# load host genus information
#---------------------------------------------------------------------------

accession_to_genus = {}
gbk_folder = (p1_outfolder / "individual_GBKs" / "full sequence")
for gbk_file in gbk_folder.glob("*.gbk"):
    try:
        record = SeqIO.read(gbk_file, "genbank")
        organism = record.annotations.get("organism", "Unknown")
        genus = organism.split()[0]
        accession_to_genus[record.id] = genus
    except Exception as e:
        print(f"Unable to read {gbk_file.name}: {e}")
print(f"\tLoaded host information for {len(accession_to_genus):,} plasmids\n")
all_genera = sorted(
    set(accession_to_genus.values()))
print(f"\tFound {len(all_genera)} host genera:\n")
for genus in all_genera:
    print("\t\t", genus)
genus_counts = {}
for genus in accession_to_genus.values():
    genus_counts[genus] = genus_counts.get(genus, 0) + 1
print("\nMost common host genera:\n")
for genus, count in sorted(
    genus_counts.items(),
    key=lambda x: x[1],
    reverse=True):
    print(f"{count:>5}  {genus}")


#---------------------------------------------------------------------------
# 3. helper functions
#---------------------------------------------------------------------------

def get_accession_from_1kb_node(node_name):
    #e.g., NZ_CP011315.1_3570_REV_30 -> NZ_CP011315.1
    return "_".join(node_name.split("_")[:-3])

def get_accession_from_collapsed_member(member):
    #e.g., NZ_CP011315.1_3570 -> NZ_CP011315.1
    return "_".join(member.split("_")[:-1])

def get_genus_from_accession(accession):
    return accession_to_genus.get(accession, "Unknown")

def get_species_from_accession(accession):
    return host_lookup.get(accession, "Unknown")

def get_collapsed_node_genus(node_name):
    genera = set()
    for member in node_name.split("|"):
        accession = get_accession_from_collapsed_member( member)
        genus = get_genus_from_accession(accession)
        genera.add(genus)
    return sorted(genera)


#---------------------------------------------------------------------------
# 4. define variables for annotations
#---------------------------------------------------------------------------
#host genus symbols
host_symbol_map = {
    "Klebsiella": 2,    # circle
    "Escherichia": 2,   # circle
    "Salmonella": 2,    # circle
    "Enterobacter": 2,  # circle
    "Shigella": 2,
    "Citrobacter": 2,
    "Edwardsiella": 2,
    "Yersinia": 2,
    "Serratia": 2,}
default_host_symbol = 1 # rectange
host_color_map = {
    "Klebsiella": "#ff7f0e",
    "Escherichia": "#1f77b4",
    "Salmonella": "#2ca02c",
    "Enterobacter": "#d62728",
    "Shigella": "#bae91e",
    "Citrobacter": "#9467bd",
    "Edwardsiella": "#8c564b",
    "Yersinia": "#e377c2",
    "Serratia": "#17becf",
    "Unknown": "#7f7f7f"}


#known functional sites
the_cer_site_accession = "NC_001371.1"
other_known_sites_accession_list = [
    "NC_001371.1",
    "NC_006881.1",
    "NC_003486.1",
    "NC_002119.1",
    "NC_002119.1",
    "NC_008053.1",
    "NC_008053.1",
    "NC_001373.1",
    "NC_002090.1",
    "NC_017330.1",
    "X53663.1"]
known_sites = set(other_known_sites_accession_list)


#clusters
MATRIX_CLUSTER_MAX_DISTANCE = 0.01
PHYLO_CLUSTER_MAX_DISTANCE = 0.02
MIN_CLUSTER_SIZE = 5
cluster_colors = [
    "#FF0000",
    "#00FF00",
    "#0000FF",
    "#FFFF00",
    "#FF00FF",
    "#00FFFF",
    "#FF8000",
    "#8000FF",
    "#0080FF",
    "#80FF00",
    "#FF0080",
    "#00FF80",
    "#800000",
    "#008000",
    "#000080",
    "#808000",
    "#800080",
    "#008080",
    "#FF4040",
    "#40FF40",
    "#4040FF",
    "#FFAA00",
    "#AA00FF",
    "#00AAFF",
    "#FF0066",
    "#66FF00",
    "#0066FF",
    "#FF6600",
    "#6600FF",
    "#00FFAA",
    "#AAFF00",
    "#FF00AA",
    "#00AA80",
    "#AA8000",
    "#8000AA",
    "#0080AA",
    "#C00000",
    "#00C000",
    "#0000C0",
    "#C0C000",
    "#C000C0",
    "#00C0C0",
    "#FF8080",
    "#80FF80",
    "#8080FF",
    "#FFC080",
    "#C080FF",
    "#80C0FF",
    "#A00040",
    "#40A000",
    "#0040A0",
    "#A06000",
    "#6000A0",
    "#00A060",
    "#202020",
    "#A28888"
]


#---------------------------------------------------------------------------
# 5. collapsed 218bp annotations - general
#---------------------------------------------------------------------------
print("\n---------------------------------------------------------------------------\nCreating collapsed 218bp annotation\n---------------------------------------------------------------------------\n")


collapsed_tree_file = (tree_folder /"central_218bp_sequences_collapsed_tree.nwk")
collapsed_tree = Phylo.read(collapsed_tree_file,"newick")


#confirmed functional sites annotation
itol_confirmed_file = (annotation_218bp_folder / "functional_sites.txt")
with open(itol_confirmed_file, "w") as file:
    file.write("TREE_COLORS\n")
    file.write("SEPARATOR TAB\n")
    file.write("LEGEND_TITLE\tConfirmed functional Xer sites\n")
    file.write("LEGEND_SHAPES\t1\t1\n")
    file.write("LEGEND_COLORS\t#0000ff\t#ff0000\n")
    file.write("LEGEND_LABELS\tColE1 cer\tOther known site\n")
    file.write("DATA\n")
    for terminal in collapsed_tree.get_terminals():
        node_name = terminal.name
        if node_name is None:
            continue
        if the_cer_site_accession in node_name:
            file.write(f"{node_name}\trange\t#0000ff\tColE1 cer\n")
        elif any(accession in node_name for accession in known_sites):
            file.write(f"{node_name}\trange\t#ff0000\tKnown Xer site\n")


#host genus symbol annotation
itol_host_file = (annotation_218bp_folder / "host_genus_symbols.txt")
with open(itol_host_file, "w") as file:
    file.write("DATASET_SYMBOL\n")
    file.write("SEPARATOR TAB\n")
    file.write("DATASET_LABEL\tHost genus\n")
    file.write("COLOR\t#000000\n")
    file.write("MAXIMUM_SIZE\t6\n")
    file.write("DATA\n")
    for terminal in collapsed_tree.get_terminals():
        node_name = terminal.name
        if node_name is None:
            continue
        genera = get_collapsed_node_genus(node_name)
        positions = [-1, -2, -3, -4, -5, -6, -7, -8, -9, -10]
        for position, genus in zip(positions, genera):
            symbol = host_symbol_map.get(genus, default_host_symbol)
            color = host_color_map.get(genus,"#7f7f7f")
            file.write(
                f"{node_name}\t"
                f"{symbol}\t"
                f"6\t"
                f"{color}\t"
                f"1\t"
                f"{position}\t"
                f"{genus}\n")


#host species popup annotation
itol_popup_file = (annotation_218bp_folder / "host_species_popup.txt")
with open(itol_popup_file, "w") as file:
    file.write("POPUP_INFO\n")
    file.write("SEPARATOR TAB\n")
    file.write("DATA\n")
    for terminal in collapsed_tree.get_terminals():
        node_name = terminal.name
        species_set = set()
        for member in node_name.split("|"):
            accession = (get_accession_from_collapsed_member(member))
            species_set.add(get_species_from_accession(accession))
        popup_text = "<br>".join(sorted(species_set))
        file.write(
            f"{node_name}\t"
            f"Host species represented\t"
            f"{popup_text}\n")
        
#---------------------------------------------------------------------------
# 6. collapsed 218bp annotations - MEME
#---------------------------------------------------------------------------
print("\tGenerating MEME\n")

#leave this section blank for now
print("\tFinished collapsed-tree annotations\n\n\n")

#---------------------------------------------------------------------------
# 7. 1kbp tree annotations - general
#---------------------------------------------------------------------------
#define both 1kb tree files
one_kb_tree_files = [tree_file
    for tree_file in tree_files_list
    if "1000bp" in tree_file.name]


#define one of them to use
xerC_1kb_tree_file = tree_folder / "XerC_and_up_1000bp_sequences_tree.nwk"
print(f"\n---------------------------------------------------------------------------\nCreating general 1kb tree annotations\n\tusing {xerC_1kb_tree_file.name}\n---------------------------------------------------------------------------")
xerC_1kb_tree = Phylo.read(xerC_1kb_tree_file, "newick")

#make an itol annotaiton for known functional sites
itol_confirmed_file = (annotation_1kb_folder / "known_functional_sites.txt")
with open(itol_confirmed_file, "w") as file:
    file.write("TREE_COLORS\n")
    file.write("SEPARATOR TAB\n")
    file.write("LEGEND_TITLE\tConfirmed functional Xer sites\n")
    file.write("LEGEND_SHAPES\t1\t1\n")
    file.write("LEGEND_COLORS\t#0000ff\t#ff0000\n")
    file.write("LEGEND_LABELS\tColE1 cer\tOther known site\n")
    file.write("DATA\n")
    for terminal in xerC_1kb_tree.get_terminals():
        node_name = terminal.name
        if node_name is None:
            continue
        if the_cer_site_accession in node_name:
            file.write(f"{node_name}\trange\t#0000ff\tColE1 cer\n")
        elif any(
            accession in node_name
            for accession in known_sites):
            file.write(f"{node_name}\trange\t#ff0000\tKnown Xer site\n")


#make an itol symbol annotaiton for host genus
itol_host_file = (annotation_1kb_folder / "host_genus_symbols.txt")
with open(itol_host_file, "w") as file:
    file.write("DATASET_SYMBOL\n")
    file.write("SEPARATOR TAB\n")
    file.write("DATASET_LABEL\tHost genus\n")
    file.write("COLOR\t#000000\n")
    file.write("MAXIMUM_SIZE\t2\n")
    file.write("DATA\n")
    for terminal in xerC_1kb_tree.get_terminals():
        node_name = terminal.name
        accession = get_accession_from_1kb_node(node_name)
        genus = get_genus_from_accession(accession)
        symbol = host_symbol_map.get(genus, default_host_symbol)
        color = host_color_map.get(genus, "#7f7f7f")
        file.write(
            f"{node_name}\t"
            f"{symbol}\t"
            f"2\t"
            f"{color}\t"
            f"1\t"
            f"-1\t"
            f"{genus}\n")

#itol popup annotation for host species
itol_popup_file = (
    annotation_1kb_folder /
    "host_species_popup.txt")
with open(itol_popup_file, "w") as file:
    file.write("POPUP_INFO\n")
    file.write("SEPARATOR TAB\n")
    file.write("DATA\n")
    for terminal in xerC_1kb_tree.get_terminals():
        node_name = terminal.name
        accession = get_accession_from_1kb_node(node_name)
        species = get_species_from_accession(accession)
        file.write(
            f"{node_name}\t"
            f"Host species\t"
            f"{species}\n")

#---------------------------------------------------------------------------
# 8. 1kbp tree cluster annotations - clusters
#---------------------------------------------------------------------------


#1 loop through both 1kb tree files
print("\n---------------------------------------------------------------------------\nMaking annotations of phylogenetic clusters\n---------------------------------------------------------------------------\n")
for tree_file in one_kb_tree_files:
    #how do I make this stop after the first itteration. it's a waste of computing power to run this more than once at this point.
    print(f"\nprocessing {tree_file.name}")
    tree = Phylo.read(tree_file, "newick")
    prefix = tree_file.stem


    #2 Afind matching distance matrix
    matrix_stem = tree_file.stem.replace("_tree", "") #define part of the name form the tree_file name
    distance_matrix_file = (distance_matrix_folder / f"{matrix_stem}_distance_matrix")
    if not distance_matrix_file.is_file():
        raise FileNotFoundError(
            f"No matching matrix file for:\n"
            f"{tree_file.name}")
    print(f"\n\tusing matrix {distance_matrix_file.name}")

 
    #3 load matrix into memory
    names = [] #names starts empty
    graph = {} #graph starts empty
    visited = set()
    accepted_clusters = []
    cluster_number = 1
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
            node_name = parts[0]
            names.append(node_name)
            graph.setdefault(node_name,set())
            distances = parts[1:]
            for previous_index, distance_text in enumerate(distances):
                distance = float(distance_text)
                previous_node = names[previous_index]
                if distance <= MATRIX_CLUSTER_MAX_DISTANCE:
                    graph[node_name].add(previous_node)
                    graph[previous_node].add(node_name)


    #4 identify connected components
    for node in graph:
        if node in visited:
            continue
        stack = [node]
        component = []
        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)
            component.append(current)
            stack.extend(
                graph[current])
        accepted_clusters.append((f"Cluster_{cluster_number}", component))
        cluster_number += 1
        accepted_clusters = [(cluster_id, nodes)
                             for cluster_id, nodes in accepted_clusters
                             if len(nodes) >= MIN_CLUSTER_SIZE]
    final_clusters = [(f"Cluster_{i+1}", members) for i, (_, members) in enumerate(accepted_clusters)]
    print(f"\nidentified {len(final_clusters)} clusters where some members are distance < {MATRIX_CLUSTER_MAX_DISTANCE}, and 5+ plasmdis in cluster")
    cluster_sizes = sorted([len(nodes) for _, nodes in final_clusters], reverse = True)
    print("\nlargest clusters:", cluster_sizes[:20])
    print(f"\nplasmids assigned to cluster: {sum(cluster_sizes)}")

#5. write a tsv file for each
    clusters_tsv_file = (annotation_1kb_folder / f"{tree_file.stem}_clusters_summary.tsv")

    with open(clusters_tsv_file, "w") as outfile:
        outfile.write("cluster_id\tnode_name\n")
        for cluster_id, node_names in final_clusters:
            for node_name in node_names:
                outfile.write(f"{cluster_id}\t" f"{node_name}\n")


#6. write an itol .txt file for each
    itol_cluster_annotation_file = (annotation_1kb_folder / f"{tree_file.stem}_clusters_annotation.txt")
    with open(itol_cluster_annotation_file, "w") as file:
        file.write("DATASET_COLORSTRIP\n")
        file.write("SEPARATOR TAB\n")
        file.write(f"DATASET_LABEL\t{tree_file.stem} Phylogenetic Clusters\n")
        file.write("COLOR\t#000000\n")
        file.write("STRIP_WIDTH\t15\n")
        file.write("DATA\n")
        for i, (cluster_id, node_names) in enumerate(final_clusters):
            cluster_color = cluster_colors[i % len(cluster_colors)]
            for node_name in node_names:
                file.write(f"{node_name}\t{cluster_color}\t{cluster_id}\n")


#---------------------------------------------------------------------------
# 9. Finish
#---------------------------------------------------------------------------
print("\n\n     ~~ Finished running script 202!! ~~\n\n")
