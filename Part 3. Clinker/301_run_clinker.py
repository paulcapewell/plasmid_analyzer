print("\n\n\n\t**Running Script: 301_run_clinker.py**\n\n")



#---------------------------------------------------------------------------
# 1. set up 
#---------------------------------------------------------------------------


#import modules
import subprocess
from pathlib import Path
import pandas as pd




#define pre-existing folders
p3_folder = Path(__file__).resolve().parent
PA_folder = p3_folder.parent
p1_folder = PA_folder / "Part 1. emboss Scan"
p1_outfolder = p1_folder / "output_folder"
indiv_gbk_folder = p1_outfolder / "individual_GBKs"
gbk_3kb_long_folder = indiv_gbk_folder / "3kb"
p2_folder = PA_folder / "Part 2. Phylogenetics"
p2_outfolder = p2_folder / "output_folder"
p2_annotations_folder = p2_outfolder / "annotations"
p2_1kb_annot_folder = p2_annotations_folder / "for_1kb_trees"

#confirm these folders exist
folders_to_check = [p1_folder, p1_outfolder, gbk_3kb_long_folder, p2_folder, p2_outfolder, p2_annotations_folder, p2_1kb_annot_folder, indiv_gbk_folder]
for folder_path in folders_to_check:
    if not folder_path.is_dir():
        print(f"\n\tfolder not found: {folder_path}\n")
        raise FileNotFoundError



#create new output folders
p3_outfolder = p3_folder / "output_folder"
p3_outfolder.mkdir(parents=True, exist_ok=True)
XerC_clinker_folder = p3_outfolder / "HTMLs from XerC 1kb clusters"
XerC_clinker_folder.mkdir(parents=True, exist_ok=True)
XerD_clinker_folder = p3_outfolder / "HTMLs from XerD 1kb clusters"
XerD_clinker_folder.mkdir(parents=True, exist_ok=True)



#define input files
XerC_1kb_clusters_file = p2_1kb_annot_folder / "XerC_and_up_1000bp_sequences_tree_clusters_summary.tsv"
XerD_1kb_clusters_file = p2_1kb_annot_folder / "XerD_and_down_1000bp_sequences_tree_clusters_summary.tsv"

cluster_files = [XerC_1kb_clusters_file, XerD_1kb_clusters_file]
for file_path in cluster_files:
    if not file_path.is_file():
        raise FileNotFoundError(f"\n\tfile not fount for:\n  {file_path}")
print("\nBoth cluster files located successfully\n")

gbk_3kb_long_files = sorted(gbk_3kb_long_folder.glob("*.gbk"))
print(f"\nnumber of 3kb GBK files found: {len(gbk_3kb_long_files)}\n")


#---------------------------------------------------------------------------
# 2. run clinker
#---------------------------------------------------------------------------

#define a subprocess to run clinker
def run_clinker_on_cluster_file(cluster_id, gbk_files_list, output_path):
    html_file = output_path / f"{cluster_id}.html"
    subprocess.run(["clinker", *[str(file) for file in gbk_files_list], "--force", "-p", str(html_file)],check=True)


#loop through XerC and XerD clusters seperately, but note which one you are working on
for cluster_file in cluster_files:
    if "XerC" in cluster_file.name:
        XerC_or_XerD = "XerC"
    elif "XerD" in cluster_file.name:
        XerC_or_XerD = "XerD" 
    else:
        raise ValueError(f"Could not determine whether file is XerC or XerD for: {cluster_file.name}")


    #read the cluster TSV file and convert it into a table and dictionary. we will use the dictionary for clinker
    cluster_table = pd.read_csv(cluster_file, sep="\t")
    cluster_dictionary = {}
    for _, row in cluster_table.iterrows():
        cluster_id = row["cluster_id"]
        plasmid_name = row["node_name"]
        if cluster_id not in cluster_dictionary:
            cluster_dictionary[cluster_id] = []
        cluster_dictionary[cluster_id].append(plasmid_name)

    #for each cluster id in the dictionary, run clinker (you need to define a few variables for each cluster to do that)
    for cluster_id in cluster_dictionary: 
        #you need to make a list of the 3kb GBK files that correspond to the plasmid names in the dictionary
        gbk_files_list = [] 
        for plasmid_name in cluster_dictionary[cluster_id]:
            gbk_path = (gbk_3kb_long_folder / f"{plasmid_name}.gbk")
            if not gbk_path.is_file():
                raise FileNotFoundError(f"\n3kb GBK file not found for {gbk_path}")
            gbk_files_list.append(gbk_path)
        #you need to define a file path
        if XerC_or_XerD == "XerC":
            output_path = XerC_clinker_folder
        elif XerC_or_XerD == "XerD":
            output_path = XerD_clinker_folder
        #run clinker
        print(f"\nRunning clinker on 1kb {XerC_or_XerD}'s {cluster_id}, which contains {len(gbk_files_list)} plasmids")
        run_clinker_on_cluster_file(cluster_id, gbk_files_list, output_path)




#---------------------------------------------------------------------------
# 3. Finished
#---------------------------------------------------------------------------
    
print("\n\n\t***Finished running script 301!!***\n\n")
