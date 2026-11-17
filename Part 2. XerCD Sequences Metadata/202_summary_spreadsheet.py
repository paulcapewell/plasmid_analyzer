#1 import modules
from pathlib import Path
print("\n\n  Running Script create_plasmid_spreadsheet.py ... \n\n")
import csv
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

#2 define file paths
p2_folder = Path(__file__).resolve().parent
plasmid_analyzer_folder = p2_folder.parent 
p2_output_folder = p2_folder / "output_folder"
p1_folder = plasmid_analyzer_folder / "Part 1. emboss Scan"
p1_input_folder = p1_folder / "input_folder"
p1_output_folder = p1_folder / "output_folder"                   

    #2.1 define unsliced gbk file paths

gbk_folder = p2_output_folder / "individual GBKs (full sequence)"
if not gbk_folder.is_dir():
    raise FileNotFoundError(f"folder not found: full sequences")
gbk_files = sorted(gbk_folder.glob("*.gbk"))
if not gbk_files:
    raise FileNotFoundError("files not found: full sequence gbk files")
print(f"\n  ~~Found {len(gbk_files)} full length GBK files~~")
for file in gbk_files[:3]:
    print("   ", file.name)


    #2.2 define sliced gbk file paths
sliced_gbk_folder = p2_output_folder / "individual GBKs (3kb)"
if not sliced_gbk_folder.is_dir():
    raise FileNotFoundError(f"folder not found: central_3kb_sequences")
sliced_gbk_files = sorted(sliced_gbk_folder.glob("*.gbk"))
if not sliced_gbk_files:
    raise FileNotFoundError("files not found: 3kb sliced gbk files")
print("    ...", f"\n\n  ~~Found {len(sliced_gbk_files)} original sliced GBK files~~")
for file in sliced_gbk_files[:3]:
    print("   ",file.name)

    #2.3 define emboss results file paths
emboss_results_folder = p1_output_folder/ "xerCD hits after ArgR"
if not emboss_results_folder.is_dir():
    raise FileNotFoundError("folder not found: xerCD hits after ArgR")
emboss_results_files = sorted(emboss_results_folder.glob("*.txt"))
if not emboss_results_files:
    raise FileNotFoundError("files not found: emboss results files, i.e. 'NEW XerCD 28 FOR.'")
print("    ...", f"\n\n  ~~Found {len(emboss_results_files)} emboss results files~~")
for file in emboss_results_files:
    print("   ",file.name)

    #2.4 define pyxer script results file paths
pyxerscript_results_file = p2_output_folder / "compiled sequences (218+ bp)" / "central_218bp_sequences.txt"
if not pyxerscript_results_file:
    raise FileNotFoundError("files not found: 218bp sequence files, i.e. 'output_NEW XerCD 28 FOR.'")


    #2.5 define outputs
output_csv = p2_output_folder / "summary_spreadsheet.csv"


#4 write the csv
with open(output_csv, "w", newline="") as csvfile:
    write = csv.writer(csvfile)

    write.writerow([
            "Plasmid Name",
            "Plasmid Accession",
            "Plasmid Length (bp)",
            "Motif Identification Percentage",
            "Coordinates of the Xer Motif",
            "Xer Motif length (bp)",
            "Strand Orientation (+/-)",
            "218bp Region Sequence",
            "XerC Binding Sequence",
            "Crossover Region Sequence",
            "XerD Binding Sequence",
            "Number of CDS Within 1.5kb (number includes +1 for the xer-site itself)",
            "Nearest Upstream Gene (within ~1.5kb)",
            "Nearest Downstream Gene (within ~1.5kb)",
            "Host Organism",
            "Host Taxonomy",
            "GC Content (%)",
        ])
    print(f"    ...\n\n  ~~Created CSV:~~\n     {output_csv} \n")
    print("  ~~headers added~~ \n")




    #4.0 loop through sliced gbk files
    for sl_gbk_file in sliced_gbk_files:
        print(f"\nAdding row for {sl_gbk_file.name}\n     ...")


            #4.1 identify a sliced GBK file and its name
        filename = sl_gbk_file.stem
            #4.2 get data from sliced GBK file name
        filename = sl_gbk_file.stem
        parts = filename.split("_")
        filename_accession = "_".join(parts[:-3])
        motif_coordinates = int(parts[-3])
        strand_orientation = parts[-2]
        motif_length = int(parts[-1])


            #4.3 identify corresponding unsliced GBK
        gbk_file = gbk_folder / f"{filename_accession}.gbk"


            #4.4 get data from unsliced files
        gbk_record = SeqIO.read(gbk_file, "genbank")
        plasmid_accession = gbk_record.id
        if filename_accession != plasmid_accession:
            raise ValueError(f"Accession mismatch... Original = {plasmid_accession}... Sliced = {filename_accession}")
        plasmid_description = gbk_record.description
        plasmid_length = len(gbk_record.seq)
        host_organism = gbk_record.annotations["organism"]
        host_taxonomy = "; ".join(gbk_record.annotations["taxonomy"])                             
        for feature in gbk_record.features:
            if feature.type == "source":
                plasmid_name = feature.qualifiers.get(
                    "plasmid",
                    ["Unknown"]
                )[0]

        gc_percent = gc_fraction(gbk_record.seq) * 100


             #4.5 get data from the sliced GBK file
        sl_gbk_record = SeqIO.read(sl_gbk_file, "genbank")
        sliced_length = len(sl_gbk_record.seq)
        if sliced_length != 3000:
            raise ValueError(f"expected 3000bp, found {sliced_length}bp")
        sliced_id = sl_gbk_record.id
        sliced_accession = sliced_id.replace("_3000bp", "")
        if sliced_accession != plasmid_accession:
            raise ValueError(f"Accession mismatch... Original = {plasmid_accession}... Sliced = {sliced_accession}")
        number_nearby_cds = 0
        for feature in sl_gbk_record.features:
            if feature.type == "CDS":
                number_nearby_cds +=1
      


            #4.6 get data from GBK files's gene features
        gene_features = []
        for feature in sl_gbk_record.features:
            if feature.type != "gene":
                continue
            gene_features.append(feature)
        xer_index = None
        for i, feature in enumerate(gene_features):
            locus_tag = feature.qualifiers.get("locus_tag", ["No_TAG"])[0]
            if locus_tag == "XM_001":
                xer_index = i
                break
        if xer_index is None:
            raise ValueError("Could not find XM_001 feature")
        xer_feature = gene_features[xer_index]
        if xer_index > 0:
            upstream_feature = gene_features[xer_index - 1]
            upstream_label = upstream_feature.qualifiers.get("gene", upstream_feature.qualifiers.get("locus_tag", ["Unknown"]))[0]
            upstream_gene = f"gene: {upstream_label}"
        else:
            upstream_gene = "None"
        if xer_index < len(gene_features) - 1:
            downstream_feature = gene_features[xer_index + 1]
            downstream_label = downstream_feature.qualifiers.get("gene", downstream_feature.qualifiers.get("locus_tag", ["Unknown"]))[0]
            downstream_gene = f"gene: {downstream_label}"
        else:
            downstream_feature = "None"
            downstream_gene = "None found"



        #4.7 identify the 218bp sequence from pyxerscript results files
        search_pyxerscript_header = (f">{plasmid_accession}_{motif_coordinates}_{strand_orientation}_{motif_length}")
        pyxer_file_header_found = False
        with open(pyxerscript_results_file) as file:
            for line in file:
                if line.strip() == search_pyxerscript_header:
                    pyxer_file_header_found = True
                    sequence_218bp = next(file).strip()
                    break
        if not pyxer_file_header_found: print ("Header not found")
        if len(sequence_218bp) != 218:
            raise ValueError(f"Expected 218bp, found {len(sequence_218bp) }bp")
        
        #4.8 identify the percentage used to identify the motif_id



        matching_emboss_file = None
        for file in emboss_results_files:
            if file.name.endswith(
                f"{motif_length} {strand_orientation}.txt"
            ):
                matching_emboss_file = file
                break
        if matching_emboss_file is None:
            raise ValueError (f"No matching emboss file found for {motif_length} {strand_orientation}")
        search_emboss_header = (f"{plasmid_accession} {motif_coordinates}")
        emboss_file_header_found = False
        with open(matching_emboss_file) as file:
            for line in file:
                if line.startswith(f"{plasmid_accession} {motif_coordinates}"):
                    emboss_file_header_found = True
                    fields = line.split()
                    motif_id_percentage = int(fields[3])
                    break
        if not emboss_file_header_found:
            raise ValueError(f"Hit not found for {plasmid_accession} {motif_coordinates}")




        #4.9  identify the xerc, crossover region, and xerd sequence
        isolate_xerc_xerd_upstream_downstream = p2_output_folder / "compiled sequences (xerC, crossover region, xerD)"

        xerc_sequence_list_file = isolate_xerc_xerd_upstream_downstream / "all_xerC_sites.txt" #xer c block
        if not xerc_sequence_list_file:
            raise ValueError("file not found: all_xerC_sites.txt")
        xerc_sequence_line_found = False
        with open(xerc_sequence_list_file) as file:
            for line in file:
                if line.startswith(f">{plasmid_accession}_{motif_coordinates}"):
                    xerc_sequence_line_found = True
                    xerc_sequence = next(file, None)
        if not xerc_sequence_line_found:
            raise ValueError(f"xerc sequence not found for {plasmid_accession}_{motif_coordinates}")

        xerd_sequence_list_file = isolate_xerc_xerd_upstream_downstream / "all_xerD_sites.txt" #xer c block
        if not xerd_sequence_list_file:
            raise ValueError("file not found: all_xerD_sites.txt")
        xerd_sequence_line_found = False
        with open(xerd_sequence_list_file) as file:
            for line in file:
                if line.startswith(f">{plasmid_accession}_{motif_coordinates}"):
                    xerd_sequence_line_found = True
                    xerd_sequence = next(file, None)
        if not xerd_sequence_line_found:
            raise ValueError(f"xerd sequence not found for {plasmid_accession}_{motif_coordinates}")




        crossover_sequence_list_file = isolate_xerc_xerd_upstream_downstream / "all_crossover_regions.txt" #crossover region
        if not crossover_sequence_list_file:
            raise ValueError("file not found: all_crossover_regions.txt")
        crossover_sequence_line_found = False
        with open(crossover_sequence_list_file) as file:
            for line in file:
                if line.startswith(f">{plasmid_accession}_{motif_coordinates}"):
                    crossover_sequence_line_found = True
                    crossover_sequence = next(file, None)
                    crossover_sequence = crossover_sequence.strip()
        if not crossover_sequence_line_found:
            raise ValueError(f"crossover sequence not found for {plasmid_accession}_{motif_coordinates}")

        
        #5 write one row, while still in the loop
        write.writerow([
            plasmid_name,
            plasmid_accession,
            plasmid_length,
            motif_id_percentage,
            motif_coordinates,
            motif_length,
            strand_orientation,
            sequence_218bp,
            xerc_sequence,
            crossover_sequence,
            xerd_sequence,
            number_nearby_cds,
            upstream_gene,
            downstream_gene,
            host_organism,
            host_taxonomy,
            gc_percent])
        print("  ~~corresponding row added~~ \n\n")

print("\n\n          Finished!!!!!\n\n")





       

