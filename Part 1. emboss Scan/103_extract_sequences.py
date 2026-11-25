#!/usr/bin/env python3

print("\nRunning Script 103_extract_sequences.py ... \n\n")
#----------------------------------------
# 1. import modules and set up email
#----------------------------------------
from pathlib import Path
import os
import gzip
import pickle
from Bio import Entrez, SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature, CompoundLocation
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

ENTREZ_EMAIL = "2738493c@student.gla.ac.uk"
Entrez.email = ENTREZ_EMAIL
print(f"\n ***Email selected: {ENTREZ_EMAIL}***\n")

#----------------------------------------
# 2. define file paths and infiles
#----------------------------------------

P1_folder = Path(__file__).resolve().parent                                 #P2 folder
plasmid_analyzer_folder = P1_folder.parent                                  #plasmid analyzer folder
P1_output_folder = P1_folder / "output_folder"                              #P1 output folder

comp_seq_folder = P1_output_folder / "compiled_sequences_near_xerCD"
comp_seq_folder.mkdir(parents=True, exist_ok=True)

indiv_GBK_folder = P1_output_folder / "individual_GBKs"
indiv_GBK_folder.mkdir(parents=True, exist_ok=True)

emboss_hits_XerCD_folder = P1_output_folder / "final_xerCD_hits"       #input files
if not emboss_hits_XerCD_folder.is_dir():
    raise FileNotFoundError(f"Input folder not found:\n{emboss_hits_XerCD_folder}")
cached_database_folder = P1_output_folder / "database_cache"       #cached plasmid database
DATABASE_CACHE = cached_database_folder / 'sequence_database_cache.bin.gz'
if not DATABASE_CACHE.is_file():
    raise RuntimeError(f"sequence_database_cache.bin.gz was not found as a file:\n{DATABASE_CACHE}")

export_folder_1 = comp_seq_folder / "xerC, crossover region, xerD"
export_folder_1.mkdir(parents=True, exist_ok=True)
export_folder_2 = comp_seq_folder / "218+ bp"
export_folder_2.mkdir(parents=True, exist_ok=True)

#----------------------------------------
# 3. xer-motif configuration
#----------------------------------------

HIT_WINDOW_LENGTH = 218
REFERENCE_MOTIF_LENGTH = 28
SLICE_FLANK = 1500
FORWARD_CENTER_OFFSET = -6
REVERSE_CENTER_OFFSET = 0
XER_FEATURE_CENTER_OFFSET = 0
XER_FEATURE_LENGTH = 6



#----------------------------------------
# 4. load cached database
#----------------------------------------

print("Loading Cached Database\n\n")
with gzip.open(DATABASE_CACHE, 'rb') as file:
    sequence_database = pickle.load(file)

print("Database imported\n\n")


#----------------------------------------
# 5. Set up logic to use later when processing xerCD hits
#----------------------------------------

print("Processing XerC/D Hits ... \n\n")

all_hits = []
used_hit_ids = set()
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    # Derive the non-motif flank so each extracted window stays HIT_WINDOW_LENGTH bases long.
window_flank = HIT_WINDOW_LENGTH - REFERENCE_MOTIF_LENGTH



all_xerc_sequences = []
all_crossover_sequences = []
all_xerd_sequences = []
all_central_3kb_sequences = []
all_218bp_sequences = []
collapsed_218bp_sequences = defaultdict(list) # create a version of all_218bp_sequences that combines exact sequence matches


desired_sequence_lengths = [250, 1000, 3000]
all_seq_xerd_down = {size: [] for size in desired_sequence_lengths}
all_seq_xerc_up = {size: [] for size in desired_sequence_lengths}

#----------------------------------------
# 6. Analyze XerCD file names, and loop through XerCD files
#----------------------------------------


for filename in sorted(os.listdir(emboss_hits_XerCD_folder)):
    print(f"analyzing {filename} ... \n\n")
    if not filename.lower().endswith('.txt'):
        continue
    sequences_218bp = []
    try:
        motif_length = int(filename.removesuffix('.txt').split('_')[-2])
    except (IndexError, ValueError):
        print(f" - Skipping input file with an unrecognized name: {filename}")
        continue

    # The filename's penultimate token is motif length; REV selects reverse-strand handling.
    is_reverse = 'REV' in filename.upper()
    direction = 'REV' if is_reverse else 'FOR'


    #----------------------------------------
    # 7. Loop through lines (individual xerCD hits)
    #----------------------------------------

    #read individual hits as seperate lines (accession and coordinate values) from the .txt file
    with open(os.path.join(emboss_hits_XerCD_folder, filename), 'r') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith('#'):
                continue


            #----------------------------------------
            # 8. Fetch corresponding plasmid from database
            #----------------------------------------

            fields = line.split()
            if len(fields) < 2:
                print(f" - Skipping malformed hit in {filename}: {line}")
                continue
            sequence_name, sequence_position_text = fields[:2]
            try:
                sequence_position = int(sequence_position_text)
            except ValueError:
                print(f" - Skipping hit with an invalid position in {filename}: {line}")
                continue
            if sequence_name not in sequence_database:
                print(f"sequence not found in data base for {filename}: {line}")
                continue
            full_sequence = sequence_database[sequence_name]
            hit_id = f"{sequence_name}_{sequence_position}_{direction}_{motif_length}"
            if not full_sequence:
                continue



            #----------------------------------------
            # 9. isolate the 218bp sequences at the xer motif and 3kb upstream/downstream sequences
            #----------------------------------------

            crossover_length = motif_length - 22
            # Tripling lets windows cross the circular sequence boundary without special cases.
            wrapped_sequence = (full_sequence * 3).upper()
            if is_reverse:
                start = sequence_position + len(full_sequence) - 1
                sequence_218bp = wrapped_sequence[start:(start + HIT_WINDOW_LENGTH)]
                sequence_218bp = "".join(complement.get(base, base) for base in reversed(sequence_218bp))

                    #also get the 250bp, 1kb, 3kb sequences
                for size in desired_sequence_lengths: 
                    raw_xerd_down = wrapped_sequence[start - size + 11 : start + 11]
                    seq_xerd_down = "".join(complement.get(base, base) for base in reversed(raw_xerd_down))
                    raw_xerc_up = wrapped_sequence[start + 11 + crossover_length: start + 11 + size + crossover_length]
                    seq_xerc_up = "".join(complement.get(base, base) for base in reversed(raw_xerc_up))

                    all_seq_xerd_down[size].append((hit_id, seq_xerd_down))
                    all_seq_xerc_up[size].append((hit_id, seq_xerc_up)) 
            else: #forward logic
                # Shift the window boundary for this motif length while keeping its total length fixed.
                start = sequence_position - window_flank + len(full_sequence) - 1
                motif_start = start - (REFERENCE_MOTIF_LENGTH - motif_length)
                sequence_218bp = wrapped_sequence[
                    start - (REFERENCE_MOTIF_LENGTH - motif_length):
                    start + window_flank + motif_length
                ]
                    #also get the 250bp, 1kb, 3kb sequences
                for size in desired_sequence_lengths:
                    xerd_start_coord = start + motif_length - 11
                    seq_xerd_down = wrapped_sequence[xerd_start_coord +190: xerd_start_coord + 190 + size]
                    seq_xerc_up = wrapped_sequence[xerd_start_coord - size + 190 - crossover_length: xerd_start_coord + 190  - crossover_length]
                   
                    all_seq_xerd_down[size].append((hit_id, seq_xerd_down))
                    all_seq_xerc_up[size].append((hit_id, seq_xerc_up)) 
            #----------------------------------------
            # 9. give each line a unique
            #----------------------------------------

            duplicate_number = 2
            while hit_id in used_hit_ids:
                hit_id = f"{sequence_name}_{sequence_position}_{direction}_{motif_length}_{duplicate_number}"
                duplicate_number += 1
            used_hit_ids.add(hit_id)
            sequences_218bp.append((f"{sequence_name}_{sequence_position}", sequence_218bp))
            all_218bp_sequences.append((hit_id, sequence_218bp))    
            # Preserve per-hit metadata; an accession can have hits in multiple files.
            all_hits.append((sequence_name, sequence_position, motif_length, is_reverse, hit_id))

            #collect the collapsed 218bp sequences
            collapsed_218bp_sequences[sequence_218bp].append(
                f"{sequence_name}_{sequence_position}")

            #----------------------------------------
            # 10. extract xerD, CR, xerC binding sequences
            #----------------------------------------
            xerd_seq = sequence_218bp[-11:]
            all_xerd_sequences.append((hit_id, xerd_seq))

            crossover_seq = sequence_218bp[-(11 + crossover_length) : -11]
            all_crossover_sequences.append((hit_id, crossover_seq))
            
            xerc_seq = sequence_218bp[-(22 + crossover_length) : -(11 + crossover_length)]
            all_xerc_sequences.append((hit_id, xerc_seq))


#----------------------------------------
# 11. Write the output file names for these
#----------------------------------------

#the rest aren't run in the loop because I only want one file
with open(os.path.join(export_folder_2, 'central_218bp_sequences.txt'), 'w') as f:
    for name, seq in all_218bp_sequences:
        f.write(f">{name}\n{seq}\n")


with open(os.path.join(export_folder_2, 'central_218bp_sequences_collapsed.txt'), 'w') as f:
    for seq, names in collapsed_218bp_sequences.items():
        combined_name = "|".join(names)
        f.write(f">{combined_name}\n{seq}\n")


with open(os.path.join(export_folder_1, 'all_xerC_sites.txt'), 'w') as f:
    for name, seq in all_xerc_sequences:
        f.write(f">{name}\n{seq}\n")

with open(os.path.join(export_folder_1, 'all_crossover_regions.txt'), 'w') as f:
    for name, seq in all_crossover_sequences:
        f.write(f">{name}\n{seq}\n")

with open(os.path.join(export_folder_1, 'all_xerD_sites.txt'), 'w') as f:
    for name, seq in all_xerd_sequences:
        f.write(f">{name}\n{seq}\n")

for size in desired_sequence_lengths:
    with open(os.path.join(export_folder_2, f"XerD_and_down_{size}bp_sequences.txt"), 'w') as f:
        for name, seq in all_seq_xerd_down[size]:
            f.write(f">{name}\n{seq}\n")

    with open(os.path.join(export_folder_2, f"XerC_and_up_{size}bp_sequences.txt"), 'w') as f:
        for name, seq in all_seq_xerc_up[size]:
            f.write(f">{name}\n{seq}\n")

#----------------------------------------
# 10. save full sequences
#----------------------------------------
os.makedirs(os.path.join(indiv_GBK_folder, 'full sequence'), exist_ok=True)
# Fetch each accession once; the resulting record is reused for all of its hits.
print(" - Saving full sequence of hit GBKs ... (may take 20+ minutes)")

sequence_names = dict.fromkeys(hit[0] for hit in all_hits)
total = len(sequence_names)
for i, sequence_name in enumerate(sequence_names, start=1):
    print(f"\r       Progress: {i:,}/{total:,}", end="", flush = True)
    gbk_path = os.path.join(indiv_GBK_folder, 'full sequence', f'{sequence_name}.gbk')
    if not os.path.exists(gbk_path):
        try:
            handle = Entrez.efetch(db="nucleotide", id=sequence_name, rettype="gbwithparts", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            with open(gbk_path, 'w') as file:
                SeqIO.write(record, file, "genbank")
        except Exception as e:
            print(f" - Failed to download GenBank file for {sequence_name}: {e}")



#----------------------------------------
# 10. set up logic to get 3kb central sequences
#----------------------------------------
os.makedirs(os.path.join(indiv_GBK_folder, '3kb'), exist_ok=True)
print("\n\n - Slicing Sequences ... (may take a few minutes)\n\n")
for sequence_name, position, motif_length, is_reverse, hit_id in all_hits:
    gbk_path = os.path.join(indiv_GBK_folder, 'full sequence', f'{sequence_name}.gbk')
    if not os.path.exists(gbk_path):
        print(f" - File not found for {sequence_name}")
        continue

    record = SeqIO.read(gbk_path, "genbank")
    original_sequence = record.seq
    sequence_length = len(original_sequence)
    if not sequence_length:
        print(f" - Skipping empty sequence for {sequence_name}")
        continue

    # Convert scanner positions to the center convention used by the 3 kb slice.
    center_offset = REVERSE_CENTER_OFFSET if is_reverse else motif_length + FORWARD_CENTER_OFFSET
    sequence_number = position + center_offset
    # Repeat enough copies to keep the full slice in bounds for short circular sequences.
    copies_each_side = (SLICE_FLANK + sequence_length - 1) // sequence_length
    copy_count = copies_each_side * 2 + 1
    middle_copy = copy_count // 2
    new_sequence = original_sequence * copy_count
    start_range = sequence_number - SLICE_FLANK + middle_copy * sequence_length
    end_range = sequence_number + SLICE_FLANK + middle_copy * sequence_length
    trimmed_sequence = new_sequence[start_range:end_range]

    if "REV" in hit_id.upper():
        central_sequence = str(trimmed_sequence.reverse_complement())
    elif "FOR" in hit_id.upper():
        central_sequence = str(trimmed_sequence)
    else:
        raise ValueError(f"expected For or Rev in file name but got: {hit_id}")

    all_central_3kb_sequences.append(
        (hit_id, central_sequence))


    #----------------------------------------
    # 11. annotate the xer site
    #----------------------------------------

    adjusted_features = []
    for copy_index in range(copy_count):
        shift = copy_index * sequence_length
        for feature in record.features:
            source_parts = feature.location.parts if isinstance(feature.location, CompoundLocation) else [feature.location]
            spans_origin = (
                isinstance(feature.location, CompoundLocation)
                and min(int(part.start) for part in source_parts) == 0
                and max(int(part.end) for part in source_parts) == sequence_length
            )
            clipped_parts = []
            for part in source_parts:
                # Move the zero-coordinate part forward one copy to preserve origin-spanning joins.
                part_shift = shift + sequence_length if spans_origin and int(part.start) == 0 else shift
                feature_start = int(part.start) + part_shift
                feature_end = int(part.end) + part_shift
                # Clip feature parts to this window, then convert them to slice-relative coordinates.
                clipped_start = max(feature_start, start_range)
                clipped_end = min(feature_end, end_range)
                if clipped_start < clipped_end:
                    clipped_parts.append(FeatureLocation(
                        clipped_start - start_range,
                        clipped_end - start_range,
                        strand=part.strand,
                    ))
            if not clipped_parts:
                continue
            if len(clipped_parts) == 1:
                new_location = clipped_parts[0]
            else:
                new_location = CompoundLocation(
                    clipped_parts,
                    operator=feature.location.operator,
                )
            adjusted_features.append(SeqFeature(
                new_location,
                type=feature.type,
                qualifiers=feature.qualifiers.copy(),
            ))
    
    #----------------------------------------
    # 12. annotate the xer site
    #----------------------------------------
    final_record = SeqRecord(
        trimmed_sequence,
        id=record.id + "_3000bp",
        features=[],
        description=f"3000bp segment around Xer site {position}",
    )
    final_record.annotations["molecule_type"] = "DNA"

    # Mark the expected Xer site; this placeholder is not a new motif-search result.
    xer_motif_start = SLICE_FLANK + XER_FEATURE_CENTER_OFFSET
    xer_motif_location = FeatureLocation(
        xer_motif_start,
        xer_motif_start + XER_FEATURE_LENGTH,
        strand=1,
    )
    xer_motif_feature = SeqFeature(xer_motif_location, type="gene", qualifiers={
        "pseudo": "",
        "locus_tag": "XM_001",
        "note": "Xer Motif placeholder"
    })
    adjusted_features.append(xer_motif_feature)
    xer_motif_feature = SeqFeature(xer_motif_location, type="CDS", qualifiers={
        "locus_tag": "XER",
        "note": "Xer Motif placeholder",
        "pseudo": "",
        "codon_start": "1",
        "transl_table": "11",
        "product": "Xer Motif placeholder",
    })

    #----------------------------------------
    # 13. save the 3kb central sequence file
    #----------------------------------------
    adjusted_features.append(xer_motif_feature)
    adjusted_features.sort(key=lambda feature: feature.location.start)
    final_record.features = adjusted_features

    with open(os.path.join(indiv_GBK_folder, '3kb', f'{hit_id}.gbk'), "w") as output_handle:
        SeqIO.write(final_record, output_handle, "genbank")

with open(os.path.join(export_folder_2, 'all_central_3kb_sequences_sorted.txt'), 'w') as f:
    for name, seq in all_central_3kb_sequences:
        f.write(f">{name}\n{seq}\n")

print("\n\n    ~~  Finished!! you have run the script 103_extract_sequences.py  ~~\n\n")
