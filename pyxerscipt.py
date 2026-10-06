#!/usr/bin/env python3

import os
import glob
import gzip
import pickle
from Bio import Entrez, SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature, CompoundLocation

from Bio.SeqRecord import SeqRecord
import subprocess

# User-configurable settings
ENTREZ_EMAIL = "paul.capewell@glasgow.ac.uk"
INPUT_DIR = 'inputfiles'
OUTPUT_DIR = 'outputfiles'
FASTA_DATABASE = 'bacterial_plasmids_curated.fa.gz'
DATABASE_CACHE = 'sequence_database.bin.gz'
HIT_WINDOW_LENGTH = 218
REFERENCE_MOTIF_LENGTH = 28
SLICE_FLANK = 1500
FORWARD_CENTER_OFFSET = -6
REVERSE_CENTER_OFFSET = 0
XER_FEATURE_CENTER_OFFSET = 0
XER_FEATURE_LENGTH = 6

Entrez.email = ENTREZ_EMAIL

# Load the FASTA database or its cached dictionary.
if not os.path.exists(DATABASE_CACHE):
    print(f'\rBuilding New Database, Progress: 0.00%', end='')
    sequence_database = {}
    with gzip.open(FASTA_DATABASE, 'rt') as file:
        total_lines = sum(1 for _ in file)
    with gzip.open(FASTA_DATABASE, 'rt') as file:
        name = ''
        sequence = ''
        for i, line in enumerate(file):
            print(f'\rBuilding New Database, Progress: {100 * (i+1) / total_lines:.2f}%', end='')
            if line.startswith('>'):
                if name:
                    sequence_database[name] = sequence
                name = line.strip()[1:]
                sequence = ''
            else:
                sequence += line.strip()
        if name:
            sequence_database[name] = sequence

    with gzip.open(DATABASE_CACHE, 'wb') as file:
        pickle.dump(sequence_database, file)
else:
    print("Loading Prebuilt Database\n")
    with gzip.open(DATABASE_CACHE, 'rb') as file:
        sequence_database = pickle.load(file)

print("Database imported\n")
os.makedirs(OUTPUT_DIR, exist_ok=True)

all_hits = []
used_hit_ids = set()
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
# Derive the non-motif flank so each extracted window stays HIT_WINDOW_LENGTH bases long.
window_flank = HIT_WINDOW_LENGTH - REFERENCE_MOTIF_LENGTH
for filename in sorted(os.listdir(INPUT_DIR)):
    if not filename.lower().endswith('.txt'):
        continue
    sequences = []
    try:
        motif_length = int(filename.split()[-2])
    except (IndexError, ValueError):
        print(f"Skipping input file with an unrecognized name: {filename}")
        continue
    # The filename's penultimate token is motif length; REV selects reverse-strand handling.
    is_reverse = 'REV' in filename.upper()
    direction = 'REV' if is_reverse else 'FOR'
    with open(os.path.join(INPUT_DIR, filename), 'r') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) < 2:
                print(f"Skipping malformed hit in {filename}: {line}")
                continue
            sequence_name, sequence_position_text = fields[:2]
            try:
                sequence_position = int(sequence_position_text)
            except ValueError:
                print(f"Skipping hit with an invalid position in {filename}: {line}")
                continue
            if sequence_name not in sequence_database:
                continue

            full_sequence = sequence_database[sequence_name]
            if not full_sequence:
                continue
            # Tripling lets windows cross the circular sequence boundary without special cases.
            wrapped_sequence = (full_sequence * 3).upper()
            if is_reverse:
                start = sequence_position + len(full_sequence) - 1
                sequence = wrapped_sequence[start:(start + HIT_WINDOW_LENGTH)]
                sequence = "".join(complement.get(base, base) for base in reversed(sequence))
            else:
                # Shift the window boundary for this motif length while keeping its total length fixed.
                start = sequence_position - window_flank + len(full_sequence) - 1
                sequence = wrapped_sequence[
                    start - (REFERENCE_MOTIF_LENGTH - motif_length):
                    start + window_flank + motif_length
                ]

            hit_id = f"{sequence_name}_{sequence_position}_{direction}_{motif_length}"
            duplicate_number = 2
            while hit_id in used_hit_ids:
                hit_id = f"{sequence_name}_{sequence_position}_{direction}_{motif_length}_{duplicate_number}"
                duplicate_number += 1
            used_hit_ids.add(hit_id)
            sequences.append((f"{sequence_name}_{sequence_position}", sequence))
            # Preserve per-hit metadata; an accession can have hits in multiple files.
            all_hits.append((sequence_name, sequence_position, motif_length, is_reverse, hit_id))
    with open(os.path.join(OUTPUT_DIR, 'output_' + filename), 'w') as file:
        for name, sequence in sequences:
            file.write(f">{name}\n{sequence}\n")

os.makedirs(os.path.join(OUTPUT_DIR, 'gbk'), exist_ok=True)
os.makedirs(os.path.join(OUTPUT_DIR, 'sliced_gbk'), exist_ok=True)

# Fetch each accession once; the resulting record is reused for all of its hits.
sequence_names = dict.fromkeys(hit[0] for hit in all_hits)
for sequence_name in sequence_names:
    gbk_path = os.path.join(OUTPUT_DIR, 'gbk', f'{sequence_name}.gbk')
    if not os.path.exists(gbk_path):
        try:
            handle = Entrez.efetch(db="nucleotide", id=sequence_name, rettype="gbwithparts", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            with open(gbk_path, 'w') as file:
                SeqIO.write(record, file, "genbank")
        except Exception as e:
            print(f"Failed to download GenBank file for {sequence_name}: {e}")

print("Slicing Sequences\n")
for sequence_name, position, motif_length, is_reverse, hit_id in all_hits:
    gbk_path = os.path.join(OUTPUT_DIR, 'gbk', f'{sequence_name}.gbk')
    if not os.path.exists(gbk_path):
        print(f"File not found for {sequence_name}")
        continue

    record = SeqIO.read(gbk_path, "genbank")
    original_sequence = record.seq
    sequence_length = len(original_sequence)
    if not sequence_length:
        print(f"Skipping empty sequence for {sequence_name}")
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
    adjusted_features.append(xer_motif_feature)
    adjusted_features.sort(key=lambda feature: feature.location.start)
    final_record.features = adjusted_features

    with open(os.path.join(OUTPUT_DIR, 'sliced_gbk', f'{hit_id}.gbk'), "w") as output_handle:
        SeqIO.write(final_record, output_handle, "genbank")

print("Running Clinker\n")
gbk_files = sorted(glob.glob(os.path.join(OUTPUT_DIR, 'sliced_gbk', '*.gbk')))
if not gbk_files:
    raise FileNotFoundError("No sliced GenBank files found for Clinker")
subprocess.run([
    'clinker',
    *gbk_files,
    '--force',
    '-p', os.path.join(OUTPUT_DIR, 'clusters.html'),
    '-o', os.path.join(OUTPUT_DIR, 'output_summary.txt'),
    '-mo', os.path.join(OUTPUT_DIR, 'matrix.csv'),
], check=True)