#!/usr/bin/env python3

import os
import gzip
import pickle
from Bio import Entrez, SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
import subprocess

Entrez.email = "paul.capewell@glasgow.ac.uk"

# if there isn't a pre-built pickle file, load in the sequence database and build a pickle binary for quick loading in future
if not os.path.exists('sequence_database.bin.gz'):
    print(f'\rBuilding New Database, Progress: 0.00%', end='') # empty progress report
    sequence_database = {} # set an emptry database
    #with open('nucleotideseq.fa', 'r') as file: # open the database file if a text file
    with gzip.open('bacterial_plasmids_curated.fa.gz', 'rt') as file:
        total_lines = sum(1 for _ in file) # get total size of database for the progress update
    with gzip.open('bacterial_plasmids_curated.fa.gz', 'rt') as file: # open the database file if a gzip file
        name = '' # var to store current sequence name from the file
        sequence = '' # var to store current sequence  from the file
            
        for i, line in enumerate(file): # go through each line until we get to the next fasta header >
            print(f'\rBuilding New Database, Progress: {100 * (i+1) / total_lines:.2f}%', end='') # Print progress
            if line.startswith('>'): # if the line contains a > character, it is a sequence name
                if name:  # save the previous name and sequence
                    sequence_database[name] = sequence # store the sequence in the database
                name = line.strip()[1:]  # remove the '>' character for storing the name
                sequence = ''  # reset the sequence var
            else:
                sequence += line.strip() # add the sequence to the sequence var with line breaks removed
        if name:  # make sure to save the last name and sequence and check it's lengths
            sequence_database[name] = sequence

    with gzip.open('sequence_database.bin.gz', 'wb') as file:     # Save the database into a pickle file
        pickle.dump(sequence_database, file)
else:
    print("Loading Prebuilt Database\n")
    with gzip.open('sequence_database.bin.gz', 'rb') as file: # if there is already a prebuilt database, load it from the pickle file
        sequence_database = pickle.load(file)

print("Database imported\n")
# get the input files ready
input_dir = 'inputfiles'
output_dir = 'outputfiles'
os.makedirs(output_dir, exist_ok=True)

all_sequences=[] # store all sequences in a list
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} # lookup table for complementing bases
for filename in os.listdir(input_dir):
    sequences = [] # reset sequences list for each file
    length = int(filename.split(' ')[-2])  # parse the filename and extract the required length
    with open(os.path.join(input_dir, filename), 'r') as file:
        for line in file:
            if not line.startswith('#'): # skip any header files
                sequence_name, sequence_position, *rest = line.split() # split the line into name and position
                sequence_position = int(sequence_position) # make sure positon is an integer
                if sequence_name in sequence_database:
                    full_sequence = sequence_database[sequence_name]
                    # Wrap around the sequence string if necessary
                    wrapped_sequence = (full_sequence + full_sequence + full_sequence).upper() # add the sequence to itself 3 times to make wrapped file. Also mke it uppercase to make sure bases consistent with lookup table
                    if 'REV' in filename:
                        start = sequence_position + len(full_sequence) - 1  # if a reverse file, start at the end of the sequence
                        sequence = wrapped_sequence[start:(start+218)] # store 190+length bases from the start position, allowing for wrap around
                        sequence = "".join(complement.get(base, base) for base in reversed(sequence)) # reverse complement the sequence if it is a reverse file
                    else:
                        start = sequence_position - 190 + len(full_sequence) - 1  # if a forward file, start at the beginning of the sequence -190
                        sequence = wrapped_sequence[start-(28-length):(start+190+length)] # store 190+length bases from the start position (correcting for specific length), allowing for wrap around
                    sequences.append((f"{sequence_name}_{sequence_position}", sequence)) # store the sequence name + length and the sequence in the single input file output array
                    all_sequences.append((f"{sequence_name}_{sequence_position}", sequence)) # store the sequence name + length and the sequence in the overall output array
    # Write the sequences to an output file (per input file)
    with open(os.path.join(output_dir, 'output_' + filename), 'w') as file:
        for name, sequence in sequences:
            file.write(f">{name}\n{sequence}\n")

# check required output files exist            
os.makedirs('outputfiles/gbk', exist_ok=True)
os.makedirs('outputfiles/sliced_gbk', exist_ok=True)

# Prepare the list of sequences to be examined (from all input files)
name_number_mapping = {} 
for sequence in all_sequences:
    parts = sequence[0].rsplit('_', 1)  # Split from the right to get the last part separately
    if len(parts) == 2:  # Ensure the sequence name contains an underscore
        sequence_name, sequence_number_str = parts
        try:
            sequence_number = int(sequence_number_str)  # Convert sequence number to integer
            name_number_mapping[sequence_name] = sequence_number # store the sequence name and motif positon so they can be looked up later
        except ValueError:
            print(f"Skipping sequence with invalid number: {sequence}")
    else:
        print(f"Skipping invalid sequence name: {sequence}")

# Download the gbk files 
for i, sequence_name in enumerate(name_number_mapping.keys()):
    if not os.path.exists(os.path.join(output_dir, 'gbk', f'{sequence_name}.gbk')): # Check if the gbk file already exists, if not download it using Biopython
        try:
            handle = Entrez.efetch(db="nucleotide", id=sequence_name, rettype="gbwithparts", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            with open(os.path.join(output_dir, 'gbk', f'{sequence_name}.gbk'), 'w') as f:
                SeqIO.write(record, f, "genbank")
        except Exception as e:
            print(f"Failed to download GenBank file for {sequence_name}: {e}") # print an rror message if download fails

#slice sequences to focus analysis on the motif
print("Slicing Sequences\n")
for sequence_name, sequence_number in name_number_mapping.items():
    if os.path.exists(os.path.join(f"outputfiles/gbk/{sequence_name}.gbk")):
        record = SeqIO.read(os.path.join(f"outputfiles/gbk/{sequence_name}.gbk"), "genbank")
        sequence = record.seq
        with open(os.path.join(f"outputfiles/gbk/{sequence_name}.gbk"), "r") as input_handle:
            record = SeqIO.read(input_handle, "genbank")
        sliced_record = record[sequence_number - 1501:sequence_number + 1500]
        # Adjust the features to match the sliced sequence
        new_features = []
        for feature in record.features:
            if feature.location.end >= sequence_number - 1500 and feature.location.start <= sequence_number + 1500:
                # Adjust the feature location
                start = sequence_number - 1500
                end = sequence_number + 1500
                new_start = max(0, feature.location.start - start)
                new_end = min(end - start, feature.location.end - start)
                new_location = FeatureLocation(new_start, new_end, strand=feature.location.strand)
                new_feature = SeqFeature(new_location, type=feature.type, qualifiers=feature.qualifiers)
                new_features.append(new_feature) 
        sliced_record.features = new_features     # Update the features of the sliced record
        
        with open(os.path.join(f"outputfiles/sliced_gbk/{sequence_name}.gbk"), "w") as output_handle:
            SeqIO.write(sliced_record, output_handle, "genbank")
        #slice_gbk(file_path,output_path,sequence_number - 1500, sequence_number + 1500)
    else:
        print(f"File not found for {sequence_name}")

# Run the clinker command
print ("Running Clinker\n")
subprocess.run(f'clinker ./outputfiles/sliced_gbk/*.gbk --force -p ./outputfiles/clusters.html -o outputfiles/output_summary.txt -mo outputfiles/matrix.csv', shell=True, check=True)