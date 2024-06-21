#!/usr/bin/env python3

import os
import gzip
import pickle
from Bio import Entrez, SeqIO
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
    print("\nLoading Prebuilt Database\n")
    with gzip.open('sequence_database.bin.gz', 'rb') as file: # if there is already a prebuilt database, load it from the pickle file
        sequence_database = pickle.load(file)

print("\nDatabase imported\n")

# get the input files ready
input_dir = 'inputfiles'
output_dir = 'outputfiles'
os.makedirs(output_dir, exist_ok=True)

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
                    sequences.append((f"{sequence_name}_{sequence_position}", sequence)) # store the sequence name + length and the sequence in the output array

    # Write the sequences to an output file
    with open(os.path.join(output_dir, 'output_' + filename), 'w') as file:
        for name, sequence in sequences:
            file.write(f">{name}\n{sequence}\n")
            
os.makedirs('outputfiles/gbk', exist_ok=True)
os.makedirs('outputfiles/gff3', exist_ok=True)

# Extract sequence names from sequences
sequence_names = [name.rsplit('_', 1)[0] for name, _ in sequences]
total_sequences = len(sequence_names)

for i, sequence_name in enumerate(sequence_names):
    # Define file paths for gbk, gff3, and fasta files
    gbk_file_path = os.path.join(output_dir, 'gbk', f'{sequence_name}.gbk')
    gff3_file_path = os.path.join(output_dir, 'gff3', f'{sequence_name}.gff3')
    fasta_file_path = os.path.join(output_dir, 'gff3', f'{sequence_name}.fasta')

    # Check if the gbk file already exists, if not download it
    if not os.path.exists(gbk_file_path):
        try:
            handle = Entrez.efetch(db="nucleotide", id=sequence_name, rettype="gbwithparts", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            with open(gbk_file_path, 'w') as f:
                SeqIO.write(record, f, "genbank")
        except Exception as e:
            print(f"Failed to download GenBank file for {sequence_name}: {e}")

    # Check if the gff3 file already exists, if not download it
    if not os.path.exists(gff3_file_path):
        try:
            handle = Entrez.efetch(db="nucleotide", id=sequence_name, rettype="gff3", retmode="text")
            data = handle.read()
            handle.close()

            with open(gff3_file_path, 'w') as f:
                f.write(data)
        except Exception as e:
            print(f"Failed to download GFF3 file for {sequence_name}: {e}")

    # Check if the fasta file already exists, if not download it
    if not os.path.exists(fasta_file_path):
        try:
            handle = Entrez.efetch(db="nucleotide", id=sequence_name, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()

            with open(fasta_file_path, 'w') as f:
                SeqIO.write(record, f, "fasta")
        except Exception as e:
            print(f"Failed to download FASTA file for {sequence_name}: {e}")

    # Print progress
    print(f'Progress: {100 * (i+1) / total_sequences:.2f}%')

# Run the clinker command
print ("Running Clinker")
subprocess.run(f'clinker {output_dir}/gff3/*.gff3 -p cluster.html', shell=True, check=True)