#!/usr/bin/env python3

import os

# Load in the sequence database
print("Reading database\n")
sequence_database = {} # set an emptry database
max_length = -1 # set max length as -1 so that we can find the longest in the loop
max_name = '' # name of the longest sequence
with open('nucleotideseq.fa', 'r') as file: # open the database file if a text file
    name = '' # var to store current sequence name from the file
    sequence = '' # var to store current sequence  from the file
    for line in file: # go through each line until we get to the next fasta header >
        if line.startswith('>'): # if the line contains a > character, it is a sequence name
            if name:  # save the previous name and sequence
                sequence_database[name] = sequence # store the sequence in the database
                if len(sequence) > max_length: # check to see if the sequence is the longest
                    max_length = len(sequence) # if it is, store it ...
                    max_name = name # ... and store the name
            name = line.strip()[1:]  # remove the '>' character for storing the name
            sequence = ''  # reset the sequence var
        else:
            sequence += line.strip() # add the sequence to the sequence var with line breaks removed
    if name:  # make sure to save the last name and sequence and check it's lengths
        sequence_database[name] = sequence
        if len(sequence) > max_length:
            max_length = len(sequence)
            max_name = name

# output some info to show database has been imported
print("Database imported\n")
print("Max size found:", max_length,"nt (",max_name,")")

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
