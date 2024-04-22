#!/usr/bin/env python3
# QC the output files to check they contain the same data as the example files.

import os

output_dir = 'outputfiles'
output_dir_example = 'outputfiles_example'
test_dir = 'validation'
os.makedirs(test_dir, exist_ok=True)

for filename in os.listdir(output_dir):
    # Only process files that exist in both directories
    if filename in os.listdir(output_dir_example):
        new_files = {}
        old_files = {}
        with open(os.path.join(output_dir, filename), 'r') as file:
            name = ''
            sequence = ''
            for line in file:
                if line.startswith('>'):
                    if name: 
                        new_files[name] = sequence
                    name = line.strip()[1:]  # remove the '>' character
                    sequence = ''  # reset the sequence
                else:
                    sequence += line.strip()
            if name:  # save the last name and sequence
                new_files[name] = sequence
        # open the complment new data file and store the sequence info
        with open(os.path.join(output_dir_example, filename), 'r') as file:
            name = ''
            sequence = ''
            for line in file:
                if line.startswith('>'):
                    if name:  # save the previous name and sequence
                        old_files[name] = sequence
                    name = line.strip()[1:]  # remove the '>' character
                    sequence = ''  # reset the sequence
                else:
                    sequence += line.strip()
            if name:  # save the last name and sequence
                old_files[name] = sequence

        with open(os.path.join(test_dir, 'test_' + filename), 'w') as file:
            for name in new_files:
                file.write(f"{name} {new_files[name] == old_files.get(name, '')}\n") # do the comparison and write the result to the file