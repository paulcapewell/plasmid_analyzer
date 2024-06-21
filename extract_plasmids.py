#!/usr/bin/env python3
import os
import gzip
import pickle


def load_fasta_into_dict(filename):
    fasta_dict = {}
    with open(filename, 'r') as file:
        for line in file:
            line = line.rstrip()
            if line.startswith('>'):
                seq_id = line[1:]
                fasta_dict[seq_id] = ''
            else:
                fasta_dict[seq_id] += line
    return fasta_dict

test_seqs = load_fasta_into_dict('outputfiles/unique_default.fa')

def rename_keys(test_seqs):
    new_test_seqs = {}
    for key in test_seqs:
        new_key = '_'.join(key.split('_')[:2])
        new_test_seqs[new_key] = test_seqs[key]
    return new_test_seqs

test_seqs = rename_keys(test_seqs)

with gzip.open('sequence_database.bin.gz', 'rb') as file: # if there is already a prebuilt database, load it from the pickle file
    sequence_database = pickle.load(file)

def export_matching_sequences(sequence_database, test_seqs, output_filename):
    with open(output_filename, 'w') as output_handle:
        for key, seq in sequence_database.items():
            if key in test_seqs:
                output_handle.write(f">{key}\n{seq}\n")

export_matching_sequences(sequence_database, test_seqs, 'reduced_database.fa')