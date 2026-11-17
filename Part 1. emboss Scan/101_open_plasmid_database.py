print("\n\n\n**Running Script: 101_open_plasmid_database.py**\n\n")
from pathlib import Path
import os
import glob
import gzip
import pickle
from Bio import Entrez, SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature, CompoundLocation
from Bio.SeqRecord import SeqRecord
import subprocess
import hashlib
from datetime import datetime







#name of plasmid database file you are analyzing
DATABASE_NAME = 'bacterial_plasmids_curated.fa.gz'








#step 2. define file paths

P1_emboss_folder = Path(__file__).resolve().parent                              #the folder for all of P1_emboss_folder
input_folder = P1_emboss_folder / "input_folder"                        #input folder
output_folder = P1_emboss_folder / "output_folder"                      #output folder
output_folder.mkdir(parents=True, exist_ok=True)

plasmid_database_folder = input_folder / "original_plasmid_database"        #plasmid database folder
ORIGINAL_DATABASE = plasmid_database_folder / DATABASE_NAME                       #plasmid database file

cached_database_folder = output_folder / "cached_plasmid_database"          #cached plasmid database folder
cached_database_folder.mkdir(parents=True, exist_ok=True)
DATABASE_CACHE = cached_database_folder / 'sequence_database_cache.bin.gz'   #cached plasmid database file
file_id_folder = output_folder / "file_identity_and_integrity"
file_id_folder.mkdir(parents=True, exist_ok=True)

#check whether original or cached plasmid database exist
if os.path.exists(DATABASE_CACHE):
    print("Cached file already exists in `output_folder / cached_database_folder` . No need to run this script")
    raise SystemExit
if not os.path.exists(ORIGINAL_DATABASE):
    print("Original database not located. Please double check the .fa.gz file's location and name")
    raise SystemExit
print(" - Original Database File Located \n\n - No Cached Database Located\n\n**Creating New Cached Database**\n       (requires interent and may take 30+ minutes)\n")


#collect some other useful data about the original database, that I'll put in a receipt.txt file later on in case other errors come up later
record_count = 0
duplicate_headers = 0
duplicate_names = []
total_bases = 0
sha256 = hashlib.sha256()
with open(ORIGINAL_DATABASE, "rb") as file:
    for chunk in iter(lambda: file.read(1024 * 1024), b""):
        sha256.update(chunk)
database_checksum = sha256.hexdigest()


# create a new cached database 
sequence_database = {}

print(" - Counting lines in FASTA database...\n")
with gzip.open(ORIGINAL_DATABASE, 'rt') as file:
    total_lines = sum(1 for _ in file)

print(" - Finished counting lines. Now building sequence dictionary...\n")
with gzip.open(ORIGINAL_DATABASE, 'rt') as file:
    name = ''
    sequence = ''
    for i, line in enumerate(file):
        print(f'\r - Building New Database, Progress: {100 * (i+1) / total_lines:.2f}%', end='')
        if line.startswith('>'):
            if name:
                record_count += 1
                total_bases += len(sequence)

                if name in sequence_database:
                    duplicate_headers += 1
                    duplicate_names.append(name)
                
                sequence_database[name] = sequence
            name = line.strip()[1:]
            sequence = ''
        else:
            sequence += line.strip()
    if name:
        record_count += 1
        total_bases += len(sequence)
        
        if name in sequence_database:
            duplicate_headers += 1
            duplicate_names.append(name)
                        
        sequence_database[name] = sequence

print("\n\n - Finished building dictionary. Now saving compressed database cache... \n       (very slow might take 25+ minutes)\n")
with gzip.open(DATABASE_CACHE, 'wb') as file:
    pickle.dump(sequence_database, file)
print(" - Finished saving database cache.\n\n**Verifying the Cache can be Opened**")


#verify the cache seems correct and make a SHA256
with gzip.open(DATABASE_CACHE, "rb") as file:
    verify_database = pickle.load(file)
if len(verify_database) != len(sequence_database):
    raise ValueError(
        f" - Verification failed. original dictionary contained {len(sequence_database)} entries but related cache contained {len(verify_database)} entries.")
print(" - Cache verification Successful\n")

cached_sha256 = hashlib.sha256()
with open(DATABASE_CACHE, "rb") as file:
    for chunk in iter(lambda: file.read(1024 * 1024), b""):
        cached_sha256.update(chunk)
cached_database_checksum = cached_sha256.hexdigest()



#make a receipt.txt file


receipt_file = file_id_folder / "receipt.txt"
receipt_text = f"""
Database Cache Creation Report
------------------------------

Creation Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

Original Database:
    file - {ORIGINAL_DATABASE.name}
    path - {ORIGINAL_DATABASE}
    Size in bytes = {ORIGINAL_DATABASE.stat().st_size:,}
    SHA256 = {database_checksum} 
       This SHA256 is unique to the original database's file contents
       it will change if any edit is made at all to the file content

Reading the Original Database:
    number of FASTA Records Encountered: {record_count:,}
    unique dictionary entries: {len(sequence_database):,}
    Total Bases Stored: {total_bases:,}
    duplicate headers: {duplicate_headers:,}
    First 20 duplicate headers, if detected:
        {chr(10).join(duplicate_names[:20])}
    
    

Cache File:
        File Was Verified to successfully reload after being written
        File's dictionary size was verified to match the origianl 
    File: {DATABASE_CACHE.name}
    Size in bytes = {DATABASE_CACHE.stat().st_size:,}
    SHA256 = {cached_database_checksum}
        this SHA256 is unique to the cached database's file contents
        it will change if any edit is made at all to the file content
"""

receipt_file.write_text(receipt_text)
print(f" - Receipt saved to:\n      {receipt_file}")
print("\n\n    ***  Script Finished!!  ***")
