print("\n\n\n** Running Script: 102_xer_matrix_scan.py **\n\n")
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
from pathlib import Path


#define file paths and locate files

P1_emboss_folder = Path(__file__).resolve().parent                      #the folder for all of P1_emboss_folder
input_folder = P1_emboss_folder / "input_folder"                        #input folder
output_folder = P1_emboss_folder / "output_folder"                      #output folder
file_id_folder = output_folder / "file_identity_and_integrity"
file_id_folder.mkdir(parents=True, exist_ok=True)

xerCD_hits_folder = output_folder / "xerCD hits before ArgR"
xerCD_hits_folder.mkdir(parents=True, exist_ok=True)
temp_folder = xerCD_hits_folder / "temporary files used in script"
temp_folder.mkdir(parents=True, exist_ok=True)

print("*locating necessary files*\n")
cached_database_folder = output_folder / "cached_plasmid_database"          #cached plasmid database folder
DATABASE_CACHE = cached_database_folder / 'sequence_database_cache.bin.gz'   #cached plasmid database file
if not DATABASE_CACHE.exists():
    print("Cached database file not found. cannot run this script")
    raise SystemExit
with gzip.open(DATABASE_CACHE, "rb") as file:
    plasmid_database = pickle.load(file)
print(f" - cached database found and loaded successfully. Number of plasmids:\n       {len(plasmid_database):,}")


xer_matrix_folder = input_folder / "Xer_frequency_matrices"                 #matrices folder
matrix_files = sorted(xer_matrix_folder.glob("*.txt"))                      #frequency matrices
print("\n - List of matrix files found:")
for file in matrix_files:
    print(f"       {file.name}")
if len(matrix_files) != 6:
    raise ValueError(f"Expected exactly 6 matrix files, found {len(matrix_files)}")
FOR_matrix_6bp = xer_matrix_folder / "Prophecy_matrix_6bp_FOR.txt"

#check that EMBOSS is working
Emboss_verify_result = subprocess.run(
    ["profit", "-help"],
    capture_output=True, text=True)
if Emboss_verify_result.returncode != 0:
    raise ValueError(" - EMBOSS PROFIT could not be executed. Double check the installation.")
print("\n - EMBOSS PROFIT is installed and can be executed")


# create SHA256 receipt for all matrix files
sha256_file = file_id_folder / "xer matrices sha256.txt"
with open(sha256_file, "w") as outfile:
    outfile.write("SHA256 Checksums for Matrix Files\n")
    outfile.write("--------------------------------\n\n")
    for matrix_file in matrix_files:
        sha256 = hashlib.sha256()
        with open(matrix_file, "rb") as file:
            for chunk in iter(lambda: file.read(1024 * 1024),b""):
                sha256.update(chunk)
        outfile.write(f"{matrix_file.name}\n")
        outfile.write(f"SHA256: {sha256.hexdigest()}\n\n")

print(f"\n*Matrix SHA256 report written to:*\n - {sha256_file.name}\n - in {file_id_folder.name}")



#loop through all 6 xer matrix files
for xer_matrix_file in matrix_files:
    print(f"\n\n\n*Running PROFIT on {xer_matrix_file.name}*")


    #for each matrix, set up a unique output file based on the unique name of the matrix file
    xer_matrix_file_name_parts = xer_matrix_file.stem.split("_")
    xer_matrix_file_length = xer_matrix_file_name_parts[2]
    modified_matrix_length = int(xer_matrix_file_length.replace("bp", "")) + 22
    xer_matrix_file_orientation = xer_matrix_file_name_parts[3]
    output_name = f"HITS XerCD {modified_matrix_length} {xer_matrix_file_orientation}.txt"
    HITS_XerCD_txt_File = xerCD_hits_folder / output_name
    print(f" - outputs will be written in this file: {output_name}")

    #write in the output.txt file
    with open(HITS_XerCD_txt_File, "w") as outfile:

        #before you loop through every 'i' in plasmid database, the first time the following variable will be true
        first_profit_file = True

        #loop through every line 'i' in plasmid_database
        for i, (header, sequence) in enumerate(plasmid_database.items()):

            #include a proress indicator
            print(f"\r - Processing plasmid {i+1:,} of {len(plasmid_database):,}", end="") #progress indicator

            # create temporary FASTA and profit_output folders for the 'i' in plasmid_database
            temporary_fasta = temp_folder / f"temp_plasmid_{i}.fasta"
            with open(temporary_fasta, "w") as file:
                file.write(f">{header}\n")
                file.write(sequence + "\n")
            profit_output = temp_folder / f"profit_output_{i}.txt"
            
            #run profit
            profit_command = [
                "profit",
                "-sequence", str(temporary_fasta),
                "-infile", str(xer_matrix_file),
                "-outfile", str(profit_output)]
            result = subprocess.run(
                profit_command,
                capture_output=True,
                text=True)

            #if profit was not run on this successfully and a file created, notify me
            if result.returncode != 0:
                print(f"\nPROFIT failed on: {header} {result.stderr}")
            if not profit_output.exists():
                print(f"\nNo PROFIT output created for {header}")
                continue

            
            #if profit was run successfully, annotate it into the output file
            with open(profit_output) as file:
                for line in file:
                    if first_profit_file: #will only be true the first time looping through every 'i' in plasmid database
                        outfile.write(line)
                    else:
                        if line.startswith("#"): #skip the annotations after the first 'i' in the plasmid database
                            continue
                        outfile.write(line)
                        

            #ensure after the first 'i' in the plasmid database that firs_profit_file becomes false:
            first_profit_file = False

            #delete temporary files
            temporary_fasta.unlink()
            profit_output.unlink()


#delete the temporary folder and any files inside it, now that all are finished
for file in temp_folder.iterdir():
    file.unlink()
temp_folder.rmdir()
       
print(f"\n\n\n\n ** Finished! **\n\n")

