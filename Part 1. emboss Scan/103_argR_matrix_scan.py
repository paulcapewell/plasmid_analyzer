# ----------------------------------------------------------------------
# 1. import modules
# ----------------------------------------------------------------------

print("\n\n\n** Running Script: 103_argR_matrix_scan.py **\n\n")
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


# ----------------------------------------------------------------------
# 2. define file paths ad verify input files
# ----------------------------------------------------------------------

print("*locating necessary files and verifying EMBOSS profit is available*")

P1_emboss_folder = Path(__file__).resolve().parent
input_folder = P1_emboss_folder / "input_folder" 

    #2.2 output folders
output_folder = P1_emboss_folder / "output_folder" 
xerCD_hits_folder = output_folder / "xerCD hits before ArgR"
if not xerCD_hits_folder.is_dir():
    raise FileNotFoundError(f"Folder not found:\n{xerCD_hits_folder}")
New_XerCD_folder = output_folder / "xerCD hits after ArgR"
New_XerCD_folder.mkdir(parents=True, exist_ok=True)
ArgR_hits_folder = output_folder / "ArgR Hits"
ArgR_hits_folder.mkdir(parents=True, exist_ok=True)
XerCD_without_ArgR_folder = output_folder / "XerCD hits failed ArgR scan"
XerCD_without_ArgR_folder.mkdir(parents=True, exist_ok=True)
temp_folder = output_folder / "temporary files used in script"
temp_folder.mkdir(parents=True, exist_ok=True)
file_id_folder = output_folder / "file_identity_and_integrity"
file_id_folder.mkdir(parents=True, exist_ok=True)


# [DEBUG] Persistent folder to review extracted 218bp sequences (easy to delete/disable later)
debug_218bp_folder = output_folder / "debug_extracted_218bp_sequences"
debug_218bp_folder.mkdir(parents=True, exist_ok=True)


    #2.3 cached plasmid database
cached_database_folder = output_folder / "cached_plasmid_database"          #cached plasmid database folder
DATABASE_CACHE = cached_database_folder / 'sequence_database_cache.bin.gz'   #cached plasmid database file
if not DATABASE_CACHE.exists():
    print("Cached database file not found. cannot run this script")
    raise SystemExit
with gzip.open(DATABASE_CACHE, "rb") as file:
    plasmid_database = pickle.load(file)
#print(f" - cached database found and loaded successfully. Number of plasmids:\n       {len(plasmid_database):,}")

    #2.4 argR matrices
argR_matrix_folder = input_folder / "ArgR_frequency_matrices"     
argR_matrix_files = sorted(argR_matrix_folder.glob("*.txt"))           
argR_matrix_FOR = argR_matrix_folder / "argr-freq-matrix-for.txt"
if not argR_matrix_FOR.exists():
    ValueError(" - ArgR Forward matrix not found. cannot run script")
argR_matrix_REV = argR_matrix_folder / "argr-freq-matrix-rev.txt"
if not argR_matrix_REV.exists():
    ValueError(" - ArgR reverse matrix not found. cannot run script")
print(" - Both ArgR frequency matrices located: FOR and REV\n")

    #2.5 XerCD hits
HITS_XerCD_FOR_files = sorted(xerCD_hits_folder.glob("*FOR.txt")) 
HITS_XerCD_REV_files = sorted(xerCD_hits_folder.glob("*REV.txt"))
HITS_XerCD_files = (HITS_XerCD_FOR_files + HITS_XerCD_REV_files)
HITS_XerCD_file_names = [file.name for file in HITS_XerCD_files]

if len(HITS_XerCD_files) != 6:
    raise ValueError(f"Expected exactly 3 `HITS XerCD FOR.txt` files, found {len(HITS_XerCD_FOR_files)}. fix and retry")
print(f" - all 6 `HITS XerCD FOR.txt` files found: {HITS_XerCD_file_names}\n")

# ----------------------------------------------------------------------
# 3. check that EMBOSS is working
# ----------------------------------------------------------------------

Emboss_verify_result = subprocess.run(
    ["profit", "-help"],
    capture_output=True, text=True)
if Emboss_verify_result.returncode != 0:
    raise ValueError(" - EMBOSS PROFIT could not be executed. Double check the installation.")
print(" - EMBOSS PROFIT is installed and can be executed\n")


# ----------------------------------------------------------------------
# 4. create SHA256 receipt for all matrix files----------------------
# ------------------------------------------------
sha256_file = file_id_folder / "argR matrices sha256.txt"
with open(sha256_file, "w") as outfile:
    outfile.write("SHA256 Checksums for Matrix Files\n")
    outfile.write("--------------------------------\n\n")
    for argR_matrix_file in argR_matrix_files:
        sha256 = hashlib.sha256()
        with open(argR_matrix_file, "rb") as file:
            for chunk in iter(lambda: file.read(1024 * 1024),b""):
                sha256.update(chunk)
        outfile.write(f"{argR_matrix_file.name}\n")
        outfile.write(f"SHA256: {sha256.hexdigest()}\n\n")
print(f"\n\n* Matrix SHA256 report written to: *\n - `{sha256_file.name}` file in `{file_id_folder.name}` folder")



# ----------------------------------------------------------------------
# 5. Set up looping logic and output files for every XerCD hit files
# ----------------------------------------------------------------------


HIT_WINDOW_LENGTH = 218
REFERENCE_MOTIF_LENGTH = 28
window_flank = HIT_WINDOW_LENGTH - REFERENCE_MOTIF_LENGTH


# 5.1 loop through each XerCD hit file
for xer_hit_file in HITS_XerCD_files:
    
    print(f"\n\n\n*Running ArgR matrix search on {xer_hit_file.name}*")

    # 5.2 for each, establish the naming logic for output files
    xer_hit_file_name_parts = xer_hit_file.stem.split(" ")
    xer_hit_file_motif_length =  xer_hit_file_name_parts[2]
    xer_hit_file_motif_orientation =  xer_hit_file_name_parts[3]

    # 5.3 select the appropriate frequency matrix to use
    if xer_hit_file_motif_orientation == "REV":
        corresponding_argR_matrix = argR_matrix_REV
    elif xer_hit_file_motif_orientation == "FOR":
        corresponding_argR_matrix = argR_matrix_FOR
    else:
        raise ValueError(f"Unexepcted text instead of FOR/REV for orientation: {xer_hit_file_motif_orientation}")
    is_reverse = (xer_hit_file_motif_orientation == "REV")

    # 5.4 name output files
    New_XerCD_file = New_XerCD_folder / f"NEW XerCD {xer_hit_file_motif_length} {xer_hit_file_motif_orientation}.txt"
    ArgR_hits_file = ArgR_hits_folder / f"HITS ArgR {xer_hit_file_motif_orientation}_{xer_hit_file_motif_length}xer_.txt"
    XerCD_without_ArgR_file = XerCD_without_ArgR_folder / f"XerCD {xer_hit_file_motif_length} {xer_hit_file_motif_orientation}_noArgR_.txt"
    # 5.5 open all output files for this XerCD input file: 
    with open(ArgR_hits_file, "w") as argR_outfile, \
        open(New_XerCD_file, "w") as new_xer_outfile, \
        open(XerCD_without_ArgR_file, "w") as no_argr_outfile:

        # 5.6 establish logic for identifying first_profit_file, which is important later for formatting output file content
        first_profit_file = True


        # ----------------------------------------------------------------------
        # 6. go through the XerCD Hit File line by line
        # ----------------------------------------------------------------------

        # 6.1 process the XerCD hit file, line by line and give a number to each line
        with open(xer_hit_file) as xer_file:
            print(f"\n - Analyzing {xer_hit_file.name} line-by-line")
            for i, xer_hit_line in enumerate(xer_file):
                print(f"\r - Lines processed so far: {i+1}", end="")
        

                # 6.2 copy annotation lines into the relevant output files
                if xer_hit_line.startswith("#"):
                     new_xer_outfile.write(xer_hit_line)
                     no_argr_outfile.write(xer_hit_line)
                     continue

            
                # ----------------------------------------------------------------------
                # 7. for any line in the XerCD Hit file that's not an annotation, do the ArgR profit search
                # ----------------------------------------------------------------------

                # 7.1 find the corresponding plasmid from the database
                xer_hit_line_data = xer_hit_line.split()
                xer_hit_line_accession = xer_hit_line_data[0]
                xer_hit_line_coordinate = xer_hit_line_data[1]
                xer_hit_line_percentage = xer_hit_line_data[3]

                header = xer_hit_line_accession
                sequence = plasmid_database.get(header)
                xer_hit_position = int(xer_hit_line_coordinate)


                # 7.2 error message if sequence was not found in the plasmid database
                if sequence is None:
                    print(f" - could not find accession for {xer_hit_line}. This line is will not appear in output files")
                    continue

                # 7.3 orhwise, if it was found in the database, identify the 250bp region you want to serach it of
                
                else:
                    motif_length = int(xer_hit_file_motif_length)
                    wrapped_sequence = (sequence * 3).upper()
                    if is_reverse:
                        start = xer_hit_position + len(sequence) - 1
                        sequence_218bp = wrapped_sequence[start:(start + HIT_WINDOW_LENGTH)]
                    else:
                        start = xer_hit_position - window_flank + len(sequence) - 1
                        motif_start = start - (REFERENCE_MOTIF_LENGTH - motif_length)
                        sequence_218bp = wrapped_sequence[
                            start - (REFERENCE_MOTIF_LENGTH - motif_length):
                            start + window_flank + motif_length]

                    #7.4 runt he ArgR Emboss search on the 250 region
                    temporary_fasta = temp_folder / f"temp_plasmid_{i}.fasta" #this block defines temporary output folders
                    with open(temporary_fasta, "w") as file:
                        file.write(f">{header}\n")
                        file.write(sequence_218bp + "\n")

                    profit_output = temp_folder / f"profit_output_{i}.txt"
                    profit_command = [ #this block defines how EMBOSS wil run
                        "profit",
                        "-sequence", str(temporary_fasta),
                        "-infile", str(corresponding_argR_matrix),
                        "-outfile", str(profit_output)]
                    result = subprocess.run(profit_command, capture_output=True, text=True) # this line actually runs EMBOSS

                # 7.5 notify about errors
                                #if profit was not run on this successfully and a file created, notify me
                if result.returncode != 0:
                    print(f"\n - PROFIT failed on: {header} {result.stderr}. this line will not be in any output file")
                if not profit_output.exists():
                    print(f"\nNo temporary PROFIT output created for {header}. this line will not be in any output file")
                    continue

                
                # ----------------------------------------------------------------------
                # 8. write the output files
                # ----------------------------------------------------------------------
                
                #8.1 make an output file with every ArgR hit
                with open(profit_output) as file:
                    for line in file:
                        if first_profit_file: #will only be true the first time looping through every 'i' in plasmid database
                            argR_outfile.write(line)
                        else:
                            if line.startswith("#"): #skip the annotations after the first one
                                continue
                            argR_outfile.write(line)
                                
                # 8.2 identify if the line from the xer_hit_file was a hit in argR screen at all, based on whether the word "Percenntage" is in the output file
                line_is_argR_hit = False
                with open(profit_output) as file:
                    if "Percentage:" in file.read():
                        line_is_argR_hit = True

                #8.3 write an output fille for any lines that were not ArgR hits
                if line_is_argR_hit == False:
                    no_argr_outfile.write(xer_hit_line)

                #8.4 write an output file for the lines that were ArgR hits 
                if line_is_argR_hit == True:
                    new_xer_outfile.write(xer_hit_line)
                                
                first_profit_file = False


# ----------------------------------------------------------------------
# 9. finish loop and clean up
# ----------------------------------------------------------------------
                temporary_fasta.unlink()
                profit_output.unlink()

        print(f"\r\n - finished {xer_hit_file.name})\n\n")

for file in temp_folder.iterdir():
    file.unlink()
temp_folder.rmdir()
       
print(f"\n\n\n\n ** Finished! **\n\n")