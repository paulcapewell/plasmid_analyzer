print("\n\n\n** Running Script: 102_emboss_matrix_scans.py **\n\n")


#------------------------------------------------------------------
# 1. imports
#------------------------------------------------------------------

from pathlib import Path
import gzip
import pickle
import subprocess
import hashlib
from datetime import datetime

#------------------------------------------------------------------
# 2. file paths
#------------------------------------------------------------------

print("*locating necessary files*\n")



# 2.1 core folders
P1_emboss_folder = Path(__file__).resolve().parent                      #the folder for all of P1_emboss_folder
input_folder = P1_emboss_folder / "input_folder"                        #input folder
output_folder = P1_emboss_folder / "output_folder"                      #output folder



# 2.2 output folders
output_subfolder = output_folder / "other_emboss_output_files"
output_subfolder.mkdir(parents=True, exist_ok=True)

xerCD_hits_folder = output_subfolder / "initial_xerCD_hits"
xerCD_hits_folder.mkdir(parents=True, exist_ok=True)

New_XerCD_folder = output_folder / "final_xerCD_hits"
New_XerCD_folder.mkdir(parents=True, exist_ok=True)

ArgR_hits_folder = output_subfolder / "argR_hits"
ArgR_hits_folder.mkdir(parents=True, exist_ok=True)

XerCD_without_ArgR_folder = output_subfolder / "failed_xerCD_hits_without_argR"
XerCD_without_ArgR_folder.mkdir(parents=True, exist_ok=True)

temp_folder = output_folder / "temporary files used in script"
temp_folder.mkdir(parents=True, exist_ok=True)

#debug_218bp_folder = output_folder / "debug_extracted_218bp_sequences"
#debug_218bp_folder.mkdir(parents=True, exist_ok=True)




# 2.3 cached plasmid database... also go ahead and get the checksum
cached_database_folder = output_folder / "database_cache"          #cached plasmid database folder
DATABASE_CACHE = cached_database_folder / 'sequence_database_cache.bin.gz'   #cached plasmid database file
if not DATABASE_CACHE.is_file():
    print("Cached database file not found. cannot run this script")
    raise FileNotFoundError(f"Cached database not found:\n{DATABASE_CACHE}")
with gzip.open(DATABASE_CACHE, "rb") as cache_file:
    plasmid_database = pickle.load(cache_file)
print(f" - cached database found and loaded successfully. Number of plasmids:\n       {len(plasmid_database):,}")

cache_sha256 = hashlib.sha256()
with open(DATABASE_CACHE, "rb") as cache_file:
    for chunk in iter( lambda: cache_file.read(1024 * 1024), b""):
        cache_sha256.update(chunk)
cached_database_checksum = (cache_sha256.hexdigest())


# 2.4 xer matrices

xer_matrix_folder = input_folder / "Xer_frequency_matrices"                 #matrices folder
xer_matrix_files = sorted(xer_matrix_folder.glob("*.txt"))                      #frequency matrices
print("\n - List of matrix files found:")
for xer_matrix_file in xer_matrix_files:
    print(f"       {xer_matrix_file.name}")
if len(xer_matrix_files) != 6:
    raise ValueError(f"Expected exactly 6 matrix files, found {len(xer_matrix_files)}")

#2.5 argR matrices
argR_matrix_folder = input_folder / "ArgR_frequency_matrices"     
argR_matrix_files = sorted(argR_matrix_folder.glob("*.txt"))           
argR_matrix_FOR = argR_matrix_folder / "argr-freq-matrix-for.txt"
if not argR_matrix_FOR.is_file():
    raise ValueError(" - ArgR Forward matrix not found. cannot run script")
argR_matrix_REV = argR_matrix_folder / "argr-freq-matrix-rev.txt"
if not argR_matrix_REV.is_file():
    raise ValueError(" - ArgR reverse matrix not found. cannot run script")
print("\n - Both ArgR frequency matrices located: FOR and REV")

#------------------------------------------------------------------
# 3. append input file information to project report
#------------------------------------------------------------------

receipt_file = (P1_emboss_folder.parent /"infiles_used_report.txt")
receipt_text = f"""
------------------------------
Matrix Files Accessed
------------------------------

Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

Cached Database Used
    File:
        {DATABASE_CACHE.name}
    SHA256:
        {cached_database_checksum}
    Dictionary Entries:
        {len(plasmid_database):,}

Xer Matrices
"""

for matrix_file in xer_matrix_files:
    sha256 = hashlib.sha256()
    with open(matrix_file, "rb") as infile:
        for chunk in iter(lambda: infile.read(1024 * 1024), b""):
            sha256.update(chunk)
    receipt_text += (f"\n    {matrix_file.name}\n        SHA256 = {sha256.hexdigest()}")



receipt_text += "\n\nArgR Matrices\n"
for matrix_file in argR_matrix_files:
    sha256 = hashlib.sha256()
    with open(matrix_file, "rb") as infile:
        for chunk in iter(lambda: infile.read(1024 * 1024), b""):
            sha256.update(chunk)
    receipt_text += (f"\n    {matrix_file.name}\n        SHA256 = {sha256.hexdigest()}")

with open(receipt_file, "a") as outfile:
    outfile.write(receipt_text)
    outfile.write("\n\n")

#------------------------------------------------------------------
# 4. define helper functions to use later for profit
#------------------------------------------------------------------

def run_profit(sequence_file, matrix_file, output_file):
    profit_command = ["profit",
                      "-sequence", str(sequence_file),
                      "-infile", str(matrix_file),
                      "-outfile", str(output_file)]
    result = subprocess.run(
        profit_command,
        capture_output=True,
        text=True)
    if result.returncode != 0:
        raise RuntimeError(f"PROFIT failed:\n{result.stderr}")
    if not output_file.is_file():
        raise RuntimeError(f"PROFIT produced no output:\n{output_file}")


#------------------------------------------------------------------
# 5. verify profit is working and write Sha256 report
#------------------------------------------------------------------
Emboss_verify_result = subprocess.run(
    ["profit", "-help"],
    capture_output=True, text=True)
if Emboss_verify_result.returncode != 0:
    raise ValueError(" - EMBOSS PROFIT could not be executed. Double check the installation.")
print("\n - EMBOSS PROFIT is installed and can be executed")



#------------------------------------------------------------------
# 6. xer scan
#------------------------------------------------------------------
print("\n - running profit scan with all 6 xerCD matrices. might take 1hr+")
#loop through all 6 xer matrix files
for xer_matrix_file in xer_matrix_files:
    print(f"\n\n\n*Running xerCD PROFIT scan with {xer_matrix_file.name}*")
    
    #for each matrix, set up a unique output file based on the unique name of the matrix file
    xer_matrix_file_name_parts = xer_matrix_file.stem.split("_")
    xer_matrix_file_length = xer_matrix_file_name_parts[2]
    modified_matrix_length = int(xer_matrix_file_length.replace("bp", "")) + 22
    xer_matrix_file_orientation = xer_matrix_file_name_parts[3]
    output_name = f"initial_xerCD_hits_{modified_matrix_length}_{xer_matrix_file_orientation}.txt"
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
            with open(temporary_fasta, "w") as fasta_file:
                fasta_file.write(f">{header}\n")
                fasta_file.write(sequence + "\n")
            profit_output = temp_folder / f"profit_output_{i}.txt"
            
            #run profit
            run_profit(
                temporary_fasta,
                xer_matrix_file,
                profit_output)
            
            #if profit was run successfully, annotate it into the output file
            with open(profit_output) as profit_file:
                for line in profit_file:
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


       
print(f"\n\n\n\n ** Finished Xer Scan. now performing ArgR scan **\n\n")



#------------------------------------------------------------------
# 7. argR scan 218bp upstream of xer hits
#------------------------------------------------------------------


HITS_XerCD_FOR_files = sorted(xerCD_hits_folder.glob("*FOR.txt")) 
HITS_XerCD_REV_files = sorted(xerCD_hits_folder.glob("*REV.txt"))
HITS_XerCD_files = (HITS_XerCD_FOR_files + HITS_XerCD_REV_files)
HITS_XerCD_file_names = [HITS_XerCD_file.name for HITS_XerCD_file in HITS_XerCD_files]
if len(HITS_XerCD_files) != 6:
    raise ValueError(f"Expected exactly 3 `HITS XerCD FOR.txt` files, found {len(HITS_XerCD_files)}. fix and retry")
print(f" - all 6 initial xerCD hits files found:\n\t\t{HITS_XerCD_file_names}")




# ----------------------------------------------------------------------
# 7.1 Set up looping logic and output files for every XerCD hit files
# ----------------------------------------------------------------------


HIT_WINDOW_LENGTH = 218
REFERENCE_MOTIF_LENGTH = 28
window_flank = HIT_WINDOW_LENGTH - REFERENCE_MOTIF_LENGTH

print("\n\n")
# 7.1.2 loop through each XerCD hit file
for xer_hit_file in HITS_XerCD_files:
    
    print(f"\n*Running ArgR matrix search on {xer_hit_file.name}*")

    # 7.1.3 for each, establish the naming logic for output files
    xer_hit_file_name_parts = xer_hit_file.stem.split("_")
    xer_hit_file_motif_length =  xer_hit_file_name_parts[3]
    xer_hit_file_motif_orientation =  xer_hit_file_name_parts[4]

    # 7.1.4 select the appropriate frequency matrix to use
    if xer_hit_file_motif_orientation == "REV":
        corresponding_argR_matrix = argR_matrix_REV
    elif xer_hit_file_motif_orientation == "FOR":
        corresponding_argR_matrix = argR_matrix_FOR
    else:
        raise ValueError(f"Unexepcted text instead of FOR/REV for orientation: {xer_hit_file_motif_orientation}")
    is_reverse = (xer_hit_file_motif_orientation == "REV")

    # 7.1.5 name output files
    New_XerCD_file = New_XerCD_folder / f"xerCD_hits_{xer_hit_file_motif_length}_{xer_hit_file_motif_orientation}.txt"
    ArgR_hits_file = ArgR_hits_folder / f"HITS ArgR {xer_hit_file_motif_orientation}_{xer_hit_file_motif_length}xer_.txt"
    XerCD_without_ArgR_file = XerCD_without_ArgR_folder / f"XerCD {xer_hit_file_motif_length} {xer_hit_file_motif_orientation}_noArgR_.txt"
    # 7.5 open all output files for this XerCD input file: 
    with open(ArgR_hits_file, "w") as argR_outfile, \
        open(New_XerCD_file, "w") as new_xer_outfile, \
        open(XerCD_without_ArgR_file, "w") as no_argr_outfile:

        # 7.1.6 establish logic for identifying first_profit_file, which is important later for formatting output file content
        first_profit_file = True


        # ----------------------------------------------------------------------
        # 7.2 go through the XerCD Hit File line by line
        # ----------------------------------------------------------------------

        # 7.2.1 process the XerCD hit file, line by line and give a number to each line
        with open(xer_hit_file) as xer_file:
            print(f" - Analyzing {xer_hit_file.name} line-by-line")
            for i, xer_hit_line in enumerate(xer_file):
                print(f"\r - Lines processed so far: {i+1}", end="")
        

                # 7.2.2 copy annotation lines into the relevant output files
                if xer_hit_line.startswith("#"):
                     new_xer_outfile.write(xer_hit_line)
                     no_argr_outfile.write(xer_hit_line)
                     continue

            
                # ----------------------------------------------------------------------
                # 7.3. for any line in the XerCD Hit file that's not an annotation, do the ArgR profit search
                # ----------------------------------------------------------------------

                # 7.3.1 find the corresponding plasmid from the database
                xer_hit_line_data = xer_hit_line.split()
                xer_hit_line_accession = xer_hit_line_data[0]
                xer_hit_line_coordinate = xer_hit_line_data[1]
                xer_hit_line_percentage = xer_hit_line_data[3]

                header = xer_hit_line_accession
                sequence = plasmid_database.get(header)
                xer_hit_position = int(xer_hit_line_coordinate)


                # 7.3.2 error message if sequence was not found in the plasmid database
                if sequence is None:
                    print(f" - could not find accession for {xer_hit_line}. This line is will not appear in output files")
                    continue

                # 7.3.3 orhwise, if it was found in the database, identify the 250bp region you want to serach it of
                
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

                    #7.3.4 runt he ArgR Emboss search on the 218bp region
                    temporary_fasta = temp_folder / f"temp_plasmid_{i}.fasta" #this block defines temporary output folders
                    with open(temporary_fasta, "w") as temp_file:
                        temp_file.write(f">{header}\n")
                        temp_file.write(sequence_218bp + "\n")
                    profit_output = temp_folder / f"profit_output_{i}.txt"
                    run_profit(
                        temporary_fasta,
                        corresponding_argR_matrix,
                        profit_output)
                
                # ----------------------------------------------------------------------
                # 7.4. write the output files
                # ----------------------------------------------------------------------
                
                #7.4.1 make an output file with every ArgR hit
                with open(profit_output) as profit_file:
                    for line in profit_file:
                        if first_profit_file: #will only be true the first time looping through every 'i' in plasmid database
                            argR_outfile.write(line)
                        else:
                            if line.startswith("#"): #skip the annotations after the first one
                                continue
                            argR_outfile.write(line)
                                
                #7.4.2 identify if the line from the xer_hit_file was a hit in argR screen at all, based on whether the word "Percenntage" is in the output file
                line_is_argR_hit = False
                with open(profit_output) as profit_file:
                    if "Percentage:" in profit_file.read():
                        line_is_argR_hit = True

                #7.4.3 write an output fille for any lines that were not ArgR hits
                if line_is_argR_hit == False:
                    no_argr_outfile.write(xer_hit_line)

                #7.4.4 write an output file for the lines that were ArgR hits 
                if line_is_argR_hit == True:
                    new_xer_outfile.write(xer_hit_line)
                                
                first_profit_file = False


# ----------------------------------------------------------------------
# 7.5. finish loop
# ----------------------------------------------------------------------
                temporary_fasta.unlink()
                profit_output.unlink()

        print(f"\r\n - finished\n\n")





# ----------------------------------------------------------------------
# 8. clean up and finish
# ----------------------------------------------------------------------

for temp_file in temp_folder.iterdir():
    temp_file.unlink()
temp_folder.rmdir()
       
print(f"\n ** Finished runnign script 102! **\n\n")


