




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


# file paths 
P1_emboss_folder = Path(__file__).resolve().parent                      #the folder for all of P1_emboss_folder
input_folder = P1_emboss_folder / "input_folder"                        #input folder
output_folder = P1_emboss_folder / "output_folder"    


import subprocess

ArgR_rev_sequences_file = # add file path for the plain text file with the ArgR reverse sequences
ArgR_prophecy_output = #add output file path for "argr-freq-matrix-rev.txt"

Emboss_verify_result = subprocess.run(
    ["prophecy", "-help"],
    capture_output=True, text=True)
if Emboss_verify_result.returncode != 0:
    raise ValueError(" - EMBOSS PROPHECY could not be executed. Double check the installation.")
print("\n - EMBOSS PROPHECY is installed and can be executed")


prophecy_cmd = ["prophecy", 
                "-sequence", str(ArgR_rev_sequences_file),  
                "-outfile", str(ArgR_prophecy_output), 
                "-threshold", "60",
                "-auto"]

prophecy_result = subprocess.run(prophecy_cmd, capture_output = True, text = True)

if prophecy_result.returncode !=0:
    print(f"Error executing prophecy. stderr: {prophecy_result.stderr}")


print("\n\n   finished\n\n")