#this was at the bottom of pyxerscript.py in Pauls's script (corresponding to 01_process_xercd_screen_hits))
import clinker
import glob
import os
import subprocess
from pathlib import Path
print("Running Clinker\n")


script_folder = Path(__file__).resolve().parent
GBK_folder = script_folder / "3Kb GBKs"
gbk_files = sorted(glob.glob(os.path.join(GBK_folder, '*.gbk')))
if not gbk_files:
    raise FileNotFoundError("No sliced GenBank files found for Clinker")
subprocess.run([
    'clinker',
    *gbk_files,
    '--force',
    '-p', os.path.join(script_folder, 'clusters.html'),
    '-o', os.path.join(script_folder, 'output_summary.txt'),
    '-mo', os.path.join(script_folder, 'matrix.csv'),
], check=True)