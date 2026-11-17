#this was at the bottom of pyxerscript.py in Pauls's script (corresponding to 01_process_xercd_screen_hits))
print("Running Clinker\n")
gbk_files = sorted(glob.glob(os.path.join(OUTPUT_DIR, 'sliced_gbk', '*.gbk')))
if not gbk_files:
    raise FileNotFoundError("No sliced GenBank files found for Clinker")
subprocess.run([
    'clinker',
    *gbk_files,
    '--force',
    '-p', os.path.join(OUTPUT_DIR, 'clusters.html'),
    '-o', os.path.join(OUTPUT_DIR, 'output_summary.txt'),
    '-mo', os.path.join(OUTPUT_DIR, 'matrix.csv'),
], check=True)