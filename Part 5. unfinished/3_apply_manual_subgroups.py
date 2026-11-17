import subprocess
from pathlib import Path







#1 establish number ranges to later parse clinker_sequenced_GBK_files.txt
    # The logic here iwll chnage for the 2 python scripts
    # Error messages TBA
    # These are 1-indexed (not 0-indexed like python)
subgroup_number_range = [
    (1, 12),
    (11, 22),
    (21, 32),
    (31, 42),
    (41, 52),
    (51, 62),
    (61, 72),
    (71, 82),
    (81, 92),
    (91, 102),
    (101, 112),
    (111, 122),
    (121, 132),
    (131, 142),
    (141, 152),
    (151, 162),
    (161, 172),
    (171, 182),
    (181, 192),
    (191, 202),
    (201, 212),
    (211, 222),
    (221, 232),
    (231, 242),
    (241, 252),
    (251, 262),
    (261, 272),
    (271, 282),
    (281, 292),
    (291, 302),
    (301, 312),
    (311, 322),
    (321, 332),
    (331, 342),
    (341, 352),
    (351, 362),
    (361, 372),
    (371, 382),
    (381, 392),
    (391, 402),
    (401, 411)
]





#defining the folder paths for this script
script_folder = Path(__file__).resolve().parent
OUTPUT_DIR = script_folder / "outputfiles"
gbk_folder = OUTPUT_DIR / "sliced_gbk"
subgroup_folder = OUTPUT_DIR / "sliced_gbk_subgroups"
subgroup_folder.mkdir(parents=True, exist_ok=True)
    # clear subgroup folder every go
for old_file in subgroup_folder.glob("*Clinker_Sequenced_GBK_Files.txt"):
    old_file.unlink()
cluster_output_folder = OUTPUT_DIR / "subgroup_clusters"
cluster_output_folder.mkdir(parents=True, exist_ok=True)



# Validate the manually entered subgroup ranges
for start, end in subgroup_number_range:
    if not isinstance(start, int) or not isinstance(end, int):
        raise ValueError("error, incluced non-integer Subgroup range values")
    if start < 1 or end < 1:
        raise ValueError("error. Subgroup range values must be 1 or greater")
    if end < start:
        raise ValueError("error. Subgroup range end cannot be smaller than its start")




        

print("Subgroup number ranges have been established \n")
print(subgroup_number_range)
print("\n")


print("Finding list of gbks in clinker sequence to reorganize \n")

#2 find the clinker sequenced file
clinker_sequence_file = OUTPUT_DIR / 'Clinker_Sequenced_GBK_Files.txt'
if not clinker_sequence_file:
    raise FileNotFoundError(
        f"Clinker sequence list not found: {clinker_sequence_file}"
    )
all_lines = clinker_sequence_file.read_text().splitlines()
if not all_lines:
    raise ValueError("Clinker sequence list is empty")


print(clinker_sequence_file)
all_lines = clinker_sequence_file.read_text().splitlines()


#3 apply the subgroup file number ranges and write the output file, in a loop
for start, end in subgroup_number_range:
    subgroup_lines = all_lines[start -1 : end]

    #trouble shooting and error messages for applying the number ranges to the the list
    if end > len(all_lines):
        raise ValueError(
            f"Range {start}to{end} exceeds the "
            f"{len(all_lines)} entries in the sequence list"
        )
    if not subgroup_lines:
        raise ValueError(f"Range {start}to{end} produced an empty subgroup")

    #create new files with the papropriate ranges
    subgroup_file_content = "\n".join(subgroup_lines)
    new_file_name = f"{start}to{end}_Clinker_Sequenced_GBK_Files.txt"
    new_file_path = subgroup_folder / new_file_name
    new_file_path.write_text(subgroup_file_content)

    print(f"Created {new_file_name}")

  
    # there will probably need to be error messages here for any thing I couldn't put an error message for in #1


#4 find the subgroup text files
subgroup_files = list(subgroup_folder.glob('*Clinker_Sequenced_GBK_Files.txt'))
if not subgroup_files:
    raise FileNotFoundError("No subgroup files found\n")


#5 Run a loop on each text file to do clinker of each subgroup
for subgroup_file in subgroup_files:
       #5a. Read one subgroup file
    gbk_names = subgroup_file.read_text().splitlines()
    print("\n")
    print("GBK file names found: \n")
    print(gbk_names)
    print("\n")



       #5b. Locate the GenBank files

    print("Loading corresponding GBK files, if found... \n")
    gbk_files = []

    for name in gbk_names:
        gbk_path = gbk_folder / name
        gbk_files.append(gbk_path)

        print(gbk_files)


    print("\n")
    print("Running Clinker...\n")


        #5c. Run clinker for this group, ensuring the-p has the proper name
    group_name = subgroup_file.name.replace(
    "_Clinker_Sequenced_GBK_Files.txt", ""
    )
    html_path = cluster_output_folder / f"{group_name}_clusters.html"
    subprocess.run([
        'clinker',
        *gbk_files,
        '--force',
        '-p', str(html_path),
    ], check=True)

print("\n")
print("finished. subgroup html files should be found in output folder named subgroup_clusters")
print("\n")
