# Plasmid Xer Analyzer

This script is a python port of Xerscript_v4.m that extracts sequence windows around Xer-site hits, fetches annotated GenBank records for the matching plasmids, writes annotated 3 kb regions, and runs Clinker to compare those regions.

Requirements

- Python 3.9 or newer
- Biopython
- Clinker, including its `clinker` command-line executable
- Internet access to fetch GenBank records from NCBI, unless the corresponding records already exist in `outputfiles/gbk`
- A gzipped FASTA database named `bacterial_plasmids_curated.fa.gz` in the project root from [figshare.com/articles/dataset/Complete_plasmids_ARG_carriage_study/19438994?file=34537640](https://figshare.com/articles/dataset/Complete_plasmids_ARG_carriage_study/19438994?file=34537640)

## Install Dependencies

From PowerShell, in the project directory:

```powershell
py -3 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install biopython clinker
```

If PowerShell prevents virtual-environment activation, invoke the environment's Python directly:

```powershell
.\.venv\Scripts\python.exe -m pip install biopython clinker
```

The `clinker` executable installed by pip must be available to the process running the script.

## Input Files

Place the hit-list text files under `inputfiles/`.

```text
NEW XerCD 28 FOR.txt
NEW XerCD 28 REV.txt
NEW XerCD 29 FOR.txt
NEW XerCD 29 REV.txt
NEW XerCD 30 FOR.txt
NEW XerCD 30 REV.txt
```

The second-to-last filename token is interpreted as the motif length (for example, `28`); a filename containing `REV` is treated as reverse strand, and other names are treated as forward strand. Use `.txt` files and retain that naming convention.

Each non-comment line must start with an accession and its sequence position. Additional columns are ignored. Lines beginning with `#` and blank lines are skipped. For example:

```text
NC_011408.1 1738 Percentage: 98
```

The accession must match a FASTA header exactly and must be fetchable through NCBI Entrez for GenBank annotation.

## Run

Run from the project root so the relative input and output paths resolve correctly:

```powershell
python pyxerscipt.py
```

The first run reads `bacterial_plasmids_curated.fa.gz` and creates `sequence_database.bin.gz`, a compressed pickle cache for faster database loading on later runs. If the FASTA database changes, delete `sequence_database.bin.gz` before rerunning so it is rebuilt from the new FASTA. NCBI asks users to provide a contact email. Before running, change the `Entrez.email` value near the top of `pyxerscipt.py` to your own email address. NCBI downloads may fail due to network problems, invalid accessions, or rate limits.

## What It Does

1. Loads the FASTA sequences into memory, using the pickle cache if available.
2. For each hit, extracts a 218-base window around the reported position. Circular sequence ends are handled by repeating the sequence; reverse hits are reverse-complemented.
3. Fetches the corresponding GenBank record from NCBI if it is not already cached in `outputfiles/gbk/`.
4. Repeats each circular plasmid and its annotated features as needed, then clips a 3,000-base region centered on the computed Xer position. It adds placeholder `gene` and `CDS` annotations for the Xer site.
5. Runs Clinker on the sliced GenBank records. Clinker is a gene-cluster comparison tool: it compares annotated genes in the 3 kb regions, groups similar genes, and shows how their order and sequence similarity relate across the input records. Here it compares regions around the reported Xer sites; it does not search for Xer motifs or establish that a gene is functionally related to Xer.

## Outputs

| Path                                    | Contents                                                                                            |
| --------------------------------------- | --------------------------------------------------------------------------------------------------- |
| `outputfiles/output_<input filename>` | FASTA windows, one output file per input hit list. Headers contain accession and reported position. |
| `outputfiles/gbk/<accession>.gbk`     | Downloaded GenBank record cache. Existing records are reused.                                       |
| `outputfiles/sliced_gbk/<hit ID>.gbk` | Annotated 3 kb GenBank slice for each hit.                                                          |
| `outputfiles/clusters.html`           | Clinker cluster visualization.                                                                      |
| `outputfiles/output_summary.txt`      | Clinker text summary.                                                                               |
| `outputfiles/matrix.csv`              | Clinker similarity matrix.                                                                          |
| `sequence_database.bin.gz`            | Local compressed pickle cache of the FASTA database.                                                |

Clinker compares the genes annotated in the input GenBank slices, clusters similar genes, and reports how gene order and sequence similarity compare between regions. `clusters.html` is the interactive visualization, `output_summary.txt` contains the clustering results, and `matrix.csv` contains pairwise cluster-similarity values. These results compare the sliced regions; Clinker does not find Xer motifs or determine whether a gene is functionally related to Xer.

The script does not clean `outputfiles/sliced_gbk/` before a run. Old GenBank slices left there can be included in the next Clinker run; remove that directory's old contents before rerunning a different analysis.

## Editable Settings

The configuration block in `pyxerscipt.py` contains the main settings:

- `ENTREZ_EMAIL`: contact email sent with NCBI Entrez requests. Set this to your own address.
- `INPUT_DIR` and `OUTPUT_DIR`: folders used for hit-list inputs and generated results. Paths are relative to the working directory used to run the script.
- `FASTA_DATABASE` and `DATABASE_CACHE`: source gzipped FASTA filename and its compressed pickle cache. Delete the cache when changing the source FASTA so it is rebuilt.
- `HIT_WINDOW_LENGTH`: total bases extracted around each hit; currently 218 bases for forward and reverse hits.
- `REFERENCE_MOTIF_LENGTH`: reference motif size used in the forward-window calculation. The remaining flank is derived as `HIT_WINDOW_LENGTH - REFERENCE_MOTIF_LENGTH` (currently 190 bases).
- `SLICE_FLANK`: bases retained on each side of the computed hit center. The resulting slice is twice this value (currently 3,000 bases total).
- `FORWARD_CENTER_OFFSET` and `REVERSE_CENTER_OFFSET`: position corrections used to place the hit center. The current formulas are `position + motif_length - 6` for forward hits and `position` for reverse hits. These depend on the hit files' coordinate convention; verify that convention before changing them.
- `XER_FEATURE_CENTER_OFFSET`: shifts the start of the placeholder Xer annotation relative to the center of the sliced region. Zero starts it at coordinate 1500 with the current `SLICE_FLANK`. Changing this will make the metadata  inaccurate.
- `XER_FEATURE_LENGTH`: span of each placeholder `gene` and `CDS` feature, currently six bases. These are annotations added for comparison, not a motif detection result.

The Clinker executable arguments and its three output filenames are configured at the end of the script.

`outputfiles/gbk/` records are reused if present. If you need to refresh an annotation, remove the corresponding cached `.gbk` file so it is downloaded again.
