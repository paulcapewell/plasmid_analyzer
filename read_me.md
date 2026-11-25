# Plasmid Analyzer

This is a four-part Python workflow, built on an initial script from Paul Capewell, that searches a plasmid database for XerCD-recognition sites that are similar in structure to the cer site on ColE1. Hits from this search are analyzed, by extracting DNA around the hit to build phylogenetic trees with iTOL annotaitons and to cluster nearby sequences if they are closely related on either xerC or xerD flanks. THese clusters are compared across the nearby DNA in clinker to highlight possible sites of xer-mediated recombination. In the last step, lots of this data is compiled into one summary spreadsheet to analyze.

The script and output folders are organized into subfolders based on function an dscripts are numbered (1xx, 2xx, 3xx, 4xx) to match these folders. Scripts **must be run in numerical order** because they read files written by earlier ones.

------------------------------------------------------------

## Workflow overview

Before you begin, you must **add the plasmid database FASTA** (ending in .fa.gz) into its appropriate folder: 
     *   `plasmid_analyzer/Part 1. emboss Scan/input_folder/original_plasmid_database`

| Part | Script | What it does | Main output |
|------|--------|--------------|-------------|
| 1 | `101_open_plasmid_database.py` | Reads the plasmid FASTA (`.fa.gz`) and saves it as a fast-loading pickled cache | `database_cache/sequence_database_cache.bin.gz` |
| 1 | `102_emboss_matrix_scans.py` | Scans every plasmid with Xer frequency matrices (using emboss profit), then keeps hits with an ArgR site inside the surrounding 218 bp window | `final_xerCD_hits/*.txt` |
| 1 | `103_extract_sequences.py` | From the final_xerCD_hits files, it extracts the 218 bp Cer-like region, XerC/crossover/XerD sub-sequences, and upstream/downstream flanks for every hit; downloads full GenBank records and cuts 3 kb GBK slices | `compiled_sequences_near_xerCD/` , `individual_GBKs/`  |
| 2 | `201_phylogenetic_trees.py` | Aligns sequence sets with MAFFT, calculates distance matrices, builds Neighbour-Joining trees | `trees/*.nwk` |
| 2 | `202_itol_annotations.py` | Writes iTOL annotation files (known sites, host genus, host species) and defines phylogenetic clusters on both xerC and xerD flanks | `annotations/` (iTOL `.txt` files + cluster `.tsv` files) |
| 3 | `301_run_clinker.py` | Runs Clinker on the 3kb GBK slices based on xerC or xerD flank clusters from script 202 | One `.html` per cluster (XerC and XerD folders) |
| 3 | `302_modify_clinker_htmls.py` | Re-renders each Clinker HTML in headless Chromium, flips reverse-strand plasmids, exports an SVG, and builds summary PDFs | `polished images of 3kb clinkers/`, `XerC/XerD_all_clusters_summary.pdf` |
| 4 | `401_summary_spreadsheet.py` | Combines important earlier results into a spreadsheet with one row per Xer hit | `summary_spreadsheet.csv` |

------------------------------------------------------------

## Requirements

**Python 3.9+** with:

- `biopython`
- `pandas`
- `playwright` (then run `playwright install chromium`) - chromium is used by 302

**Command-line tools** (must be on your `PATH`):

- [EMBOSS](https://emboss.sourceforge.net/) (specifically `profit`), used by 102
- [MAFFT](https://mafft.cbrc.jp/alignment/software/), used by 201
- [Clinker](https://github.com/gamcil/clinker), used by 301

**Other:**

- Internet access for script 103 (downloads GenBank records from NCBI via Entrez). Set your own email in the `ENTREZ_EMAIL` variable near the top of the script for NCBI
- [iTOL](https://itol.embl.de/) (web tool). The annotation files from 202 are uploaded there manually along with the `.nwk` trees; no script runs iTOL.

------------------------------------------------------------


## Folder structure
```
Plasmid Analyzer/                                 <- project folder
├── infiles_used_report.txt                       <- created automatically. Let's you detect changes in input files used across versions (previously an issue)
├── Part 1. emboss Scan/               
│   ├── 101_open_plasmid_database.py
│   ├── 102_emboss_matrix_scans.py
│   ├── 103_extract_sequences.py
│   ├── input_folder/
│   │   ├── original_plasmid_database/
│   │   ├── Xer_frequency_matrices/       
│   │   └── ArgR_frequency_matrices/
│   └── output_folder/                 
├── Part 2. Phylogenetics/             
│   ├── 201_phylogenetic_trees.py                 <- all 3 can be uploaded to iTOL as they are
│   ├── 202_itol_annotations.py                   
│   │   ├── for_1kb_trees/                        <- all .txt files here can be uploaded to iTOL for both 1kb trees
│   │   └── for_collapsed_218bp_tree/             <- all files here can be uploaded to iTOL for the 218bp tree
│   └── output_folder/
├── Part 3. Clinker/
│   ├── 301_run_clinker.py
│   ├── 302_modify_clinker_htmls.py
│   └── output_folder/
│       └── [...]_all_clusters_summary.pdf        <- 2 very helpful pdfs, containing all clinkers
├── Part 4. Sumamary Spreadsheet/
│   ├── 401_summary_spreadsheet.py
│   └── output_folder/
│        └── summary_spreadsheet.csv              <- very helpful output file with lots of data on cer-like sites
├── reade_me.md
└── Reference. original sequences/                <- for reference (not essential to the scipt)
                                                    lists sequences used by Dr Emanuele Conte to generate the XerCD and ArgR freq matrices
```

------------------------------------------------------------
## Running the workflow
------------------------------------------------------------

Run each script from its own folder (or by full path), in this order:
     (will take several hours to run on a laptop)

```
python 101_open_plasmid_database.py
python 102_emboss_matrix_scans.py
python 103_extract_sequences.py
python 201_phylogenetic_trees.py
python 202_itol_annotations.py
python 301_run_clinker.py
python 302_modify_clinker_htmls.py
python 401_summary_spreadsheet.py
```




------------------------------------------------------------
### Part 1. EMBOSS scan
------------------------------------------------------------
**101: open plasmid database.** Set `DATABASE_NAME` at the top if your file has a different name. The script reads the gzipped FASTA, builds a `{header: sequence}` dictionary, pickles it, and checks the cache can be reloaded. It records SHA256 checksums of the original and cached files in `infiles_used_report.txt`. If the cache already exists, the script exits without doing anything, so delete the cache to force a rebuild. Expect 30+ minutes on a large database.

**102: EMBOSS matrix scans.**
1. Scans every plasmid with each of the 6 Xer matrices (28bp/29bp/30bp × FOR/REV strands). Matrix filenames are parsed by splitting on `_`: the third part must be the length (e.g. `6bp`) and the fourth the orientation (`FOR`/`REV`). The motif length written to the hits file is the matrix length + 22 bp (the XerC and XerD binding sites). Results go to `other_emboss_output_files/initial_xerCD_hits/`.
2. For every hit, takes the 218 bp window to the xerD direction and scans it with the matching ArgR matrix (`FOR`/`REV`).
3. Hits **with** an ArgR site are written to `final_xerCD_hits/`; hits without one go to `other_emboss_output_files/failed_xerCD_hits_without_argR/`. All ArgR hits go to `other_emboss_output_files/argR_hits/`.

Plasmids are treated as circular, so windows that cross the origin wrap around. Be prepared for this to possibly take a couple of hours.

**103: extract sequences.** For every hit in `final_xerCD_hits/`, writes (in `compiled_sequences_near_xerCD/`):

- `218+ bp/central_218bp_sequences.txt` (one 218 bp window per hit) and `central_218bp_sequences_collapsed.txt` (identical 218 bp sequences merged into one entry, with member names joined by `|`)
- `218+ bp/XerC_and_up_{250,1000,3000}bp_sequences.txt` and `XerD_and_down_{250,1000,3000}bp_sequences.txt` (flanking sequence sets, oriented consistently regardless of strand)
- `xerC, crossover region, xerD/all_xerC_sites.txt`, `all_crossover_regions.txt`, `all_xerD_sites.txt`

It then downloads the full GenBank record for each accession (`individual_GBKs/full sequence/`, skipped if already present) and cuts a 3,000 bp slice centred on each site (`individual_GBKs/3kb/`). Each slice gets a placeholder feature (locus tag `XM_001`) marking the Xer site so it shows up in Clinker and in script 401.

Hit IDs follow the pattern `{accession}_{position}_{FOR|REV}_{motif length}`. This ID is used in sequence headers, tree node names, and GBK filenames throughout the rest of the workflow.




------------------------------------------------------------
### Part 2. Phylogenetics
------------------------------------------------------------
**201: phylogenetic trees.** Processes three sequence sets: the collapsed 218 bp sequences, the XerC-and-up 1000 bp sequences, and the XerD-and-down 1000 bp sequences. For each: MAFFT alignment (`--auto`), pairwise "identity" distance matrix (Biopython), and a Neighbour-Joining tree saved as Newick. This can take a long time (20-60 minutes per tree).


**202: iTOL annotations.** Reads the trees and distance matrices, looks up host organism/genus from the full-length GBK files, and writes to `annotations/`:

Clusters are defined as connected components of plasmids whose pairwise distance is ≤ `MATRIX_CLUSTER_MAX_DISTANCE` (0.01), keeping only clusters of at least `MIN_CLUSTER_SIZE` (5) members. Both settings are near the top of the script. Clusters are defined in `for_1kb_trees/*_summary.tsv` file for later use

The accession lists for known functional sites (e.g. ColE1 *cer*, `NC_001371.1`) are also defined there.

- `for_collapsed_218bp_tree/`: 
     *    known functional sites (blue ColE1, red others) helps orient the phylogenetic tree
     *    host genus symbols (different coloured circles around the perimeter) points out when the same 218bp sequence is hosted by multiple bacteria genus
     *    host species popup, showign the species whenever your mouse goes over the region inside the tree
- `for_1kb_trees/`: 
     *    the same three annotations (built from the XerC 1 kb tree), but here aren no duplicates on these trees so how genus are less helpful
     *    plus for each 1 kb tree a `*_clusters_annotation.txt` (iTOL colour strip for related branches)

Upload each `.nwk` tree to iTOL, then upload in the matching annotation files.





------------------------------------------------------------
### Part 3. Clinker
------------------------------------------------------------

**301: run Clinker.** Reads the two cluster `.tsv` files from 202 and, for each cluster, runs `clinker` on the corresponding 3 kb GBK slices, writing one interactive HTML per cluster into `HTMLs from XerC 1kb clusters/` and `HTMLs from XerD 1kb clusters/` depending on which 1kb flank the clusters are baed on

**302: modify Clinker HTMLs.** Opens each HTML in headless Chromium (via Playwright) and patches a temporary copy so that:

- plasmids with `_REV_` in their name are flipped, so all sites read in the same direction
- plot size, font sizes and gene shapes are set to fixed values (`clinker_settings` near the top of the script)

It then exports a clean SVG for each cluster to `polished images of 3kb clinkers/` and combines them into `XerC_all_clusters_summary.pdf` and `XerD_all_clusters_summary.pdf`. The original Clinker HTMLs are never modified, and the temporary copies are deleted at the end. The script patches specific strings in Clinker's HTML output, so if a Clinker update changes that HTML, affected files are skipped and listed in the "failed" summary at the end.






------------------------------------------------------------
### Part 4. Summary spreadsheet
------------------------------------------------------------
**401: summary spreadsheet.** Writes `summary_spreadsheet.csv` with one row per 3 kb GBK slice (i.e. per Xer hit). Columns:
- Plasmid name
- plasmid accession
- plasmid length (bp)
- host organism
- host taxonomy
- GC content
- Full-length GBK
- Motif identification percentage (from xerCD emboss profit scan)
- xer motif coordinates
- xer motif length
- xer motif strand
- 218bp region sequence
- XerC binding sequence
- Crossover region sequence
- XerD binding sequence
- Number of CDS within 1.5kb (and nearest upsream/downstream genes)
- is it a known functional site? (cer-site / other known functional site / matches a known 218bp sequence exactly / no)
- xerC cluster number (or N/A if not in a cluster)
- xerD cluster number (or N/A if not in a cluster)
- MOB type (not yet implemented)











(claude was used to help write and script 302 and to proofread my scripts and help me write this read_me, because Microsoft CoPilot Struggled to help at those parts. Copilot was used elsewhere, but all scripts were tested laboriously and largely thought up by Anthony)