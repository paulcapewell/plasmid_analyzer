

------------
overall
-------------
This plasmid_analyze goes through a full data analysis workflow to examine xer sites in a plasmid database
* it reads a big plasmid database and identifies plasmids that match a frequency matrix for xerC/D bidning sites and ArgR binding sites
* then it pulls the relevant xerC/D bindign sequences + metadata for further analysis by phylogenetic trees and clinker maps of recombination


The `workflow_diagram.txt` may be very helpful to understand the workflow to use this
* you can either manually run the constituent scripts of P1, then P2, then P3, then P4,
* or you can run 001_clear_outfiles.py and 002_run_all_scripts.py to do the full analysis from scratch, but that would take a very long time


For technical information, look at the read_me files within P1, P2, P3, and P4

-----------------------
Output Files Overview
-----------------------
* P1 EMBOSS Scan outputs include data from EMBOSS profit scans for XerCD and ArgR (plasmid accession, position, % match)

* P2 XerCD Sequences Metadata outputs include a useful summary spreadsheet and trimed sequences such as xerC and xerD site isolated seperatly

* P3 Phylogenetics has phylogenetic trees

* P4 Clinker has schematic diagrams of sequence changes acrosss populations
