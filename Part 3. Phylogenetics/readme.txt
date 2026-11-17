#Phylogenetic Analysis of Xer Sites (upstream vs downstream)

*Input folder*
 
 sequences for phylogenetic analysis are copied from this directory '1. Analyze EMBOSS hits of XerCD sites /output_folder/isolate_xerc_xerd_upstream_downstream'
 In this case, 3kb sequences upstream and downstream from the xer crossover region are analyzed: XerC and upstream, XerD and downstream


*workflow*

Part 1) a phylogenetic tree is generated for each input file, using Jalview
 Using Jalview 2.11.5.2 on Mac, for each file a tree is claculated using neighbor joining and DNA
 This is saved in output_folder / `phylogenetic tree file`


Part 2) these phylogenetic trees are  uploaded to iTOL for analysis
	