# _Metarhizium acridum_ pangenome
This repository hosts data and R code for Parker DM., Wilson AM., Nogueira C., Nielsen KN., Hansen LH., Nielsen TK., Habig M., De Fine Licht HH. Genome compartmentalization in a host-specific fungal insect pathogen reveals a putative mating type locus on an accessory chromosome. bioRxiv preprint 2025.08.27.672719. <br />

This is the ERDA (data repository host by the University of Copenhagen) sharelink that contains all scripts (also included on this Github) and additional processed data from outputs of specified programs for figure generation. 
[ERDA data repository] (https://sid.erda.dk/cgi-sid/ls.py?share_id=G4eUxX0bO8) <br />

## DATA & FILE OVERVIEW

**Figure 1**

1. **r01_map_pangenome_Fig1a.R** <br />
   R script to generate map showcasing the distribution of isolates used in the pangenome study of _M. acridum_, Fig 1a. <br />
   The data required to run this script is:
   - d01a_acridum_map_pg_metadata.xlsx 

2. **r02_table_formattable_Fig1b.R** <br />
   R script to generate table showcasing details of genome assemblies and other relevant genome information, Fig 1b. <br />
   The data required to run this script is:
   - d02a_table_genome_info.txt


**Figure 2**

1. **r01_link_to_genespace_Fig2a.txt** <br />
   Link to evernote describing the methods used to generate the GENESPACE gene synteny file, Fig 2a. <br />
   The data required to run this script is:
   - d01_genespace.zip

2. **r02_finding_orthogroups_2026.R** <br />
   R script to generate circular plot of core and accessory genes of all _Metarhizium_ isolates, Fig 2b. <br />
   The data required to run this script is:
   - d02_Orthogroups.tsv
   - d03_Orthogroups_UnassignedGenes.tsv
   - d04_ortho_detail_2.txt
   - d05_orthogroup_info_2.txt
   - d06_summary_table_s8.txt
     

**Figure 3**
1. **r01_fig3.Rmd** <br />
   R script to generate distribution plots of different gene subsets according to whether they were categorized as core, accessory, or singletons, Fig 3a-f. <br />
   The data required to run this script is:
   - d02_Orthogroups.tsv (from Figure_2 datafiles)
   - d03_Orthogroups_UnassignedGenes.tsv (from Figure_2 datafiles)
   - d01_KVL_data_subsets.zip (contains all annotation information required for running the script, must be unzipped)
   - d02_MAC_length.txt (contains length information on genes, not included in this figure)




**Figure 4**

**Figure 5**

**Figure 6**
