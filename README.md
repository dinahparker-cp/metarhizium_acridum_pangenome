# _Metarhizium acridum_ pangenome
This repository hosts data and R code for Parker DM., Wilson AM., Nogueira C., Nielsen KN., Hansen LH., Nielsen TK., Habig M., De Fine Licht HH. Genome compartmentalization in a host-specific fungal insect pathogen reveals a putative mating type locus on an accessory chromosome. bioRxiv preprint 2025.08.27.672719. <br />

All scripts used for figure generation are included in Github. All data can be found on the linked ERDA (data repository host by the University of Copenhagen) sharelink that contains all scripts (also included on this Github) and additional processed data from outputs of specified programs for figure generation, details of which can be found in the manuscript. 
[ERDA data repository] (https://sid.erda.dk/cgi-sid/ls.py?share_id=G4eUxX0bO8) <br />

## DATA & FILE OVERVIEW

### **Figure 1**

1. **r01_map_pangenome_Fig1a.R** <br />
   R script to generate map showcasing the distribution of isolates used in the pangenome study of _M. acridum_, Fig 1a. <br />
   The data required to run this script is:
   - d01a_acridum_map_pg_metadata.xlsx 

2. **r02_table_formattable_Fig1b.R** <br />
   R script to generate table showcasing details of genome assemblies and other relevant genome information, Fig 1b. <br />
   The data required to run this script is:
   - d02a_table_genome_info.txt


### **Figure 2**

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
     

### **Figure 3**
1. **r01_fig3.Rmd** <br />
   R script to generate distribution plots of different gene subsets according to whether they were categorized as core, accessory, or singletons, Fig 3a-f. <br />
   The data required to run this script is:
   - d02_Orthogroups.tsv (from Figure_2 datafiles)
   - d03_Orthogroups_UnassignedGenes.tsv (from Figure_2 datafiles)
   - d01_KVL_data_subsets.zip (contains all annotation information required for running the script, must be unzipped)
   - d02_MAC_length.txt (contains length information on genes, not included in this figure)

### **Figure 4**

1. **r01_distribution_plots.R** <br />
   R script to generate distribution plots of of where different gene groups are found across the genome, by distance, Fig 4a-f. <br />
   The data required to run this script is:
   - distribution_data.zip (contains all intergenic distance information for each genome assembly, required for plot generation, must be unzipped)

### **Figure 5**

1. **script_shinyCircos_140426.R** <br />
   R script to generate circos plot shocasing mapping of reads to indicate accessory chromosome present on ARSEF 3391, Fig 5a. Additional information can be found in README_circos_plots.txt regarding generation of input files. <br />
   The data required to run this script is:
   - ARSEF_3391_GenomeCoord_10contigs.csv
   - KVL1801_on_KVL1801.regions_10contigs_logtrans.bed
   - KVL1802_on_KVL1801.regions_10contig_logtrans.bed
   - KVL1803_on_KVL1801.regions_10contigs_logtrans.bed.txt
   - KVL1806_on_KVL1801.regions_10contigs_logtrans.bed.txt
   - KVL1807_on_KVL1801.regions_10contigs_logtrans.bed.txt
   - KVL1808_on_KVL1801.regions_10contigs_logtrans.bed.txt
   - TE_windows_ARSEF3391.bed.txt
  
     
2. **R_script_karyotype_plot_170325.R** <br />
   R script to generate karyotype information regarding gene density vs. repeat density for ARSEF 3391, Fig 5d.  <br />
   The data required to run this script is:
   - GeneVsRepeat_plot.xlsx
   - Ma3391_geneLenght.txt
   - Macridum.filteredRepeats.bed
   - Macridum.filteredRepeats.gff
   
### **Figure 6**

1. **r01_fig_6.Rmd** <br />
   R script to generate accessory chromosome distributoin plots, to understand the functions that exist on the accessory chromosome. <br />
   The data required to run this script is:
   - d02_Orthogroups.tsv (from Figure_2 datafiles)
   - d03_Orthogroups_UnassignedGenes.tsv (from Figure_2 datafiles)
   - d01_Metarhizium_acridum_KVL_1801_H324.gff3
  
## METHODOLOGICAL INFORMATION

1. Methods for processing the data: R. Only processed data downstream of programs specified in the article is included in this Github repository, to explain figure generation.
2. People involved with data formatting and analysis: Dinah Parker, Henrik H. De Fine Licht (please contact dinahmparker@gmail.com for questions regarding data or scripts). 
