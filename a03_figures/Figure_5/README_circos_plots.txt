# README file for circos plots on M. acridum ARSEF 3391

# First generate a genome coordinates file with the contig lengths that will form the base of the plot
# the file has three columns: name, start, end and looks like this (name=ARSEF_3391_genomeCoord_10contigs.csv
scaffold-,start,end
scaffold-1,1,9909219
scaffold-2,1,9046696
scaffold-3,1,7782126
scaffold-4,1,5820841
scaffold-5,1,4615087
scaffold-6,1,3579740
scaffold-7,1,3148061
scaffold-8,1,903127
scaffold-9,1,526792
scaffold-10,1,90960


# Next generate read coverage files from mapped BAM files generated with minimap
# First index the bamfile
samtools indext *.bam

# Use mosdepth to calculate coverage in 10.000 bp regions across the genome and output as bed file
mosdepth -b 10000 -t 10 -n Input.bam.file

# Use AWK to convert coverage reads to Log10 of reads (nb. log in AWK is the natural log, hence this calculation)
# ex. command - run for each
awk -F"\t" '{a = log(1+$4)/log(10); print $1, $2, $3, a, "\n"}' KVL1803_on_KVL1801.regions.bed > KVL1803_on_KVL1801.regions_logtrans.bed

# -F species tab delimiter, a = calculates the log transformed  values, and print prints the different columns, substituting the original coverages with log10 transformed values

############### TE's ########################
# From the output earlGrey use the bed file output with the annotated TE's

# Inside Rstudio run the following to calculate the percentage TE's in 10000 bp windows across the genome
# use the function: genomicDensity (from circlize package)
test <- genomicDensity(Macridum.filteredRepeats_modified4column, window.size = 10000, overlap = FALSE, count_by = "percent")

# write the table to a file
write.table(test,"TE_windows_ARSEF3391.bed", sep = "  ", row.names = FALSE)

############### Circos Plot #################
# inside Rstudio the following is run to start the Rshiny app
install.packages("shiny")  
install.packages("circlize")  
install.packages("RColorBrewer")
install.packages("data.table")
install.packages("RLumShiny")  

source("https://bioconductor.org/biocLite.R")  
biocLite("GenomicRanges")

# This starts the shiny app
shiny::runGitHub("shinyCircos", "venyao")  


# in the shiny app choose:
bar plots for all tracks
Chromosome height = 0.04
unidirectional bar direction
Custom for data with multi-column (color) = #5C85FF
Leave baselines color and Background colors blank

############################ OUTPUT SCRIPT FROM shinyCircos  ##############################

# Output script from ShinyCircos, i.e. the actual code/commands used are stored in script_shinyCircos_140425.R
