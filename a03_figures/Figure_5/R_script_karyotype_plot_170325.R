#Load required packages
library(karyoploteR)
library(readr)
library(dplyr)
library(readr)



#Plot an ideogram using a custom genome

custom.genome <- toGRanges(data.frame(chr=c("chr_1", "chr_2", "chr_3", "chr_4", "chr_5", "chr_6", "chr_7", "chr_8", "chr_9", "chr_10"), start=c(1, 1), end=c(9909219, 9046696, 7782126, 5820841, 4615087, 3579740, 3148061, 903127, 526792, 90960))) 

ma_ig <- plotKaryotype(genome=custom.genome)


##ONLY GENES

#Load gff file for genes
kvl3391_genes_gff <- read_table("C:/Users/psb126/Documents_Cdrev/Henrik_temp/geneDensityPlot_Mb_Dinah/ARSEF_3391_genedensity/Metarhizium_acridum_ARSEF_3391.gff3", col_names = FALSE, comment = "#")
#Filter gff and convert to dataframe
kvl3391_genes <- kvl3391_genes_gff %>% filter(X3 == "exon") %>% select(X1, X4, X5)
kvl3391_genes <- as.data.frame(kvl3391_genes)

#Convert to GRanges object
kvl3391_genes_range <- toGRanges(kvl3391_genes) 

#Prepare custom genome
ma_ig <- plotKaryotype(genome=custom.genome, plot.type = 2, ideogram.plotter = NULL)


#Plot the genes on the positive side of the plot, adjust parameters as needed
kpPlotDensity(ma_ig, data=kvl3391_genes_range, r0=0, r1=2, window.size = 25000, col = "darkblue") ##4c83c3")


##GENES and TE

#Load gff file for TE
kvl3391_TE_gff <- read_table("C:/Users/psb126/Documents_Cdrev/Henrik_temp/geneDensityPlot_Mb_Dinah/ARSEF_3391_genedensity/Macridum.filteredRepeats.gff", col_names = FALSE, comment = "#")

#Filter gff and convert to dataframe
kvl3391_TE <- kvl3391_TE_gff %>% select(X1, X4, X5)
kvl3391_TE <- as.data.frame(kvl3391_TE)

#Convert to GRanges object
kvl3391_TE_range <- toGRanges(kvl3391_TE)

#Plot the repetitive elements on the negative side of the plot, adjust parameters as needed
kpPlotDensity(ma_ig, data=kvl3391_TE_range, data.panel = 1, r0=0, r1=-2, window.size = 25000, col = "royalblue") #" #CCFFAA")

#kpAxis(ma_ig, ymin = 0, ymax=1, data.panel= 2, side = 2)
#kpAxis(ma_ig, ymin = 0, ymax=1, data.panel = 1, side = 2)

test <- ma_ig$latest.plot$computed.values$max.density

# plot a rectangle highlighting specific positions
# Mating type locus on Chr_2 in this case
kpRect(ma_ig, chr="chr_2", x0=3473000, x1=3580000, y0=0.2, y1=0.8)

# Mating type locus on Chr_8 in this case
kpRect(ma_ig, chr="chr_8", x0=305000, x1=405000, y0=0.2, y1=0.8)
