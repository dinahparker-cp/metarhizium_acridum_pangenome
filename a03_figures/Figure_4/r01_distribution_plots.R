# Load required libraries
library(dplyr)
library(hexbin)
library(RColorBrewer)
library(tidyverse)
library(cowplot)
library(gmodels)
library(DescTools)
library(ggplot2)
library(qqplotr)
library(grid)
library(forcats)
library(scales)
library(grid)
library(gridtext)
library(gridExtra)
library(rcompanion)

# File paths and genome IDs
file_paths <- list.files(
  path = "/Users/mbs625/Downloads/distribution_data/p0_gff",
  pattern = "Metarhizium_acridum_KVL_.*_H324.gff3",
  full.names = TRUE,
  recursive = TRUE
)

genome_ids <- c("MAC1", "MAC2", "MAC3", "MAC6", "MAC7", "MAC8") # Updated to match desired format

# Function to read and preprocess a single GFF3 file
read_gff <- function(file_path, genome_id) {
  read.delim(file_path, header = FALSE, comment.char = "#") %>%
    mutate(genome = genome_id, V9 = gsub("FUN", genome_id, V9)) # Correct replacement for "FUN"
}

# Read all files and combine into a single data frame
MAC_gff3 <- map2_dfr(file_paths, genome_ids, read_gff)

# Filter for "gene" features
MAC_gff3 <- MAC_gff3 %>%
  filter(V3 == "gene") %>%
  mutate(
    V9 = gsub("ID=", "", V9), # Remove "ID=" prefix
    V9 = gsub(";.*", "", V9)  # Remove everything after the first semicolon
  ) %>%
  select(V1, V4, V5, V9, genome)

# Remove single-gene contigs
MAC_gff3 <- MAC_gff3 %>%
  group_by(V1) %>%
  filter(n() != 1) %>%
  ungroup()

# Transform variables into numeric
MAC_gff3$V4 <- as.numeric(MAC_gff3$V4)
MAC_gff3$V5 <- as.numeric(MAC_gff3$V5)

# Calculate 3' and 5' lagging values
MAC_gff3 <- MAC_gff3 %>%
  group_by(V5) %>% 
  group_by(V1) %>%
  mutate(
    ig3 = (V4 - lag(V5)) / 1000, # 3' intergenic region
    ig5 = lead(ig3)              # 5' intergenic region
  ) %>%
  ungroup()

# Remove negative and NA values
MAC_35p_gff <- MAC_gff3 %>%
  filter(!is.na(ig3), !is.na(ig5), ig3 >= 0, ig5 >= 0)


# File paths for effector genes and corresponding genome IDs
effector_files <- list.files(
  path = "/Users/mbs625/Downloads/distribution_data/p1_effectors",
  pattern = "eff_F_.*\\.txt",
  full.names = TRUE,
  recursive = TRUE
)

genome_ids <- c("MAC1", "MAC2", "MAC3", "MAC6", "MAC7", "MAC8")

# Function to read and preprocess effector files
read_effector <- function(file_path, genome_id) {
  read.delim(file_path, header = FALSE, comment.char = "#") %>%
    filter(!grepl("^Y", V4)) %>% # Remove rows with "Y.*" in column V4
    select(V1) %>%              # Keep only column V1
    mutate(
      V1 = gsub("-T1", "", V1),
      V1 = gsub("FUN", genome_id, V1)
    )
}

# Read and preprocess all effector files
effector_data <- map2_dfr(effector_files, genome_ids, read_effector)

# Merge effector data with the GFF3 dataset
eff_MAC <- effector_data %>%
  inner_join(MAC_35p_gff, by = c("V1" = "V9"))

# Add a new column to MAC_35p_gff to check effector presence
MAC_35p_gff$check <- eff_MAC$V1[match(MAC_35p_gff$V9, eff_MAC$V1)]
MAC_35p_gff$check[is.na(MAC_35p_gff$check)] <- 0
MAC_35p_gff$check[MAC_35p_gff$check != 0] <- 1


#PLOTTING EFFECTORS
##gg_hex
p_grob <- grobTree(richtext_grob("p<0.0001", x=0.05, y=0.95, hjust=0,
                                 gp=gpar(col="black", fontsize=13, fontface="bold"), 
                                 box_gp = gpar(col = "black"),
                                 padding = unit(c(6, 6, 4, 6), "pt")))
hex_MAC <- ggplot(MAC_35p_gff, aes(ig3, ig5)) +
  geom_hex(bins=40, show.legend = TRUE) +
  scale_fill_continuous(low = "grey80", high = "black") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  geom_point(data=eff_MAC, aes(x=ig3, y=ig5), 
             color = "#0096FF", size =1.5, alpha = 0.5, inherit.aes = FALSE) +
  theme_minimal() +
  labs(x ="5' Flanking intergenic region (kbp)", y = "3' Flanking intergenic region (kbp)") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position ="none",  axis.text.x = element_text(color="black", size=25), axis.text.y = element_text(color="black", size=25)) +
  geom_hline(yintercept=1, color = "black", linetype="dashed", alpha=0.8) + 
  geom_vline(xintercept=1, color = "black", linetype="dashed", alpha=0.8) 

#log trans reverse 
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

#make density plots for ig3 and ig5
n1 <- ggplot(data=MAC_35p_gff, aes(ig3, colour = check)) +
  geom_line(show.legend = FALSE, stat = "density", alpha = 0.5) +
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(name = "Density", breaks = c(0,1.0), limits = c(0,1.25)) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.line.x = element_line("grey"), axis.line.y = element_line("grey"), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(color="black", size=20, face="bold")) +
  scale_color_manual(values=c("darkgrey", "#0096FF"))

n2 <- ggplot(data=MAC_35p_gff, aes(x=ig5, colour = check)) +
  geom_line(show.legend = FALSE, stat = "density", alpha = 0.5) +
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(name = "Density", breaks = c(0,1.0), limits = c(0,1.2)) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(color="black", size=20, face="bold"), axis.line.x = element_line("grey"), axis.line.y = element_line("grey"), axis.ticks.y = element_blank()) +
  scale_color_manual(values=c("darkgrey", "#0096FF")) + coord_flip()

#stack hex and density plots 
grid.newpage()
n1 <- n1 + theme(plot.margin = unit(c(0, 0, 0, 2.7), "lines"))
hex_MAC <- hex_MAC +
  theme(plot.margin = unit(c(0, 0.5, 0, 0), "lines"))
n2 <- n2 + theme(plot.margin = unit(c(0, 0, 1.2, 0), "lines"))
grid.arrange(n1, hex_MAC, n2,
             ncol = 2,
             widths = c(10,1),
             heights = c(1.2,10),
             layout_matrix = rbind(c(1, NA),
                                   c(2, 3)))



#Effector descriptive stats

MAC_35p_gff$check <- as.factor(MAC_35p_gff$check)

# Descriptive statistics by group
MAC_35p_gff %>%
  dplyr::select(check, ig3) %>%
  group_by(check) %>%
  summarise(
    n = n(),
    mean = mean(ig3, na.rm = TRUE),
    sd = sd(ig3, na.rm = TRUE),
    stderr = sd/sqrt(n),
    LCL = mean - qt(1 - (0.05 / 2), n - 1) * stderr,
    UCL = mean + qt(1 - (0.05 / 2), n - 1) * stderr,
    median = median(ig3, na.rm = TRUE),
    min = min(ig3, na.rm = TRUE),
    max = max(ig3, na.rm = TRUE),
    IQR = IQR(ig3, na.rm = TRUE),
    LCLmed = MedianCI(ig3, na.rm=TRUE)[2],
    UCLmed = MedianCI(ig3, na.rm=TRUE)[3]
  )

# Mann-Whitney U tests
m1 <- wilcox.test(ig3 ~ check, data = MAC_35p_gff,
                  exact = FALSE, conf.int = TRUE)
print(m1)

m2 <- wilcox.test(ig5 ~ check, data = MAC_35p_gff,
                  exact = FALSE, conf.int = TRUE)
print(m2)

# Rank biserial correlation
# For ig3
r_rb_ig3 <- wilcoxonR(x = MAC_35p_gff$ig3, g = MAC_35p_gff$check)
r_rb_ig3

# For ig5
r_rb_ig5 <- wilcoxonR(x = MAC_35p_gff$ig5, g = MAC_35p_gff$check)
r_rb_ig5

# SIGNIFICANT

pvals <- c(m1$p.value, m2$p.value)
pvals_adj <- p.adjust(pvals, method = "BH")
pvals_adj



#SECRETOME
# Define file paths and genome IDs for secretome files
sec_file_paths <- list.files(
  path = "/Users/mbs625/Downloads/distribution_data/p2_secretome",
  pattern = "MAC_.*_secf.txt",
  full.names = TRUE,
  recursive = TRUE
)

sec_genome_ids <- c("MAC1", "MAC2", "MAC3", "MAC6", "MAC7", "MAC8")

# Function to process each secretome file
process_secretome <- function(file_path, genome_id) {
  read.delim(file_path, header = FALSE, comment.char = "#") %>%
    mutate(
      V1 = gsub("FUN", genome_id, V1),
      V1 = gsub("-T1", "", V1)
    )
}

# Process and merge all secretome files
MAC_sec <- map2_dfr(sec_file_paths, sec_genome_ids, process_secretome)

#merge gff3 and secretome
secretome_MAC <- MAC_sec %>% 
  inner_join(MAC_35p_gff, by = c("V1" = "V9"))

#create new column
MAC_35p_gff$check_sec <- secretome_MAC$V1[match(MAC_35p_gff$V9, secretome_MAC$V1)]
MAC_35p_gff$check_sec[is.na(MAC_35p_gff$check_sec)] <- 0
MAC_35p_gff$check_sec[MAC_35p_gff$check_sec != 0] <- 1

#plotting secretome
##gg_hex
p_grob <- grobTree(richtext_grob("p<0.0001", x=0.05, y=0.95, hjust=0,
                                 gp=gpar(col="black", fontsize=13, fontface="bold"), 
                                 box_gp = gpar(col = "black"),
                                 padding = unit(c(6, 6, 4, 6), "pt")))
sec_MAC_plot <- ggplot(MAC_35p_gff, aes(ig3, ig5)) +
  geom_hex(bins=40, show.legend = TRUE) +
  scale_fill_continuous(low = "grey80", high = "black") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  geom_point(data=secretome_MAC, aes(x=ig3, y=ig5), 
             color = "#0096FF", size =1.5, alpha = 0.3, inherit.aes = TRUE) +
  theme_minimal() +
  labs(x ="5' Flanking intergenic region (kbp)", y = "3' Flanking intergenic region (kbp)") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position ="none",  axis.text.x = element_text(color="black", size=25), axis.text.y = element_text(color="black", size=25)) +
  geom_hline(yintercept=1, color = "black", linetype="dashed", alpha=0.8) + 
  geom_vline(xintercept=1, color = "black", linetype="dashed", alpha=0.8) 

#make density plots for ig3 and ig5
n3 <- ggplot(data=MAC_35p_gff, aes(ig3, colour = check_sec)) +
  geom_line(show.legend = FALSE, stat = "density", alpha = 0.5) +
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(name = "Density", breaks = c(0,1.0), limits = c(0,1.25)) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.line.x = element_line("grey"), axis.line.y = element_line("grey"), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(color="black", size=20, face="bold")) +
  scale_color_manual(values=c("darkgrey", "#0096FF"))

n4 <- ggplot(data=MAC_35p_gff, aes(x=ig5, colour = check_sec)) +
  geom_line(show.legend = FALSE, stat = "density", alpha = 0.5) +
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(name = "Density", breaks = c(0,1.0), limits = c(0,1.2)) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(color="black", size=20, face="bold"), axis.line.x = element_line("grey"), axis.line.y = element_line("grey"), axis.ticks.y = element_blank()) +
  scale_color_manual(values=c("darkgrey", "#0096FF")) + coord_flip()

#stack hex and density plots 
grid.newpage()
n3 <- n3 + theme(plot.margin = unit(c(0, 0, 0, 2.7), "lines"))
sec_MAC_plot <- sec_MAC_plot + theme(plot.margin = unit(c(0, 0.5, 0, 0), "lines"))
n4 <- n4 + theme(plot.margin = unit(c(0, 0, 1.2, 0), "lines"))
grid.arrange(n3, sec_MAC_plot, n4,
             ncol = 2,
             widths = c(10,1),
             heights = c(1.2,10),
             layout_matrix = rbind(c(1, NA),
                                   c(2, 3)))



#Secretome descriptive stats

MAC_35p_gff$check_sec <- as.factor(MAC_35p_gff$check_sec)

# Descriptive statistics by group
MAC_35p_gff %>%
  dplyr::select(check_sec, ig3) %>%
  group_by(check_sec) %>%
  summarise(
    n = n(),
    mean = mean(ig3, na.rm = TRUE),
    sd = sd(ig3, na.rm = TRUE),
    stderr = sd/sqrt(n),
    LCL = mean - qt(1 - (0.05 / 2), n - 1) * stderr,
    UCL = mean + qt(1 - (0.05 / 2), n - 1) * stderr,
    median = median(ig3, na.rm = TRUE),
    min = min(ig3, na.rm = TRUE),
    max = max(ig3, na.rm = TRUE),
    IQR = IQR(ig3, na.rm = TRUE),
    LCLmed = MedianCI(ig3, na.rm=TRUE)[2],
    UCLmed = MedianCI(ig3, na.rm=TRUE)[3]
  )



# Mann-Whitney U tests
m1_sec <- wilcox.test(ig3 ~ check_sec, data = MAC_35p_gff,
                  exact = FALSE, conf.int = TRUE)
print(m1_sec)

m2_sec <- wilcox.test(ig5 ~ check_sec, data = MAC_35p_gff,
                  exact = FALSE, conf.int = TRUE)
print(m2_sec)

# NOT SIGNIFICANT

pvals_sec <- c(m1_sec$p.value, m2_sec$p.value)
pvals_adj_sec <- p.adjust(pvals_sec, method = "BH")
pvals_adj_sec






#CAZYME
# Define file paths and genome IDs for CAZyme files
caz_file_paths <- list.files(
  path = "/Users/mbs625/Downloads/distribution_data/p3_cazymes",
  pattern = "annotations.dbCAN.txt",
  full.names = TRUE,
  recursive = TRUE
)

caz_genome_ids <- c("MAC1", "MAC2", "MAC3", "MAC6", "MAC7", "MAC8")

# Function to process each CAZyme file
process_cazyme <- function(file_path) {
  read.delim(file_path, header = FALSE, comment.char = "#") %>%
    select(V1) %>% # Remove unwanted columns
    mutate(V1 = gsub("-T1", "", V1))
}

# Process and merge all CAZyme files
MAC_caz <- map_dfr(caz_file_paths, process_cazyme)

#merge gff3 and secretome
cazyme_MAC <- MAC_caz %>% 
  inner_join(MAC_35p_gff, by = c("V1" = "V9"))

#create new column
MAC_35p_gff$check_caz <- cazyme_MAC$V1[match(MAC_35p_gff$V9, cazyme_MAC$V1)]
MAC_35p_gff$check_caz[is.na(MAC_35p_gff$check_caz)] <- 0
MAC_35p_gff$check_caz[MAC_35p_gff$check_caz != 0] <- 1

#plotting cazymes
##gg_hex
p_grob <- grobTree(richtext_grob("p<0.0001", x=0.05, y=0.95, hjust=0,
                                 gp=gpar(col="black", fontsize=13, fontface="bold"), 
                                 box_gp = gpar(col = "black"),
                                 padding = unit(c(6, 6, 4, 6), "pt")))
caz_MAC_plot <- ggplot(MAC_35p_gff, aes(ig3, ig5)) +
  geom_hex(bins=40, show.legend = TRUE) +
  scale_fill_continuous(low = "grey80", high = "black") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  geom_point(data=cazyme_MAC, aes(x=ig3, y=ig5), 
             color = "#0096FF", size =1.5, alpha = 0.6, inherit.aes = TRUE) +
  theme_minimal() +
  labs(x ="5' Flanking intergenic region (kbp)", y = "3' Flanking intergenic region (kbp)") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position ="none",  axis.text.x = element_text(color="black", size=25), axis.text.y = element_text(color="black", size=25)) +
  geom_hline(yintercept=1, color = "black", linetype="dashed", alpha=0.8) + 
  geom_vline(xintercept=1, color = "black", linetype="dashed", alpha=0.8) 


#make density plots for ig3 and ig5
n5 <- ggplot(data=MAC_35p_gff, aes(ig3, colour = check_caz)) +
  geom_line(show.legend = FALSE, stat = "density", alpha = 0.5) +
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(name = "Density", breaks = c(0,1.0), limits = c(0,1.25)) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.line.x = element_line("grey"), axis.line.y = element_line("grey"), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(color="black", size=20, face="bold")) +
  scale_color_manual(values=c("darkgrey", "#0096FF"))

n6 <- ggplot(data=MAC_35p_gff, aes(x=ig5, colour = check_caz)) +
  geom_line(show.legend = FALSE, stat = "density", alpha = 0.5) +
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(name = "Density", breaks = c(0,1.0), limits = c(0,1.2)) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(color="black", size=20, face="bold"), axis.line.x = element_line("grey"), axis.line.y = element_line("grey"), axis.ticks.y = element_blank()) +
  scale_color_manual(values=c("darkgrey", "#0096FF")) + coord_flip()

#stack hex and density plots 
grid.newpage()
n5 <- n5 + theme(plot.margin = unit(c(0, 0, 0, 2.7), "lines"))
sec_MAC_plot <- caz_MAC_plot + theme(plot.margin = unit(c(0, 0.5, 0, 0), "lines"))
n6 <- n6  + theme(plot.margin = unit(c(0, 0, 1.2, 0), "lines"))
grid.arrange(n5, sec_MAC_plot, n6,
             ncol = 2,
             widths = c(10,1),
             heights = c(1.2,10),
             layout_matrix = rbind(c(1, NA),
                                   c(2, 3)))



#CAZYME descriptive stats

MAC_35p_gff$check_caz <- as.factor(MAC_35p_gff$check_caz)

# Descriptive statistics by group
MAC_35p_gff %>%
  dplyr::select(check_caz, ig3) %>%
  group_by(check_caz) %>%
  summarise(
    n = n(),
    mean = mean(ig3, na.rm = TRUE),
    sd = sd(ig3, na.rm = TRUE),
    stderr = sd/sqrt(n),
    LCL = mean - qt(1 - (0.05 / 2), n - 1) * stderr,
    UCL = mean + qt(1 - (0.05 / 2), n - 1) * stderr,
    median = median(ig3, na.rm = TRUE),
    min = min(ig3, na.rm = TRUE),
    max = max(ig3, na.rm = TRUE),
    IQR = IQR(ig3, na.rm = TRUE),
    LCLmed = MedianCI(ig3, na.rm=TRUE)[2],
    UCLmed = MedianCI(ig3, na.rm=TRUE)[3]
  )

# Mann-Whitney U tests
m1_caz <- wilcox.test(ig3 ~ check_caz, data = MAC_35p_gff,
                  exact = FALSE, conf.int = TRUE)
print(m1_caz)

m2_caz <- wilcox.test(ig5 ~ check_caz, data = MAC_35p_gff,
                  exact = FALSE, conf.int = TRUE)
print(m2_caz)



# SIGNIFICANT
pvals_caz <- c(m1_caz$p.value, m2_caz$p.value)
pvals_adj_caz <- p.adjust(pvals_caz, method = "BH")
pvals_adj_caz


#PEPTIDASE
# Define file paths for Peptidase_S8 files
pep_file_paths <- list.files(
  path = "/Users/mbs625/Downloads/distribution_data/p5_peptidase",
  pattern = "peptidaseS8.txt",
  full.names = TRUE,
  recursive = TRUE
)

# Function to process each Peptidase_S8 file
process_peptidase <- function(file_path) {
  read.delim(file_path, header = FALSE, comment.char = "#") %>%
    mutate(V1 = gsub("-T1", "", V1))
}

# Process and merge all Peptidase_S8 files
MAC_pep <- map_dfr(pep_file_paths, process_peptidase)

#merge gff3 and secretome
peptidase_MAC <- MAC_pep %>% 
  inner_join(MAC_35p_gff, by = c("V1" = "V9"))

#create new column
MAC_35p_gff$check_pep <- peptidase_MAC$V1[match(MAC_35p_gff$V9, peptidase_MAC$V1)]
MAC_35p_gff$check_pep[is.na(MAC_35p_gff$check_pep)] <- 0
MAC_35p_gff$check_pep[MAC_35p_gff$check_pep != 0] <- 1

#plotting peptidases
##gg_hex
p_grob <- grobTree(richtext_grob("p<0.0001", x=0.05, y=0.95, hjust=0,
                                 gp=gpar(col="black", fontsize=13, fontface="bold"), 
                                 box_gp = gpar(col = "black"),
                                 padding = unit(c(6, 6, 4, 6), "pt")))
pep_MAC_plot <- ggplot(MAC_35p_gff, aes(ig3, ig5)) +
  geom_hex(bins=40, show.legend = TRUE) +
  scale_fill_continuous(low = "grey80", high = "black") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  geom_point(data=peptidase_MAC, aes(x=ig3, y=ig5), 
             color = "#0096FF", size = 2.5, alpha = 0.8, inherit.aes = TRUE) +
  theme_minimal() +
  labs(x ="5' Flanking intergenic region (kbp)", y = "3' Flanking intergenic region (kbp)") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position ="none",  axis.text.x = element_text(color="black", size=25), axis.text.y = element_text(color="black", size=25)) +
  geom_hline(yintercept=1, color = "black", linetype="dashed", alpha=0.8) + 
  geom_vline(xintercept=1, color = "black", linetype="dashed", alpha=0.8) 

#make density plots for ig3 and ig5
n7 <- ggplot(data=MAC_35p_gff, aes(ig3, colour = check_pep)) +
  geom_line(show.legend = FALSE, stat = "density", alpha = 0.5) +
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(name = "Density", breaks = c(0,1.0), limits = c(0,1.25)) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.line.x = element_line("grey"), axis.line.y = element_line("grey"), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(color="black", size=20, face="bold")) +
  scale_color_manual(values=c("darkgrey", "#0096FF"))

n8 <- ggplot(data=MAC_35p_gff, aes(x=ig5, colour = check_pep)) +
  geom_line(show.legend = FALSE, stat = "density", alpha = 0.5) +
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(name = "Density", breaks = c(0,1.0), limits = c(0,1.2)) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(color="black", size=20, face="bold"), axis.line.x = element_line("grey"), axis.line.y = element_line("grey"), axis.ticks.y = element_blank()) +
  scale_color_manual(values=c("darkgrey", "#0096FF")) + coord_flip()


#stack hex and density plots 
grid.newpage()
n7 <- n7 + theme(plot.margin = unit(c(0, 0, 0, 2.7), "lines"))
pep_MAC_plot <- pep_MAC_plot + theme(plot.margin = unit(c(0, 0.5, 0, 0), "lines"))
n8 <- n8 + theme(plot.margin = unit(c(0, 0, 1.2, 0), "lines"))
grid.arrange(n7, pep_MAC_plot, n8,
             ncol = 2,
             widths = c(10,1),
             heights = c(1.2,10),
             layout_matrix = rbind(c(1, NA),
                                   c(2, 3)))

#Peptidase descriptive stats

MAC_35p_gff$check_pep <- as.factor(MAC_35p_gff$check_pep)

# Descriptive statistics by group
MAC_35p_gff %>%
  dplyr::select(check_pep, ig3) %>%
  group_by(check_pep) %>%
  summarise(
    n = n(),
    mean = mean(ig3, na.rm = TRUE),
    sd = sd(ig3, na.rm = TRUE),
    stderr = sd/sqrt(n),
    LCL = mean - qt(1 - (0.05 / 2), n - 1) * stderr,
    UCL = mean + qt(1 - (0.05 / 2), n - 1) * stderr,
    median = median(ig3, na.rm = TRUE),
    min = min(ig3, na.rm = TRUE),
    max = max(ig3, na.rm = TRUE),
    IQR = IQR(ig3, na.rm = TRUE),
    LCLmed = MedianCI(ig3, na.rm=TRUE)[2],
    UCLmed = MedianCI(ig3, na.rm=TRUE)[3]
  )

# Mann-Whitney U tests
m1_pep <- wilcox.test(ig3 ~ check_pep, data = MAC_35p_gff,
                      exact = FALSE, conf.int = TRUE)
print(m1_pep)

m2_pep <- wilcox.test(ig5 ~ check_pep, data = MAC_35p_gff,
                      exact = FALSE, conf.int = TRUE)
print(m2_pep)



# SIGNIFICANT
pvals_pep <- c(m1_pep$p.value, m2_pep$p.value)
pvals_adj_pep <- p.adjust(pvals_pep, method = "BH")
pvals_adj_pep




#Chitinases
#Define file paths for Chitinase files
chit_file_paths <- list.files(
  path = "/Users/mbs625/Downloads/distribution_data/p4_chitinase",
  pattern = "chitinase.txt",
  full.names = TRUE,
  recursive = TRUE
)

# Function to process each Chitinase file
process_chitinase <- function(file_path) {
  read.delim(file_path, header = FALSE, comment.char = "#") %>%
    mutate(V1 = gsub("-T1", "", V1))
}

# Process and merge all Chitinase files
MAC_chit <- map_dfr(chit_file_paths, process_chitinase)

#merge gff3 and chitins
chitinase_MAC <- MAC_chit %>% 
  inner_join(MAC_35p_gff, by = c("V1" = "V9"))

#create new column
MAC_35p_gff$check_chit <- chitinase_MAC$V1[match(MAC_35p_gff$V9, chitinase_MAC$V1)]
MAC_35p_gff$check_chit[is.na(MAC_35p_gff$check_chit)] <- 0
MAC_35p_gff$check_chit[MAC_35p_gff$check_chit != 0] <- 1


#plotting chitinases
##gg_hex
p_grob <- grobTree(richtext_grob("p<0.0001", x=0.05, y=0.95, hjust=0,
                                 gp=gpar(col="black", fontsize=13, fontface="bold"), 
                                 box_gp = gpar(col = "black"),
                                 padding = unit(c(6, 6, 4, 6), "pt")))
chit_MAC_plot <- ggplot(MAC_35p_gff, aes(ig3, ig5)) +
  geom_hex(bins=40, show.legend = TRUE) +
  scale_fill_continuous(low = "grey80", high = "black") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  geom_point(data=chitinase_MAC, aes(x=ig3, y=ig5), 
             color = "#0096FF", size =2.5, alpha = 0.8, inherit.aes = TRUE) +
  theme_minimal()  +
  labs(x ="5' Flanking intergenic region (kbp)", y = "3' Flanking intergenic region (kbp)") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position ="none",  axis.text.x = element_text(color="black", size=25), axis.text.y = element_text(color="black", size=25)) +
  geom_hline(yintercept=1, color = "black", linetype="dashed", alpha=0.8) + 
  geom_vline(xintercept=1, color = "black", linetype="dashed", alpha=0.8) 

#make density plots for ig3 and ig5
n9 <- ggplot(data=MAC_35p_gff, aes(ig3, colour = check_chit)) +
  geom_line(show.legend = FALSE, stat = "density", alpha = 0.5) +
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(name = "Density", breaks = c(0,1.0), limits = c(0,1.25)) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.line.x = element_line("grey"), axis.line.y = element_line("grey"), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(color="black", size=20, face="bold")) +
  scale_color_manual(values=c("darkgrey", "#0096FF"))

n10 <- ggplot(data=MAC_35p_gff, aes(x=ig5, colour = check_pep)) +
  geom_line(show.legend = FALSE, stat = "density", alpha = 0.5) +
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(name = "Density", breaks = c(0,1.0), limits = c(0,1.2)) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(color="black", size=20, face="bold"), axis.line.x = element_line("grey"), axis.line.y = element_line("grey"), axis.ticks.y = element_blank()) +
  scale_color_manual(values=c("darkgrey", "#0096FF")) + coord_flip()


#stack hex and density plots 
grid.newpage()
n9 <- n9 + theme(plot.margin = unit(c(0, 0, 0, 2.7), "lines"))
chit_MAC_plot <- chit_MAC_plot + theme(plot.margin = unit(c(0, 0.5, 0, 0), "lines"))
n10 <- n10 + theme(plot.margin = unit(c(0, 0, 1.2, 0), "lines"))
grid.arrange(n9, chit_MAC_plot, n10,
             ncol = 2,
             widths = c(10,1),
             heights = c(1.2,10),
             layout_matrix = rbind(c(1, NA),
                                   c(2, 3)))


#Chitinase descriptive stats

MAC_35p_gff$check_chit <- as.factor(MAC_35p_gff$check_chit)

# Descriptive statistics by group
MAC_35p_gff %>%
  dplyr::select(check_chit, ig3) %>%
  group_by(check_chit) %>%
  summarise(
    n = n(),
    mean = mean(ig3, na.rm = TRUE),
    sd = sd(ig3, na.rm = TRUE),
    stderr = sd/sqrt(n),
    LCL = mean - qt(1 - (0.05 / 2), n - 1) * stderr,
    UCL = mean + qt(1 - (0.05 / 2), n - 1) * stderr,
    median = median(ig3, na.rm = TRUE),
    min = min(ig3, na.rm = TRUE),
    max = max(ig3, na.rm = TRUE),
    IQR = IQR(ig3, na.rm = TRUE),
    LCLmed = MedianCI(ig3, na.rm=TRUE)[2],
    UCLmed = MedianCI(ig3, na.rm=TRUE)[3]
  )

# Mann-Whitney U tests
m1_chit <- wilcox.test(ig3 ~ check_chit, data = MAC_35p_gff,
                      exact = FALSE, conf.int = TRUE)
print(m1_chit)

m2_chit <- wilcox.test(ig5 ~ check_chit, data = MAC_35p_gff,
                      exact = FALSE, conf.int = TRUE)
print(m2_chit)



# SIGNIFICANT
pvals_chit <- c(m1_chit$p.value, m2_chit$p.value)
pvals_adj_chit <- p.adjust(pvals_chit, method = "BH")
pvals_adj_chit



#ANTISMASH
# Define file paths for secondary metabolites (antiSMASH) files
anti_file_paths <- list.files(
  path = "/Users/mbs625/Downloads/distribution_data/p6_antismash",
  pattern = "annotations.antismash.txt",
  full.names = TRUE,
  recursive = TRUE
)

# Function to process each antiSMASH file
process_antismash <- function(file_path) {
  read.delim(file_path, header = FALSE, comment.char = "#") %>%
    select(-V3) %>% 
    mutate(V1 = gsub("-T1", "", V1))
}

# Process and merge all antiSMASH files
MAC_anti <- map_dfr(anti_file_paths, process_antismash)


#merge gff3 and secretome
antismash_MAC <- MAC_anti %>% 
  inner_join(MAC_35p_gff, by = c("V1" = "V9"))

#create new column
MAC_35p_gff$check_anti <- antismash_MAC$V1[match(MAC_35p_gff$V9, antismash_MAC$V1)]
MAC_35p_gff$check_anti[is.na(MAC_35p_gff$check_anti)] <- 0
MAC_35p_gff$check_anti[MAC_35p_gff$check_anti != 0] <- 1

#plotting antismash
##gg_hex
p_grob <- grobTree(richtext_grob("p<0.0001", x=0.05, y=0.95, hjust=0,
                                 gp=gpar(col="black", fontsize=13, fontface="bold"), 
                                 box_gp = gpar(col = "black"),
                                 padding = unit(c(6, 6, 4, 6), "pt")))
anti_MAC_plot <- ggplot(MAC_35p_gff, aes(ig3, ig5)) +
  geom_hex(bins=40, show.legend = TRUE) +
  scale_fill_continuous(low = "grey80", high = "black") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  geom_point(data=antismash_MAC, aes(x=ig3, y=ig5), 
             color = "#0096FF", size =2.5, alpha = 0.8, inherit.aes = TRUE) +
  theme_minimal()  +
  labs(x ="5' Flanking intergenic region (kbp)", y = "3' Flanking intergenic region (kbp)") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position ="none",  axis.text.x = element_text(color="black", size=25), axis.text.y = element_text(color="black", size=25)) +
  geom_hline(yintercept=1, color = "black", linetype="dashed", alpha=0.8) + 
  geom_vline(xintercept=1, color = "black", linetype="dashed", alpha=0.8) 

#make density plots for ig3 and ig5
n9 <- ggplot(data=MAC_35p_gff, aes(ig3, colour = check_anti)) +
  geom_line(show.legend = FALSE, stat = "density", alpha = 0.5) +
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(name = "Density", breaks = c(0,1.0), limits = c(0,1.25)) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.line.x = element_line("grey"), axis.line.y = element_line("grey"), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(color="black", size=20, face="bold")) +
  scale_color_manual(values=c("darkgrey", "#0096FF"))

n10 <- ggplot(data=MAC_35p_gff, aes(x=ig5, colour = check_anti)) +
  geom_line(show.legend = FALSE, stat = "density", alpha = 0.5) +
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(name = "Density", breaks = c(0,1.0), limits = c(0,1.2)) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(color="black", size=20, face="bold"), axis.line.x = element_line("grey"), axis.line.y = element_line("grey"), axis.ticks.y = element_blank()) +
  scale_color_manual(values=c("darkgrey", "#0096FF")) + coord_flip()


#stack hex and density plots 
grid.newpage()
n9 <- n9 + theme(plot.margin = unit(c(0, 0, 0, 2.7), "lines"))
anti_MAC_plot <- anti_MAC_plot + theme(plot.margin = unit(c(0, 0.5, 0, 0), "lines"))
n10 <- n10 + theme(plot.margin = unit(c(0, 0, 1.2, 0), "lines"))
grid.arrange(n9, anti_MAC_plot, n10,
             ncol = 2,
             widths = c(10,1),
             heights = c(1.2,10),
             layout_matrix = rbind(c(1, NA),
                                   c(2, 3)))

#Antismash descriptive stats

MAC_35p_gff$check_anti <- as.factor(MAC_35p_gff$check_anti)

# Descriptive statistics by group
MAC_35p_gff %>%
  dplyr::select(check_anti, ig3) %>%
  group_by(check_anti) %>%
  summarise(
    n = n(),
    mean = mean(ig3, na.rm = TRUE),
    sd = sd(ig3, na.rm = TRUE),
    stderr = sd/sqrt(n),
    LCL = mean - qt(1 - (0.05 / 2), n - 1) * stderr,
    UCL = mean + qt(1 - (0.05 / 2), n - 1) * stderr,
    median = median(ig3, na.rm = TRUE),
    min = min(ig3, na.rm = TRUE),
    max = max(ig3, na.rm = TRUE),
    IQR = IQR(ig3, na.rm = TRUE),
    LCLmed = MedianCI(ig3, na.rm=TRUE)[2],
    UCLmed = MedianCI(ig3, na.rm=TRUE)[3]
  )

# Mann-Whitney U tests
m1_anti <- wilcox.test(ig3 ~ check_anti, data = MAC_35p_gff,
                       exact = FALSE, conf.int = TRUE)
print(m1_anti)

m2_anti <- wilcox.test(ig5 ~ check_anti, data = MAC_35p_gff,
                       exact = FALSE, conf.int = TRUE)
print(m2_anti)



# NOT SIGNIFICANT
pvals_anti <- c(m1_anti$p.value, m2_anti$p.value)
pvals_adj_anti <- p.adjust(pvals_anti, method = "BH")
pvals_adj_anti
