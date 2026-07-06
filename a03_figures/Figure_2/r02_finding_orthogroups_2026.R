library(dplyr)
library(stringr)

# --- Load input data ---
ortho <- read.delim("d02_Orthogroups.tsv")
ortho_unassigned <- read.delim("d03_Orthogroups_UnassignedGenes.tsv")
ortho_total <- rbind(ortho, ortho_unassigned)

ortho_detail <- read.delim("/Users/mbs625/Downloads/orthofinder_pangenome/fig_2/d04__ortho_detail_2.txt")

# --- Create clean mapping table (genome ID <-> species) ---
ortho_mapping_clean <- ortho_detail %>%
  mutate(
    orthofinder_name = case_when(
      str_detect(accessions, "GCA_") ~ str_extract(accessions, "GCA_\\d+\\.\\d+"),
      str_detect(accessions, "KVL") ~ paste0(str_extract(accessions, "KVL\\d+"), "_proteins"),
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(V1 = orthofinder_name, V2 = species) %>%
  distinct() %>%
  filter(!is.na(V1)) %>%
  mutate(V1 = gsub("^KVL(\\d+)_proteins$", "KVL_\\1", V1),
         V1 = trimws(V1))

# --- Convert orthogroup table to binary presence/absence ---
ortho_binary <- ortho_total
ortho_binary[,-1] <- ifelse(ortho_binary[,-1] != "", 1, 0)

# --- Prepare lookup lists ---
all_genomes  <- intersect(colnames(ortho_binary)[-1], ortho_mapping_clean$V1)
species_list <- unique(ortho_mapping_clean$V2)

sp_cols_list <- lapply(
  setNames(nm = species_list),
  function(sp) intersect(ortho_mapping_clean$V1[ortho_mapping_clean$V2 == sp], all_genomes)
)

og_ids <- ortho_binary$Orthogroup

# --- Define orthogroup categories ---
global_core <- og_ids[rowSums(ortho_binary[,-1]) == ncol(ortho_binary[,-1])]

species_specific_sets <- lapply(sp_cols_list, function(cols){
  if(length(cols) == 0) return(character(0))
  sp_all   <- rowSums(ortho_binary[, cols, drop=FALSE]) == length(cols)
  others   <- setdiff(all_genomes, cols)
  others_0 <- if(length(others)==0) TRUE else rowSums(ortho_binary[, others, drop=FALSE]) == 0
  og_ids[sp_all & others_0]
})

species_core_nonunique_sets <- lapply(names(sp_cols_list), function(sp){
  cols <- sp_cols_list[[sp]]
  if(length(cols) == 0) return(character(0))
  core_sp <- og_ids[rowSums(ortho_binary[, cols, drop=FALSE]) == length(cols)]
  setdiff(core_sp, c(global_core, species_specific_sets[[sp]]))
})
names(species_core_nonunique_sets) <- names(sp_cols_list)

singleton_set <- og_ids[rowSums(ortho_binary[,-1]) == 1]

# species presence matrix
species_presence <- sapply(names(sp_cols_list), function(sp){
  cols <- sp_cols_list[[sp]]
  if(length(cols) == 0) return(rep(FALSE, nrow(ortho_binary)))
  rowSums(ortho_binary[, cols, drop=FALSE]) > 0
})
species_count_per_og <- rowSums(as.matrix(species_presence))

species_exclusive_accessory_sets <- lapply(names(sp_cols_list), function(sp){
  og_ids[
    species_presence[, sp] &
      species_count_per_og == 1 &
      !(og_ids %in% c(global_core,
                      unlist(species_specific_sets, use.names = FALSE),
                      unlist(species_core_nonunique_sets, use.names = FALSE),
                      singleton_set))
  ]
})
names(species_exclusive_accessory_sets) <- names(sp_cols_list)

# --- Build per-genome results ---
results <- data.frame(
  Genome = all_genomes,
  Species = ortho_mapping_clean$V2[match(all_genomes, ortho_mapping_clean$V1)],
  Species_core_nonunique = 0L,
  Species_specific_core  = 0L,
  Other_shared           = 0L,
  Other_shared_species_exclusive = 0L,
  Singletons             = 0L,
  stringsAsFactors = FALSE
)

for(g in all_genomes){
  sp <- results$Species[results$Genome == g]
  present_ogs <- og_ids[ortho_binary[[g]] == 1]
  
  excl <- unique(c(global_core,
                   species_specific_sets[[sp]],
                   species_core_nonunique_sets[[sp]],
                   singleton_set))
  
  results$Species_core_nonunique[results$Genome == g] <- sum(present_ogs %in% species_core_nonunique_sets[[sp]])
  results$Species_specific_core [results$Genome == g] <- sum(present_ogs %in% species_specific_sets[[sp]])
  results$Singletons            [results$Genome == g] <- sum(present_ogs %in% singleton_set)
  
  other_shared_set <- setdiff(present_ogs, excl)
  results$Other_shared[results$Genome == g] <- length(other_shared_set)
  results$Other_shared_species_exclusive[results$Genome == g] <-
    sum(other_shared_set %in% species_exclusive_accessory_sets[[sp]])
}


# 1) Unique count of "Other_shared_species_exclusive" per species
other_shared_excl_unique_per_species <- data.frame(
  Species = names(species_exclusive_accessory_sets),
  other_shared_excl_unique = sapply(
    species_exclusive_accessory_sets,
    function(ogs) length(unique(ogs))
  ),
  row.names = NULL
)

# 2) Total singletons per species (safe to sum per-genome; each singleton belongs to one genome)
singletons_per_species <- results %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(total_singletons = sum(Singletons, na.rm = TRUE), .groups = "drop")

# 3) Species specific core
# Count genomes per species
species_n_genomes <- results %>%
  dplyr::count(Species, name = "n_genomes")

# Extract one Species_specific_core value per species
species_core_per_species <- results %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    species_specific_core = first(Species_specific_core),
    .groups = "drop"
  ) %>%
  left_join(species_n_genomes, by = "Species") %>%
  dplyr::mutate(
    species_specific_core = ifelse(
      n_genomes >= 2,
      species_specific_core,
      0
    )
  ) %>%
  dplyr::select(-n_genomes)

# 4) Final per-species summary + combined total (no double counting)
species_summary <- other_shared_excl_unique_per_species %>%
  left_join(singletons_per_species, by = "Species") %>%
  left_join(species_core_per_species, by = "Species") %>%
  mutate(
    total_combined =
      other_shared_excl_unique +
      total_singletons +
      species_specific_core
  )

species_summary

species_summary <- as.data.frame(species_summary)

# Outputs:
results          # per-genome breakdown
results <- as.data.frame(results)

# Add a per-genome "non-core total" column
results <- results %>%
  mutate(
    Non_core_total = Other_shared + Other_shared_species_exclusive + Singletons
  )


# Select only the columns we want
results_long <- results %>%
  dplyr::select(Genome, Species, Species_core_nonunique, Species_specific_core, Non_core_total) %>%
  pivot_longer(
    cols = c(Species_core_nonunique, Species_specific_core, Non_core_total),
    names_to = "Category",
    values_to = "Count"
  )

# Modify ortho_detail
ortho_detail <- ortho_detail %>%
  mutate(
    orthofinder_name = case_when(
      str_detect(accessions, "GCA_") ~ str_extract(accessions, "GCA_\\d+\\.\\d+"),
      str_detect(accessions, "KVL") ~ paste0(str_extract(accessions, "KVL\\d+"), "_proteins"),
      TRUE ~ NA_character_   # For 'space' or any other non-genome rows
    )
  )
ortho_detail <- ortho_detail %>%
  mutate(
    orthofinder_name = ifelse(
      orthofinder_name == "" | is.na(orthofinder_name),
      accessions,
      orthofinder_name
    )
  )

results_long <- results_long %>%
  dplyr::rename(
    orthofinder_name = Genome,
    species          = Species,
    key              = Category,
    value            = Count
  )


# Modify table externally 
#write.table(results_long, "/Users/mbs625/Downloads/orthofinder_pangenome/orthogroup_info_2.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# Read in modified table
orthogroup_info <- read.delim("d05_orthogroup_info_2.txt")

# Plot the new information
orthogroup_info$fill <- as.factor(orthogroup_info$fill)
orthogroup_info$fill <- factor(orthogroup_info$fill, levels = c("a", "b", "c", "d", "e", "f", "none "))
ortho_levels <- unique(orthogroup_info$orthofinder_name)
orthogroup_info$orthofinder_name <- factor(orthogroup_info$orthofinder_name, levels = ortho_levels)


# Get the name and the y position of each label
nObsType=nlevels(as.factor(orthogroup_info$key))
orthogroup_info$id=rep( seq(1, nrow(orthogroup_info)/nObsType) , each=nObsType)
label_data= orthogroup_info %>% group_by(id, rank) %>% summarize(tot=sum(value))
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5)/number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust<-ifelse(angle < -90, 1.2, 0)
label_data$angle<-ifelse(angle < -90, angle+180, angle)

cols <- c("darkred", "#ff9d9d", "#ff9d9d", "black", "#e3e3e3", "#e3e3e3", "#CECECE")

orthogroup_info$fill <- factor(
  orthogroup_info$fill,
  levels = c("c", "b", "a", "f", "e", "d")
)

test <- ggplot(orthogroup_info) +
  # Make custom panel grid
  geom_hline(
    aes(yintercept = y), 
    data.frame(y = c(110, 1110, 2110, 3110, 4110, 5110)),
    color = "lightgrey", linetype=2) +
  geom_bar(aes(x=orthofinder_name, y=value, fill=fill), stat="identity", alpha=1) +
  scale_fill_manual(values = cols) +
  theme_minimal() +
  theme(
    # Remove axis ticks and text
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(), 
    panel.grid = element_blank(), 
    axis.text.y = element_blank(),
    # Use gray text for the region names
    # Move the legend to the bottom
  ) +
  coord_polar() +
  # Labels are difficult to orient altogether using hjust
  #geom_text(data=label_data, aes(x=id, y=tot+10, label=species, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=4, angle=label_data$angle, inherit.aes = FALSE ) +
  ylim(-4000,max(label_data$tot + 2000 , na.rm=T)) 

test


# Single copy orthogroups M acridum
m_acridum_genomes <- ortho_mapping_clean$V1[
  ortho_mapping_clean$V2 == "M.acridum"
]
unique(ortho_mapping_clean$V2)

copy_counts <- ortho_total

copy_counts[m_acridum_genomes] <-
  lapply(copy_counts[m_acridum_genomes], function(x) {
    ifelse(
      x == "",
      0,
      stringr::str_count(x, ",") + 1
    )
  })

m_acridum_single_copy <- copy_counts$Orthogroup[
  rowSums(copy_counts[, m_acridum_genomes] == 1) ==
    length(m_acridum_genomes)
]


# Total orthogroups for each genome
summary <- data.frame(
  Genome = names(colSums(ortho_binary[,-1])),
  total_orthogroups = as.numeric(colSums(ortho_binary[,-1]))
)
summary$accessory_og <- summary$total_orthogroups - length(global_core)

summary$singletons <- colSums(ortho_binary[og_ids %in% singleton_set, -1, drop = FALSE])

gene_counts <- ortho_total
gene_counts[,-1] <- lapply(gene_counts[,-1], function(x)
  ifelse(x == "", 0, lengths(strsplit(x, " ")))
)

single_copy_per_genome <- colSums(gene_counts[,-1] == 1)

summary$single_copy_orthogroups <- single_copy_per_genome[summary$Genome]

paralog_orthogroups_per_genome <- colSums(gene_counts[,-1] > 1)

summary$orthogroups_with_paralogs <- paralog_orthogroups_per_genome[summary$Genome]

#write.table(summary, "d06_summary_table_s8.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
