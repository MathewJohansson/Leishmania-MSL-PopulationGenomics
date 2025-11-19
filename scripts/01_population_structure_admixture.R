

# ADMIXTURE POPULATION STRUCTURE ANALYSIS
# Leishmania infantum Miltefosine Resistance Study


# Author: Mathew Johansson
# Institution: University of York, Jeffares Lab
# Date: November 2025

# Description:
# This script processes and visualises ADMIXTURE results for three complementary
# datasets to investigate population structure in L. infantum:

#   1. All L. infantum + L. donovani outgroup (K=12)
#      Purpose: Establish species-level separation and major geographic structure

#   2. All L. infantum only (K=13)
#      Purpose: Fine-scale cryptic population structure without outgroup influence

#   3. Americas-only L. infantum (K=13)
#      Purpose: High-resolution structure in New World where MSL deletions predominate

# Key Outputs:
#   - High-confidence population assignments (≥90% ancestry threshold)
#   - Cross-validation error plots justifying optimal K selection
#   - Stacked ancestry barplots showing population composition
#   - Population summary tables with geographic distribution

# Requirements:
#   - ADMIXTURE output files (.Q matrices, .CV files)
#   - Sample metadata with geographic and MSL coverage information

 



# 1. SETUP AND CONFIGURATION ---------------------------------------------------

# 1.1. Load Required Libraries -------------------------------------------------

library(tidyverse)
library(readr)
library(gtools)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(colorspace)
library(forcats)


# 1.2. Define Project Directories ----------------------------------------------

# Set base directory - ONLY CHANGE THIS LINE to run on different machines
BASE_DIR <- "~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Leishmania-MSL-Population-Genomics"

# Define subdirectories
METADATA_DIR <- file.path(BASE_DIR, "data/metadata")
DATA_DIR <- file.path(BASE_DIR, "data/population_structure")
FIGURES_DIR <- file.path(BASE_DIR, "figures/admixture")
RESULTS_DIR <- file.path(BASE_DIR, "results/population_structure")

# Create output directories if they don't exist
dir.create(file.path(FIGURES_DIR, "All_Linfantum_With_Outgroup"), 
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIGURES_DIR, "Linfantum_Only"), 
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIGURES_DIR, "Americas_Only"), 
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULTS_DIR, "High_Confidence_Samples"), 
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULTS_DIR, "Population_Assignments"), 
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULTS_DIR, "Population_Labels"), 
           recursive = TRUE, showWarnings = FALSE)


# 1.3. Load Metadata -----------------------------------------------------------

metadata <- read_tsv(file.path(METADATA_DIR, "vcf_samples29_location_and_MSL_data_2025-06-19_corrected.tsv"))


# 1.4. Define Colour Palettes --------------------------------------------------

# Master palette for All L. infantum + Outgroup (K=12)
master_palette_ordered <- setNames(
  c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
    "#A65628", "#F781BF", "#999999", "#66C2A5", "#B8860B", "#8DA0CB"),
  paste0("V", 1:12)
)

# Alphabetically sorted version for legends
master_palette_alpha <- master_palette_ordered[sort(names(master_palette_ordered))]

# Generic master palette (using same colours)
master_palette <- master_palette_ordered

# Population palette for L. infantum only (K=13) 
population_palette_base <- setNames(
  c(brewer.pal(12, "Paired"), "#8dd3c7"),
  c("Brazil, Honduras & Colombia Major Cluster", 
    "Brazil, Honduras & Panama Major Cluster",  
    "North Africa & Europe Mixed Cluster (Tunisia, France, Israel)", 
    "Western Europe & North Africa Cluster (France, Portugal, Tunisia)", 
    "Americas Minor Mixed Cluster (Honduras, Panama, France)", 
    "Southern Cone Major Cluster (Brazil, Uruguay)", 
    "Europe & Americas Mixed Cluster (Spain, France, Brazil, Morocco)", 
    "Middle East Cluster (Israel, Yemen)", 
    "Brazil Minor Population Cluster",  
    "Brazil-Paraguay Major Cluster", 
    "Southern South America Cluster (Brazil, Uruguay, Paraguay)",
    "Brazil, Honduras & Panama Major Cluster 2",
    "Central & East Asia with Middle East Mix (Uzbekistan, China, Israel)")
)

# Americas palette (K=13) - bright colours to match other datasets
americas_palette_base <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
  "#A65628", "#F781BF", "#66C2A5", "#FC8D62", "#8DA0CB", "#B15928", "#FDBF6F"
)


cat("\n=== SETUP COMPLETE ===\n")
cat("Base directory:", BASE_DIR, "\n")
cat("Metadata loaded:", nrow(metadata), "samples\n\n")




# 2. All L.infantum + Outgroup -------------------------------------------------


# 2.1. Confidence Filtering and Plotting Clusters ------------------------------

# Load the K=12 ADMIXTURE results for all L.infantum samples + outgroup
all_admix_data <- read_tsv(file.path(DATA_DIR, "population_membership_table.allsamples.k12.tsv"))

# Join metadata
all_admix_data <- all_admix_data %>%
  left_join(metadata, by = c("sample" = "sra_run_accession"))

# Filter by ≥90% (high-confidence) assignment
high_conf_all <- all_admix_data %>%
  filter(population_proportion >= 0.9) %>%
  mutate(population = as.numeric(sub("V", "", population)))

# Clean Host_Species column
high_conf_all <- high_conf_all %>%
  mutate(host_species_clean = case_when(
    grepl("Canis familiaris|Dog|Nyctomys procyonoides", host_species, ignore.case = TRUE) ~ "Dog",
    grepl("Homo sapiens|Human", host_species, ignore.case = TRUE) ~ "Human",
    TRUE ~ NA_character_
  ))


# Define population labels
all_linfantum_plus_outgroup_pop_labels <- tibble(
  population = paste0("V", 1:12),
  population_label = c(
    "Mediterranean & Middle East Cluster (Albania, Egypt, Greece, Italy, Turkey, Yemen)",
    "Brazil Minor Population Cluster",
    "Central & East Asia with Middle East Mix (Uzbekistan, China, Israel)",
    "Western Europe & North Africa Mixed Cluster (France, Portugal, Tunisia)",
    "Europe & Americas Mixed Minor Cluster (Spain, France, Honduras, Morocco, Panama)",
    "Brazil Major Population Cluster A",
    "Brazil Major Population Cluster B",
    "North Africa & Middle East Mixed Cluster (Tunisia, Israel)",
    "Brazil Minor Mixed Mix Cluster",
    "Brazil Major Population Cluster C",
    "Southern South America Cluster (Brazil, Uruguay, Paraguay)",
    "East African Cluster (Ethiopia, Sudan)"
  )
)


# Summarise per population
table_data_clean_numeric_all_samples <- high_conf_all %>%
  group_by(population) %>%
  summarise(
    N_Samples = n(),
    Countries_Represented = paste(unique(source_country), collapse = ", "),
    Host_Species = paste(unique(na.omit(host_species_clean)), collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(population = paste0("V", population)) %>%
  left_join(all_linfantum_plus_outgroup_pop_labels, by = "population") %>%
  mutate(Notes = "") %>%
  select(population, population_label, N_Samples, Countries_Represented, Host_Species, Notes)

table_data_clean_numeric_all_samples

# Save population table with cols: population, population_label (geo region), N_Samples, Countries_Represented, Host_Species, Notes
write_tsv(table_data_clean_numeric_all_samples,
          file.path(RESULTS_DIR, "all_linfantum_plus_outgroup_population_summary_k12.tsv")
)


# Join population labels back to high_conf_all for plotting
high_conf_all <- high_conf_all %>%
  mutate(population = paste0("V", population)) %>% 
  select(-starts_with("population_label")) %>%
  left_join(all_linfantum_plus_outgroup_pop_labels, by = "population")


# Save high-confidence samples
write_tsv(high_conf_all, 
          file.path(RESULTS_DIR, "High_Confidence_Samples", 
                    "all_linfantum_plus_outgroup_high_confidence_samples_k12.tsv"))


# Ensure factor levels are alphabetically sorted for legend
all_pop_order <- sort(all_linfantum_plus_outgroup_pop_labels$population_label)
high_conf_all <- high_conf_all %>%
  mutate(population_label = factor(population_label, levels = all_pop_order))


# Create palette properly mapped to alphabetically sorted labels
# First get the V codes in alphabetical order of their labels
sorted_mapping <- all_linfantum_plus_outgroup_pop_labels %>%
  arrange(population_label)

# Create palette with labels as names, in alphabetical order
master_palette_for_plot <- setNames(
  master_palette_ordered[sorted_mapping$population],
  sorted_mapping$population_label
)


# Assignment barplot with forced legend order
ggplot(high_conf_all, aes(x = population_label, fill = population_label)) +
  geom_bar() +
  scale_fill_manual(values = master_palette_for_plot, drop = FALSE) +
  theme_minimal() +
  labs(
    title = "ADMIXTURE-Inferred Population Assignments (K=12)\nL.infantum Samples with L.donovani Outgroup",
    x = NULL,
    y = "Number of Samples",
    fill = "Population"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


# Save the plot as a jpeg
ggsave(file.path(FIGURES_DIR, "All_Linfantum_With_Outgroup", 
                 "all_linfantum_with_outgroup_population_assignments_k12.jpeg"), 
       width = 10, height = 6)

# Save the plot as a pdf
ggsave(file.path(FIGURES_DIR, "All_Linfantum_With_Outgroup", 
                 "all_linfantum_with_outgroup_population_assignments_k12.pdf"), 
       width = 10, height = 6)


# Save population labels lookup
write_tsv(
  all_linfantum_plus_outgroup_pop_labels, 
  file.path(RESULTS_DIR, "Population_Labels", 
            "all_linfantum_plus_outgroup_population_labels_lookup.tsv")
)


# The above has: 

# 1. Read the ADMIXTURE K=12 population table for L.infantum + outgroup.
# 2. Filtered samples with ≥90% confidence ancestry to a single population.
# 3. Saved the samples to this file: all_linfantum_plus_outgroup_high_confidence_samples_k12.tsv
# 4. Counts how many samples are in each confident population. 
# 5. Generates and saves a bar plot showing how many samples are in each ADMIXTURE-defined population. Saved as all_linfantum_with_outgroup_population_assignments_k12.jpeg/pdf.




# 2.2. K=12 Justification ------------------------------------------------------

# Despite K=12 being chosen and the ADMIXTURE runs already completed before this analysis began,
#    justification of K being 12 must still be stated. 

# The following code loads the cross-validation data for K=12 and plots the cross-validation error curve to justify the choice of K. 
# The lowest CV-Error value on the plot indicates the optimal number of populations (K). 


# Load the cross-validation data for K=12
cv_data_linfantum_with_outgroup <- read_table(file.path(DATA_DIR, "Linfantum_samples_with_outgroup_2025-03-12.CV.txt"),
                                              col_names = c("K", "CV_Error"))
cv_data_linfantum_with_outgroup


# Plot the cross-validation error curve
ggplot(cv_data_linfantum_with_outgroup,
       aes(x = K,
           y = CV_Error)) +
  geom_line(color = "steelblue",
            size = 1.2) + 
  geom_point(color = "darkred",
             size = 2) +
  theme_minimal() + 
  labs(title = "ADMIXTURE Cross-Validation Error vs. K",
       subtitle = "L.infantum Samples with L.donovani Outgroup",
       x = "Number of Populations (K)",
       y = "Cross-Validation Error") +
  geom_vline(xintercept = 12,
             linetype = "dashed",
             color = "darkgreen") +
  annotate("text",
           x = 12.5,
           y = min(cv_data_linfantum_with_outgroup$CV_Error),
           label = "K = 12",
           color = "darkgreen",
           hjust = 0)

# Save the plot as a jpeg
ggsave(file.path(FIGURES_DIR, "All_Linfantum_With_Outgroup", 
                 "linfantum_samples_with_outgroup_cv_error_plot_k1_to_k20.jpeg"), 
       width = 10, height = 6)

# Save the plot as a pdf
ggsave(file.path(FIGURES_DIR, "All_Linfantum_With_Outgroup", 
                 "linfantum_samples_with_outgroup_cv_error_plot_k1_to_k20.pdf"), 
       width = 10, height = 6)




# 2.3. Stacked Ancestry Bar Plot -----------------------------------------------

# These allow for visualisation of each sample's proportional ancestry components from ADMIXTURE for selected K values.


# Read Q-matrix
q_file <- file.path(DATA_DIR, "Linfantum_samples_with_outgroup_2025-03-12.q30.missing0.8.sampMissing0.25.biallelic.depth.annotate.pruned0.5_renamed.12.Q")
admix_wide <- read_table(q_file, col_names = FALSE)

# Add sample IDs (ensure order matches metadata)
admix_wide <- admix_wide %>%
  mutate(sample = metadata$sra_run_accession) %>%
  relocate(sample)

# Name columns as V1..V12
colnames(admix_wide)[-1] <- paste0("V", 1:(ncol(admix_wide) - 1))

# Pivot to long format
admix_long <- admix_wide %>%
  pivot_longer(cols = -sample,
               names_to = "population",
               values_to = "population_proportion")

# Convert population to factor with correct ordering
ordered_pops <- paste0("V", 1:12)
admix_long <- admix_long %>%
  mutate(population = factor(population, levels = ordered_pops))

# Verify the factor levels are correct
print("Population factor levels:")
print(levels(admix_long$population))

# Join population labels (optional, can be useful for labels or ordering)
admix_long <- admix_long %>%
  left_join(all_linfantum_plus_outgroup_pop_labels, by = "population")

# Assign each sample to its dominant population (for ordering)
dominant_pop <- admix_long %>%
  group_by(sample) %>%
  slice_max(population_proportion, n = 1, with_ties = FALSE) %>%
  select(sample, dominant_population = population_label)

admix_long <- admix_long %>%
  left_join(dominant_pop, by = "sample")

# Reorder samples: group by dominant population, then by max ancestry proportion
sample_order <- admix_long %>%
  group_by(sample) %>%
  summarise(dominant_population = first(dominant_population),
            max_ancestry = max(population_proportion)) %>%
  arrange(dominant_population, desc(max_ancestry)) %>%
  pull(sample)

admix_long <- admix_long %>%
  mutate(sample = factor(sample, levels = sample_order))

# Plot using properly mapped palette
ggplot(admix_long, aes(x = sample, y = population_proportion, fill = population_label)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = master_palette_for_plot, drop = FALSE) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  labs(title = "All L.infantum With Outgroup ADMIXTURE Ancestry Proportions (K=12)",
       x = "Samples ordered by population",
       y = "Ancestry Proportion",
       fill = "Population")


# Save the barplot as a jpeg 
ggsave(file.path(FIGURES_DIR, "All_Linfantum_With_Outgroup", 
                 "all_linfantum_with_outgroup_ancestry_barplot.jpeg"), 
       width = 10, height = 6)

# Save the barplot as a pdf 
ggsave(file.path(FIGURES_DIR, "All_Linfantum_With_Outgroup", 
                 "all_linfantum_with_outgroup_ancestry_barplot.pdf"), 
       width = 10, height = 6)


# The above code has:

# 1. Loaded ADMIXTURE Q‑matrix (12 populations, K=12) from the .Q file for the All L.infantum with Outgroup data.
# 2. Added sample IDs by matching Q‑matrix rows to the sample metadata.
# 3. Pivoted data to long format so each row = (sample, population, ancestry proportion).
# 4. Joined population labels (human‑readable cluster names) from your lookup table.
# 5. Determined plotting order for samples by their dominant ancestry cluster.
# 6. Plotted a stacked ancestry barplot.






# 3. All L.infantum Only -------------------------------------------------------



# 3.1. Confidence Filtering and Plotting Clusters ------------------------------

# K=13 was determined to be the optimal number of populations for L.infantum only.


# Define population order (consistent across plots)
population_order <- c(
  "Brazil, Honduras & Colombia Major Cluster", 
  "Brazil, Honduras & Panama Major Cluster",  
  "North Africa & Europe Mixed Cluster (Tunisia, France, Israel)", 
  "Western Europe & North Africa Cluster (France, Portugal, Tunisia)", 
  "Americas Minor Mixed Cluster (Honduras, Panama, France)", 
  "Southern Cone Major Cluster (Brazil, Uruguay)", 
  "Europe & Americas Mixed Cluster (Spain, France, Brazil, Morocco)", 
  "Middle East Cluster (Israel, Yemen)", 
  "Brazil Minor Population Cluster",  
  "Brazil-Paraguay Major Cluster", 
  "Southern South America Cluster (Brazil, Uruguay, Paraguay)",
  "Brazil, Honduras & Panama Major Cluster 2",
  "Central & East Asia with Middle East Mix (Uzbekistan, China, Israel)" 
)

# Apply the same palette in the same order
population_palette <- setNames(
  c(brewer.pal(12, "Paired"), "#8dd3c7"),
  population_order
)

# Alphabetical order for legend
pop_order_alpha <- sort(population_order)

# Reorder the palette to match alphabetical order
population_palette_alpha <- population_palette[pop_order_alpha]



# Load the K=13 ADMIXTURE results for all L.infantum only
linfantum_k13_admix_data <- read_tsv(file.path(DATA_DIR, "population_membership_table.infantum_only.k13.tsv"))

# Filter by ≥90% high-confidence and join metadata
high_conf_linfantum_k13 <- linfantum_k13_admix_data %>%
  filter(population_proportion >= 0.9) %>%
  left_join(metadata, by = c("sample" = "sra_run_accession")) %>%
  mutate(population = as.numeric(sub("V", "", population)))

# Summarise per population
high_conf_linfantum_k13_summary <- high_conf_linfantum_k13 %>%
  group_by(population) %>%
  summarise(
    N_Samples = n(),
    Countries_Represented = paste(unique(source_country), collapse = ", "),
    Host_Species = paste(unique(host_species), collapse = ", "),
    .groups = "drop"
  )
high_conf_linfantum_k13_summary

# Save high-confidence samples to a new file
write_tsv(high_conf_linfantum_k13, 
          file.path(RESULTS_DIR, "High_Confidence_Samples", 
                    "linfantum_only_high_confidence_samples_k13.tsv"))


# Clean Host_Species column
high_conf_linfantum_k13_clean <- high_conf_linfantum_k13 %>%
  mutate(host_species_clean = case_when(
    grepl("Canis familiaris|Dog", host_species, ignore.case = TRUE) ~ "Dog",
    grepl("Homo sapiens|Human", host_species, ignore.case = TRUE) ~ "Human",
    grepl("Nyctomys procyonoides", host_species, ignore.case = TRUE) ~ "Dog", 
    TRUE ~ NA_character_
  ))


# Create a tibble mapping population codes (V1–V13) to meaningful labels
linfantum_only_pop_labels <- tibble(
  population = 1:13,  
  population_label = c(
    "Brazil, Honduras & Colombia Major Cluster", 
    "Brazil, Honduras & Panama Major Cluster",  
    "North Africa & Europe Mixed Cluster (Tunisia, France, Israel)", 
    "Western Europe & North Africa Cluster (France, Portugal, Tunisia)", 
    "Americas Minor Mixed Cluster (Honduras, Panama, France)", 
    "Southern Cone Major Cluster (Brazil, Uruguay)", 
    "Europe & Americas Mixed Cluster (Spain, France, Brazil, Morocco)", 
    "Middle East Cluster (Israel, Yemen)", 
    "Brazil Minor Population Cluster",  
    "Brazil-Paraguay Major Cluster", 
    "Southern South America Cluster (Brazil, Uruguay, Paraguay)",
    "Brazil, Honduras & Panama Major Cluster 2",
    "Central & East Asia with Middle East Mix (Uzbekistan, China, Israel)" 
  )
)


# Create summary table 
table_data_clean_numeric_all_linfantum_k13 <- high_conf_linfantum_k13_clean %>%
  group_by(population) %>%
  summarise(
    N_Samples = n(),
    Countries_Represented = paste(unique(source_country), collapse = ", "),
    Host_Species = paste(unique(na.omit(host_species_clean)), collapse = ", "),
    .groups = "drop"
  ) %>%
  left_join(linfantum_only_pop_labels, by = "population") %>% 
  mutate(Notes = "") %>%
  select(population, population_label, N_Samples, Countries_Represented, Host_Species, Notes)
table_data_clean_numeric_all_linfantum_k13


# Save population table with cols: population, population_label (geo region), N_Samples, Countries_Represented, Host_Species, Notes
write_tsv(table_data_clean_numeric_all_linfantum_k13,
          file.path(RESULTS_DIR, "linfantum_only_population_summary_k13.tsv")
)


# Alphabetical order for legend
pop_order_alpha <- sort(linfantum_only_pop_labels$population_label)

# Join labels to high_conf data
high_conf_linfantum_k13 <- high_conf_linfantum_k13 %>%
  left_join(linfantum_only_pop_labels, by = "population") %>%
  mutate(population_label = factor(population_label, levels = pop_order_alpha))


ggplot(high_conf_linfantum_k13, aes(x = population_label, fill = population_label)) +
  geom_bar() +
  scale_fill_manual(values = population_palette_alpha, drop = FALSE) +
  guides(fill = guide_legend(reverse = FALSE)) +
  theme_minimal() +
  labs(
    title = "ADMIXTURE-Inferred Population Assignments (K=13)\nL.infantum Samples Only",
    x = NULL, 
    y = "Number of Samples",
    fill = "Population"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )



# Save the plot as a jpeg
ggsave(file.path(FIGURES_DIR, "Linfantum_Only", 
                 "linfantum_only_population_assignments_k13.jpeg"), 
       width = 10, height = 6)

# Save the plot as a pdf
ggsave(file.path(FIGURES_DIR, "Linfantum_Only", 
                 "linfantum_only_population_assignments_k13.pdf"), 
       width = 10, height = 6)



# The above code has: 

# 1. Read the ADMIXTURE K=13 population table for L.infantum only.
# 2. Filtered samples with ≥90% confidence ancestry to a single population.
# 3. Saved the samples to: linfantum_only_high_confidence_samples_k13.tsv 
# 4. Counts how many samples are in each confident population. 
# 5. Generates and saves bar plot showing how many samples are in each ADMIXTURE-defined population. Saved as linfantum_only_population_assignments_k13.jpeg/pdf.




# 3.2. K=13 Justification ------------------------------------------------------

# Load the cross-validation data
cv_data_linfantum_only <- read_table(file.path(DATA_DIR, "infantum_only_n463.CV.txt"),
                                     col_names = c("K", "CV_Error"))
cv_data_linfantum_only

# Plot the cross-validation error curve highlighting K=13
ggplot(cv_data_linfantum_only,
       aes(x = K, y = CV_Error)) +
  geom_line(color = "steelblue", size = 1.2) + 
  geom_point(color = "darkred", size = 2) +
  theme_minimal() + 
  labs(title = "ADMIXTURE Cross-Validation Error vs. K",
       subtitle = "All L.infantum Only",
       x = "Number of Populations (K)",
       y = "Cross-Validation Error") +
  geom_vline(xintercept = 13, linetype = "dashed", color = "purple") +
  annotate("text",
           x = 13.5,
           y = min(cv_data_linfantum_only$CV_Error),
           label = "K = 13",
           color = "purple",
           hjust = 0)

# Save the plot as a jpeg
ggsave(file.path(FIGURES_DIR, "Linfantum_Only", 
                 "linfantum_only_cv_error_plot_k1_to_k20_for_k13_only.jpeg"), 
       width = 10, height = 6)

# Save the plot as a pdf
ggsave(file.path(FIGURES_DIR, "Linfantum_Only", 
                 "linfantum_only_cv_error_plot_k1_to_k20_for_k13_only.pdf"), 
       width = 10, height = 6)




# K=13 was chosen as the optimal number of populations for L.infantum only. 




# 3.3. Define Populations -----------------------------------------------------

# Load the high-confidence assignment table for K=13
high_conf_linfantum_k13


# Summarise metadata by population

# 1. Number of samples per population, descending order
high_conf_linfantum_k13 %>%
  group_by(population) %>%
  summarise(num_samples = n()) %>%
  arrange(desc(num_samples))


# 2. Geographic distribution

# Join geographic region info to high-confidence population assignments
high_conf_linfantum_k13 <- high_conf_linfantum_k13 %>%
  left_join(metadata %>% select(sra_run_accession, source_geographic_region),
            by = c("sample" = "sra_run_accession"))

# Remove duplicated columns from previous joins
high_conf_linfantum_k13 <- high_conf_linfantum_k13 %>%
  mutate(source_geographic_region = coalesce(source_geographic_region.x, source_geographic_region.y)) %>%
  select(-ends_with(".x"), -ends_with(".y"))

# Summarise the number of samples per population and region
high_conf_linfantum_k13 %>%
  group_by(population, source_geographic_region) %>%
  summarise(num_samples = n(), .groups = "drop") %>%
  arrange(desc(num_samples))


# Save the lookup table for future reference
write_tsv(
  linfantum_only_pop_labels, 
  file.path(RESULTS_DIR, "Population_Labels", 
            "linfantum_only_population_labels_lookup.tsv")
)


# Merge descriptive population labels into the high-confidence assignment table
high_conf_linfantum_k13 <- high_conf_linfantum_k13 %>%
  select(-starts_with("population_label")) %>%  
  left_join(linfantum_only_pop_labels, by = "population")

# Final output table: sample ID, population code, and descriptive label
final_population_assignments_linfantum_only <- high_conf_linfantum_k13 %>%
  select(sample, population, population_label)
final_population_assignments_linfantum_only

# Save the final assignments
write_tsv(
  final_population_assignments_linfantum_only, 
  file.path(RESULTS_DIR, "Defined_Populations", 
            "final_population_assignments_linfantum_only.tsv")
)



# The above has: 

# 1. Inspected the high-confidence population assignment table.
# 2. Summarised sample counts per population.
# 3. Integrated geographic metadata into the assignment table.
# 4. Cleaned up duplicated columns created during joins.
# 5. Summarised sample counts per population and geographic region.
# 6. Created and saved a lookup table mapping population codes to descriptive labels.
# 7. Merged high-confidence assignments with descriptive labels to generate a final table 
#    of strain IDs with population group and population label, then saved it for downstream use.




# 3.4. Stacked Ancestry Bar Plot -----------------------------------------------


# Create lookup for ancestry components -> geographic names
ancestry_labels <- tibble(
  ancestry_component = paste0("X", 1:13),
  population_label = linfantum_only_pop_labels$population_label
)


# Load the Q matrix for K=13
Q_file_k13 <- file.path(DATA_DIR, "infantum_only_n463.q30.missing0.8.sampMissing0.25.biallelic.depth.annotate.pruned0.5_renamed.13.Q")
Qmat_k13 <- read_table2(Q_file_k13, col_names = FALSE)

# Load membership table with population assignments for K=13
pop_membership_k13 <- read_tsv(file.path(DATA_DIR, "population_membership_table.infantum_only.k13.tsv"))


# Make sure sample counts match
stopifnot(nrow(Qmat_k13) == nrow(pop_membership_k13))

# Combine sample names and population info with Q matrix
Qmat_k13_labeled <- Qmat_k13 %>%
  mutate(sample = pop_membership_k13$sample,
         population = pop_membership_k13$population) %>%
  relocate(sample, population)

# Order samples by population, then by highest ancestry proportion within population
Qmat_k13_labeled <- Qmat_k13_labeled %>%
  rowwise() %>%
  mutate(max_ancestry = max(c_across(starts_with("X")))) %>%
  ungroup() %>%
  arrange(population, desc(max_ancestry))

# Create a factor for sample ordered by population + ancestry (this orders x-axis)
Qmat_k13_labeled$sample <- factor(Qmat_k13_labeled$sample, levels = Qmat_k13_labeled$sample)

# Convert to long format for plotting
Qmat_k13_long <- Qmat_k13_labeled %>%
  pivot_longer(cols = starts_with("X"),
               names_to = "ancestry_component",
               values_to = "proportion") %>%
  mutate(ancestry_component = factor(ancestry_component, levels = paste0("X", 1:13)))

Qmat_k13_long <- Qmat_k13_long %>%
  left_join(ancestry_labels, by = "ancestry_component") %>%
  mutate(population_label = factor(population_label, 
                                   levels = linfantum_only_pop_labels$population_label))



Qmat_k13_long <- Qmat_k13_long %>%
  mutate(population_label = factor(population_label, levels = pop_order_alpha))

# Use the alphabetical palette for the fill
ggplot(Qmat_k13_long, 
       aes(x = sample, y = proportion, fill = population_label)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = population_palette_alpha, drop = FALSE) +
  theme_minimal() +
  labs(title = "L. infantum ADMIXTURE Ancestry Proportions (K=13)",
       x = "Samples ordered by population",
       y = "Ancestry Proportion", 
       fill = "Population") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


# Save the plot as a jpeg
ggsave(file.path(FIGURES_DIR, "Linfantum_Only", 
                 "linfantum_only_stacked_ancestry_barplot_k13.jpeg"), 
       width = 15, height = 8)

# Save the plot as a pdf
ggsave(file.path(FIGURES_DIR, "Linfantum_Only", 
                 "linfantum_only_stacked_ancestry_barplot_k13.pdf"), 
       width = 15, height = 8)


# The above code has:

# 1. Loaded ADMIXTURE Q‑matrix for K=13 for the All L.infantum Only data. 
# 2. Added sample IDs by matching Q‑matrix rows to the sample metadata.
# 3. Pivoted data to long format so each row = (sample, population, ancestry proportion).
# 4. Joined population labels (human‑readable cluster names) from the lookup table.
# 5. Determined plotting order for samples by their dominant ancestry cluster.
# 6. Plotted a stacked ancestry barplot.






# 4. Americas Only -------------------------------------------------------------

# 4.1. Confidence Filtering and Plotting Clusters ------------------------------

# Load the K=13 ADMIXTURE results for Americas Only
americas_only_k13_admix_data <- read_tsv(file.path(DATA_DIR, "population_membership_table.Americas_only.k13.tsv"))

# Convert population to numeric
americas_only_k13_admix_data <- americas_only_k13_admix_data %>%
  mutate(population = as.numeric(sub("V", "", population)))

# Filter by ≥90% high-confidence and join metadata
high_conf_americas_k13 <- americas_only_k13_admix_data %>%
  filter(population_proportion >= 0.9) %>%
  left_join(
    metadata %>% select(sra_run_accession, source_country, source_state_or_region, host_species),
    by = c("sample" = "sra_run_accession")
  )

# Clean host species
high_conf_americas_k13_clean <- high_conf_americas_k13 %>%
  mutate(host_species_clean = case_when(
    grepl("Canis familiaris|Dog", host_species, ignore.case = TRUE) ~ "Dog",
    grepl("Homo sapiens|Human", host_species, ignore.case = TRUE) ~ "Human",
    grepl("Nyctomys procyonoides", host_species, ignore.case = TRUE) ~ "Dog", 
    TRUE ~ NA_character_
  ))

# Create accurate labels based on actual geographic distribution from your data
americas_only_pop_labels <- tibble(
  population = 1:13,
  population_label = c(
    "Brazil Southeast (Espírito Santo, Minas Gerais, Piauí)",
    "Brazil Northeast A",
    "Brazil Mixed Central-South (Multi-state)",
    "Brazil Central-West A",
    "Brazil Central-West B",
    "Brazil Northeast-Central (Multi-state)",
    "Brazil Northeast B",
    "Brazil North-Northeast (Maranhão, Pará, Piauí)",
    "Brazil South (Santa Catarina)",
    "Brazil Mixed (Mato Grosso, Rio Grande do Norte, São Paulo)",
    "Brazil-Paraguay-Uruguay Cluster",
    "Brazil-Honduras Mixed (Rio de Janeiro + Honduras)",
    "Central America (Honduras, Panama)"
  )
)

# Join labels to high-confidence data
high_conf_americas_k13_clean <- high_conf_americas_k13_clean %>%
  left_join(americas_only_pop_labels, by = "population") 

# Create population summary table
table_data_clean_numeric_americas_only <- high_conf_americas_k13_clean %>%
  group_by(population, population_label) %>%
  summarise(
    N_Samples = n(),
    Countries_Represented = paste(sort(unique(source_country)), collapse = ", "),
    Regions_Represented   = paste(sort(unique(na.omit(source_state_or_region))), collapse = ", "),
    Host_Species = paste(sort(unique(na.omit(host_species_clean))), collapse = ", "),
    Notes = "",
    .groups = "drop"
  ) %>%
  arrange(population_label) 

# Save population table
write_tsv(
  table_data_clean_numeric_americas_only,
  file.path(RESULTS_DIR, "americas_only_population_summary_k13.tsv")
)

# Force population labels into alphabetical order
high_conf_americas_k13_clean$population_label <- factor(
  high_conf_americas_k13_clean$population_label,
  levels = sort(unique(high_conf_americas_k13_clean$population_label))
)

# Define 13-colour palette and match to population labels
americas_palette <- americas_palette_base
names(americas_palette) <- sort(americas_only_pop_labels$population_label)

# Create barplot (x-axis labels hidden, legend/key visible)
ggplot(high_conf_americas_k13_clean, aes(x = population_label, fill = population_label)) +
  geom_bar() +
  scale_fill_manual(values = americas_palette) +
  labs(
    title = "ADMIXTURE-Inferred Population Assignments (K=13)\nAmericas-Only L. infantum Samples",
    x = NULL, 
    y = "Number of Samples",
    fill = "Population"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

# Save the plot as a jpeg
ggsave(file.path(FIGURES_DIR, "Americas_Only", 
                 "americas_only_linfantum_population_assignments_k13.jpeg"), 
       width = 10, height = 6)

# Save the plot as a pdf
ggsave(file.path(FIGURES_DIR, "Americas_Only", 
                 "americas_only_linfantum_population_assignments_k13.pdf"), 
       width = 10, height = 6)



# The above code has: 

# 1. Read the ADMIXTURE K=13 population table for the Americas Only L.infantum.
# 2. Filtered samples with ≥90% confidence ancestry to a single population.
# 3. Saved the samples to this file: americas_only_high_confidence_samples_k13.tsv
# 4. Counts how many samples are in each confident population. 
# 5. Generates and saves a bar plot showing how many samples are in each ADMIXTURE-defined population. Saved as americas_only_linfantum_population_assignments_k13.jpeg/pdf.




# 4.2. K=13 Justification ------------------------------------------------------

# Despite K=13 being chosen and the ADMIXTURE runs already completed before this analysis began,
#    justification of K being 13 must still be stated. 

# Load the cross-validation data for K=13
cv_data_americas_only <- read_table(file.path(DATA_DIR, "Americas_samples_n388.CV.txt"),
                                    col_names = c("K", "CV_Error"))
cv_data_americas_only

# Plot the cross-validation error curve
ggplot(cv_data_americas_only,
       aes(x = K,
           y = CV_Error)) +
  geom_line(color = "steelblue",
            size = 1.2) + 
  geom_point(color = "darkred",
             size = 2) +
  theme_minimal() + 
  labs(title = "ADMIXTURE Cross-Validation Error vs. K",
       subtitle = "L.infantum Samples, Americas Only",
       x = "Number of Populations (K)",
       y = "Cross-Validation Error") +
  geom_vline(xintercept = 13,
             linetype = "dashed",
             color = "darkgreen") +
  annotate("text",
           x = 13.5,
           y = min(cv_data_americas_only$CV_Error),
           label = "K = 13",
           color = "darkgreen",
           hjust = 0)

# Save the plot as a jpeg
ggsave(file.path(FIGURES_DIR, "Americas_Only", 
                 "americas_only_cv_error_plot_k1_to_k20.jpeg"), 
       width = 10, height = 6)

# Save the plot as a pdf
ggsave(file.path(FIGURES_DIR, "Americas_Only", 
                 "americas_only_cv_error_plot_k1_to_k20.pdf"), 
       width = 10, height = 6)




# 4.3. Define Populations ------------------------------------------------------

# Load the high-confidence assignment table for K=13
americas_only_k13_admix_data


# Summarise metadata by population

# 1. Number of samples per population, descending order
americas_only_k13_admix_data %>%
  group_by(population) %>%
  summarise(num_samples = n()) %>%
  arrange(desc(num_samples))


# 2. Geographic distribution

# Remove any existing source_country columns before join to avoid suffixes
high_conf_americas_k13 <- high_conf_americas_k13 %>%
  select(-starts_with("source_country"))

# Join metadata's source_country to high_conf_americas_k13 by sample ID
high_conf_americas_k13 <- high_conf_americas_k13 %>%
  left_join(
    select(metadata, sra_run_accession, source_country),
    by = c("sample" = "sra_run_accession")
  )

# Now summarise by population and source_country
high_conf_americas_k13 %>%
  group_by(population, source_country) %>%
  summarise(num_samples = n(), .groups = "drop") %>%
  arrange(population, desc(num_samples))


# Save the lookup table for future reference
write_tsv(
  americas_only_pop_labels, 
  file.path(RESULTS_DIR, "Population_Labels", 
            "americas_only_population_labels_lookup.tsv")
)


# Merge descriptive population labels into the high-confidence assignment table
high_conf_americas_k13 <- high_conf_americas_k13 %>%
  select(-starts_with("population_label")) %>% 
  left_join(americas_only_pop_labels, by = "population")

# Final output table: sample ID, population code, and descriptive label
final_population_assignments_americas_only <- high_conf_americas_k13 %>%
  select(sample, population, population_label)
final_population_assignments_americas_only

# Save the final assignments
write_tsv(final_population_assignments_americas_only, 
          file.path(RESULTS_DIR, "Defined_Populations", 
                    "final_population_assignments_americas_only.tsv")
)



# The above has: 

# 1. Inspected the high-confidence assignment table.
# 2. Summarised the number of samples per population, by geographic distribution, and by sample name.
# 3. Assigned descriptive labels to population codes.
# 4. Created and saved a lookup table for future reference.
# 5. Merged the ADMIXTURE outputs with the population labels to get full sample metadata and population label.




# 4.4. Stacked Ancestry Bar Plot -----------------------------------------------

# Read Q-matrix for Americas only K=13
q_file_americas <- file.path(DATA_DIR, "Americas_samples_n388.q30.missing0.8.sampMissing0.25.biallelic.depth.annotate.pruned0.5_renamed.13.Q")
admix_wide_americas <- read_table(q_file_americas, col_names = FALSE)

# Use sample vector from population membership table
samples_americas <- americas_only_k13_admix_data$sample

# Add sample names to the Q-matrix
admix_wide_americas$sample <- samples_americas

# Create column names for populations (V1 to V13)
colnames(admix_wide_americas)[1:13] <- paste0("V", 1:13)

# Convert wide format to long format
admix_long_americas <- admix_wide_americas %>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "population", 
               values_to = "population_proportion") %>%
  # Convert population factor to numeric for the join
  mutate(population = as.numeric(sub("V", "", population)))

# Join population labels
admix_long_americas <- admix_long_americas %>%
  left_join(americas_only_pop_labels, by = "population")

# Assign dominant population for ordering
dominant_pop_americas <- admix_long_americas %>%
  group_by(sample) %>%
  slice_max(population_proportion, n = 1, with_ties = FALSE) %>%
  select(sample, dominant_population = population_label)

admix_long_americas <- admix_long_americas %>%
  left_join(dominant_pop_americas, by = "sample")

# Reorder samples by dominant population and max ancestry proportion
sample_order_americas <- admix_long_americas %>%
  group_by(sample) %>%
  summarise(dominant_population = first(dominant_population),
            max_ancestry = max(population_proportion),
            .groups = "drop") %>%
  arrange(dominant_population, desc(max_ancestry)) %>%
  pull(sample)

admix_long_americas <- admix_long_americas %>%
  mutate(sample = factor(sample, levels = sample_order_americas))

# Ensure population labels are factors with consistent ordering for legend
admix_long_americas <- admix_long_americas %>%
  mutate(population_label = factor(population_label, 
                                   levels = americas_only_pop_labels$population_label))

# Reuse americas_palette from above
names(americas_palette) <- americas_only_pop_labels$population_label

# Create the stacked barplot
ancestry_plot_americas <- ggplot(admix_long_americas, 
                                 aes(x = sample, y = population_proportion, fill = population_label)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = americas_palette, 
                    breaks = americas_only_pop_labels$population_label,
                    name = "Population") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", size = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 11),
    plot.title = element_text(size = 14, hjust = 0.5)
  ) +
  labs(title = "Americas-Only L.infantum ADMIXTURE Ancestry Proportions (K=13)",
       x = "Samples ordered by population",
       y = "Ancestry Proportion") +
  guides(fill = guide_legend(ncol = 1))

# Display the plot
print(ancestry_plot_americas)

# Save the plot as a jpeg
ggsave(file.path(FIGURES_DIR, "Americas_Only", 
                 "americas_only_stacked_ancestry_barplot_k13.jpeg"), 
       width = 10, height = 6)

# Save the plot as a pdf
ggsave(file.path(FIGURES_DIR, "Americas_Only", 
                 "americas_only_stacked_ancestry_barplot_k13.pdf"), 
       width = 10, height = 6)




# SCRIPT COMPLETE --------------------------------------------------------------

cat("\n")
cat("====================================================================\n")
cat("ADMIXTURE POPULATION STRUCTURE ANALYSIS COMPLETE\n")
cat("====================================================================\n")
cat("Output directories:\n")
cat("  Figures:", FIGURES_DIR, "\n")
cat("  Analysis outputs:", RESULTS_DIR, "\n")
cat("====================================================================\n\n")

