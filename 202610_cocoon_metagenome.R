##-----------------------------------------------------------Cocoon metagenomics

library(phyloseq)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

cocoon_meta_abundance <- read.csv("cocoon_abundance.csv")
cocoon_kegg_path      <- read.csv("coocon_egg_prot_annotations.csv", sep = ",")

# RA - part de la protéine dans la couverture totale (les 177 383)
meta_ab <- cocoon_meta_abundance %>%
  mutate(RA = meta_1 / sum(meta_1, na.rm = TRUE)) %>%
  arrange(desc(RA))

# seulement les protéines présentes dans eggNOG 
# leur RA a été calculé AVANT, donc sum(RA) après merge ≈ 0.685
# = 68,5 % de la couverture a une annotation eggNOG
merge_ab_an <- merge(meta_ab, cocoon_kegg_path, by.x = "Contig", by.y = "query")

# ---------------------------------------------------------------------------
# COG (C énergie, M membrane, S inconnu…)
# ---------------------------------------------------------------------------
functional_analysis <- merge_ab_an %>%
  mutate(COG_category_final = case_when(is.na(COG_category) | COG_category == "" | COG_category == "-" ~ "Not assigned", TRUE ~ COG_category
    ),
    Functional_category = case_when(
      COG_category_final == "Not assigned" ~ "No annotation found",  
      COG_category_final == "S" ~ "COG S - Function unknown",
      str_detect(Description, "DUF") ~ "Domain of unknown function",
      str_detect(Description, "unknown function|Unknown function|hypothetical protein") ~ "Protein of unknown function",
      is.na(Description) | Description == "-" | Description == "" ~ "Protein of unknown function",
      TRUE ~ "Known function"
    )
  )

nrow(cocoon_meta_abundance)  # 177383  protéines Prodigal
nrow(cocoon_kegg_path)       # 118153  hit eggNOG
sum(merge_ab_an$RA)          # 0.685   couverture annotée

# % = fraction de toute la couverture (177k), pas des seules annotées
cog_summary <- functional_analysis %>%
  group_by(COG_category_final) %>%
  summarise(
    total_contigs = n(),             
    total_RA = sum(RA, na.rm = TRUE),
    percentage = round(total_RA * 100, 2),
    .groups = "drop"
  ) %>%
  arrange(desc(total_RA))

functional_summary <- functional_analysis %>%
  group_by(Functional_category) %>%
  summarise(
    total_contigs = n(),
    total_RA = sum(RA, na.rm = TRUE),
    percentage = round(total_RA * 100, 2),
    .groups = "drop"
  ) %>%
  arrange(desc(total_RA))

# ---------------------------------------------------------------------------
merged_kegg <- merge_ab_an %>%
  filter(!is.na(KEGG_Pathway), KEGG_Pathway != "-")

merged_kegg_expanded <- merged_kegg %>%
  separate_rows(KEGG_Pathway, sep = ",")

kegg_pathway_summary <- merged_kegg_expanded %>%
  group_by(KEGG_Pathway) %>%
  summarise(
    total_RA = sum(RA, na.rm = TRUE),  # somme des RA des protéines dans la map
    n_contigs = n(),                    # nb de protéines dans cette map
    .groups = "drop"
  ) %>%
  arrange(desc(total_RA))

# seulement les "map" attention à ne pas ajouter les doublons, et donc 426
kegg_clean_summary <- kegg_pathway_summary %>%
  filter(str_detect(KEGG_Pathway, "^map")) %>%
  arrange(desc(total_RA))

# ATTENTION : somme des RA de 426 maps (chevauchements)
# = 83,3 %  ≠  "83,3 % des protéines annotées ont une voie"
clean_coverage <- sum(kegg_clean_summary$total_RA)

kegg_clean_summary %>%
  top_n(15, total_RA) %>%
  ggplot(aes(x = reorder(KEGG_Pathway, total_RA), y = total_RA)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 15 KEGG maps",
    x = "KEGG Pathway",
    y = "Total Relative Abundance"
  ) +
  theme_minimal()

# prendre les KO et les insérer dans kegg mapper
library(dplyr)
library(stringr)

# Extraire tous les KO de la colonne KEGG_ko
all_ko <- merged_kegg_expanded %>%
  select(KEGG_ko) %>%
  filter(!is.na(KEGG_ko) & KEGG_ko != "") %>%
  pull(KEGG_ko)

all_ko_split <- str_split(all_ko, ",") %>%
  unlist()

unique_ko_clean <- all_ko_split %>%
  str_trim() %>%
  str_replace("^ko:", "") %>%
  unique() %>%
  sort()
writeLines(unique_ko_clean, "kegg_unique_ko_list.txt")

# Analyser les enzymes par classes EC
total_proteins <- nrow(merge_ab_an)
total_RA <- sum(merge_ab_an$RA, na.rm = TRUE)

# Protéines avec EC 
enzymes <- merge_ab_an %>%
  filter(!is.na(EC) & EC != "-")

enzyme_count <- nrow(enzymes)
enzyme_RA <- sum(enzymes$RA, na.rm = TRUE)
#[1] 0.2061559

proportion_enzymes_count <- round((enzyme_count / total_proteins) * 100, 2)
proportion_enzymes_RA <- round((enzyme_RA / total_RA) * 100, 2)

ec_analysis <- merge_ab_an %>%
  filter(!is.na(EC) & EC != "-") %>%
  mutate(
    EC_class = case_when(
      str_detect(EC, "^1\\.") ~ "EC1 - Oxidoreductases",
      str_detect(EC, "^2\\.") ~ "EC2 - Transferases", 
      str_detect(EC, "^3\\.") ~ "EC3 - Hydrolases",
      str_detect(EC, "^4\\.") ~ "EC4 - Lyases",
      str_detect(EC, "^5\\.") ~ "EC5 - Isomerases",
      str_detect(EC, "^6\\.") ~ "EC6 - Ligases",
      TRUE ~ "EC - Autres"
    )
  ) %>%
  group_by(EC_class) %>%
  summarise(
    protein_count = n(),
    total_RA = sum(RA, na.rm = TRUE),
    percentage = round(total_RA * 100, 3),
    .groups = 'drop'
  ) %>%
  arrange(desc(total_RA))

print(ec_analysis)

# # A tibble: 6 × 4
# EC_class              protein_count total_RA percentage
# <chr>                         <int>    <dbl>      <dbl>
# 1 EC2 - Transferases            11725   0.0670       6.70
# 2 EC3 - Hydrolases               7807   0.0450       4.50
# 3 EC1 - Oxidoreductases          7944   0.0441       4.41
# 4 EC6 - Ligases                  3777   0.0215       2.15
# 5 EC4 - Lyases                   3001   0.0167       1.67
# 6 EC5 - Isomerases               2080   0.0118       1.18

##-------------------------------------------------------------------- CAZymes
# --- CAZy classes --- #
cazy_analysis <- merge_ab_an %>%
  filter(!is.na(CAZy), CAZy != "-") %>%
  mutate(
    CAZy_class = case_when(
      str_detect(CAZy, "^GH") ~ "GH",
      str_detect(CAZy, "^GT") ~ "GT",
      str_detect(CAZy, "^PL") ~ "PL",
      str_detect(CAZy, "^CE") ~ "CE",
      str_detect(CAZy, "^AA") ~ "AA",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(CAZy_class) %>%
  summarise(
    n = n(),
    RA = sum(RA, na.rm = TRUE),
    pct = round(sum(RA) * 100, 3),
    .groups = "drop"
  ) %>%
  arrange(desc(RA))
cazy_analysis


# --- familles GH attention à corriger --- #
gh_analysis <- merge_ab_an %>%
  filter(str_detect(CAZy, "GH[0-9]+")) %>%
  mutate(GH = str_extract(CAZy, "GH[0-9]+")) %>%
  group_by(GH) %>%
  summarise(
    n = n(),
    RA = sum(RA, na.rm = TRUE),
    pct = round(sum(RA) * 100, 3),
    .groups = "drop"
  ) %>%
  arrange(desc(pct))
gh_analysis
