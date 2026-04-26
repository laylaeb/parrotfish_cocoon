##-----------------------------------------------------------Cocoon metagenomics

library(phyloseq)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)

# loading data
cocoon_meta_abundance <- read.csv("cocoon_abundance.csv")
cocoon_kegg_path <- read.csv("coocon_egg_prot_annotations.csv", sep=",")

## what are the most abundant ones? (calculate the relative abundance per contig)
meta_ab <- cocoon_meta_abundance %>%
  mutate(RA = meta_1 / sum(meta_1)) %>%
  select(Contig, everything()) %>%
  arrange(desc(RA))

merge_ab_an <- merge(meta_ab, cocoon_kegg_path, by.x = "Contig", by.y = "query")

metagenome_functional_analysis <- merge_ab_an %>%
  mutate(
    COG_category_final = case_when(
      is.na(COG_category) | COG_category == "" | COG_category == "-" ~ "Not assigned",
      TRUE ~ COG_category
    ),
    
    Functional_category = case_when(
      COG_category_final == "Not assigned" ~ "No annotation found",
      COG_category_final == "S" ~ "COG S - Function unknown",  # ← CORRECTION ICI
      str_detect(Description, "DUF") ~ "Domain of unknown function",
      str_detect(Description, "unknown function|Unknown function|hypothetical protein") ~ "Protein of unknown function",
      is.na(Description) | Description == "-" | Description == "" ~ "Protein of unknown function",
      TRUE ~ "Known function"
    )
  )

# analyse par catégorie COG seulement
metagenome_cog_summary <- metagenome_functional_analysis %>%
  group_by(COG_category_final) %>%
  summarise(
    total_contigs = n(),
    total_RA = sum(RA, na.rm = TRUE),
    percentage = round(total_RA * 100, 2),
    .groups = 'drop'
  ) %>%
  arrange(desc(total_RA))

print("=== Top 10 COG categories - Métagénome ===")
head(metagenome_cog_summary, 10)

# analyse par statut fonctionnel
metagenome_functional_summary <- metagenome_functional_analysis %>%
  group_by(Functional_category) %>%
  summarise(
    total_contigs = n(),
    total_RA = sum(RA, na.rm = TRUE),
    percentage = round(total_RA * 100, 2),
    .groups = 'drop'
  ) %>%
  arrange(desc(total_RA))

print("=== Résumé par statut fonctionnel - Métagénome ===")
metagenome_functional_summary


metagenome_clear_summary <- metagenome_functional_analysis %>%
  mutate(
    Final_category = case_when(
      Functional_category %in% c("No annotation found", 
                                 "COG S - Function unknown", 
                                 "Domain of unknown function", 
                                 "Protein of unknown function") ~ "UNKNOWN FUNCTIONS",
      TRUE ~ paste0("COG ", COG_category_final, " - Known functions")
    )
  ) %>%
  group_by(Final_category) %>%
  summarise(
    total_contigs = n(),
    total_RA = sum(RA, na.rm = TRUE),
    percentage = round(total_RA * 100, 2),
    .groups = 'drop'
  ) %>%
  arrange(desc(total_RA))

print("=== RÉSUMÉ CLAIR MÉTAGÉNOME : Connu vs Inconnu ===")
metagenome_clear_summary

## Abundant KEGG pathways
merged_kegg <- merge_ab_an %>%
  filter(!is.na(KEGG_Pathway), KEGG_Pathway != "-")

merged_kegg_expanded <- merged_kegg %>%
  separate_rows(KEGG_Pathway, sep = ",")

kegg_pathway_summary <- merged_kegg_expanded %>%
  group_by(KEGG_Pathway) %>%
  summarise(
    total_RA = sum(RA, na.rm = TRUE),
    n_contigs = n()
  ) %>%
  arrange(desc(total_RA))

head(kegg_pathway_summary, 10)

# KEGG_Pathway total_RA n_contigs
# <chr>           <dbl>     <int>
# 1 ko01100        0.123      21854
# 2 map01100       0.123      21854
# 3 ko01110        0.0539      9558
# 4 map01110       0.0539      9558
# 5 ko01130        0.0427      7609
# 6 map01130       0.0427      7609
# 7 ko01120        0.0398      7197
# 8 map01120       0.0398      7197
# 9 ko01230        0.0232      4161
# 10 map01230      0.0232      4161

kegg_clean_summary <- kegg_pathway_summary %>%
  filter(str_detect(KEGG_Pathway, "^map")) %>%
  arrange(desc(total_RA))

# Check coverage 
clean_coverage <- sum(kegg_clean_summary$total_RA)
cat("Clean KEGG coverage:", round(clean_coverage * 100, 2), "%\n")

head(kegg_clean_summary, 10)

kegg_clean_summary %>%
  top_n(15, total_RA) %>%
  ggplot(aes(x = reorder(KEGG_Pathway, total_RA), y = total_RA)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 10 KEGG Pathways by Relative Abundance",
    x = "KEGG Pathway",
    y = "Total Relative Abundance"
  )

kegg_clean_summary %>%
  top_n(15, total_RA) %>%
  ggplot(aes(x = reorder(KEGG_Pathway, total_RA), y = total_RA)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 15 KEGG Pathways by Relative Abundance",
    x = "KEGG Pathway",
    y = "Total Relative Abundance"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5)
  )

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

# Analyser les sous-classes EC2 (Transferases)
ec2_analysis <- merge_ab_an %>%
  filter(str_detect(EC, "^2\\.")) %>%
  mutate(
    EC2_subclass = case_when(
      str_detect(EC, "^2\\.1\\.") ~ "EC2.1 - Transferases (one-carbon groups)",
      str_detect(EC, "^2\\.2\\.") ~ "EC2.2 - Transferases (aldehyde or ketone groups)",
      str_detect(EC, "^2\\.3\\.") ~ "EC2.3 - Acyltransferases",
      str_detect(EC, "^2\\.4\\.") ~ "EC2.4 - Glycosyltransferases",
      str_detect(EC, "^2\\.5\\.") ~ "EC2.5 - Transferases (alkyl or aryl groups)",
      str_detect(EC, "^2\\.6\\.") ~ "EC2.6 - Transferases (nitrogenous groups)",
      str_detect(EC, "^2\\.7\\.") ~ "EC2.7 - Transferases (phosphorus-containing groups)",
      str_detect(EC, "^2\\.8\\.") ~ "EC2.8 - Transferases (sulfur-containing groups)",
      TRUE ~ "EC2 - Other transferases"
    )
  ) %>%
  group_by(EC2_subclass) %>%
  summarise(
    protein_count = n(),
    total_RA = sum(RA, na.rm = TRUE),
    percentage = round(total_RA * 100, 3),
    .groups = 'drop'
  ) %>%
  arrange(desc(total_RA))

print(ec2_analysis)

# EC2_subclass                                        protein_count total_RA percentage
# <chr>                                                       <int>    <dbl>      <dbl>
# 1 EC2.7 - Transferases (phosphorus-containing groups)          5013 0.0292        2.92 
# 2 EC2.1 - Transferases (one-carbon groups)                     1775 0.0103        1.03 
# 3 EC2.3 - Acyltransferases                                     1439 0.00791       0.791
# 4 EC2.4 - Glycosyltransferases                                  968 0.00561       0.561
# 5 EC2.5 - Transferases (alkyl or aryl groups)                   998 0.00552       0.552
# 6 EC2.6 - Transferases (nitrogenous groups)                     643 0.00352       0.352
# 7 EC2.8 - Transferases (sulfur-containing groups)               502 0.00278       0.278
# 8 EC2.2 - Transferases (aldehyde or ketone groups)              326 0.00180       0.18 
# 9 EC2 - Other transferases                                       61 0.000332      0.033

# Analyser les sous-classes EC3 (Hydrolases)
ec3_analysis <- merge_ab_an %>%
  filter(str_detect(EC, "^3\\.")) %>%
  mutate(
    EC3_subclass = case_when(
      str_detect(EC, "^3\\.1\\.") ~ "EC3.1 - Hydrolases (ester bonds)",
      str_detect(EC, "^3\\.2\\.") ~ "EC3.2 - Glycosidases",
      str_detect(EC, "^3\\.3\\.") ~ "EC3.3 - Hydrolases (ether bonds)",
      str_detect(EC, "^3\\.4\\.") ~ "EC3.4 - Peptidases",
      str_detect(EC, "^3\\.5\\.") ~ "EC3.5 - Hydrolases (carbon-nitrogen bonds)",
      str_detect(EC, "^3\\.6\\.") ~ "EC3.6 - Hydrolases (acid anhydrides)",
      str_detect(EC, "^3\\.7\\.") ~ "EC3.7 - Hydrolases (carbon-carbon bonds)",
      str_detect(EC, "^3\\.8\\.") ~ "EC3.8 - Hydrolases (halide bonds)",
      str_detect(EC, "^3\\.9\\.") ~ "EC3.9 - Hydrolases (phosphorus-nitrogen bonds)",
      str_detect(EC, "^3\\.10\\.") ~ "EC3.10 - Hydrolases (sulfur-nitrogen bonds)",
      str_detect(EC, "^3\\.11\\.") ~ "EC3.11 - Hydrolases (carbon-phosphorus bonds)",
      str_detect(EC, "^3\\.12\\.") ~ "EC3.12 - Hydrolases (sulfur-sulfur bonds)",
      str_detect(EC, "^3\\.13\\.") ~ "EC3.13 - Hydrolases (carbon-sulfur bonds)",
      TRUE ~ "EC3 - Other hydrolases"
    )
  ) %>%
  group_by(EC3_subclass) %>%
  summarise(
    protein_count = n(),
    total_RA = sum(RA, na.rm = TRUE),
    percentage = round(total_RA * 100, 3),
    .groups = 'drop'
  ) %>%
  arrange(desc(total_RA))

# print(ec3_analysis)
# 
# EC3_subclass                                  protein_count  total_RA percentage
# <chr>                                                 <int>     <dbl>      <dbl>
# 1 EC3.1 - Hydrolases (ester bonds)                       2065 0.0119         1.19 
# 2 EC3.6 - Hydrolases (acid anhydrides)                   1687 0.00979        0.979
# 3 EC3.4 - Peptidases                                     1432 0.00813        0.813
# 4 EC3.2 - Glycosidases                                   1200 0.00744        0.744
# 5 EC3.5 - Hydrolases (carbon-nitrogen bonds)             1246 0.00676        0.676
# 6 EC3.3 - Hydrolases (ether bonds)                         75 0.000428       0.043
# 7 EC3.8 - Hydrolases (halide bonds)                        39 0.000207       0.021
# 8 EC3.7 - Hydrolases (carbon-carbon bonds)                 40 0.000194       0.019
# 9 EC3.10 - Hydrolases (sulfur-nitrogen bonds)              10 0.0000546      0.005
# 10 EC3.13 - Hydrolases (carbon-sulfur bonds)                 6 0.0000443      0.004
# 11 EC3.11 - Hydrolases (carbon-phosphorus bonds)             7 0.0000313      0.003

# Analyser les sous-classes EC1 (Oxidoreductases)
ec1_analysis <- merge_ab_an %>%
  filter(str_detect(EC, "^1\\.")) %>%
  mutate(
    EC1_subclass = case_when(
      str_detect(EC, "^1\\.1\\.") ~ "EC1.1 - Oxidoreductases (CH-OH group)",
      str_detect(EC, "^1\\.2\\.") ~ "EC1.2 - Oxidoreductases (aldehyde/oxo group)",
      str_detect(EC, "^1\\.3\\.") ~ "EC1.3 - Oxidoreductases (CH-CH group)",
      str_detect(EC, "^1\\.4\\.") ~ "EC1.4 - Oxidoreductases (CH-NH2 group)",
      str_detect(EC, "^1\\.5\\.") ~ "EC1.5 - Oxidoreductases (CH-NH group)",
      str_detect(EC, "^1\\.6\\.") ~ "EC1.6 - Oxidoreductases (NADH/NADPH)",
      str_detect(EC, "^1\\.7\\.") ~ "EC1.7 - Oxidoreductases (other nitrogenous compounds)",
      str_detect(EC, "^1\\.8\\.") ~ "EC1.8 - Oxidoreductases (sulfur group)",
      str_detect(EC, "^1\\.9\\.") ~ "EC1.9 - Oxidoreductases (heme group)",
      str_detect(EC, "^1\\.10\\.") ~ "EC1.10 - Oxidoreductases (diphenols)",
      str_detect(EC, "^1\\.11\\.") ~ "EC1.11 - Oxidoreductases (peroxide)",
      str_detect(EC, "^1\\.12\\.") ~ "EC1.12 - Oxidoreductases (hydrogen)",
      str_detect(EC, "^1\\.13\\.") ~ "EC1.13 - Oxidoreductases (single oxygen)",
      str_detect(EC, "^1\\.14\\.") ~ "EC1.14 - Oxidoreductases (paired donors)",
      str_detect(EC, "^1\\.15\\.") ~ "EC1.15 - Oxidoreductases (superoxide)",
      str_detect(EC, "^1\\.16\\.") ~ "EC1.16 - Oxidoreductases (oxidizing metal ions)",
      str_detect(EC, "^1\\.17\\.") ~ "EC1.17 - Oxidoreductases (CH or CH2 groups)",
      str_detect(EC, "^1\\.18\\.") ~ "EC1.18 - Oxidoreductases (iron-sulfur proteins)",
      str_detect(EC, "^1\\.19\\.") ~ "EC1.19 - Oxidoreductases (reducing flavodoxin)",
      str_detect(EC, "^1\\.20\\.") ~ "EC1.20 - Oxidoreductases (phosphorus/arsenic)",
      str_detect(EC, "^1\\.21\\.") ~ "EC1.21 - Oxidoreductases (X-H and Y-H to X-Y)",
      str_detect(EC, "^1\\.22\\.") ~ "EC1.22 - Oxidoreductases (H2O2 as acceptor)",
      str_detect(EC, "^1\\.23\\.") ~ "EC1.23 - Oxidoreductases (reducing O2 to H2O2)",
      str_detect(EC, "^1\\.97\\.") ~ "EC1.97 - Other oxidoreductases",
      TRUE ~ "EC1 - Other oxidoreductases"
    )
  ) %>%
  group_by(EC1_subclass) %>%
  summarise(
    protein_count = n(),
    total_RA = sum(RA, na.rm = TRUE),
    percentage = round(total_RA * 100, 3),
    .groups = 'drop'
  ) %>%
  arrange(desc(total_RA))

print(ec1_analysis)

# EC1_subclass                                 protein_count total_RA percentage
# <chr>                                                <int>    <dbl>      <dbl>
#   1 EC1.1 - Oxidoreductases (CH-OH group)                 1867  0.0105       1.04 
# 2 EC1.2 - Oxidoreductases (aldehyde/oxo group)          1045  0.00575      0.575
# 3 EC1.3 - Oxidoreductases (CH-CH group)                  753  0.00414      0.414
# 4 EC1.6 - Oxidoreductases (NADH/NADPH)                   771  0.00410      0.41 
# 5 EC1.4 - Oxidoreductases (CH-NH2 group)                 530  0.00306      0.306
# 6 EC1.17 - Oxidoreductases (CH or CH2 groups)            500  0.00283      0.283
# 7 EC1.8 - Oxidoreductases (sulfur group)                 463  0.00270      0.27 
# 8 EC1.14 - Oxidoreductases (paired donors)               451  0.00234      0.234
# 9 EC1.5 - Oxidoreductases (CH-NH group)                  356  0.00205      0.205
# 10 EC1.9 - Oxidoreductases (heme group)                   264  0.00153      0.153
# # ℹ 11 more rows

##-------------------------------------------------------------------- CAZymes

cazy_analysis <- merge_ab_an %>%
  filter(!is.na(CAZy) & CAZy != "-") %>%
  mutate(
    CAZy_family = case_when(
      str_detect(CAZy, "^GH") ~ "GH - Glycoside Hydrolases",
      str_detect(CAZy, "^GT") ~ "GT - Glycosyl Transferases",
      str_detect(CAZy, "^PL") ~ "PL - Polysaccharide Lyases",
      str_detect(CAZy, "^CE") ~ "CE - Carbohydrate Esterases",
      str_detect(CAZy, "^AA") ~ "AA - Auxiliary Activities",
      TRUE ~ "Other CAZy"
    )
  ) %>%
  group_by(CAZy_family) %>%
  summarise(
    protein_count = n(),
    total_RA = sum(RA, na.rm = TRUE),
    percentage = round(total_RA * 100, 3),
    .groups = 'drop'
  ) %>%
  arrange(desc(total_RA))

print(cazy_analysis)

# Analyser les familles GH spécifiques
gh_analysis <- merge_ab_an %>%
  filter(str_detect(CAZy, "^GH")) %>%
  mutate(
    GH_family = case_when(
      str_detect(CAZy, "GH1") ~ "GH1 - β-glucosidases",
      str_detect(CAZy, "GH2") ~ "GH2 - β-galactosidases",
      str_detect(CAZy, "GH3") ~ "GH3 - β-glucosidases",
      str_detect(CAZy, "GH4") ~ "GH4 - α-glucosidases",
      str_detect(CAZy, "GH5") ~ "GH5 - Cellulases",
      str_detect(CAZy, "GH6") ~ "GH6 - Cellulases",
      str_detect(CAZy, "GH7") ~ "GH7 - Cellulases",
      str_detect(CAZy, "GH8") ~ "GH8 - Cellulases",
      str_detect(CAZy, "GH9") ~ "GH9 - Cellulases",
      str_detect(CAZy, "GH10") ~ "GH10 - Xylanases",
      str_detect(CAZy, "GH11") ~ "GH11 - Xylanases",
      str_detect(CAZy, "GH12") ~ "GH12 - Cellulases",
      str_detect(CAZy, "GH13") ~ "GH13 - α-amylases",
      str_detect(CAZy, "GH14") ~ "GH14 - β-amylases",
      str_detect(CAZy, "GH15") ~ "GH15 - Glucoamylases",
      str_detect(CAZy, "GH16") ~ "GH16 - β-glucanases",
      str_detect(CAZy, "GH17") ~ "GH17 - β-glucanases",
      str_detect(CAZy, "GH18") ~ "GH18 - Chitinases",
      str_detect(CAZy, "GH19") ~ "GH19 - Chitinases",
      str_detect(CAZy, "GH20") ~ "GH20 - β-hexosaminidases",
      str_detect(CAZy, "GH23") ~ "GH23 - Lysozymes",
      str_detect(CAZy, "GH24") ~ "GH24 - Lysozymes",
      str_detect(CAZy, "GH25") ~ "GH25 - Lysozymes",
      str_detect(CAZy, "GH26") ~ "GH26 - β-mannanases",
      str_detect(CAZy, "GH27") ~ "GH27 - α-galactosidases",
      str_detect(CAZy, "GH28") ~ "GH28 - Polygalacturonases",
      str_detect(CAZy, "GH29") ~ "GH29 - α-fucosidases",
      str_detect(CAZy, "GH30") ~ "GH30 - β-xylosidases",
      TRUE ~ "Other GH families"
    )
  ) %>%
  group_by(GH_family) %>%
  summarise(
    count = n(),
    total_RA = sum(RA, na.rm = TRUE),
    percentage = round(total_RA * 100, 3),
    .groups = 'drop'
  ) %>%
  arrange(desc(total_RA))

print(gh_analysis)

# GH_family              count  total_RA percentage
# <chr>                  <int>     <dbl>      <dbl>
# 1 GH3 - β-glucosidases     290 0.00186        0.186
# 2 GH2 - β-galactosidases   223 0.00154        0.154
# 3 GH5 - Cellulases         129 0.000777       0.078
# 4 GH1 - β-glucosidases     126 0.000693       0.069
# 5 GH6 - Cellulases          24 0.000140       0.014
# 6 GH7 - Cellulases          25 0.000136       0.014
# 7 GH9 - Cellulases          15 0.0000916      0.009
# 8 GH8 - Cellulases           9 0.0000562      0.006
# 9 GH4 - α-glucosidases      10 0.0000536      0.005
