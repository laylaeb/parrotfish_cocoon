library(dplyr)
library(stringr)
library(tidyr)

bin_meta_abundance <- read.csv("depth.csv")
bin_eggnog <- read.csv("bin1_eggnog_annontations.csv")

bin_eggnog$contigName <- sub("_[0-9]+$", "", bin_eggnog$X.query)
merge_ab_an_mag <- merge(bin_meta_abundance, bin_eggnog,
                         by.x = "contigName", by.y = "contigName")

total_proteins_mag <- nrow(merge_ab_an_mag)
total_RA_mag <- sum(merge_ab_an_mag$Relative_abundance, na.rm = TRUE)

# --- KEGG maps ---
kegg_clean_summary_bin <- merge_ab_an_mag %>%
  filter(!is.na(KEGG_Pathway), KEGG_Pathway != "-") %>%
  separate_rows(KEGG_Pathway, sep = ",") %>%
  filter(str_detect(KEGG_Pathway, "^map")) %>%
  group_by(KEGG_Pathway) %>%
  summarise(
    n = n(),
    RA = sum(Relative_abundance, na.rm = TRUE),
    pct = round(sum(Relative_abundance) / total_RA_mag * 100, 2),
    .groups = "drop"
  ) %>%
  arrange(desc(RA))

kegg_clean_summary_bin

# KO pour KEGG Mapper
merge_ab_an_mag %>%
  filter(!is.na(KEGG_ko), KEGG_ko != "-", KEGG_ko != "") %>%
  pull(KEGG_ko) %>%
  str_split(",") %>%
  unlist() %>%
  str_trim() %>%
  str_replace("^ko:", "") %>%
  unique() %>%
  sort() %>%
  writeLines("kegg_unique_ko_list_bin.txt")

# --- EC Classification ---
enzymes_mag <- merge_ab_an_mag %>%
  filter(!is.na(EC), EC != "-")

cat("EC count:", nrow(enzymes_mag), "/", total_proteins_mag,
    "=", round(nrow(enzymes_mag) / total_proteins_mag * 100, 1), "%\n")

# --- CAZy avec le GH à inclure  ---
cazy_analysis_mag <- merge_ab_an_mag %>%
  filter(!is.na(CAZy), CAZy != "-") %>%
  mutate(
    CAZy_class = case_when(
      str_detect(CAZy, "^GH") ~ "GH",
      str_detect(CAZy, "^GT") ~ "GT",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(CAZy_class) %>%
  summarise(
    n = n(),
    pct = round(sum(Relative_abundance) / total_RA_mag * 100, 3),
    .groups = "drop"
  ) %>%
  arrange(desc(pct))

cazy_analysis_mag

gh_analysis_mag <- merge_ab_an_mag %>%
  filter(str_detect(CAZy, "GH[0-9]+")) %>%
  mutate(GH = str_extract(CAZy, "GH[0-9]+")) %>%
  group_by(GH) %>%
  summarise(
    n = n(),
    pct = round(sum(Relative_abundance) / total_RA_mag * 100, 3),
    .groups = "drop"
  ) %>%
  arrange(desc(pct))

gh_analysis_mag

as.data.frame(table(lenti_micro$microtrait_result$domains_detected))

# ------------ Microtraits ----------

library(dplyr)
library(stringr)
library(ggplot2)

#lenti_micro <- readRDS("~/USER/lenti_micro.rds")

# granularity for figure 4
main_traits <- lenti_micro$microtrait_result$trait_counts_atgranularity1 %>%
  filter(`microtrait_trait-value` > 0) %>%
  arrange(desc(`microtrait_trait-value`))

main_traits

ggplot(main_traits,
       aes(x = reorder(`microtrait_trait-displaynameshort`, `microtrait_trait-value`),
           y = `microtrait_trait-value`)) +
  geom_col() +
  coord_flip() +
  labs(title = "MicroTrait functional traits (MAG)",
       x = NULL, y = "Number of genes") +
  theme_minimal()

main_traits %>%
  filter(str_detect(`microtrait_trait-name`, "degradation")) %>%
  dplyr::select(`microtrait_trait-displaynameshort`, `microtrait_trait-value`)

# degradation complémentaire
degrad_g3 <- lenti_micro$microtrait_result$trait_counts_atgranularity3 %>%
  filter(`microtrait_trait-value` > 0,
         str_detect(`microtrait_trait-name`, "degradation|depolymer")) %>%
  dplyr::select(`microtrait_trait-displaynameshort`, `microtrait_trait-value`) %>%
  arrange(desc(`microtrait_trait-value`))

degrad_g3 # avec détails

# stress traits
stress_traits <- main_traits %>%
  filter(str_detect(`microtrait_trait-displaynameshort`,
                    regex("stress|tolerance|temperature|pH|desiccation", ignore_case = TRUE)))

stress_traits

