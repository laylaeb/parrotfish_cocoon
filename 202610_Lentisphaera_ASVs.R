library(phyloseq)
library(dplyr)
library(Biostrings)

## RA sur tout cocoon_r
ps_rel <- transform_sample_counts(cocoon_r, function(x) x / sum(x))

tt <- data.frame(tax_table(ps_rel), check.names = FALSE)
tt$ASV <- rownames(tt)
lenti_ids <- tt$ASV[
  grepl("Lentisphaer", tt$Order,  ignore.case = TRUE) |
  grepl("Lentisphaer", tt$Family, ignore.case = TRUE) |
  grepl("Lentisphaera", tt$Genus, ignore.case = TRUE)
]
length(lenti_ids)

df <- psmelt(ps_rel) %>%
  filter(OTU %in% lenti_ids)

by_sample_type <- df %>%
  group_by(Sample, what) %>%
  summarise(lenti_ra = sum(Abundance), .groups = "drop") %>%
  group_by(what) %>%
  summarise(
    mean_ra = mean(lenti_ra),sd_ra   = sd(lenti_ra),n = n(),.groups = "drop") %>%
  arrange(desc(mean_ra))
by_sample_type

# occurence dans les différents groupes d'échantillon
occ <- df %>%
  group_by(OTU, what) %>%
  summarise(present = any(Abundance > 0), .groups = "drop") %>%
  group_by(OTU) %>%
  summarise(in_envelope = any(what == "envelope" & present),in_other    = any(what != "envelope" & present),.groups = "drop")

sum(occ$in_envelope)                   # dans le cocoon (mais peuvent être détectés ailleurs)
sum(occ$in_envelope & !occ$in_other)   # dans le cocoon seulement
sum(occ$in_envelope & occ$in_other)    # ailleurs seulement
sum(!occ$in_envelope)                  # absent du cocoon

env_ids <- occ$OTU[occ$in_envelope]

env_df <- df %>% filter(what == "envelope")

env_df %>% # l'ASV le plus abondant dans le cocon à checker sa provenance de microbeatlas
  filter(OTU == "ASV_00001") %>%
  summarise(
    mean_ra    = mean(Abundance),
    sd_ra      = sd(Abundance),
    prevalence = mean(Abundance > 0),
    n          = n()
  )

env_df %>%
  filter(OTU == "ASV_04020") %>%
  summarise(
    mean_ra    = mean(Abundance),
    sd_ra      = sd(Abundance),
    prevalence = mean(Abundance > 0),
    n          = n()
  )

# ## FASTA des 41 détectés dans l'envelope (BLAST)
# fa_all <- readDNAStringSet("2025_cocoon.fasta")
# writeXStringSet(
#   fa_all[names(fa_all) %in% env_ids],
#   "Lentisphaera_41_ASV_to_BLAST.fasta"
# )
