###### ---- Beta diversity cocoon ---- ######

library(phyloseqCompanion)
library(phyloseq)

cocoon_metadata <- sample.data.frame(cocoon_r)
asv_table_cocoon <- as.matrix((otu_table(cocoon_r, taxa_are_rows=T)))
asv_table_cocoon <- asv_table_cocoon[rowSums(asv_table_cocoon[])>0,]

meta_cocoon_ait <- sample.data.frame(prune_cocoon)

physeq_clr<- microbiome::transform(prune_cocoon, "clr")

asv_fish_clr <- otu_table(physeq_clr, taxa_are_rows=T)

fish_aitchison <- phyloseq::distance(t(asv_fish_clr), method = "euclidean")
aitchison_pcoa <- ecodist::pco(fish_aitchison)

aitchison_fish_poca_df<- data.frame(pcoa1 = aitchison_pcoa$vectors[,1],
                                    pcoa2 = aitchison_pcoa$vectors[,2],
                                    pcoa3 = aitchison_pcoa$vectors[,3] )

aitchison_fish_poca_df <- cbind(aitchison_fish_poca_df,
                                meta_fish = meta_cocoon_ait$what, meta_shape=meta_cocoon_ait$site)

pcoa_aitch_fish_plot <- ggplot(data = aitchison_fish_poca_df, aes(x=pcoa1, y=pcoa3,color=meta_fish)) +
  geom_point() +
  labs(x = "PC1",
       y = "PC3",
       title = "Atichison PCoA") +
  theme(title = element_text(size = 10)) 

pcoa_aitch_fish_plot

palette_colors <- fish(option = "Oncorhynchus_kisutch", 5, direction = -1)
palette_colors[c(1, length(palette_colors))] <- palette_colors[c(length(palette_colors), 1)]

ellipse_aitchison <- ggplot(data = aitchison_fish_poca_df, aes(x = pcoa1, y = pcoa2, color = meta_fish, shape = meta_shape)) +
  scale_color_manual(values = palette_colors) +
  scale_fill_manual(values = palette_colors) +
  stat_ellipse(aes(group = meta_fish),
               color = "black",
               fill = "white",
               alpha = 0.2,
               level = 0.95,
               type = "norm",
               geom = "polygon") +
  geom_point(size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

print(ellipse_aitchison)
ellipse_aitchison <- ellipse_aitchison +
  labs(
    x = paste0("PCoA1 (", round(variance_proportion_ai[1] * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(variance_proportion_ai[2] * 100, 1), "%)")
  )

variance_proportion_ai <-aitchison_pcoa$values / sum(aitchison_pcoa$values)

library(vegan)
stopifnot(identical(attr(fish_aitchison, "Labels"), rownames(meta_cocoon_ait)))
sites <- unique(meta_cocoon_ait$site)
bdisp_by_site <- list()

for (s in sites) {
  cat("\n=== Site :", s, "===\n")
 
  idx <- meta_cocoon_ait$site == s
  
    dist_site_mat <- as.matrix(fish_aitchison)[idx, idx]
  dist_site     <- as.dist(dist_site_mat)
  
  what_site <- factor(meta_cocoon_ait$what[idx])         

  bdisp_site <- betadisper(dist_site, group = what_site, type = "centroid")
  bdisp_by_site[[as.character(s)]] <- bdisp_site
  
  print(anova(bdisp_site))
  print(permutest(bdisp_site, pairwise = TRUE, permutations = 999))
  
}

## ---- PERMANOVA ---- ##
cocoon_ado <- adonis2(
  fish_aitchison ~ what + site,
  data = meta_cocoon_ait,
  permutations = 999,
  by = "margin"
)
cocoon_ado
summary(cocoon_ado)

pairwise_cocoon <- pairwise.adonis2(fish_aitchison ~ what,
                                    data = meta_cocoon_ait)


interaction <- adonis2(fish_aitchison ~ what * site,
                    data = meta_cocoon_ait, permutations = 999, by = "margin")

meta_cocoon_ait$what_site <- interaction(meta_cocoon_ait$what, meta_cocoon_ait$site, sep = "_")
table(meta_cocoon_ait$what_site)

common <- intersect(rownames(meta_cocoon_ait), rownames(fish_aitchison))
meta2  <- meta_cocoon_ait[common, ]
dist2  <- fish_aitchison[common, common]

library(pairwiseAdonis)

pair_all <- pairwise.adonis(fish_aitchison,
  factors = meta_cocoon_ait$what_site,
  p.adjust.m="fdr",
  
)
pair_all

## ---- Indicator species ---- ##
cocoon_r_bt <- subset_samples(cocoon_r, site %in% "BT")
cocoon_r_g <- subset_samples(cocoon_r, site %in% "G")

tax_genus_bt <- tax_glom(cocoon_r_bt, taxrank=rank_names(cocoon_r_bt)[6], NArm=F)
cocoon_new_physeq_bt <- prune_taxa(taxa_sums(tax_genus_bt) > 5, tax_genus_bt)  # par ex. >10 reads au total
tax_genus_relabund_bt <- transform_sample_counts(cocoon_new_physeq_bt, function(x) x / sum(x))
tax_table_all <- tax_table(tax_genus_bt)
write.csv(tax_table_all, "tax_table_all.csv")

tax_genus_g <- tax_glom(cocoon_r_g, taxrank=rank_names(cocoon_r_g)[6], NArm=F)
cocoon_new_physeq_g <- prune_taxa(taxa_sums(tax_genus_g) > 5, tax_genus_g)  # par ex. >10 reads au total
tax_genus_relabund_g <- transform_sample_counts(cocoon_new_physeq_g, function(x) x / sum(x))


otu_bt <- as.data.frame(t(otu_table(tax_genus_relabund_bt)))
group_bt <- sample_data(tax_genus_relabund_bt)$what

otu_g <- as.data.frame(t(otu_table(tax_genus_relabund_g)))
group_g <- sample_data(tax_genus_relabund_g)$what

library(indicspecies)
set.seed(123)
indval_res_g <- multipatt(otu_g, group_g, func = "IndVal.g", duleg = T, control = how(nperm=999))
indval_res_bt <- multipatt(otu_bt, group_bt, func = "IndVal.g", duleg = T, control = how(nperm=999))

#Fidélité : fréquence de l’ASV dans les échantillons du groupe.
#Spécificité : abondance relative dans le groupe vs les autres.
#duleg=F signifie qu il va chercher les taxa qui sont indicateurs pour un groupe en particulier 

res_g <- summary(indval_res_g, indvalcomp = TRUE, alpha = 0.05)
res_bt <- summary(indval_res_bt, indvalcomp = TRUE, alpha = 0.05)

tax_genus <- tax_glom(cocoon_r, taxrank=rank_names(cocoon_r)[6], NArm=F)
cocoon_new_physeq<- prune_taxa(taxa_sums(tax_genus) > 10, tax_genus)  # par ex. >10 reads au total
tax_genus_relabund<- transform_sample_counts(cocoon_new_physeq, function(x) x / sum(x))
tax_table_all <- tax_table(tax_genus_relabund)
write.csv(tax_table_all, "update_tax_table_all.csv")

group <- sample_data(tax_genus)$what
otu <- as.data.frame(t(otu_table(tax_genus_relabund)))

ctrl <- how(nperm = 999)
setBlocks(ctrl) <- cocoon_metadata$site  

indval_all <- multipatt(otu, group, func = "IndVal.g", 
                        duleg = TRUE, control = ctrl) 
res_all <-  summary(indval_all, indvalcomp = TRUE, alpha = 0.05) 

res_df <- as.data.frame(indval_res$sign)

res_df$ASV <- rownames(res_df)
taxa_df <- as.data.frame(tax_table(cocoon_new_physeq ))
taxa_df$ASV <- rownames(taxa_df)

final <- merge(res_df, taxa_df, by = "ASV")

final_sig <- final[final$p.value <= 0.05, ]
head(final_sig)

write.csv(final_sig,"202505_indicspecies_ASV.csv")
write.csv(final_sig,"202505_indicspecies_genus.csv")

############################################################################
# Relative Abundance
############################################################################
tax_family_all <- tax_glom(cocoon_r, taxrank=rank_names(cocoon_r)[5], NArm=F)

s_envelope <- subset_samples(cocoon_r, what == "envelope")
s_skin <- subset_samples(cocoon_r, what == "fish_skin")
s_mouth <- subset_samples(cocoon_r, what == "mouth")
s_sediment <- subset_samples(cocoon_r, what == "sediment")
s_water <- subset_samples(cocoon_r, what == "water")

tax_family_env <- tax_glom(s_envelope, taxrank=rank_names(s_envelope)[6], NArm=F)
tax_family_env_RA <- transform_sample_counts(tax_family_env, function(x) x / sum(x))
asv_family_en  <- otu_table(tax_family_env_RA, taxa_are_rows=T)
tax_family_en <- tax_table(tax_family_env)
write.csv(asv_family_en, "cocoon_genus.csv")
write.csv(tax_family_en, "2025_cocoon_genus.csv")

tax_family_skin <- tax_glom(s_skin, taxrank=rank_names(s_skin)[6], NArm=F)
tax_family_skin_RA <- transform_sample_counts(tax_family_skin, function(x) x / sum(x))
asv_family_skin  <- otu_table(tax_family_skin_RA, taxa_are_rows=T)
tax_family_skin <- tax_table(tax_family_skin)
write.csv(asv_family_skin, "skin_family.csv")
write.csv(tax_family_skin, "2025_skin_taxa.csv")

tax_family_mouth <- tax_glom(s_mouth, taxrank=rank_names(s_mouth)[6], NArm=F)
tax_family_mouth_RA <- transform_sample_counts(tax_family_mouth, function(x) x / sum(x))
asv_family_mouth  <- otu_table(tax_family_mouth_RA, taxa_are_rows=T)
tax_family_mouth <- tax_table(tax_family_mouth)
write.csv(asv_family_mouth, "mouth_family.csv")
write.csv(tax_family_mouth, "2025_mouth_taxa.csv")

tax_family_sed <- tax_glom(s_sediment, taxrank=rank_names(s_sediment)[6], NArm=F)
tax_family_sed_RA <- transform_sample_counts(tax_family_sed, function(x) x / sum(x))
asv_family_sed  <- otu_table(tax_family_sed_RA, taxa_are_rows=T)
tax_family_sed <- tax_table(tax_family_sed)
write.csv(asv_family_sed, "sed_family.csv")
write.csv(tax_family_sed, "2025_sed_taxa.csv")

tax_family_wat <- tax_glom(s_water, taxrank=rank_names(s_water)[6], NArm=F)
tax_family_wat_RA <- transform_sample_counts(tax_family_wat, function(x) x / sum(x))
asv_family_wat  <- otu_table(tax_family_wat_RA, taxa_are_rows=T)
tax_family_wat <- tax_table(tax_family_wat)
write.csv(asv_family_wat, "wat_family.csv")
write.csv(tax_family_wat, "2025_wat_taxa.csv")

## ---- STACKPLOT  ---- ##
cocoon_all<- sample.data.frame(tax_family_all)
cocoon_Family = transform_sample_counts(tax_family_all, function(x) x / sum(x))
TopASV_p <- names(sort(taxa_sums(cocoon_Family), TRUE)[1:10])
top10_popfish_p <- prune_species(TopASV_p, cocoon_Family )
top10_popfish_p <- prune_taxa(taxa_sums(top10_popfish_p)>0, top10_popfish_p)
top_Family<-as.data.frame(tax_table(top10_popfish_p))

popfish_merge_by = merge_samples(tax_family_all, "what")
blank_RA = transform_sample_counts(popfish_merge_by, function(x) x / sum(x))
data_blank <- psmelt(blank_RA) # create dataframe from phyloseq object
data_blank$Family <-as.character(data_blank$Family) #convert to character

sumtot_blank <-
  data_blank %>% group_by(Family) %>% summarize(sum = sum(Abundance)) %>%
  filter(Family %in% top_Family$Family) %>% filter(!(Family %in% c(""," f__uncultured",NA,"Unknown Family")))

data_blank$Family[!(data_blank$Family %in% sumtot_blank$Family)] <- "Other"

data_blank$blank <- "blank"
data_blank$totalAbundance <- sum(data_blank$Abundance)

data_blank_mod <- data_blank%>%
  group_by(Family, Sample)%>%
  summarise(abundance = sum(Abundance))%>%
  distinct()

#data_blank_mod$Sample<- factor(data_blank_mod$Sample, levels = c("water", "sed","turf","fish"))
#data_blank_mod$Sample <- factor(data_blank_mod$Sample, levels = c("I", "P","H"))

barplot_for_all <- ggplot(data=data_blank_mod, aes(x=Sample, y=abundance, fill=Family))+
  geom_bar(aes(), stat="identity", position="stack")+
  scale_fill_manual(values=c(
    "#FF6B6B", "#6B6BFF", "#FFD06B",
    "#6BD0FF", "#D36B6B",
    "#FF966B",
    "#D396FF",
    "#96D3D3", "#D3D396", "#6BD396",
    "#F7CAC9", "#FFE4B5", "#FFA07A", 
    "#FFD700", "#98FB98", "#F0E68C", 
    "#FF69B4", "#90EE90", "#DA70D6", "#32CD32"
  )) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank()
  )

barplot_for_all + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


## ---- Stackplot including variabilities across samples ---- ## 

type_levels <- c("envelope", "fish_skin", "mouth", "sediment", "water")
type_labels <- c(envelope = "Envelope", fish_skin  = "Skin", mouth = "Mouth", sediment = "Sediment", water= "Water")

phylo <- subset_samples(tax_family_all, what %in% type_levels)
phylo <- prune_taxa(taxa_sums(phylo) > 0, phylo)
phylo_ra <- transform_sample_counts(phylo, function(x) x / sum(x))

data_all <- psmelt(phylo_ra) %>%
  mutate(Family = as.character(Family), what = as.character(what), site = as.character(site), site_lab = recode(site, BT = "S1", G = "S2"))

n_top <- 12 # top 12 families
top_families <- data_all %>%
  mutate(Family_clean = case_when(is.na(Family) | Family %in% c("", "Unknown Family") ~ NA_character_, grepl("uncultured", Family, ignore.case = TRUE) ~ NA_character_, TRUE ~ Family)) %>%
  filter(!is.na(Family_clean)) %>%
  group_by(Family_clean) %>%
  summarise(tot = sum(Abundance), .groups = "drop") %>%
  arrange(desc(tot)) %>%
  slice_head(n = n_top) %>%
  pull(Family_clean)

data_per_family <- data_all %>%
  mutate(Family = ifelse(is.na(Family) | Family %in% c("", "Unknown Family") | grepl("uncultured", Family, ignore.case = TRUE) |!(Family %in% top_families),"Other",Family)) %>%
  group_by(Sample, what, site, site_lab, Family) %>%
  summarise(abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(abundance = abundance / sum(abundance)) %>%
  ungroup() %>%
  mutate(what = factor(what, levels = type_levels, labels = type_labels),site_lab = factor(site_lab, levels = c("S1", "S2"))
  ) %>%
  arrange(what, site_lab, Sample)

data_per_family$Sample <- factor(data_per_family$Sample, levels = unique(as.character(data_per_family$Sample)))

x_labels <- data_per_family %>%
  distinct(Sample, site_lab) %>%
  { setNames(as.character(.$site_lab), .$Sample) }

fam_lev <- c(top_families, "Other")
data_per_family$Family <- factor(data_per_family$Family, levels = fam_lev)

# palette pour les différentes familles

pal_named <- c(
  Flavobacteriaceae     = "#3B6C8F",
  Lentisphaeraceae      = "#C44E00",
  Vibrionaceae          = "#3D8B6E",
  Rhodobacteraceae      = "#B5658A",
  Rubritaleaceae        = "#A67C2A",
  Desulfosarcinaceae    = "#4A6B7C",
  Pirellulaceae         = "#D4A017",
  Cyclobacteriaceae     = "#6B4C7A",
  Woeseiaceae           = "#6E8B5A",
  Kiloniellaceae        = "#8B4E8B",
  Fusobacteriaceae      = "#A64B3C",
  Endozoicomonadaceae   = "#2A9D8F",
  Other                 = "#E6DDC8"
)

pal <- pal_named[fam_lev]
if (any(is.na(pal))) {
  extra <- setdiff(fam_lev, names(pal_named))
  pal[extra] <- c("#7B8FA1", "#C47B6A", "#5B7C99", "#9A7B4F")[seq_along(extra)]
  pal["Other"] <- "#E6DDC8"
}

barplot_st <- ggplot(data_per_family, aes(x = Sample, y = abundance, fill = Family)) +
  geom_col(position = "stack", width = 0.9, colour = NA) +
  facet_grid(~ what, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = x_labels) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = pal) +
  labs(x = NULL, y = "Relative abundance", fill = "Family") +
  theme_classic(base_size = 11) +
  theme(
    axis.line.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x  = element_text(size = 7, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8),
    panel.spacing = unit(0.4, "lines")
  )

barplot_st

# Relative abundance across sample types for Lenti ASVs (Kruskal wallis and dunn test)
library(phyloseq)
library(dplyr)
library(emmeans)

type_levels <- c("envelope", "fish_skin", "mouth", "sediment", "water")

ps <- subset_samples(cocoon_r, what %in% type_levels)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps_ra <- transform_sample_counts(ps, function(x) x / sum(x))

lenti <- subset_taxa(
  ps_ra,
  (!is.na(Order)  & grepl("Lentisphaer", Order,  ignore.case = TRUE)) |
    (!is.na(Family) & grepl("Lentisphaer", Family, ignore.case = TRUE)) |
    (!is.na(Genus)  & grepl("Lentisphaera", Genus, ignore.case = TRUE))
)

otu <- as(otu_table(lenti), "matrix")
if (!taxa_are_rows(lenti)) otu <- t(otu)
lenti_ra <- colSums(otu)

dat <- data.frame(sample_data(ps_ra), check.names = FALSE)
dat$lenti_ra <- lenti_ra[rownames(dat)]
dat$what <- factor(dat$what, levels = type_levels)

kruskal.test(lenti_ra ~ what, data = dat)
dunn.test::dunn.test(dat$lenti_ra, dat$what, method = "bh")

dat$y <- asin(log10(pmin(pmax(dat$lenti_ra, 0), 1)))
m <- glm(y ~ what + site, data = dat)
summary(m)
pairs(emmeans(m, ~ what), adjust = "BH")

library(DHARMa)
library(performance)

simulationOutput <- simulateResiduals(fittedModel = m, plot = FALSE)

plot(simulationOutput)
testUniformity(simulationOutput)
testDispersion(simulationOutput)
testOutliers(simulationOutput)
testZeroInflation(simulationOutput)

plotResiduals(simulationOutput, dat$what)
plotResiduals(simulationOutput, dat$site)

check_model(m)

# BOF 
dat$lenti_ra <- pmin(pmax(dat$lenti_ra, 1e-6), 1 - 1e-6)
m2 <- glmmTMB(lenti_ra ~ what + site, family = beta_family(), data = dat)
plot(simulateResiduals(m2))
pairs(emmeans(m2, ~ what), adjust = "BH")
