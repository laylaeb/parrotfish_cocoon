library(phyloseq)
library(dplyr)
library(cowplot)
library(fishualize)
library(MASS)
library(performance)
library(dunn.test)
library(MiscMetabar)
library(hillR)

## Alpha diversity - cocoon 
alphadt_cocoon<- microbiome::alpha(
  cocoon_r,
  index = c("observed",
            "diversity_shannon", "evenness_pielou")
)

alphadt_indices_cocoon <- phyloseq::sample_data(alphadt_cocoon)
cocoon_physeq_rar <- phyloseq::merge_phyloseq(cocoon_r, alphadt_indices_cocoon)
meta_cocoon<- sample.data.frame(cocoon_r)

ncolors <- lengths(levels(alphadt_cocoon$what))
primary <- 'what'
primary_color_list <- fish(option="Oncorhynchus_kisutch", 5, direction=-1)

observed_plot <- ggplot(cocoon_metadata, aes(x=what, y=observed, fill=what))+
  geom_boxplot()+
  ylab('Observed')+
  xlab(primary)+
  #facet_wrap(facets = ~ cocoon_metadata$site)+
  scale_fill_manual(name=primary, values=primary_color_list)+
  labs(fill=primary)+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=45,hjust=1), axis.title.x=element_blank())

observed_plot

richness_cocoon <- cocoon_metadata %>% 
  group_by(what)%>%summarise(median=median(observed), x = quantile(observed, c(0.25, 0.5, 0.75),type=2))

# Stats
b <- boxcox(lm(cocoon_metadata$observed ~ 1))
lambda <- b$x[which.max(b$y)]
lambda
new_x_exact <- (cocoon_metadata$observed ^ lambda - 1) / lambda
cocoon_metadata$observed_box <- new_x_exact

model_fixed <- glm(observed ~ what + site, 
                   family = Gamma(link = "log"), 
                   data = cocoon_metadata)

check_model(model_fixed)
testDispersion(model_fixed)
simulationOutput <- simulateResiduals(fittedModel = model_fixed, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
plotQQunif(simulationOutput) 
plotResiduals(simulationOutput) 

emmeans_what <- emmeans(model_fixed, ~ what)
emmeans_what
pairs(emmeans_what, adjust = "tukey")  

model_fixed <- glm(observed ~ what + site, 
                     family = Gamma(link = "log"), 
                     data = cocoon_metadata)

mean_observed_by_type <- cocoon_metadata %>%
  group_by(what) %>%
  summarise(
    mean_observed = mean(observed, na.rm = TRUE),
    sd_observed = sd(observed, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(mean_observed))

print("Mean observed richness by sample type:")
print(mean_observed_by_type)

cocoon_mean <- mean_observed_by_type$mean_observed[mean_observed_by_type$what == "envelope"]
mouth_mean <- mean_observed_by_type$mean_observed[mean_observed_by_type$what == "mouth"]
skin_mean <- mean_observed_by_type$mean_observed[mean_observed_by_type$what == "fish_skin"]

cocoon_mouth_ratio <- cocoon_mean / mouth_mean
cocoon_skin_ratio <- cocoon_mean / skin_mean

print(paste("Cocoon/Mouth ratio:", round(cocoon_mouth_ratio, 2)))
print(paste("Cocoon/Skin ratio:", round(cocoon_skin_ratio, 2)))

# Calcuate hillR q=1 
cocoon_asv <- t(otu_table(cocoon_r, taxa_are_rows=T))
hill_shannon <- hill_taxa(cocoon_asv, q = 1, MARGIN = 1, base = exp(1))

hill_df <- tibble::tibble(
  SampleID = names(hill_shannon),
  Hill_Shannon = as.numeric(hill_shannon)
)
cocoon_alpha <- merge(cocoon_metadata, hill_df, by.x="Sample", by.y="SampleID")

HShannon_plot <- ggplot(cocoon_alpha, aes(x=what, y=Hill_Shannon, fill=what))+
  geom_boxplot()+
  ylab('Exp_Shannon_Entropy')+
  xlab(primary)+
  #facet_wrap(facets = ~ cocoon_metadata$site)+
  scale_fill_manual(name=primary, values=primary_color_list)+
  labs(fill=primary)+
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=45,hjust=1), axis.title.x=element_blank())

HShannon_plot 

glm_hill1 <- glm (Hill_Shannon ~ what + site, family=gaussian(link="log"), data=cocoon_alpha)
check_model(glm_hill1)
testDispersion(glm_hill1)
simulationOutput <- simulateResiduals(fittedModel = glm_hill1, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
plotQQunif(simulationOutput) 
plotResiduals(simulationOutput) 

emmeans_what <- emmeans(glm_hill1, ~ what)
emmeans_what
pairs(emmeans_what, adjust = "tukey")  

mean_shannon_by_type <- cocoon_alpha %>%
  group_by(what) %>%
  summarise(
    mean_shannon = mean(Hill_Shannon, na.rm = TRUE),
    sd_shannon = sd(Hill_Shannon, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(mean_shannon))

ggarrange(observed_plot, HShannon_plot)

## Statistical Tests to assess differences in relative abundance
# Get Lentisphaera (ASV_00001) abundance data across all sample types
lentisphaera_data <- all_relative_abundance_data %>%
  filter(ASV == "ASV_00001") %>%
  select(ASV, Sample, what, relative_abundance)

kruskal_site <- kruskal.test(relative_abundance ~ site, data = lentisphaera_data)
kruskal_result <- kruskal.test(relative_abundance ~ what, data = lentisphaera_data)

if(kruskal_result$p.value < 0.05) {
  library(dunn.test)
  dunn_result <- dunn.test(lentisphaera_data$relative_abundance, 
                           lentisphaera_data$what, 
                           method = "bonferroni")
  
  print("Dunn's post-hoc test results:")
  print(dunn_result)
}

## ------------------------------------------ Extraire les ASVs liés à Lentisphaera, Lentisphaeraceae et Lentisphaerales et Blaster
ps_rel_all <- transform_sample_counts(only_cocoon, function(x) x / sum(x))

sub_Lenti_rel <- subset_taxa(
  ps_rel_all,
  (!is.na(Order)  & grepl("Lentisphaer", Order,  ignore.case = TRUE)) |
    (!is.na(Family) & grepl("Lentisphaer", Family, ignore.case = TRUE)) |
    (!is.na(Genus)  & grepl("Lentisphaera", Genus, ignore.case = TRUE))
)


tt <- as.data.frame(tax_table(ps_rel_all))
check <- tt %>%
  mutate(
    match_order  = !is.na(Order)  & grepl("Lentisphaer", Order,  ignore.case = TRUE),
    match_family = !is.na(Family) & grepl("Lentisphaer", Family, ignore.case = TRUE),
    match_genus  = !is.na(Genus)  & grepl("Lentisphaera", Genus, ignore.case = TRUE),
    selected = match_order | match_family | match_genus
  ) %>%
  filter(selected)
table(check$match_order, check$match_family, check$match_genus)

df <- psmelt(sub_Lenti_rel)

summary_lenti <- df %>%
  group_by(OTU) %>%
  summarise(
    mean_rel_abund_total = mean(Abundance, na.rm = TRUE),  # proportion du total
    prevalence = mean(Abundance > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_rel_abund_total))

asv_id <- "ASV_00001"

stats_asv <- psmelt(ps_rel_all) %>%
  filter(OTU == asv_id) %>%
  summarise(
    mean_rel_abund = mean(Abundance, na.rm = TRUE),
    sd_rel_abund   = sd(Abundance, na.rm = TRUE),
    prevalence     = mean(Abundance > 0, na.rm = TRUE),
    n_samples      = n()
  )

stats_asv

# Lentisphaera ASVs
lenti_ids <- taxa_names(sub_Lenti_rel)   
length(lenti_ids)

fa_all <- readDNAStringSet("2025_cocoon.fasta")

sum(names(fa_all) %in% lenti_ids)

fa_lenti <- fa_all[names(fa_all) %in% lenti_ids]
writeXStringSet(fa_lenti, "Lentisphaera_like_ASV_53_for_BLAST.fasta")

ids <- c("ASV_13827", "ASV_04020")

abund_2 <- psmelt(ps_rel_all) %>%
  filter(OTU %in% ids) %>%
  group_by(OTU) %>%
  summarise(
    mean_rel_abund = mean(Abundance, na.rm = TRUE),
    sd_rel_abund   = sd(Abundance, na.rm = TRUE),
    prevalence     = mean(Abundance > 0, na.rm = TRUE),
    .groups = "drop"
  )

abund_2


ps_rel <- transform_sample_counts(cocoon_r, function(x) x / sum(x))

tt <- as.data.frame(tax_table(ps_rel))
lenti_ids <- rownames(tt %>% 
                        mutate(
                          keep = (!is.na(Order)  & grepl("Lentisphaer", Order,  ignore.case = TRUE)) |
                            (!is.na(Family) & grepl("Lentisphaer", Family, ignore.case = TRUE)) |
                            (!is.na(Genus)  & grepl("Lentisphaera", Genus, ignore.case = TRUE))
                        ) %>% filter(keep))

length(lenti_ids)   

df <- psmelt(ps_rel) %>%
  filter(OTU %in% lenti_ids)

by_type <- df %>%
  group_by(Sample, what) %>%
  summarise(lenti_rel_abund = sum(Abundance), .groups = "drop") %>%
  group_by(what) %>%
  summarise(
    mean_rel_abund = mean(lenti_rel_abund),
    sd_rel_abund   = sd(lenti_rel_abund),
    median_rel_abund = median(lenti_rel_abund),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_rel_abund))

by_type
