##------------------------------------------------------------ Indicator species 

library(phyloseq)
library(indicspecies)

# ============================================================================
# ANALYSES PAR SITE
# ============================================================================

cocoon_r_bt <- subset_samples(cocoon_r, site %in% "BT")
cocoon_r_g <- subset_samples(cocoon_r, site %in% "G")

cocoon_new_physeq <- prune_taxa(taxa_sums(cocoon_r) > 10, cocoon_r)
tax_genus_relabund <- transform_sample_counts(cocoon_new_physeq, function(x) x / sum(x))

group <- sample_data(tax_genus_relabund)$what
otu <- as.data.frame(t(otu_table(tax_genus_relabund)))

indval_global <- multipatt(otu, group, func = "IndVal.g", 
                           duleg = TRUE, control = ctrl_global)

summary(indval_global, indvalcomp = TRUE, alpha = 0.05)
write.csv(res_global, "202602_indicator_cocoon_ASV.csv")

# Site BT
tax_genus_bt <- tax_glom(cocoon_r_bt, taxrank=rank_names(cocoon_r_bt)[6], NArm=F)
cocoon_new_physeq_bt <- prune_taxa(taxa_sums(tax_genus_bt) > 5, tax_genus_bt)
tax_genus_relabund_bt <- transform_sample_counts(cocoon_new_physeq_bt, function(x) x / sum(x))

otu_bt <- as.data.frame(t(otu_table(tax_genus_relabund_bt)))
group_bt <- sample_data(tax_genus_relabund_bt)$what

# Site G
tax_genus_g <- tax_glom(cocoon_r_g, taxrank=rank_names(cocoon_r_g)[6], NArm=F)
cocoon_new_physeq_g <- prune_taxa(taxa_sums(tax_genus_g) > 5, tax_genus_g)
tax_genus_relabund_g <- transform_sample_counts(cocoon_new_physeq_g, function(x) x / sum(x))

otu_g <- as.data.frame(t(otu_table(tax_genus_relabund_g)))
group_g <- sample_data(tax_genus_relabund_g)$what

set.seed(123)
indval_res_g <- multipatt(otu_g, group_g, func = "IndVal.g", duleg = T, control = how(nperm=999))
indval_res_bt <- multipatt(otu_bt, group_bt, func = "IndVal.g", duleg = T, control = how(nperm=999))

res_g <- summary(indval_res_g, indvalcomp = TRUE, alpha = 0.05)
res_bt <- summary(indval_res_bt, indvalcomp = TRUE, alpha = 0.05)

# ============================================================================
# ANALYSE AVEC BLOCKS 
# ============================================================================

tax_genus <- tax_glom(cocoon_r, taxrank=rank_names(cocoon_r)[6], NArm=F)
cocoon_new_physeq <- prune_taxa(taxa_sums(tax_genus) > 10, tax_genus)
tax_genus_relabund <- transform_sample_counts(cocoon_new_physeq, function(x) x / sum(x))

group <- sample_data(tax_genus_relabund)$what
otu <- as.data.frame(t(otu_table(tax_genus_relabund)))

ctrl_controlled <- how(nperm = 999)
setBlocks(ctrl_controlled) <- cocoon_metadata$site  # Permutations restreintes par site

indval_controlled <- multipatt(otu, group, func = "IndVal.g", 
                               duleg = TRUE, control = ctrl_controlled)

res_controlled <- summary(indval_controlled, indvalcomp = TRUE, alpha = 0.05)

# ============================================================================
# ANALYSE GLOBALE SANS BLOCKS 
# ============================================================================

ctrl_global <- how(nperm = 999)
# PAS de setBlocks() ici !

indval_global <- multipatt(otu, group, func = "IndVal.g", 
                           duleg = TRUE, control = ctrl_global)

res_global <- summary(indval_global, indvalcomp = TRUE, alpha = 0.05)

taxa_df <- as.data.frame(tax_table(cocoon_new_physeq))
taxa_df$ASV <- rownames(taxa_df)

indicators_both <- intersect(final_sig_controlled$ASV, final_sig_global$ASV)
cat("Indicateurs dans les DEUX analyses (robustes):", length(indicators_both), "\n")

indicators_only_global <- setdiff(final_sig_global$ASV, final_sig_controlled$ASV)
cat("Indicateurs seulement dans l'analyse globale:", length(indicators_only_global), "\n")

indicators_only_controlled <- setdiff(final_sig_controlled$ASV, final_sig_global$ASV)
cat("Indicateurs seulement dans l'analyse contrôlée:", length(indicators_only_controlled), "\n\n")

##-------------------------------------------------Rank abundance curves according to sample types
comm<-as(otu_table(cocoon_r), "matrix")
env<-as(sample.data.frame(cocoon_r), "matrix")

comm<-as.data.frame(comm)
env<-as.data.frame(env)

comm_melt <- reshape2::melt(as.matrix(comm), na.rm=TRUE)
colnames(comm_melt)<-c("Var1","sample","Abundance")
head(comm_melt)
comm_melt<-as.data.frame(comm_melt)
comm_melt_scaling<-comm_melt[comm_melt$Abundance!=0,] 
head(comm_melt_scaling)

comm_melt_RA<-comm_melt_scaling %>%
  group_by(sample) %>%
  mutate(Scale_RA=Abundance/sum(Abundance)) 
head(comm_melt_RA)

##Sanity check
comm_melt_RA %>%
  group_by(sample)%>% summarize(somme=sum(Scale_RA))

comm_melt_RANK<-comm_melt_RA %>%
  group_by(sample) %>%
  mutate(my_ranks = order(order(Scale_RA,decreasing=T))) #Ranks

Rank_merge<-merge(comm_melt_RANK, env, by.x="sample",by.y="Sample")

Rank_merge_mod <- expand.grid(sample_type = unique(Rank_merge$what),
                              Rank = 1:max(Rank_merge$my_ranks))

Rank_merge_mod$Scale_RA = vapply(1:nrow(Rank_merge_mod), 
                                 function(i) mean(Rank_merge$Scale_RA[(Rank_merge$my_ranks == Rank_merge_mod$Rank[i]) & 
                                                                        (Rank_merge$what == Rank_merge_mod$sample_type[i])]), FUN.VALUE = numeric(1))

library(fishualize)

palette_colors <- fish(option = "Oncorhynchus_kisutch", 5, direction = -1)
names(palette_colors) <- c("envelope", "fish_skin", "mouth", "sediment", "water")

palette_colors[c(1, 5)] <- palette_colors[c(5, 1)]

rare_curve <- ggplot() +
  geom_line(data = Rank_merge, aes(x = my_ranks, y = log(Scale_RA), group = sample), size = 0.2, colour = "black", alpha = 0.2) +
  geom_line(data = Rank_merge_mod, aes(x = Rank, y = log(Scale_RA), colour = sample_type)) +
  scale_color_manual(values = palette_colors) +
  xlab("Rank") +
  facet_wrap(~ sample_type, nrow = 5, ncol = 2) +
  theme(panel.background = element_rect(fill = 'white', color = "grey60"))
