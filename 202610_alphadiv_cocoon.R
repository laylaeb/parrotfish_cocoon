library(phyloseq)
library(dplyr)
library(cowplot)
library(fishualize)
library(MASS)
library(performance)
library(dunn.test)
library(MiscMetabar)
library(hillR)

## ---- Alpha diversity ---- ##
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

## ---- Stats ---- ##
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
