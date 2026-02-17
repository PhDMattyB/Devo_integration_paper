##############################
## MUTLIVATIATE REACTION NORM MADNESS
##
## Matt Brachmann (PhDMattyB)
##
## 11.02.2026
##
##############################

setwd("~/Parsons_postdoc/Stickleback_Morphometric_data/Updated Landmarks/")

library(tidyverse)
library(lme4)
library(lmerTest)
library(brms)
library(paran)


read_csv("F2_Original_univariate_traits_FIXED_11.02.2026.csv") %>% 
  dplyr::select(Lake_morph, 
                Ecotype_pair, 
                lake_morph_Pair_Full_Temp, 
                Morph, 
                Offspring_temp, 
                Parent_temp, 
                jaw_length:caudal2_15_17) %>% View()
  write_csv('F2_Orig_traits_scaled_formatted.csv', 
            col_names = F)
F2_orig_traits = read_csv('F2_Orig_traits_scaled_formatted.csv', 
         col_names = c('Lake_morph', 
                       'Ecotype_pair', 
                       'lake_morph_pair_full_temp', 
                       'morph', 
                       'offspring_temp', 
                       'parent_temp', 
                       paste0("trait_", 1:33))) %>% 
  separate(col = lake_morph_pair_full_temp, 
           into = c('trash', 
                    'full_fac'), 
           sep = '_') %>%
  pivot_longer(
    cols = starts_with("trait_"),
    names_to = "trait",
    values_to = "value"
  ) %>%
  mutate(
    ecotype = factor(morph, levels = c("Cold", "Warm")),  # adjust labels
    offspring_temp = factor(offspring_temp, levels = c(12, 18)),
    population = factor(Ecotype_pair)
  )

F2_orig_traits$offspring_temp = as.character(F2_orig_traits$offspring_temp)
F2_orig_traits$parent_temp = as.character(F2_orig_traits$parent_temp)

## F2 plasticity
ggplot(F2_orig_traits, 
       aes(x = offspring_temp,
               y = value,
               color = morph,
               group = morph)) +
  
  stat_summary(fun = mean, geom = "line", linewidth = 0.6, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 1.5) +
  
  facet_wrap(~ trait) +
  
  scale_color_manual(values = c("Cold" = "#1b9e77",
                                "Warm" = "#d95f02")) +
  
  labs(x = "Temperature (°C)",
       y = "Standardized trait value (z-score)",
       color = "Ecotype") +
  
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "top"
  )

# ggplot(F2_orig_traits,
#        aes(x = offspring_temp,
#            y = value,
#            color = ecotype,
#            group = interaction(Lake_morph, trait))) +
#   
#   geom_line(alpha = 0.15) +   # individual plasticity
#   
#   stat_summary(aes(group = interaction(Lake_morph)),
#                fun = mean, geom = "line", linewidth = 0.9) +
#   
#   facet_wrap(~ trait, scales = "free_y", ncol = 6) +
#   theme_bw()

## F1 plasticity

ggplot(F2_orig_traits, 
       aes(x = parent_temp,
           y = value,
           color = morph,
           group = morph)) +
  
  stat_summary(fun = mean, geom = "line", linewidth = 0.6, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 1.5) +
  
  facet_wrap(~ trait) +
  
  scale_color_manual(values = c("Cold" = "#1b9e77",
                                "Warm" = "#d95f02")) +
  
  labs(x = "Temperature (°C)",
       y = "Standardized trait value (z-score)",
       color = "Ecotype") +
  
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "top"
  )

ggplot(F2_orig_traits,
       aes(x = full_fac,
           y = value,
           color = morph,
           group = morph)) +
  
  stat_summary(fun = mean, geom = "line", linewidth = 0.6, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 1.5) +
  
  facet_wrap(~ trait) +
  
  scale_color_manual(values = c("Cold" = "#1b9e77",
                                "Warm" = "#d95f02")) +
  
  labs(x = "Temperature (°C)",
       y = "Standardized trait value (z-score)",
       color = "Ecotype") +
  
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "top"
  )

# WGP reaction norm plasticity --------------------------------------------
read_csv("WGP_TRAITS_SCALED_FIXED_11.02.2026.csv")%>% 
  dplyr::select(Lake_morph, 
                Ecotype_pair, 
                lake_morph_Pair_Full_Temp, 
                Morph, 
                Offspring_temp, 
                Parent_temp, 
                jaw_length:caudal2_15_17) %>% 
write_csv('WGP_traits_scaled_formatted.csv', 
          col_names = F)

WGP_traits = read_csv('WGP_traits_scaled_formatted.csv', 
                          col_names = c('Lake_morph', 
                                        'Ecotype_pair', 
                                        'lake_morph_pair_full_temp', 
                                        'morph', 
                                        'offspring_temp', 
                                        'parent_temp', 
                                        paste0("trait_", 1:33))) %>%
  separate(col = lake_morph_pair_full_temp, 
           into = c('trash', 
                    'full_fac'), 
           sep = '_') %>%   
    pivot_longer(
    cols = starts_with("trait_"),
    names_to = "trait",
    values_to = "value"
  ) %>%
  mutate(
    ecotype = factor(morph, levels = c("Cold", "Warm")),  # adjust labels
    offspring_temp = factor(offspring_temp, levels = c(12, 18)),
    population = factor(Ecotype_pair), 
    full_fac = factor(full_fac)
  )

WGP_traits$offspring_temp = as.character(WGP_traits$offspring_temp)
WGP_traits$parent_temp = as.character(WGP_traits$parent_temp)

## F2 plasticity
ggplot(WGP_traits, 
       aes(x = offspring_temp,
           y = value,
           color = morph,
           group = morph)) +
  
  stat_summary(fun = mean, geom = "line", linewidth = 0.6, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 1.5) +
  
  facet_wrap(~ trait) +
  
  scale_color_manual(values = c("Cold" = "#1b9e77",
                                "Warm" = "#d95f02")) +
  
  labs(x = "Temperature (°C)",
       y = "Standardized trait value (z-score)",
       color = "Ecotype") +
  
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "top"
  )


## F1 plasticity

ggplot(WGP_traits, 
       aes(x = parent_temp,
           y = value,
           color = morph,
           group = morph)) +
  
  stat_summary(fun = mean, geom = "line", linewidth = 0.6, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 1.5) +
  
  facet_wrap(~ trait) +
  
  scale_color_manual(values = c("Cold" = "#1b9e77",
                                "Warm" = "#d95f02")) +
  
  labs(x = "Temperature (°C)",
       y = "Standardized trait value (z-score)",
       color = "Ecotype") +
  
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "top"
  )


## full fac plasticity

ggplot(WGP_traits, 
       aes(x = full_fac,
           y = value,
           color = morph,
           group = morph)) +
  
  stat_summary(fun = mean, geom = "line", linewidth = 0.6, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 1.5) +
  
  facet_wrap(~ trait) +
  
  scale_color_manual(values = c("Cold" = "#1b9e77",
                                "Warm" = "#d95f02")) +
  
  labs(x = "Temperature (°C)",
       y = "Standardized trait value (z-score)",
       color = "Ecotype") +
  
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "top"
  )




# TGP reaction norm plasticity --------------------------------------------


read_csv('TGP_TRAITS_SCALED_FIXED_11.02.2026')%>% 
  dplyr::select(Lake_morph, 
                Ecotype_pair, 
                lake_morph_Pair_Full_Temp, 
                Morph, 
                Offspring_temp, 
                Parent_temp, 
                jaw_length:caudal2_15_17) %>% 
  write_csv('TGP_traits_scaled_formatted.csv', 
            col_names = F)

TGP_traits = read_csv('TGP_traits_scaled_formatted.csv', 
                      col_names = c('Lake_morph', 
                                    'Ecotype_pair', 
                                    'lake_morph_pair_full_temp', 
                                    'morph', 
                                    'offspring_temp', 
                                    'parent_temp', 
                                    paste0("trait_", 1:33))) %>%
  separate(col = lake_morph_pair_full_temp, 
           into = c('trash', 
                    'full_fac'), 
           sep = '_') %>%   
  pivot_longer(
    cols = starts_with("trait_"),
    names_to = "trait",
    values_to = "value"
  ) %>%
  mutate(
    ecotype = factor(morph, levels = c("Cold", "Warm")),  # adjust labels
    offspring_temp = factor(offspring_temp, levels = c(12, 18)),
    population = factor(Ecotype_pair), 
    full_fac = factor(full_fac)
  )

TGP_traits$offspring_temp = as.character(TGP_traits$offspring_temp)
TGP_traits$parent_temp = as.character(TGP_traits$parent_temp)

## F2 plasticity
ggplot(TGP_traits, 
       aes(x = offspring_temp,
           y = value,
           color = morph,
           group = morph)) +
  
  stat_summary(fun = mean, geom = "line", linewidth = 0.6, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 1.5) +
  
  facet_wrap(~ trait) +
  
  scale_color_manual(values = c("Cold" = "#1b9e77",
                                "Warm" = "#d95f02")) +
  
  labs(x = "Temperature (°C)",
       y = "Standardized trait value (z-score)",
       color = "Ecotype") +
  
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "top"
  )


## F1 plasticity

ggplot(TGP_traits, 
       aes(x = parent_temp,
           y = value,
           color = morph,
           group = morph)) +
  
  stat_summary(fun = mean, geom = "line", linewidth = 0.6, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 1.5) +
  
  facet_wrap(~ trait) +
  
  scale_color_manual(values = c("Cold" = "#1b9e77",
                                "Warm" = "#d95f02")) +
  
  labs(x = "Temperature (°C)",
       y = "Standardized trait value (z-score)",
       color = "Ecotype") +
  
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "top"
  )


## full fac plasticity

ggplot(TGP_traits, 
       aes(x = full_fac,
           y = value,
           color = morph,
           group = morph)) +
  
  stat_summary(fun = mean, geom = "line", linewidth = 0.6, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 1.5) +
  
  facet_wrap(~ trait) +
  
  scale_color_manual(values = c("Cold" = "#1b9e77",
                                "Warm" = "#d95f02")) +
  
  labs(x = "Temperature (°C)",
       y = "Standardized trait value (z-score)",
       color = "Ecotype") +
  
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "top"
  )



# PCA approach ------------------------------------------------------------


# PCA F2 orginal traits z-scores ----------------------------------------------

F2_data = read_csv("F2_Original_univariate_traits_FIXED_11.02.2026.csv") %>% 
  dplyr::select(Lake_morph, 
                Ecotype_pair, 
                lake_morph_Pair_Full_Temp, 
                Morph, 
                Offspring_temp, 
                Parent_temp, 
                jaw_length:caudal2_15_17)


F2_orig_trait_mat <- F2_data %>%
  dplyr::select(jaw_length:caudal2_15_17) %>% 
  # dplyr::select(starts_with('value'))
  # dplyr::select(all_of(trait_names)) %>%
  as.matrix()

pca <- prcomp(F2_orig_trait_mat, 
              center = F, 
              scale. = F)

paran(F2_orig_trait_mat, iterations = 1000, graph = TRUE)

scores <- as.data.frame(pca$x[, 1:7])
F2_PCA <- bind_cols(F2_data, scores)

F2_off_means <- F2_PCA %>%
  group_by(Morph,
           Offspring_temp) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    .groups = "drop"
  )


PC_strong_loadings = pca$rotation %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  rename(traits = rowname) %>% 
  as_tibble() %>%
  dplyr::select(traits, 
                PC1, 
                PC2)%>% 
  pivot_longer(cols = c(PC1, PC2),
               names_to = "PC",
               values_to = "loading")

load_cols = c('#005f73', 
              '#ca6702')
trait_loading_rank = ggplot(PC_strong_loadings,
       aes(x = reorder(traits, loading),
           y = loading,
           fill = PC)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = load_cols)+
  facet_wrap(~ PC, scales = "free_y") +
  coord_flip() +
  theme_bw() +
  labs(x = "Trait", y = "Loading", 
       title = 'PC loadings')+
  theme(strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'))

PC_trait_loadings = pca$rotation %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  rename(traits = rowname) %>% 
  as_tibble() %>%
  dplyr::select(traits, 
                PC1, 
                PC2)

trait_loading_gg = ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(data = PC_trait_loadings, 
             aes(PC1, PC2), 
             size = 3)+
  geom_text_repel(data = PC_trait_loadings, 
                  aes(x = PC1, 
                      y = PC2, 
                      label = traits), 
                  size = 4, 
                  max.overlaps = 50)+
  theme_bw() +
  labs(x = "PC1 loading", y = "PC2 loading", 
       title = 'Trait loadings')


F2_off_means$Offspring_temp = as.character(F2_off_means$Offspring_temp)

F2_off_vectors <- F2_off_means %>%
  pivot_wider(
    names_from = Offspring_temp,
    values_from = c(PC1, PC2)
  )

F2_off_dashed_vec = F2_off_means %>%
  pivot_wider(
    names_from = Morph,
    values_from = c(PC1, PC2)
  )

F2_parent_means <- F2_PCA %>%
  group_by(Morph,
           Parent_temp) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    .groups = "drop"
  )
F2_parent_means$Parent_temp = as.character(F2_parent_means$Parent_temp)

F2_par_vectors <- F2_parent_means %>%
  pivot_wider(
    names_from = Parent_temp,
    values_from = c(PC1, PC2)
  )
F2_parent_dashed_vec = F2_parent_means %>%
  pivot_wider(
    names_from = Morph,
    values_from = c(PC1, PC2)
  )


normy_cols = c('#3d348b',
               '#ff006e')


multivar_reaction_norm = ggplot() +
  geom_point(data = F2_off_means,
             aes(PC1,
                 PC2,
                 color = Morph,
                 shape = Offspring_temp),
             size = 4) +
  geom_segment(data = F2_off_vectors,
               aes(x = PC1_12, y = PC2_12,
                   xend = PC1_18, yend = PC2_18,
                   color = Morph),
               arrow = arrow(length = unit(0.25, "cm")),
               linewidth = 1.2) +
  xlim(-2, 2)+
  ylim(-2, 2)+
  geom_hline(yintercept = 0, 
             linetype = 'dashed')+
  geom_vline(xintercept = 0, 
             linetype = 'dashed')+
  labs(x = "PC1",
       y = "PC2", 
       title = 'Offspring temperature') +
  scale_color_manual(values = normy_cols)+
  theme_bw() +
  theme(legend.position = 'none')+

ggplot() +
  geom_point(data = F2_parent_means,
             aes(PC1,
                 PC2,
                 color = Morph,
                 shape = Parent_temp),
             size = 4) +
  geom_segment(data = F2_par_vectors,
               aes(x = PC1_12, y = PC2_12,
                   xend = PC1_18, yend = PC2_18,
                   color = Morph),
               arrow = arrow(length = unit(0.25, "cm")),
               linewidth = 1.2) +
  xlim(-2, 2)+
  ylim(-2, 2)+
  geom_hline(yintercept = 0, 
             linetype = 'dashed')+
  geom_vline(xintercept = 0, 
             linetype = 'dashed')+
  # 
  labs(x = "PC1",
       y = "PC2",
       color = "Ecotype", 
       title = 'Parental temperature') +
  scale_colour_manual(values = normy_cols)+
  theme_bw() +
  theme(legend.position = 'none')


Multivariate_reaction_plot = trait_loading_rank + trait_loading_gg / multivar_reaction_norm

# ggsave('Multivariate_reaction_norm_graph_12.02.2026.tiff', 
#        plot = Multivariate_reaction_plot, 
#        dpi = 'retina', 
#        units = 'cm', 
#        width = 50, 
#        height = 20)  
# 
# 
# ggsave('Multivariate_reaction_norm_graph_12.02.2026.svg', 
#        plot = Multivariate_reaction_plot, 
#        dpi = 'retina', 
#        units = 'cm', 
#        width = 50, 
#        height = 20)  

# F2 orginal plasticity model ---------------------------------------------

F2_pc1_mod = lmer(PC1 ~ Morph * Offspring_temp * Parent_temp + (1 | Ecotype_pair),
                  data = F2_PCA)

summary(F2_pc1_mod)


F2_pc2_mod = lmer(PC2 ~ Morph * Offspring_temp * Parent_temp + (1 | Ecotype_pair), 
                  data = F2_PCA)

summary(F2_pc2_mod)


F2_brms_fit <- brm(
  mvbind(PC1, PC2, PC3, PC4, PC5, PC6, PC7) ~ Morph * Offspring_temp * Parent_temp + (1 | Ecotype_pair),
  data = F2_PCA, 
  iter = 10000
)

manova_fit <- manova(cbind(PC1, PC2) ~ Morph * Offspring_temp * Parent_temp * Ecotype_pair, 
                     data = F2_PCA)
summary(manova_fit, test = "Pillai")


# WGP multivariate reaction norms -----------------------------------------


WGP_data = read_csv("WGP_TRAITS_SCALED_FIXED_11.02.2026.csv") %>% 
  dplyr::select(Lake_morph, 
                Ecotype_pair, 
                lake_morph_Pair_Full_Temp, 
                Morph, 
                Offspring_temp, 
                Parent_temp, 
                jaw_length:caudal2_15_17)


WGP_orig_trait_mat <- WGP_data %>%
  dplyr::select(jaw_length:caudal2_15_17) %>% 
  # dplyr::select(starts_with('value'))
  # dplyr::select(all_of(trait_names)) %>%
  as.matrix()

WGP_pca <- prcomp(WGP_orig_trait_mat, 
              center = F, 
              scale. = F)

paran(WGP_orig_trait_mat, iterations = 1000, graph = TRUE)

WGP_scores <- as.data.frame(WGP_pca$x[, 1:7])
WGP_PCA <- bind_cols(WGP_data, WGP_scores)

WGP_off_means <- WGP_PCA %>%
  group_by(Morph,
           Offspring_temp) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    .groups = "drop"
  )


WGP_PC_strong_loadings = WGP_pca$rotation %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  rename(traits = rowname) %>% 
  as_tibble() %>%
  dplyr::select(traits, 
                PC1, 
                PC2)%>% 
  pivot_longer(cols = c(PC1, PC2),
               names_to = "PC",
               values_to = "loading")

load_cols = c('#005f73', 
              '#ca6702')
WGP_trait_loading_rank = ggplot(WGP_PC_strong_loadings,
                            aes(x = reorder(traits, loading),
                                y = loading,
                                fill = PC)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = load_cols)+
  facet_wrap(~ PC, scales = "free_y") +
  coord_flip() +
  theme_bw() +
  labs(x = "Trait", y = "Loading", 
       title = 'PC loadings')+
  theme(strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'))

WGP_PC_trait_loadings = WGP_pca$rotation %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  rename(traits = rowname) %>% 
  as_tibble() %>%
  dplyr::select(traits, 
                PC1, 
                PC2)

WGP_trait_loading_gg = ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(data = WGP_PC_trait_loadings, 
             aes(PC1, PC2), 
             size = 3)+
  geom_text_repel(data = WGP_PC_trait_loadings, 
                  aes(x = PC1, 
                      y = PC2, 
                      label = traits), 
                  size = 4, 
                  max.overlaps = 50)+
  theme_bw() +
  labs(x = "PC1 loading", y = "PC2 loading", 
       title = 'Trait loadings')


WGP_off_means$Offspring_temp = as.character(WGP_off_means$Offspring_temp)

WGP_off_vectors <- WGP_off_means %>%
  pivot_wider(
    names_from = Offspring_temp,
    values_from = c(PC1, PC2)
  )

WGP_off_dashed_vec = WGP_off_means %>%
  pivot_wider(
    names_from = Morph,
    values_from = c(PC1, PC2)
  )

WGP_parent_means <- WGP_PCA %>%
  group_by(Morph,
           Parent_temp) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    .groups = "drop"
  )
WGP_parent_means$Parent_temp = as.character(WGP_parent_means$Parent_temp)

WGP_par_vectors <- WGP_parent_means %>%
  pivot_wider(
    names_from = Parent_temp,
    values_from = c(PC1, PC2)
  )
WGP_parent_dashed_vec = WGP_parent_means %>%
  pivot_wider(
    names_from = Morph,
    values_from = c(PC1, PC2)
  )


normy_cols = c('#3d348b',
               '#ff006e')


WGP_multivar_reaction_norm = ggplot() +
  geom_point(data = WGP_off_means,
             aes(PC1,
                 PC2,
                 color = Morph,
                 shape = Offspring_temp),
             size = 4) +
  geom_segment(data = WGP_off_vectors,
               aes(x = PC1_12, y = PC2_12,
                   xend = PC1_18, yend = PC2_18,
                   color = Morph),
               arrow = arrow(length = unit(0.25, "cm")),
               linewidth = 1.2) +
  xlim(-2, 2)+
  ylim(-2, 2)+
  geom_hline(yintercept = 0, 
             linetype = 'dashed')+
  geom_vline(xintercept = 0, 
             linetype = 'dashed')+
  labs(x = "PC1",
       y = "PC2", 
       title = 'Offspring temperature') +
  scale_color_manual(values = normy_cols)+
  theme_bw() +
  theme(legend.position = 'none')+
  
  ggplot() +
  geom_point(data = WGP_parent_means,
             aes(PC1,
                 PC2,
                 color = Morph,
                 shape = Parent_temp),
             size = 4) +
  geom_segment(data = WGP_par_vectors,
               aes(x = PC1_12, y = PC2_12,
                   xend = PC1_18, yend = PC2_18,
                   color = Morph),
               arrow = arrow(length = unit(0.25, "cm")),
               linewidth = 1.2) +
  xlim(-2, 2)+
  ylim(-2, 2)+
  geom_hline(yintercept = 0, 
             linetype = 'dashed')+
  geom_vline(xintercept = 0, 
             linetype = 'dashed')+
  # 
  labs(x = "PC1",
       y = "PC2",
       color = "Ecotype", 
       title = 'Parental temperature') +
  scale_colour_manual(values = normy_cols)+
  theme_bw() +
  theme(legend.position = 'none')


WGP_Multivariate_reaction_plot = WGP_trait_loading_rank + WGP_trait_loading_gg / WGP_multivar_reaction_norm

ggsave('WGP_Multivariate_reaction_norm_graph_12.02.2026.tiff', 
       plot = WGP_Multivariate_reaction_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 50, 
       height = 20)  


ggsave('WGP_Multivariate_reaction_norm_graph_12.02.2026.svg', 
       plot = WGP_Multivariate_reaction_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 50, 
       height = 20)  

# WGP plasticity model ---------------------------------------------

WGP_pc1_mod = lmer(PC1 ~ Morph * Offspring_temp * Parent_temp + (1 | Ecotype_pair),
                  data = WGP_PCA)

summary(WGP_pc1_mod)


WGP_pc2_mod = lmer(PC2 ~ Morph * Offspring_temp * Parent_temp + (1 | Ecotype_pair), 
                  data = WGP_PCA)

summary(WGP_pc2_mod)


WGP_brms_fit <- brm(
  mvbind(PC1, PC2, PC3, PC4, PC5, PC6, PC7) ~ Morph * Offspring_temp * Parent_temp + (1 | Ecotype_pair),
  data = WGP_PCA, 
  iter = 10000
)

manova_fit <- manova(cbind(PC1, PC2) ~ Morph * Offspring_temp * Parent_temp * Ecotype_pair, 
                     data = WGP_PCA)
summary(manova_fit, test = "Pillai")


# TGP multivariate reaction norms -----------------------------------------


TGP_data = read_csv("TGP_TRAITS_SCALED_FIXED_11.02.2026.csv") %>% 
  dplyr::select(Lake_morph, 
                Ecotype_pair, 
                lake_morph_Pair_Full_Temp, 
                Morph, 
                Offspring_temp, 
                Parent_temp, 
                jaw_length:caudal2_15_17)


TGP_orig_trait_mat <- TGP_data %>%
  dplyr::select(jaw_length:caudal2_15_17) %>% 
  # dplyr::select(starts_with('value'))
  # dplyr::select(all_of(trait_names)) %>%
  as.matrix()

TGP_pca <- prcomp(TGP_orig_trait_mat, 
              center = F, 
              scale. = F)

paran(TGP_orig_trait_mat, iterations = 1000, graph = TRUE)

TGP_scores <- as.data.frame(TGP_pca$x[, 1:7])
TGP_PCA <- bind_cols(TGP_data, TGP_scores)

TGP_off_means <- TGP_PCA %>%
  group_by(Morph,
           Offspring_temp) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    .groups = "drop"
  )


TGP_PC_strong_loadings = TGP_pca$rotation %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  rename(traits = rowname) %>% 
  as_tibble() %>%
  dplyr::select(traits, 
                PC1, 
                PC2)%>% 
  pivot_longer(cols = c(PC1, PC2),
               names_to = "PC",
               values_to = "loading")

load_cols = c('#005f73', 
              '#ca6702')
TGP_trait_loading_rank = ggplot(TGP_PC_strong_loadings,
                            aes(x = reorder(traits, loading),
                                y = loading,
                                fill = PC)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = load_cols)+
  facet_wrap(~ PC, scales = "free_y") +
  coord_flip() +
  theme_bw() +
  labs(x = "Trait", y = "Loading", 
       title = 'PC loadings')+
  theme(strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'))

TGP_PC_trait_loadings = TGP_pca$rotation %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  rename(traits = rowname) %>% 
  as_tibble() %>%
  dplyr::select(traits, 
                PC1, 
                PC2)

TGP_trait_loading_gg = ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(data = TGP_PC_trait_loadings, 
             aes(PC1, PC2), 
             size = 3)+
  geom_text_repel(data = TGP_PC_trait_loadings, 
                  aes(x = PC1, 
                      y = PC2, 
                      label = traits), 
                  size = 4, 
                  max.overlaps = 50)+
  theme_bw() +
  labs(x = "PC1 loading", y = "PC2 loading", 
       title = 'Trait loadings')


TGP_off_means$Offspring_temp = as.character(TGP_off_means$Offspring_temp)

TGP_off_vectors <- TGP_off_means %>%
  pivot_wider(
    names_from = Offspring_temp,
    values_from = c(PC1, PC2)
  )

TGP_off_dashed_vec = TGP_off_means %>%
  pivot_wider(
    names_from = Morph,
    values_from = c(PC1, PC2)
  )

TGP_parent_means <- TGP_PCA %>%
  group_by(Morph,
           Parent_temp) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    .groups = "drop"
  )
TGP_parent_means$Parent_temp = as.character(TGP_parent_means$Parent_temp)

TGP_par_vectors <- TGP_parent_means %>%
  pivot_wider(
    names_from = Parent_temp,
    values_from = c(PC1, PC2)
  )

# TGP_PCA_means = TGP_PCA %>%
#   group_by(Morph,
#            Offspring_temp,
#            Parent_temp) %>%
#   summarise(
#     PC1 = mean(PC1),
#     PC2 = mean(PC2),
#     .groups = "drop"
#   ) %>%
#   unite(col = 'temp_rearing',
#         c('Offspring_temp',
#           'Parent_temp'),
#         sep = '_',
#         remove = F)
# 
# vectors <- TGP_PCA_means %>%
#   pivot_wider(
#     names_from = temp_rearing,
#     values_from = c(PC1, PC2)
#   )

normy_cols = c('#3d348b',
               '#ff006e')


TGP_multivar_reaction_norm = ggplot() +
  geom_point(data = TGP_off_means,
             aes(PC1,
                 PC2,
                 color = Morph,
                 shape = Offspring_temp),
             size = 4) +
  geom_segment(data = TGP_off_vectors,
               aes(x = PC1_12, y = PC2_12,
                   xend = PC1_18, yend = PC2_18,
                   color = Morph),
               arrow = arrow(length = unit(0.25, "cm")),
               linewidth = 1.2) +
  xlim(-2, 2)+
  ylim(-2, 2)+
  geom_hline(yintercept = 0, 
             linetype = 'dashed')+
  geom_vline(xintercept = 0, 
             linetype = 'dashed')+
  labs(x = "PC1",
       y = "PC2", 
       title = 'Offspring temperature') +
  scale_color_manual(values = normy_cols)+
  theme_bw() +
  theme(legend.position = 'none')+
  
  ggplot() +
  geom_point(data = TGP_parent_means,
             aes(PC1,
                 PC2,
                 color = Morph,
                 shape = Parent_temp),
             size = 4) +
  geom_segment(data = TGP_par_vectors,
               aes(x = PC1_12, y = PC2_12,
                   xend = PC1_18, yend = PC2_18,
                   color = Morph),
               arrow = arrow(length = unit(0.25, "cm")),
               linewidth = 1.2) +
  xlim(-2, 2)+
  ylim(-2, 2)+
  geom_hline(yintercept = 0, 
             linetype = 'dashed')+
  geom_vline(xintercept = 0, 
             linetype = 'dashed')+
  # 
  labs(x = "PC1",
       y = "PC2",
       color = "Ecotype", 
       title = 'Parental temperature') +
  scale_colour_manual(values = normy_cols)+
  theme_bw() +
  theme(legend.position = 'none')


TGP_Multivariate_reaction_plot = TGP_trait_loading_rank + TGP_trait_loading_gg / TGP_multivar_reaction_norm

ggsave('TGP_Multivariate_reaction_norm_graph_12.02.2026.tiff', 
       plot = TGP_Multivariate_reaction_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 50, 
       height = 20)  


ggsave('TGP_Multivariate_reaction_norm_graph_12.02.2026.svg', 
       plot = TGP_Multivariate_reaction_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 50, 
       height = 20)  



# TGP plasticity model ---------------------------------------------

TGP_pc1_mod = lmer(PC1 ~ Morph * Offspring_temp * Parent_temp + (1 | Ecotype_pair),
                   data = TGP_PCA)

summary(TGP_pc1_mod)


TGP_pc2_mod = lmer(PC2 ~ Morph * Offspring_temp * Parent_temp + (1 | Ecotype_pair), 
                   data = TGP_PCA)

summary(TGP_pc2_mod)


TGP_brms_fit <- brm(
  mvbind(PC1, PC2, PC3, PC4, PC5, PC6, PC7) ~ Morph * Offspring_temp * Parent_temp + (1 | Ecotype_pair),
  data = TGP_PCA, 
  iter = 10000
)

manova_fit <- manova(cbind(PC1, PC2) ~ Morph * Offspring_temp * Parent_temp * Ecotype_pair, 
                     data = TGP_PCA)
summary(manova_fit, test = "Pillai")


# BRMS plots --------------------------------------------------------------
library(broom.mixed)

F2_brms_fit_sum = tidy(F2_brms_fit, 
                       effects = c('fixed', 
                                   'ran_pars', 
                                   'ran_vals'), 
                       conf.int = T) %>% 
  mutate(dataset = 'F2 original')
WGP_brms_fit_sum = tidy(WGP_brms_fit, 
                        effects = c('fixed', 
                                    'ran_pars', 
                                    'ran_vals'), 
                        conf.int = T) %>% 
  mutate(dataset = 'Within-generational plasticity')
TGP_brms_fit_sum = tidy(TGP_brms_fit, 
                        effects = c('fixed', 
                                    'ran_pars', 
                                    'ran_vals'), 
                        conf.int = T) %>% 
  mutate(dataset = 'Trans-generational plasticity')

brms_fits_full = bind_rows(F2_brms_fit_sum, 
                           WGP_brms_fit_sum, 
                           TGP_brms_fit_sum)

# brms_fits_full %>% 
#   write_csv('Multivariate_plast_BRMS_results.csv')


