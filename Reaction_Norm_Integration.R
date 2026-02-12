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

scores <- as.data.frame(pca$x[, 1:2])
F2_PCA <- bind_cols(F2_data, scores)

F2_off_means <- F2_PCA %>%
  group_by(Morph,
           Offspring_temp) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    .groups = "drop"
  )

# PC1_strong_loading = pca$rotation %>% 
#   as.data.frame() %>% 
#   rownames_to_column() %>% 
#   rename(traits = rowname) %>% 
#   as_tibble() %>%
#   dplyr::select(traits, 
#                 PC1, 
#                 PC2)%>%
#   mutate(strong_PC1 = abs(PC1) > 0.2,
#          strong_PC2 = abs(PC2) > 0.2) %>%
#   filter(strong_PC1 == 'TRUE') 
# 
# PC2_strong_loading = pca$rotation %>% 
#   as.data.frame() %>% 
#   rownames_to_column() %>% 
#   rename(traits = rowname) %>% 
#   as_tibble() %>%
#   dplyr::select(traits, 
#                 PC1, 
#                 PC2)%>%
#   mutate(strong_PC1 = abs(PC1) > 0.2,
#          strong_PC2 = abs(PC2) > 0.2) %>%
#   filter(strong_PC2 == 'TRUE')
# 
# 
# PC_strong_loadings = bind_rows(PC1_strong_loading, 
#                                PC2_strong_loading) %>% 
#   pivot_longer(cols = c(PC1, PC2),
#                names_to = "PC",
#                values_to = "loading")

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
  # xlim(-2, 2)+
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = load_cols)+
  facet_wrap(~ PC, scales = "free_y") +
  # xlim(-2, 2)+
  # ylim(-2, 2)+
  coord_flip() +
  theme_bw() +
  labs(x = "Trait", y = "Loading", 
       title = 'A) PC loadings')+
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
# %>% 
#   pivot_longer(cols = c(PC1, PC2),
#                names_to = "PC",
#                values_to = "loading")



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
  # geom_point(data = PC1_strong_loading,
  #            aes(PC1, PC2),
  #            size = 2) +
  # geom_text_repel(data = PC1_strong_loading, 
  #                 aes(x = PC1, 
  #                     y = PC2, 
  #                     label = traits), size = 3, 
  #                 max.overlaps = 50) +
  # geom_point(data = PC2_strong_loading, 
  #            aes(PC1, PC2), 
  #            size = 2)+
  # geom_text_repel(data = PC2_strong_loading, 
  #                 aes(x = PC1, 
  #                     y = PC2, 
  #                     label = traits), size = 3, 
  #                 max.overlaps = 50) +
  # xlim(-2, 2)+
  # ylim(-2, 2)+
  theme_bw() +
  labs(x = "PC1 loading", y = "PC2 loading", 
       title = 'B) Trait loadings')


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

# F2_PCA_means = F2_PCA %>%
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
# vectors <- F2_PCA_means %>%
#   pivot_wider(
#     names_from = temp_rearing,
#     values_from = c(PC1, PC2)
#   )

normy_cols = c('#3d348b',
               '#ff006e')

# normy_cols2 = c('#3d348b',
#                '#ff006e')




multivar_reaction_norm = ggplot() +
  # geom_point(data = F2_PCA_means, 
  #            aes(PC1, 
  #                PC2, 
  #                color = temp_rearing, 
  #                shape = Morph), 
  #            size = 4)+
  geom_point(data = F2_off_means,
             aes(PC1,
                 PC2,
                 color = Morph,
                 shape = Offspring_temp),
             size = 4) +
  # geom_point(data = F2_parent_means, 
  #            aes(PC1, 
  #                PC2, 
  #                color = Morph, 
  #                shape = Parent_temp), 
  #            size = 2)
  
  geom_segment(data = F2_off_vectors,
               aes(x = PC1_12, y = PC2_12,
                   xend = PC1_18, yend = PC2_18,
                   color = Morph),
               arrow = arrow(length = unit(0.25, "cm")),
               linewidth = 1.2) +
  # geom_segment(data = F2_off_dashed_vec,
  #              aes(x = PC1_Cold, y = PC2_Cold,
  #                  xend = PC1_Warm, yend = PC2_Warm),
  #              linetype = 'dashed',
  #              arrow = arrow(length = unit(0.25, "cm")),
  #              linewidth = 1.2) +
  
  xlim(-2, 2)+
  ylim(-2, 2)+
  geom_hline(yintercept = 0, 
             linetype = 'dashed')+
  geom_vline(xintercept = 0, 
             linetype = 'dashed')+
  labs(x = "PC1",
       y = "PC2", 
       title = 'C) Offspring temperature') +
  scale_color_manual(values = normy_cols)+
  theme_bw() +
  theme(legend.position = 'none')+
  # guides(fill ='none')+
  

ggplot() +
  # geom_point(data = F2_PCA_means, 
  #            aes(PC1, 
  #                PC2, 
  #                color = temp_rearing, 
  #                shape = Morph), 
  #            size = 4)+
  geom_point(data = F2_parent_means,
             aes(PC1,
                 PC2,
                 color = Morph,
                 shape = Parent_temp),
             size = 4) +
  # geom_point(data = F2_parent_means, 
  #            aes(PC1, 
  #                PC2, 
  #                color = Morph, 
  #                shape = Parent_temp), 
  #            size = 2)
  
  geom_segment(data = F2_par_vectors,
               aes(x = PC1_12, y = PC2_12,
                   xend = PC1_18, yend = PC2_18,
                   color = Morph),
               arrow = arrow(length = unit(0.25, "cm")),
               linewidth = 1.2) +
  # geom_segment(data = F2_parent_dashed_vec,
  #              aes(x = PC1_Cold, y = PC2_Cold,
  #                  xend = PC1_Warm, yend = PC2_Warm),
  #              linetype = 'dashed',
  #              arrow = arrow(length = unit(0.25, "cm")),
  #              linewidth = 1.2) +
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
       title = 'D) Parental temperature') +
  scale_colour_manual(values = normy_cols)+
  theme_bw() +
  theme(legend.position = 'none')


Multivariate_reaction_plot = trait_loading_rank + trait_loading_gg / multivar_reaction_norm

ggsave('Multivariate_reaction_norm_graph_12.02.2026.tiff', 
       plot = Multivariate_reaction_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 50, 
       height = 40)  


ggsave('Multivariate_reaction_norm_graph_12.02.2026.svg', 
       plot = Multivariate_reaction_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 50, 
       height = 40)  
