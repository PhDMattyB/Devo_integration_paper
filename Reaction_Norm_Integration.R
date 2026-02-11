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
                jaw_length:caudal2_15_17) %>%
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


WGP_traits = read_csv("WGP_TRAITS_SCALED_FIXED_11.02.2026.csv")
TGP_traits = read_csv('TGP_TRAITS_SCALED_FIXED_11.02.2026')

## F2 plasticity
ggplot(F2_orig_traits, 
       aes(x = offspring_temp,
               y = value,
               color = morph,
               group = interaction(morph, Ecotype_pair))) +
  
  stat_summary(fun = mean, geom = "line", linewidth = 0.6, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 1.5) +
  
  facet_wrap(~ trait, scales = "free_y") +
  
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
           group = interaction(morph, Ecotype_pair))) +
  
  stat_summary(fun = mean, geom = "line", linewidth = 0.6, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 1.5) +
  
  facet_wrap(~ trait, scales = "free_y") +
  
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
