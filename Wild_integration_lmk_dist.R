##############################
## wild vrel interlandmark dist
##
## Matt Brachmann (PhDMattyB)
##
## 18.11.2024
##
##############################


# functions ---------------------------------------------------------------
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}


# start -------------------------------------------------------------------


setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/Updated Landmarks/')

library(tidyverse)
library(geomorph)
library(reshape2)
library(viridis)
library(hrbrthemes)

# wild traits mag integration ---------------------------------------------

wild_lmk_dist = read_csv('Wild_Univariate_traits.csv')

wild_dist = wild_lmk_dist %>% 
  select(2:29)

# lmk_dist = geomorph.data.frame(lmk_dist)
wild_lmk_matrix = as.matrix(wild_dist)
wild_lmk_array = arrayspecs(wild_lmk_matrix, 14, 2)

wild_lmk_sub = coords.subset(wild_lmk_array, 
              wild_lmk_dist$Lake_morph)


vrel_wild_lmkdist = Map(function(x) integration.Vrel(x), 
                       wild_lmk_sub)

ASHN_compare_wild = compare.ZVrel(vrel_wild_lmkdist$ASHNC, 
                             vrel_wild_lmkdist$ASHNW)
MYV_compare_wild = compare.ZVrel(vrel_wild_lmkdist$MYVC, 
                            vrel_wild_lmkdist$MYVW)

SKR_compare_wild = compare.ZVrel(vrel_wild_lmkdist$SKRC, 
                            vrel_wild_lmkdist$SKRW)

GTS_CSWY_compare_wild = compare.ZVrel(vrel_wild_lmkdist$CSWY, 
                                 vrel_wild_lmkdist$GTS)

Wild_vrel_compare = compare.ZVrel(vrel_wild_lmkdist$ASHNC, 
              vrel_wild_lmkdist$ASHNW, 
              vrel_wild_lmkdist$MYVC, 
              vrel_wild_lmkdist$MYVW, 
              vrel_wild_lmkdist$SKRC, 
              vrel_wild_lmkdist$SKRW, 
              vrel_wild_lmkdist$CSWY, 
              vrel_wild_lmkdist$GTS)

# Wild_vrel_compare$pairwise.z

wild_zscore = Wild_vrel_compare$pairwise.z %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  # as_tibble() %>% 
  melt(id.vars = c('rowname')) %>% 
  as_tibble()


Wild_pval = Wild_vrel_compare$pairwise.P %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  # as_tibble() %>% 
  melt(id.vars = c('rowname')) %>% 
  as_tibble()

Wild_int_data = bind_cols(wild_zscore,
                          Wild_pval) %>% 
  select(1:3, 
         6) %>% 
  rename(Ecotype1 = 1, 
         Ecotype2 = 2, 
         zscore = 3, 
         pvalue = 4) %>% 
  separate(col = Ecotype1, 
           into = c('trash', 
                    'Ecotype1'), 
           sep = "[$]") %>% 
  separate(col = Ecotype2, 
           into = c('trash2', 
                    'Ecotype2'), 
           sep = '[$]') %>% 
  select(-trash, 
         -trash2)%>% 
  mutate(across(where(is.numeric),
                ~ round(., 3))) 

# Wild_int_data %>% 
#   write_csv('Wild_integration_metric.csv')


# F2 Uncorrected integration ----------------------------------------------


F2_raw_lmk_dist = read_csv('F2_Original_univariate_traits.csv')
F2_raw_dist = F2_raw_lmk_dist %>% 
  select(2:29)

# lmk_dist = geomorph.data.frame(lmk_dist)
F2_raw_lmk_matrix = as.matrix(F2_raw_dist)
F2_raw_lmk_array = arrayspecs(F2_raw_lmk_matrix, 14, 2)

F2_raw_lmk_sub = coords.subset(F2_raw_lmk_array, 
                           F2_raw_lmk_dist$Lake_morph)


vrel_F2_raw_lmkdist = Map(function(x) integration.Vrel(x), 
                      F2_raw_lmk_sub)

# ASHN_compare_raw = compare.ZVrel(vrel_F2_raw_lmkdist$ASHNC, 
#                              vrel_F2_raw_lmkdist$ASHNW)
# MYV_compare_raw = compare.ZVrel(vrel_F2_raw_lmkdist$MYVC, 
#                             vrel_F2_raw_lmkdist$MYVW)
# 
# SKR_compare_raw = compare.ZVrel(vrel_F2_raw_lmkdist$SKRC, 
#                             vrel_F2_raw_lmkdist$SKRW)
# 
# GTS_CSWY_compare_raw = compare.ZVrel(vrel_F2_raw_lmkdist$CSWY, 
#                                  vrel_F2_raw_lmkdist$GTS)

F2_raw_vrel_compare = compare.ZVrel(vrel_F2_raw_lmkdist$ASHNC, 
                                    vrel_F2_raw_lmkdist$ASHNW, 
                                    vrel_F2_raw_lmkdist$MYVC, 
                                    vrel_F2_raw_lmkdist$MYVW, 
                                    vrel_F2_raw_lmkdist$SKRC, 
                                    vrel_F2_raw_lmkdist$SKRW, 
                                    vrel_F2_raw_lmkdist$CSWY, 
                                    vrel_F2_raw_lmkdist$GTS)


F2_raw_zscore = F2_raw_vrel_compare$pairwise.z %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  # as_tibble() %>% 
  melt(id.vars = c('rowname')) %>% 
  as_tibble()


F2_raw_pval = F2_raw_vrel_compare$pairwise.P %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  # as_tibble() %>% 
  melt(id.vars = c('rowname')) %>% 
  as_tibble()

F2_raw_int_data = bind_cols(F2_raw_zscore,
                          F2_raw_pval) %>% 
  select(1:3, 
         6) %>% 
  rename(Ecotype1 = 1, 
         Ecotype2 = 2, 
         zscore = 3, 
         pvalue = 4) %>% 
  separate(col = Ecotype1, 
           into = c('trash', 
                    'Ecotype1'), 
           sep = "[$]") %>% 
  separate(col = Ecotype2, 
           into = c('trash2', 
                    'Ecotype2'), 
           sep = '[$]') %>% 
  select(-trash, 
         -trash2)%>% 
  mutate(across(where(is.numeric),
                ~ round(., 3))) 

F2_raw_int_data %>%
  write_csv('F2_uncorrected_integration_metric.csv')

##
# F1 effect mag integration -----------------------------------------------
F1_lmk_dist = read_csv('F1_Plasticity_Corrected.csv')
F1_dist = F1_lmk_dist %>% 
  select(2:29)

# lmk_dist = geomorph.data.frame(lmk_dist)
F1_lmk_matrix = as.matrix(F1_dist)
F1_lmk_array = arrayspecs(F1_lmk_matrix, 14, 2)

F1_lmk_sub = coords.subset(F1_lmk_array, 
                             F1_lmk_dist$Lake_morph)


vrel_F1_lmkdist = Map(function(x) integration.Vrel(x), 
                        F1_lmk_sub)

# ASHN_compare_f1 = compare.ZVrel(vrel_F1_lmkdist$ASHNC, 
#                              vrel_F1_lmkdist$ASHNW)
# MYV_compare_f1 = compare.ZVrel(vrel_F1_lmkdist$MYVC, 
#                             vrel_F1_lmkdist$MYVW)
# 
# SKR_compare_f1 = compare.ZVrel(vrel_F1_lmkdist$SKRC, 
#                             vrel_F1_lmkdist$SKRW)
# 
# GTS_CSWY_compare_f1 = compare.ZVrel(vrel_F1_lmkdist$CSWY, 
#                                  vrel_F1_lmkdist$GTS)

F1_effect_vrel_compare = compare.ZVrel(vrel_F1_lmkdist$ASHNC, 
                                       vrel_F1_lmkdist$ASHNW, 
                                       vrel_F1_lmkdist$MYVC, 
                                       vrel_F1_lmkdist$MYVW, 
                                       vrel_F1_lmkdist$SKRC, 
                                       vrel_F1_lmkdist$SKRW, 
                                       vrel_F1_lmkdist$CSWY, 
                                       vrel_F1_lmkdist$GTS)


F1_effect_zscore = F1_effect_vrel_compare$pairwise.z %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  # as_tibble() %>% 
  melt(id.vars = c('rowname')) %>% 
  as_tibble()


F1_effect_pval = F1_effect_vrel_compare$pairwise.P %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  # as_tibble() %>% 
  melt(id.vars = c('rowname')) %>% 
  as_tibble()

F1_effect_int_data = bind_cols(F1_effect_zscore,
                            F1_effect_pval) %>% 
  select(1:3, 
         6) %>% 
  rename(Ecotype1 = 1, 
         Ecotype2 = 2, 
         zscore = 3, 
         pvalue = 4) %>% 
  separate(col = Ecotype1, 
           into = c('trash', 
                    'Ecotype1'), 
           sep = "[$]") %>% 
  separate(col = Ecotype2, 
           into = c('trash2', 
                    'Ecotype2'), 
           sep = '[$]') %>% 
  select(-trash, 
         -trash2)%>% 
  mutate(across(where(is.numeric),
                ~ round(., 3))) 

F1_effect_int_data %>%
  write_csv('F1_Effect_integration_metric.csv')


# F2 effect mag integration -----------------------------------------------

F2_lmk_dist = read_csv('F2_Corrected_F2_temp_only.csv')
F2_dist = F2_lmk_dist %>% 
  select(2:29)

# lmk_dist = geomorph.data.frame(lmk_dist)
F2_lmk_matrix = as.matrix(F2_dist)
F2_lmk_array = arrayspecs(F2_lmk_matrix, 14, 2)

F2_lmk_sub = coords.subset(F2_lmk_array, 
                           F2_lmk_dist$Lake_morph)


vrel_F2_lmkdist = Map(function(x) integration.Vrel(x), 
                      F2_lmk_sub)

# ASHN_compare_f2 = compare.ZVrel(vrel_F2_lmkdist$ASHNC, 
#                              vrel_F2_lmkdist$ASHNW)
# MYV_compare_f2 = compare.ZVrel(vrel_F2_lmkdist$MYVC, 
#                             vrel_F2_lmkdist$MYVW)
# 
# SKR_compare_f2 = compare.ZVrel(vrel_F2_lmkdist$SKRC, 
#                             vrel_F2_lmkdist$SKRW)
# 
# GTS_CSWY_compare_f2 = compare.ZVrel(vrel_F2_lmkdist$CSWY, 
#                                  vrel_F2_lmkdist$GTS)
# 

F2_effect_vrel_compare = compare.ZVrel(vrel_F2_lmkdist$ASHNC, 
                                       vrel_F2_lmkdist$ASHNW, 
                                       vrel_F2_lmkdist$MYVC, 
                                       vrel_F2_lmkdist$MYVW, 
                                       vrel_F2_lmkdist$SKRC, 
                                       vrel_F2_lmkdist$SKRW, 
                                       vrel_F2_lmkdist$CSWY, 
                                       vrel_F2_lmkdist$GTS)


F2_effect_zscore = F2_effect_vrel_compare$pairwise.z %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  # as_tibble() %>% 
  melt(id.vars = c('rowname')) %>% 
  as_tibble()


F2_effect_pval = F2_effect_vrel_compare$pairwise.P %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  # as_tibble() %>% 
  melt(id.vars = c('rowname')) %>% 
  as_tibble()

F2_effect_int_data = bind_cols(F2_effect_zscore,
                               F2_effect_pval) %>% 
  select(1:3, 
         6) %>% 
  rename(Ecotype1 = 1, 
         Ecotype2 = 2, 
         zscore = 3, 
         pvalue = 4) %>% 
  separate(col = Ecotype1, 
           into = c('trash', 
                    'Ecotype1'), 
           sep = "[$]") %>% 
  separate(col = Ecotype2, 
           into = c('trash2', 
                    'Ecotype2'), 
           sep = '[$]') %>% 
  select(-trash, 
         -trash2)%>% 
  mutate(across(where(is.numeric),
                ~ round(., 3))) 

F2_effect_int_data %>%
  write_csv('F2_Effect_integration_metric.csv')




# Magnitude integration graphs --------------------------------------------

## Wild integration plots
Wild_int_data = read_csv('Wild_integration_metric.csv')
Wild_int_data$Ecotype1 = factor(Wild_int_data$Ecotype1, 
                                levels = c('ASHNC', 
                                           'ASHNW', 
                                           'MYVC', 
                                           'MYVW', 
                                           'SKRC', 
                                           'SKRW', 
                                           'CSWY', 
                                           'GTS'))
Wild_int_data$Ecotype2 = factor(Wild_int_data$Ecotype2, 
                                levels = c('ASHNC', 
                                           'ASHNW', 
                                           'MYVC', 
                                           'MYVW', 
                                           'SKRC', 
                                           'SKRW', 
                                           'CSWY', 
                                           'GTS'))
Wild_int_data$stars = cut(Wild_int_data$pvalue, 
                          breaks = c(-Inf, 
                                     0.001, 
                                     0.01,
                                     0.05, 
                                     Inf), 
                          label=c("***", "**", "*", ""))

Wild_int_plot = ggplot(Wild_int_data, 
                       aes(Ecotype1, 
                           Ecotype2, 
                           fill= zscore)) + 
  geom_tile(col = 'white') +
  # geom_text(aes(label = zscore),
  #           color = "black",
  #           size = 2,
  #           fontface = 'bold')+
  geom_text(aes(label=stars), 
            color="black", 
            size=5) + 
  scale_fill_viridis(discrete=FALSE, 
                     # direction = -1, 
                     option = 'D') +
  labs(title = 'A) Grandpartental (wild) generation', 
       fill = 'Pairwise z-score')+
  # theme_ipsum()+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.ticks = element_blank(), 
        panel.border = element_blank(), 
        # legend.justification = c('left', 'top'), 
        legend.position = c(.05, 0.98),
        legend.justification = c("left", "top"),
        legend.box.just = "left")

## F2 uncorrected integration plots
F2_raw_int_data = read_csv('F2_uncorrected_integration_metric.csv')
F2_raw_int_data$Ecotype1 = factor(F2_raw_int_data$Ecotype1, 
                                levels = c('ASHNC', 
                                           'ASHNW', 
                                           'MYVC', 
                                           'MYVW', 
                                           'SKRC', 
                                           'SKRW', 
                                           'CSWY', 
                                           'GTS'))
F2_raw_int_data$Ecotype2 = factor(F2_raw_int_data$Ecotype2, 
                                levels = c('ASHNC', 
                                           'ASHNW', 
                                           'MYVC', 
                                           'MYVW', 
                                           'SKRC', 
                                           'SKRW', 
                                           'CSWY', 
                                           'GTS'))
F2_raw_int_data$stars = cut(F2_raw_int_data$pvalue, 
                          breaks = c(-Inf, 
                                     0.001, 
                                     0.01,
                                     0.05, 
                                     Inf), 
                          label=c("***", "**", "*", ""))

F2_raw_int_plot = ggplot(F2_raw_int_data, 
                       aes(Ecotype1, 
                           Ecotype2, 
                           fill= zscore)) + 
  geom_tile(col = 'white') +
  # geom_text(aes(label = zscore),
  #           color = "black",
  #           size = 2,
  #           fontface = 'bold')+
  geom_text(aes(label=stars), 
            color="black", 
            size=5) + 
  scale_fill_viridis(discrete=FALSE, 
                     # direction = -1, 
                     option = 'D') +
  labs(title = 'B) F2 generation (uncorrected)', 
       fill = 'Pairwise z-score')+
  # theme_ipsum()+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.ticks = element_blank(), 
        panel.border = element_blank(), 
        # legend.justification = c('left', 'top'), 
        legend.position = c(.05, 0.98),
        legend.justification = c("left", "top"),
        legend.box.just = "left")


## F1 effect magnitude integration plot
F1_effect_int_data = read_csv('F1_Effect_integration_metric.csv')
F1_effect_int_data$Ecotype1 = factor(F1_effect_int_data$Ecotype1, 
                                  levels = c('ASHNC', 
                                             'ASHNW', 
                                             'MYVC', 
                                             'MYVW', 
                                             'SKRC', 
                                             'SKRW', 
                                             'CSWY', 
                                             'GTS'))
F1_effect_int_data$Ecotype2 = factor(F1_effect_int_data$Ecotype2, 
                                  levels = c('ASHNC', 
                                             'ASHNW', 
                                             'MYVC', 
                                             'MYVW', 
                                             'SKRC', 
                                             'SKRW', 
                                             'CSWY', 
                                             'GTS'))
F1_effect_int_data$stars = cut(F1_effect_int_data$pvalue, 
                            breaks = c(-Inf, 
                                       0.001, 
                                       0.01,
                                       0.05, 
                                       Inf), 
                            label=c("***", "**", "*", ""))

F1_effect_int_plot = ggplot(F1_effect_int_data, 
                         aes(Ecotype1, 
                             Ecotype2, 
                             fill= zscore)) + 
  geom_tile(col = 'white') +
  # geom_text(aes(label = zscore),
  #           color = "black",
  #           size = 2,
  #           fontface = 'bold')+
  geom_text(aes(label=stars), 
            color="black", 
            size=5) + 
  scale_fill_viridis(discrete=FALSE, 
                     # direction = -1, 
                     option = 'D') +
  labs(title = 'B) F2 generation (uncorrected)', 
       fill = 'Pairwise z-score')+
  # theme_ipsum()+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.ticks = element_blank(), 
        panel.border = element_blank(), 
        # legend.justification = c('left', 'top'), 
        legend.position = c(.05, 0.98),
        legend.justification = c("left", "top"),
        legend.box.just = "left")

## F2 effect magnitude integration plot
F2_effect_int_data = read_csv('F2_effect_integration_metric.csv')
F2_effect_int_data$Ecotype1 = factor(F2_effect_int_data$Ecotype1, 
                                     levels = c('ASHNC', 
                                                'ASHNW', 
                                                'MYVC', 
                                                'MYVW', 
                                                'SKRC', 
                                                'SKRW', 
                                                'CSWY', 
                                                'GTS'))
F2_effect_int_data$Ecotype2 = factor(F2_effect_int_data$Ecotype2, 
                                     levels = c('ASHNC', 
                                                'ASHNW', 
                                                'MYVC', 
                                                'MYVW', 
                                                'SKRC', 
                                                'SKRW', 
                                                'CSWY', 
                                                'GTS'))
F2_effect_int_data$stars = cut(F2_effect_int_data$pvalue, 
                               breaks = c(-Inf, 
                                          0.001, 
                                          0.01,
                                          0.05, 
                                          Inf), 
                               label=c("***", "**", "*", ""))

F2_effect_int_plot = ggplot(F2_effect_int_data, 
                            aes(Ecotype1, 
                                Ecotype2, 
                                fill= zscore)) + 
  geom_tile(col = 'white') +
  # geom_text(aes(label = zscore),
  #           color = "black",
  #           size = 2,
  #           fontface = 'bold')+
  geom_text(aes(label=stars), 
            color="black", 
            size=5) + 
  scale_fill_viridis(discrete=FALSE, 
                     # direction = -1, 
                     option = 'D') +
  labs(title = 'B) F2 generation (uncorrected)', 
       fill = 'Pairwise z-score')+
  # theme_ipsum()+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.ticks = element_blank(), 
        panel.border = element_blank(), 
        # legend.justification = c('left', 'top'), 
        legend.position = c(.05, 0.98),
        legend.justification = c("left", "top"),
        legend.box.just = "left")
