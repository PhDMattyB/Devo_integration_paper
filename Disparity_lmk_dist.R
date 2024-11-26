##############################
## Disparity analysis linear distances
##
## Matt Brachmann (PhDMattyB)
##
## 18.11.2024
##
##############################

setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/Updated Landmarks/')

library(tidyverse)
library(geomorph)


# wild disparity ----------------------------------------------------------

wild_lmk_dist = read_csv('Wild_Univariate_traits.csv') %>% 
  filter(Lake %in% c('ASHN', 
                     'MYV', 
                     'SKR', 
                     'GTS', 
                     'CSWY')) %>% 
  mutate(Lake_morph= as.factor(case_when(
    Lake_morph == 'GTS' ~ 'GTSW',
    Lake_morph == 'CSWY' ~ 'CSWYC',
    Lake_morph == 'ASHNW' ~ 'ASHNW',
    Lake_morph == 'ASHNC' ~ 'ASHNC',
    Lake_morph == 'MYVW' ~ 'MYVW',
    Lake_morph == 'MYVC' ~ 'MYVC',
    Lake_morph == 'SKRW' ~ 'SKRW',
    Lake_morph == 'SKRC' ~ 'SKRC')))

wild_dist = wild_lmk_dist %>% 
  dplyr::select(2:29)

wild_mat = as.matrix(wild_dist)

Wild_disparity = morphol.disparity(wild_mat ~ 1,
                                       groups = ~Lake_morph,
                                       data = wild_lmk_dist,
                                       iter = 999)

wild_pval = Wild_disparity$PV.dist.Pval
wild_disparity_dist = Wild_disparity$PV.dist
wild_proc_var = Wild_disparity$Procrustes.var

wild_disp_pval = wild_pval %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  # as_tibble() %>% 
  melt(id.vars = c('rowname')) %>% 
  as_tibble()


wild_disp_dist = wild_disparity_dist %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  # as_tibble() %>% 
  melt(id.vars = c('rowname')) %>% 
  as_tibble()

wild_disp_data = bind_cols(wild_disp_dist,
                           wild_disp_pval) %>% 
  dplyr::select(1:3, 
         6) %>% 
  rename(Ecotype1 = 1, 
         Ecotype2 = 2, 
         zscore = 3, 
         pvalue = 4) %>% 
  mutate(across(where(is.numeric),
                ~ round(., 3))) %>% 
  mutate(label = 'Grandparental (wild) generation')

wild_disp_data %>% 
  write_csv('Wild_disparity_data.csv')

# F2 uncorrected data -----------------------------------------------------

F2_raw_lmk_dist = read_csv('F2_Original_univariate_traits.csv')
F2_raw_dist = F2_raw_lmk_dist %>% 
  dplyr::select(2:29)

F2_raw_mat = as.matrix(F2_raw_dist)

F2_raw_disparity = morphol.disparity(F2_raw_mat ~ 1,
                                   groups = ~Lake_morph,
                                   data = F2_raw_lmk_dist,
                                   iter = 999)

F2_raw_pval = F2_raw_disparity$PV.dist.Pval
F2_raw_disparity_dist = F2_raw_disparity$PV.dist
F2_raw_proc_var = F2_raw_disparity$Procrustes.var

F2_raw_disp_pval = F2_raw_pval %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  # as_tibble() %>% 
  melt(id.vars = c('rowname')) %>% 
  as_tibble()


F2_raw_disp_dist = F2_raw_disparity_dist %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  # as_tibble() %>% 
  melt(id.vars = c('rowname')) %>% 
  as_tibble()

F2_raw_disp_data = bind_cols(F2_raw_disp_dist,
                           F2_raw_disp_pval) %>% 
  dplyr::select(1:3, 
                6) %>% 
  rename(Ecotype1 = 1, 
         Ecotype2 = 2, 
         zscore = 3, 
         pvalue = 4) %>% 
  mutate(across(where(is.numeric),
                ~ round(., 3))) %>% 
  mutate(label = 'F2 generation')

F2_raw_disp_data %>% 
  write_csv('F2_raw_disparity_data.csv')

# TGP effect on disparity --------------------------------------------------
F1_lmk_dist = read_csv('F1_Plasticity_Corrected.csv')
F1_dist = F1_lmk_dist %>% 
  dplyr::select(2:29)

TGP_mat = as.matrix(F1_dist)

TGP_disparity = morphol.disparity(TGP_mat ~ 1,
                                     groups = ~Lake_morph,
                                     data = F2_raw_lmk_dist,
                                     iter = 999)

TGP_pval = TGP_disparity$PV.dist.Pval
TGP_disparity_dist = TGP_disparity$PV.dist
TGP_proc_var = TGP_disparity$Procrustes.var

TGP_disp_pval = TGP_pval %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  # as_tibble() %>% 
  melt(id.vars = c('rowname')) %>% 
  as_tibble()


TGP_disp_dist = TGP_disparity_dist %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  # as_tibble() %>% 
  melt(id.vars = c('rowname')) %>% 
  as_tibble()

TGP_disp_data = bind_cols(TGP_disp_dist,
                           TGP_disp_pval) %>% 
  dplyr::select(1:3, 
                6) %>% 
  rename(Ecotype1 = 1, 
         Ecotype2 = 2, 
         zscore = 3, 
         pvalue = 4) %>% 
  mutate(across(where(is.numeric),
                ~ round(., 3))) %>% 
  mutate(label = 'Trans-generational plasticity')

TGP_disp_data %>% 
  write_csv('TGP_disparity_data.csv')

# WGP effect on disparity -------------------------------------------------
F2_lmk_dist = read_csv('F2_Corrected_F2_temp_only.csv')
F2_dist = F2_lmk_dist %>% 
  dplyr::select(2:29)

WGP_mat = as.matrix(F2_dist)

WGP_disparity = morphol.disparity(WGP_mat ~ 1,
                                  groups = ~Lake_morph,
                                  data = F2_raw_lmk_dist,
                                  iter = 999)

WGP_pval = WGP_disparity$PV.dist.Pval
WGP_disparity_dist = WGP_disparity$PV.dist
WGP_proc_var = WGP_disparity$Procrustes.var

WGP_disp_pval = WGP_pval %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  # as_tibble() %>% 
  melt(id.vars = c('rowname')) %>% 
  as_tibble()


WGP_disp_dist = WGP_disparity_dist %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  # as_tibble() %>% 
  melt(id.vars = c('rowname')) %>% 
  as_tibble()

WGP_disp_data = bind_cols(WGP_disp_dist,
                           WGP_disp_pval) %>% 
  dplyr::select(1:3, 
                6) %>% 
  rename(Ecotype1 = 1, 
         Ecotype2 = 2, 
         zscore = 3, 
         pvalue = 4) %>% 
  mutate(across(where(is.numeric),
                ~ round(., 3))) %>% 
  mutate(label = 'Within-generational plasticity')

WGP_disp_data %>% 
  write_csv('WGP_disparity_data.csv')


# Graph combined disparity data -------------------------------------------

wild_disp_data = read_csv('Wild_disparity_data.csv')
F2_disp_data = read_csv('F2_raw_disparity_data.csv')
TGP_disp_data = read_csv('TGP_disparity_data.csv')
WGP_disp_data = read_csv('WGP_disparity_data.csv')


full_disp_data = bind_rows(wild_disp_data, 
                           F2_disp_data, 
                           TGP_disp_data, 
                           WGP_disp_data)

full_disp_data$label = factor(full_disp_data$label, 
                              levels = c('Grandparental (wild) generation', 
                                         'F2 generation', 
                                         'Trans-generational plasticity', 
                                         'Within-generational plasticity'))
full_disp_data$stars = cut(full_disp_data$pvalue, 
                          breaks = c(-Inf, 
                                     0.001, 
                                     0.01,
                                     0.05, 
                                     Inf), 
                          label=c("***", "**", "*", ""))


Disparity_plot = ggplot(full_disp_data, 
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
  # geom_rect(aes(fill = Eco_pair), 
  #           colour = "#fb6f92")+
  facet_wrap(~label)+
  scale_fill_viridis(discrete=FALSE, 
                     # direction = -1, 
                     option = 'D') +
  labs(fill = 'Disparity distance')+
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
        legend.box.just = "left", 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold', 
                                  size = 14))


ggsave('Full_disparity_figure.tiff', 
       plot = Disparity_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 15)
