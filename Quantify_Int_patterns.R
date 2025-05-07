##############################
## Quantifying patterns of integration
##
## Matt Brachmann (PhDMattyB)
##
## 
## 07.05.2025
##############################

setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/Updated Landmarks/')


# WILD --------------------------------------------------------------------

ASHN_wild_cor_0.3 = read_csv("ASHN_wild_Pattern_integration_Ecotype_diffs_0.3cutoff.csv") %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_', 
        remove = F)

MYV_wild_cor_0.3 = read_csv('MYV_wild_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_', 
        remove = F)

SKR_wild_cor_0.3 = read_csv('SKR_wild_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_', 
        remove = F)

GTS_CSWY_wild_cor_0.3 = read_csv('GTS_CSWY_wild_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_', 
        remove = F)

inner_join(ASHN_wild_cor_0.3, 
           MYV_wild_cor_0.3, 
           by = 'Integrated_traits') %>% 
  arrange(Integrated_traits) 

inner_join(ASHN_wild_cor_0.3, 
           SKR_wild_cor_0.3, 
           by = 'Integrated_traits')

inner_join(ASHN_wild_cor_0.3, 
           GTS_CSWY_wild_cor_0.3, 
           by = 'Integrated_traits')

inner_join(MYV_wild_cor_0.3, 
           SKR_wild_cor_0.3, 
           by = 'Integrated_traits')

inner_join(MYV_wild_cor_0.3, 
           GTS_CSWY_wild_cor_0.3, 
           by = 'Integrated_traits')

inner_join(SKR_wild_cor_0.3, 
           GTS_CSWY_wild_cor_0.3, 
           by = 'Integrated_traits')


inner_join(ASHN_wild_cor_0.3, 
           MYV_wild_cor_0.3, 
           by = 'Integrated_traits') %>% 
  inner_join(., 
             SKR_wild_cor_0.3, 
             by = 'Integrated_traits') %>% 
  inner_join(., 
             GTS_CSWY_wild_cor_0.3, 
             by = 'Integrated_traits') %>% 
  arrange(Integrated_traits) %>% 
  write_csv('WILD_Parallel_Pattern_Integration_cor0.3.csv')
  


# F2 Unaltered traits -----------------------------------------------------

ASHN_F2_cor_0.3 = read_csv("ASHN_F2Orig_Pattern_integration_Ecotype_diffs_0.3cutoff.csv") %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_', 
        remove = F)

MYV_F2_cor_0.3 = read_csv('MYV_F2Orig_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_', 
        remove = F)

SKR_F2_cor_0.3 = read_csv('SKR_F2Orig_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_', 
        remove = F)

GTS_CSWY_F2_cor_0.3 = read_csv('GTS_CSWY_F2Orig_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_', 
        remove = F)

inner_join(ASHN_F2_cor_0.3, 
           MYV_F2_cor_0.3, 
           by = 'Integrated_traits') 

inner_join(ASHN_F2_cor_0.3, 
           SKR_F2_cor_0.3, 
           by = 'Integrated_traits')

inner_join(ASHN_F2_cor_0.3, 
           GTS_CSWY_F2_cor_0.3, 
           by = 'Integrated_traits')

inner_join(MYV_F2_cor_0.3, 
           SKR_F2_cor_0.3, 
           by = 'Integrated_traits')

inner_join(MYV_F2_cor_0.3, 
           GTS_CSWY_F2_cor_0.3, 
           by = 'Integrated_traits')

inner_join(SKR_F2_cor_0.3, 
           GTS_CSWY_F2_cor_0.3, 
           by = 'Integrated_traits')


inner_join(ASHN_F2_cor_0.3, 
           MYV_F2_cor_0.3, 
           by = 'Integrated_traits') %>% 
  inner_join(., 
             SKR_F2_cor_0.3, 
             by = 'Integrated_traits') %>% 
  inner_join(., 
             GTS_CSWY_F2_cor_0.3, 
             by = 'Integrated_traits') %>% 
  arrange(Integrated_traits) %>% 
  write_csv('F2_Parallel_Pattern_Integration_cor0.3.csv')



# WGP ---------------------------------------------------------------------

ASHN_WGP_cor_0.3 = read_csv("ASHN_WGP_Pattern_integration_Ecotype_diffs_0.3cutoff.csv") %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_', 
        remove = F)

MYV_WGP_cor_0.3 = read_csv('MYV_WGP_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_', 
        remove = F)

SKR_WGP_cor_0.3 = read_csv('SKR_WGP_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_', 
        remove = F)

GTS_CSWY_WGP_cor_0.3 = read_csv('GTS_CSWY_WGP_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_', 
        remove = F)

inner_join(ASHN_WGP_cor_0.3, 
           MYV_WGP_cor_0.3, 
           by = 'Integrated_traits') 

inner_join(ASHN_WGP_cor_0.3, 
           SKR_WGP_cor_0.3, 
           by = 'Integrated_traits')

inner_join(ASHN_WGP_cor_0.3, 
           GTS_CSWY_WGP_cor_0.3, 
           by = 'Integrated_traits')

inner_join(MYV_WGP_cor_0.3, 
           SKR_WGP_cor_0.3, 
           by = 'Integrated_traits')

inner_join(MYV_WGP_cor_0.3, 
           GTS_CSWY_WGP_cor_0.3, 
           by = 'Integrated_traits')

inner_join(SKR_WGP_cor_0.3, 
           GTS_CSWY_WGP_cor_0.3, 
           by = 'Integrated_traits')


inner_join(ASHN_WGP_cor_0.3, 
           MYV_WGP_cor_0.3, 
           by = 'Integrated_traits') %>% 
  inner_join(., 
             SKR_WGP_cor_0.3, 
             by = 'Integrated_traits') %>% 
  inner_join(., 
             GTS_CSWY_WGP_cor_0.3, 
             by = 'Integrated_traits') %>% 
  arrange(Integrated_traits) %>% 
  write_csv('WGP_Parallel_Pattern_Integration_cor0.3.csv')


# TGP ---------------------------------------------------------------------

ASHN_TGP_cor_0.3 = read_csv("ASHN_TGP_Pattern_integration_Ecotype_diffs_0.3cutoff.csv") %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_', 
        remove = F)

MYV_TGP_cor_0.3 = read_csv('MYV_TGP_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_', 
        remove = F)

SKR_TGP_cor_0.3 = read_csv('SKR_TGP_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_', 
        remove = F)

GTS_CSWY_TGP_cor_0.3 = read_csv('GTS_CSWY_TGP_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_', 
        remove = F)

inner_join(ASHN_TGP_cor_0.3, 
           MYV_TGP_cor_0.3, 
           by = 'Integrated_traits') 

inner_join(ASHN_TGP_cor_0.3, 
           SKR_TGP_cor_0.3, 
           by = 'Integrated_traits')

inner_join(ASHN_TGP_cor_0.3, 
           GTS_CSWY_TGP_cor_0.3, 
           by = 'Integrated_traits')

inner_join(MYV_TGP_cor_0.3, 
           SKR_TGP_cor_0.3, 
           by = 'Integrated_traits')

inner_join(MYV_TGP_cor_0.3, 
           GTS_CSWY_TGP_cor_0.3, 
           by = 'Integrated_traits')

inner_join(SKR_TGP_cor_0.3, 
           GTS_CSWY_TGP_cor_0.3, 
           by = 'Integrated_traits')


inner_join(ASHN_TGP_cor_0.3, 
           MYV_TGP_cor_0.3, 
           by = 'Integrated_traits') %>% 
  inner_join(., 
             SKR_TGP_cor_0.3, 
             by = 'Integrated_traits') %>% 
  inner_join(., 
             GTS_CSWY_TGP_cor_0.3, 
             by = 'Integrated_traits') %>% 
  arrange(Integrated_traits) %>% 
  write_csv('TGP_Parallel_Pattern_Integration_cor0.3.csv')





# WILD VS F2 original -----------------------------------------------------

Wild_parallel = read_csv('WILD_Parallel_Pattern_Integration_cor0.3.csv')
F2_parallel = read_csv("F2_Parallel_Pattern_Integration_cor0.3.csv")
WGP_parallel = read_csv('WGP_Parallel_Pattern_Integration_cor0.3.csv')
TGP_parallel = read_csv("TGP_Parallel_Pattern_Integration_cor0.3.csv")

inner_join(Wild_parallel, 
           F2_parallel, 
           by = 'Integrated_traits') 

inner_join(Wild_parallel, 
          WGP_parallel, 
          by = 'Integrated_traits')

inner_join(Wild_parallel, 
           TGP_parallel, 
           by = 'Integrated_traits')

inner_join(F2_parallel, 
           WGP_parallel, 
           by = 'Integrated_traits')

inner_join(F2_parallel, 
           TGP_parallel, 
           by = 'Integrated_traits')

inner_join(WGP_parallel, 
           TGP_parallel, 
           by = 'Integrated_traits')

inner_join(Wild_parallel, 
           F2_parallel, 
           by = 'Integrated_traits') %>% 
  inner_join(.,
             WGP_parallel, 
             by = 'Integrated_traits') %>% 
  inner_join(., 
             TGP_parallel, 
             by = 'Integrated_traits') %>% 
  arrange(Integrated_traits) %>% 
  write_csv('Super_parallel_Integrated_traits.csv')


inner_join(F2_parallel, 
           WGP_parallel, 
           by = 'Integrated_traits') %>% 
  inner_join(.,
             TGP_parallel, 
             by = 'Integrated_traits') %>% 
  arrange(Integrated_traits) 



# Summarize the data ------------------------------------------------------

Wild_mean_parallel = Wild_parallel%>% 
  mutate(mean_cor_value = rowMeans(across(starts_with('value')))) %>%
  dplyr::select(Integrated_traits, 
                Var1.x, 
                Var2.x,
                mean_cor_value) %>% 
  mutate(Parallel_value = 'Wild parallel')

F2_mean_parallel = F2_parallel%>% 
  mutate(mean_cor_value = rowMeans(across(starts_with('value')))) %>%
  dplyr::select(Integrated_traits,
                Var1.x, 
                Var2.x,
                mean_cor_value)%>% 
  mutate(Parallel_value = 'F2 parallel')

WGP_mean_parallel = WGP_parallel%>% 
  mutate(mean_cor_value = rowMeans(across(starts_with('value')))) %>%
  dplyr::select(Integrated_traits, 
                Var1.x, 
                Var2.x,
                mean_cor_value) %>% 
  mutate(Parallel_value = 'WGP parallel')

TGP_mean_parallel = TGP_parallel%>% 
  mutate(mean_cor_value = rowMeans(across(starts_with('value')))) %>%
  dplyr::select(Integrated_traits, 
                Var1.x, 
                Var2.x,
                mean_cor_value) %>% 
  mutate(Parallel_value = 'TGP parallel')

super_mean_parallel = read_csv('Super_parallel_Integrated_traits.csv') %>% 
  mutate(mean_cor_value = rowMeans(across(starts_with('value')))) %>%
  dplyr::select(Integrated_traits,
                Var1.x.x.x, 
                Var2.x.x.x,
                mean_cor_value) %>% 
  mutate(Parallel_value = 'Super parallel')


background_traits = read_csv('ASHN_wild_cor_mat.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_', 
        remove = F) %>% 
  arrange(Var1) %>% 
  mutate(mean_cor_value = 1) %>% 
  dplyr::select(-value) %>% 
  mutate(Parallel_value = 'Not Parallel') %>% 
  rename(Var1.x = Var1, 
         Var2.x = Var2)


super_background_traits = read_csv('ASHN_wild_cor_mat.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_', 
        remove = F) %>% 
  arrange(Var1) %>% 
  mutate(mean_cor_value = 1) %>% 
  dplyr::select(-value) %>% 
  mutate(Parallel_value = 'Not Parallel') %>% 
  rename(Var1.x.x.x = Var1, 
         Var2.x.x.x = Var2)

## Need the background traits that are not parallel

wild_background_traits = anti_join(background_traits, 
          Wild_mean_parallel, 
          by = 'Integrated_traits')

F2_background_traits = anti_join(background_traits, 
                                   F2_mean_parallel, 
                                   by = 'Integrated_traits')

WGP_background_traits = anti_join(background_traits, 
                                  WGP_mean_parallel, 
                                  by = 'Integrated_traits')

TGP_background_traits = anti_join(background_traits, 
                                  TGP_mean_parallel, 
                                  by = 'Integrated_traits')

super_background_traits = anti_join(super_background_traits, 
                                  super_mean_parallel, 
                                  by = 'Integrated_traits')


## Graph time
## Wild parallel traits

wild_parallel_graph_data = bind_rows(Wild_mean_parallel, 
          wild_background_traits)


wild_parallel_graph = ggplot(wild_parallel_graph_data,
       aes(x = Var1.x,
           y = Var2.x,
           fill = mean_cor_value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ffe5d9",
                       mid = "#ff006e",
                       high = "#ffe5d9") +
  labs(title = 'A) Wild parallel traits')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        plot.title = element_text(size = 22),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        legend.position = 'none',
        # axis.text.x = element_blank(), 
        axis.text.x = element_text(angle = 90),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill = 'transparent', 
                                       color = NA))



### F2 parallel graph

F2_parallel_graph_data = bind_rows(F2_mean_parallel, 
                                     F2_background_traits)


F2_parallel_graph = ggplot(F2_parallel_graph_data,
                             aes(x = Var1.x,
                                 y = Var2.x,
                                 fill = mean_cor_value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ffe5d9",
                       mid = "#ff006e",
                       high = "#ffe5d9") +
  labs(title = 'B) F2 parallel traits')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        plot.title = element_text(size = 22),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        legend.position = 'none',
        # axis.text.x = element_blank(), 
        axis.text.x = element_text(angle = 90),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill = 'transparent', 
                                       color = NA))


### WGP parallel graph

WGP_parallel_graph_data = bind_rows(WGP_mean_parallel, 
                                   WGP_background_traits)


WGP_parallel_graph = ggplot(WGP_parallel_graph_data,
                           aes(x = Var1.x,
                               y = Var2.x,
                               fill = mean_cor_value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ffe5d9",
                       mid = "#ff006e",
                       high = "#ffe5d9") +
  labs(title = 'C) WGP parallel traits')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        plot.title = element_text(size = 22),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        legend.position = 'none',
        # axis.text.x = element_blank(), 
        axis.text.x = element_text(angle = 90),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill = 'transparent', 
                                       color = NA))

### TGP parallel graph

TGP_parallel_graph_data = bind_rows(TGP_mean_parallel, 
                                    TGP_background_traits)


TGP_parallel_graph = ggplot(TGP_parallel_graph_data,
                            aes(x = Var1.x,
                                y = Var2.x,
                                fill = mean_cor_value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ffe5d9",
                       mid = "#ff006e",
                       high = "#ffe5d9") +
  labs(title = 'D) TGP parallel traits')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        plot.title = element_text(size = 22),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        legend.position = 'none',
        # axis.text.x = element_blank(), 
        axis.text.x = element_text(angle = 90),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill = 'transparent', 
                                       color = NA))


### SUPER parallel graph

Super_parallel_graph_data = bind_rows(super_mean_parallel, 
                                    super_background_traits)


Super_parallel_graph = ggplot(Super_parallel_graph_data,
                            aes(x = Var1.x.x.x,
                                y = Var2.x.x.x,
                                fill = mean_cor_value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ffe5d9",
                       mid = "#ff006e",
                       high = "#ffe5d9") +
  labs(title = 'D) TGP parallel traits')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        plot.title = element_text(size = 22),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        legend.position = 'none',
        # axis.text.x = element_blank(), 
        axis.text.x = element_text(angle = 90),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill = 'transparent', 
                                       color = NA))
