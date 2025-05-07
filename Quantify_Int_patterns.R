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
        sep = '_')

MYV_wild_cor_0.3 = read_csv('MYV_wild_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_')

SKR_wild_cor_0.3 = read_csv('SKR_wild_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_')

GTS_CSWY_wild_cor_0.3 = read_csv('GTS_CSWY_wild_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_')

inner_join(ASHN_wild_cor_0.3, 
           MYV_wild_cor_0.3, 
           by = 'Integrated_traits') %>% 
  arrange(Integrated_traits) %>% View()

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
        sep = '_')

MYV_F2_cor_0.3 = read_csv('MYV_F2Orig_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_')

SKR_F2_cor_0.3 = read_csv('SKR_F2Orig_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_')

GTS_CSWY_F2_cor_0.3 = read_csv('GTS_CSWY_F2Orig_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_')

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


# WILD VS F2 original -----------------------------------------------------

Wild_parallel = read_csv('WILD_Parallel_Pattern_Integration_cor0.3.csv')
F2_parallel = read_csv("F2_Parallel_Pattern_Integration_cor0.3.csv")

inner_join(Wild_parallel, 
           F2_parallel, 
           by = 'Integrated_traits') %>% 
  arrange(Integrated_traits) %>% 
  View()

