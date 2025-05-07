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



# WGP ---------------------------------------------------------------------

ASHN_WGP_cor_0.3 = read_csv("ASHN_WGP_Pattern_integration_Ecotype_diffs_0.3cutoff.csv") %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_')

MYV_WGP_cor_0.3 = read_csv('MYV_WGP_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_')

SKR_WGP_cor_0.3 = read_csv('SKR_WGP_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_')

GTS_CSWY_WGP_cor_0.3 = read_csv('GTS_CSWY_WGP_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_')

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
        sep = '_')

MYV_TGP_cor_0.3 = read_csv('MYV_TGP_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_')

SKR_TGP_cor_0.3 = read_csv('SKR_TGP_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_')

GTS_CSWY_TGP_cor_0.3 = read_csv('GTS_CSWY_TGP_Pattern_integration_Ecotype_diffs_0.3cutoff.csv') %>% 
  unite(col = 'Integrated_traits', 
        c('Var1', 
          'Var2'), 
        sep = '_')

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
  View()


