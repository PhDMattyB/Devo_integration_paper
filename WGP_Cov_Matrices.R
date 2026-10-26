setwd("~/Parsons_postdoc/Stickleback_Morphometric_data/Updated Landmarks/")

library(tidyverse)
library(lme4)
library(lmerTest)
library(brms)
library(paran)
library(broom.mixed)
library(MatrixCorrelation)
library(vegan)
library(evolqg)


WGP_data = read_csv("WGP_TRAITS_SCALED_FIXED_11.02.2026.csv") %>% 
  dplyr::select(Lake_morph, 
                Ecotype_pair, 
                lake_morph_Pair_Full_Temp, 
                Morph, 
                Offspring_temp, 
                Parent_temp, 
                jaw_length:caudal2_15_17)%>% 
  unite(col = Exp_treats, 
        c('Parent_temp', 
          'Offspring_temp'), 
        sep = '_')

# WGP_orig data - KrzCov and procrustes ---------------------------------------
# ASHN --------------------------------------------------------------------

ashnc_WGP_12_12_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'ASHNC', 
         Exp_treats == '12_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
ashnc_WGP_12_18_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'ASHNC', 
         Exp_treats == '12_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
ashnc_WGP_18_12_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'ASHNC', 
         Exp_treats == '18_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
ashnc_WGP_18_18_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'ASHNC', 
         Exp_treats == '18_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)


ashnw_WGP_12_12_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'ASHNW', 
         Exp_treats == '12_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
ashnw_WGP_12_18_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'ASHNW', 
         Exp_treats == '12_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
ashnw_WGP_18_12_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'ASHNW', 
         Exp_treats == '18_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
ashnw_WGP_18_18_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'ASHNW', 
         Exp_treats == '18_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

ashnc_WGP_12_12_cov_mat = cov(ashnc_WGP_12_12_trait_mat, 
                             use = 'complete.obs')
ashnc_WGP_12_18_cov_mat = cov(ashnc_WGP_12_18_trait_mat, 
                             use = 'complete.obs')
ashnc_WGP_18_12_cov_mat = cov(ashnc_WGP_18_12_trait_mat, 
                             use = 'complete.obs')
ashnc_WGP_18_18_cov_mat = cov(ashnc_WGP_18_18_trait_mat, 
                             use = 'complete.obs')

ashnw_WGP_12_12_cov_mat = cov(ashnw_WGP_12_12_trait_mat, 
                             use = 'complete.obs')
ashnw_WGP_12_18_cov_mat = cov(ashnw_WGP_12_18_trait_mat, 
                             use = 'complete.obs')
ashnw_WGP_18_12_cov_mat = cov(ashnw_WGP_18_12_trait_mat, 
                             use = 'complete.obs')
ashnw_WGP_18_18_cov_mat = cov(ashnw_WGP_18_18_trait_mat, 
                             use = 'complete.obs')

ashn_WGP_12_12_krzcov = evolqg::KrzCor(ashnc_WGP_12_12_cov_mat, 
                                      ashnw_WGP_12_12_cov_mat)

ashn_WGP_12_18_krzcov = evolqg::KrzCor(ashnc_WGP_12_18_cov_mat, 
                                      ashnw_WGP_12_18_cov_mat)

ashn_WGP_18_12_krzcov = evolqg::KrzCor(ashnc_WGP_18_12_cov_mat, 
                                      ashnw_WGP_18_12_cov_mat)

ashn_WGP_18_18_krzcov = evolqg::KrzCor(ashnc_WGP_18_18_cov_mat, 
                                      ashnw_WGP_18_18_cov_mat)


ashn_WGP_12_12_proc_res = procrustes(ashnc_WGP_12_12_cov_mat, 
                                    ashnw_WGP_12_12_cov_mat)
ashn_WGP_12_18_proc_res = procrustes(ashnc_WGP_12_18_cov_mat, 
                                    ashnw_WGP_12_18_cov_mat)
ashn_WGP_18_12_proc_res = procrustes(ashnc_WGP_18_12_cov_mat, 
                                    ashnw_WGP_18_12_cov_mat)
ashn_WGP_18_18_proc_res = procrustes(ashnc_WGP_18_18_cov_mat, 
                                    ashnw_WGP_18_18_cov_mat)

summary(ashn_WGP_12_12_proc_res)
summary(ashn_WGP_12_18_proc_res)
summary(ashn_WGP_18_12_proc_res)
summary(ashn_WGP_18_18_proc_res)

ashn_WGP_12_12_proc_test = protest(ashnc_WGP_12_12_cov_mat, 
                                  ashnw_WGP_12_12_cov_mat, 
                                  permutations = 999)
ashn_WGP_12_18_proc_test = protest(ashnc_WGP_12_18_cov_mat, 
                                  ashnw_WGP_12_18_cov_mat, 
                                  permutations = 999)
ashn_WGP_18_12_proc_test = protest(ashnc_WGP_18_12_cov_mat, 
                                  ashnw_WGP_18_12_cov_mat, 
                                  permutations = 999)
ashn_WGP_18_18_proc_test = protest(ashnc_WGP_18_18_cov_mat, 
                                  ashnw_WGP_18_18_cov_mat, 
                                  permutations = 999)


ashn_WGP_12_12_residuals = residuals(ashn_WGP_12_12_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

ashn_WGP_12_18_residuals = residuals(ashn_WGP_12_18_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

ashn_WGP_18_12_residuals = residuals(ashn_WGP_18_12_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

ashn_WGP_18_18_residuals = residuals(ashn_WGP_18_18_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)


ashn_WGP_12_12_resid_plot = ggplot(ashn_WGP_12_12_residuals, 
                                  aes(x = factor(trait), 
                                      y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

ashn_WGP_12_18_resid_plot = ggplot(ashn_WGP_12_18_residuals, 
                                  aes(x = factor(trait), 
                                      y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

ashn_WGP_18_12_resid_plot = ggplot(ashn_WGP_18_12_residuals, 
                                  aes(x = factor(trait), 
                                      y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

ashn_WGP_18_18_resid_plot = ggplot(ashn_WGP_18_18_residuals, 
                                  aes(x = factor(trait), 
                                      y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

ashn_exp_treat_combo = (ashn_WGP_12_12_resid_plot+ashn_WGP_12_18_resid_plot)/(ashn_WGP_18_12_resid_plot+ashn_WGP_18_18_resid_plot)

ashn_WGP_12_18_resid_plot+ashn_WGP_18_12_resid_plot+ashn_WGP_18_18_resid_plot

# MYV ---------------------------------------------------------------------
myvc_WGP_12_12_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'MYVC', 
         Exp_treats == '12_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
myvc_WGP_12_18_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'MYVC', 
         Exp_treats == '12_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
myvc_WGP_18_12_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'MYVC', 
         Exp_treats == '18_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
myvc_WGP_18_18_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'MYVC', 
         Exp_treats == '18_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)


myvw_WGP_12_12_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'MYVW', 
         Exp_treats == '12_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
myvw_WGP_12_18_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'MYVW', 
         Exp_treats == '12_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
myvw_WGP_18_12_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'MYVW', 
         Exp_treats == '18_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
myvw_WGP_18_18_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'MYVW', 
         Exp_treats == '18_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

myvc_WGP_12_12_cov_mat = cov(myvc_WGP_12_12_trait_mat, 
                            use = 'complete.obs')
myvc_WGP_12_18_cov_mat = cov(myvc_WGP_12_18_trait_mat, 
                            use = 'complete.obs')
myvc_WGP_18_12_cov_mat = cov(myvc_WGP_18_12_trait_mat, 
                            use = 'complete.obs')
myvc_WGP_18_18_cov_mat = cov(myvc_WGP_18_18_trait_mat, 
                            use = 'complete.obs')

myvw_WGP_12_12_cov_mat = cov(myvw_WGP_12_12_trait_mat, 
                            use = 'complete.obs')
myvw_WGP_12_18_cov_mat = cov(myvw_WGP_12_18_trait_mat, 
                            use = 'complete.obs')
myvw_WGP_18_12_cov_mat = cov(myvw_WGP_18_12_trait_mat, 
                            use = 'complete.obs')
myvw_WGP_18_18_cov_mat = cov(myvw_WGP_18_18_trait_mat, 
                            use = 'complete.obs')

myv_WGP_12_12_krzcov = evolqg::KrzCor(myvc_WGP_12_12_cov_mat, 
                                     myvw_WGP_12_12_cov_mat)

myv_WGP_12_18_krzcov = evolqg::KrzCor(myvc_WGP_12_18_cov_mat, 
                                     myvw_WGP_12_18_cov_mat)

myv_WGP_18_12_krzcov = evolqg::KrzCor(myvc_WGP_18_12_cov_mat, 
                                     myvw_WGP_18_12_cov_mat)

myv_WGP_18_18_krzcov = evolqg::KrzCor(myvc_WGP_18_18_cov_mat, 
                                     myvw_WGP_18_18_cov_mat)


myv_WGP_12_12_proc_res = procrustes(myvc_WGP_12_12_cov_mat, 
                                   myvw_WGP_12_12_cov_mat)
myv_WGP_12_18_proc_res = procrustes(myvc_WGP_12_18_cov_mat, 
                                   myvw_WGP_12_18_cov_mat)
myv_WGP_18_12_proc_res = procrustes(myvc_WGP_18_12_cov_mat, 
                                   myvw_WGP_18_12_cov_mat)
myv_WGP_18_18_proc_res = procrustes(myvc_WGP_18_18_cov_mat, 
                                   myvw_WGP_18_18_cov_mat)

summary(myv_WGP_12_12_proc_res)
summary(myv_WGP_12_18_proc_res)
summary(myv_WGP_18_12_proc_res)
summary(myv_WGP_18_18_proc_res)

myv_WGP_12_12_proc_test = protest(myvc_WGP_12_12_cov_mat, 
                                 myvw_WGP_12_12_cov_mat, 
                                 permutations = 999)
myv_WGP_12_18_proc_test = protest(myvc_WGP_12_18_cov_mat, 
                                 myvw_WGP_12_18_cov_mat, 
                                 permutations = 999)
myv_WGP_18_12_proc_test = protest(myvc_WGP_18_12_cov_mat, 
                                 myvw_WGP_18_12_cov_mat, 
                                 permutations = 999)
myv_WGP_18_18_proc_test = protest(myvc_WGP_18_18_cov_mat, 
                                 myvw_WGP_18_18_cov_mat, 
                                 permutations = 999)


myv_WGP_12_12_residuals = residuals(myv_WGP_12_12_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

myv_WGP_12_18_residuals = residuals(myv_WGP_12_18_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

myv_WGP_18_12_residuals = residuals(myv_WGP_18_12_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

myv_WGP_18_18_residuals = residuals(myv_WGP_18_18_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)


myv_WGP_12_12_resid_plot = ggplot(myv_WGP_12_12_residuals, 
                                 aes(x = factor(trait), 
                                     y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

myv_WGP_12_18_resid_plot = ggplot(myv_WGP_12_18_residuals, 
                                 aes(x = factor(trait), 
                                     y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

myv_WGP_18_12_resid_plot = ggplot(myv_WGP_18_12_residuals, 
                                 aes(x = factor(trait), 
                                     y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

myv_WGP_18_18_resid_plot = ggplot(myv_WGP_18_18_residuals, 
                                 aes(x = factor(trait), 
                                     y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

myv_exp_treat_combo = (myv_WGP_12_12_resid_plot+myv_WGP_12_18_resid_plot)/(myv_WGP_18_12_resid_plot+myv_WGP_18_18_resid_plot)

# SKR ---------------------------------------------------------------------
skrc_WGP_12_12_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'SKRC', 
         Exp_treats == '12_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
skrc_WGP_12_18_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'SKRC', 
         Exp_treats == '12_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
skrc_WGP_18_12_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'SKRC', 
         Exp_treats == '18_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
skrc_WGP_18_18_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'SKRC', 
         Exp_treats == '18_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)


skrw_WGP_12_12_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'SKRW', 
         Exp_treats == '12_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
skrw_WGP_12_18_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'SKRW', 
         Exp_treats == '12_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
skrw_WGP_18_12_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'SKRW', 
         Exp_treats == '18_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
skrw_WGP_18_18_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'SKRW', 
         Exp_treats == '18_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

skrc_WGP_12_12_cov_mat = cov(skrc_WGP_12_12_trait_mat, 
                            use = 'complete.obs')
skrc_WGP_12_18_cov_mat = cov(skrc_WGP_12_18_trait_mat, 
                            use = 'complete.obs')
skrc_WGP_18_12_cov_mat = cov(skrc_WGP_18_12_trait_mat, 
                            use = 'complete.obs')
skrc_WGP_18_18_cov_mat = cov(skrc_WGP_18_18_trait_mat, 
                            use = 'complete.obs')

skrw_WGP_12_12_cov_mat = cov(skrw_WGP_12_12_trait_mat, 
                            use = 'complete.obs')
skrw_WGP_12_18_cov_mat = cov(skrw_WGP_12_18_trait_mat, 
                            use = 'complete.obs')
skrw_WGP_18_12_cov_mat = cov(skrw_WGP_18_12_trait_mat, 
                            use = 'complete.obs')
skrw_WGP_18_18_cov_mat = cov(skrw_WGP_18_18_trait_mat, 
                            use = 'complete.obs')

skr_WGP_12_12_krzcov = evolqg::KrzCor(skrc_WGP_12_12_cov_mat, 
                                     skrw_WGP_12_12_cov_mat)

skr_WGP_12_18_krzcov = evolqg::KrzCor(skrc_WGP_12_18_cov_mat, 
                                     skrw_WGP_12_18_cov_mat)

skr_WGP_18_12_krzcov = evolqg::KrzCor(skrc_WGP_18_12_cov_mat, 
                                     skrw_WGP_18_12_cov_mat)

skr_WGP_18_18_krzcov = evolqg::KrzCor(skrc_WGP_18_18_cov_mat, 
                                     skrw_WGP_18_18_cov_mat)


skr_WGP_12_12_proc_res = procrustes(skrc_WGP_12_12_cov_mat, 
                                   skrw_WGP_12_12_cov_mat)
skr_WGP_12_18_proc_res = procrustes(skrc_WGP_12_18_cov_mat, 
                                   skrw_WGP_12_18_cov_mat)
skr_WGP_18_12_proc_res = procrustes(skrc_WGP_18_12_cov_mat, 
                                   skrw_WGP_18_12_cov_mat)
skr_WGP_18_18_proc_res = procrustes(skrc_WGP_18_18_cov_mat, 
                                   skrw_WGP_18_18_cov_mat)

summary(skr_WGP_12_12_proc_res)
summary(skr_WGP_12_18_proc_res)
summary(skr_WGP_18_12_proc_res)
summary(skr_WGP_18_18_proc_res)

skr_WGP_12_12_proc_test = protest(skrc_WGP_12_12_cov_mat, 
                                 skrw_WGP_12_12_cov_mat, 
                                 permutations = 999)
skr_WGP_12_18_proc_test = protest(skrc_WGP_12_18_cov_mat, 
                                 skrw_WGP_12_18_cov_mat, 
                                 permutations = 999)
skr_WGP_18_12_proc_test = protest(skrc_WGP_18_12_cov_mat, 
                                 skrw_WGP_18_12_cov_mat, 
                                 permutations = 999)
skr_WGP_18_18_proc_test = protest(skrc_WGP_18_18_cov_mat, 
                                 skrw_WGP_18_18_cov_mat, 
                                 permutations = 999)


skr_WGP_12_12_residuals = residuals(skr_WGP_12_12_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

skr_WGP_12_18_residuals = residuals(skr_WGP_12_18_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

skr_WGP_18_12_residuals = residuals(skr_WGP_18_12_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

skr_WGP_18_18_residuals = residuals(skr_WGP_18_18_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)


skr_WGP_12_12_resid_plot = ggplot(skr_WGP_12_12_residuals, 
                                 aes(x = factor(trait), 
                                     y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

skr_WGP_12_18_resid_plot = ggplot(skr_WGP_12_18_residuals, 
                                 aes(x = factor(trait), 
                                     y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

skr_WGP_18_12_resid_plot = ggplot(skr_WGP_18_12_residuals, 
                                 aes(x = factor(trait), 
                                     y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

skr_WGP_18_18_resid_plot = ggplot(skr_WGP_18_18_residuals, 
                                 aes(x = factor(trait), 
                                     y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

skr_exp_treat_combo = (skr_WGP_12_12_resid_plot+skr_WGP_12_18_resid_plot)/(skr_WGP_18_12_resid_plot+skr_WGP_18_18_resid_plot)



# GAR/GTS -----------------------------------------------------------------
CSWYC_WGP_12_12_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'CSWYC', 
         Exp_treats == '12_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
CSWYC_WGP_12_18_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'CSWYC', 
         Exp_treats == '12_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
CSWYC_WGP_18_12_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'CSWYC', 
         Exp_treats == '18_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
CSWYC_WGP_18_18_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'CSWYC', 
         Exp_treats == '18_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)


GTSW_WGP_12_12_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'GTSW', 
         Exp_treats == '12_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
GTSW_WGP_12_18_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'GTSW', 
         Exp_treats == '12_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
GTSW_WGP_18_12_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'GTSW', 
         Exp_treats == '18_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
GTSW_WGP_18_18_trait_mat = WGP_data %>% 
  filter(Lake_morph == 'GTSW', 
         Exp_treats == '18_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

CSWYC_WGP_12_12_cov_mat = cov(CSWYC_WGP_12_12_trait_mat, 
                             use = 'complete.obs')
CSWYC_WGP_12_18_cov_mat = cov(CSWYC_WGP_12_18_trait_mat, 
                             use = 'complete.obs')
CSWYC_WGP_18_12_cov_mat = cov(CSWYC_WGP_18_12_trait_mat, 
                             use = 'complete.obs')
CSWYC_WGP_18_18_cov_mat = cov(CSWYC_WGP_18_18_trait_mat, 
                             use = 'complete.obs')

GTSW_WGP_12_12_cov_mat = cov(GTSW_WGP_12_12_trait_mat, 
                            use = 'complete.obs')
GTSW_WGP_12_18_cov_mat = cov(GTSW_WGP_12_18_trait_mat, 
                            use = 'complete.obs')
GTSW_WGP_18_12_cov_mat = cov(GTSW_WGP_18_12_trait_mat, 
                            use = 'complete.obs')
GTSW_WGP_18_18_cov_mat = cov(GTSW_WGP_18_18_trait_mat, 
                            use = 'complete.obs')

gts_gar_WGP_12_12_krzcov = evolqg::KrzCor(CSWYC_WGP_12_12_cov_mat, 
                                         GTSW_WGP_12_12_cov_mat)

gts_gar_WGP_12_18_krzcov = evolqg::KrzCor(CSWYC_WGP_12_18_cov_mat, 
                                         GTSW_WGP_12_18_cov_mat)

gts_gar_WGP_18_12_krzcov = evolqg::KrzCor(CSWYC_WGP_18_12_cov_mat, 
                                         GTSW_WGP_18_12_cov_mat)

gts_gar_WGP_18_18_krzcov = evolqg::KrzCor(CSWYC_WGP_18_18_cov_mat, 
                                         GTSW_WGP_18_18_cov_mat)


gts_gar_WGP_12_12_proc_res = procrustes(CSWYC_WGP_12_12_cov_mat, 
                                       GTSW_WGP_12_12_cov_mat)
gts_gar_WGP_12_18_proc_res = procrustes(CSWYC_WGP_12_18_cov_mat, 
                                       GTSW_WGP_12_18_cov_mat)
gts_gar_WGP_18_12_proc_res = procrustes(CSWYC_WGP_18_12_cov_mat, 
                                       GTSW_WGP_18_12_cov_mat)
gts_gar_WGP_18_18_proc_res = procrustes(CSWYC_WGP_18_18_cov_mat, 
                                       GTSW_WGP_18_18_cov_mat)

summary(gts_gar_WGP_12_12_proc_res)
summary(gts_gar_WGP_12_18_proc_res)
summary(gts_gar_WGP_18_12_proc_res)
summary(gts_gar_WGP_18_18_proc_res)

gts_gar_WGP_12_12_proc_test = protest(CSWYC_WGP_12_12_cov_mat, 
                                     GTSW_WGP_12_12_cov_mat, 
                                     permutations = 999)
gts_gar_WGP_12_18_proc_test = protest(CSWYC_WGP_12_18_cov_mat, 
                                     GTSW_WGP_12_18_cov_mat, 
                                     permutations = 999)
gts_gar_WGP_18_12_proc_test = protest(CSWYC_WGP_18_12_cov_mat, 
                                     GTSW_WGP_18_12_cov_mat, 
                                     permutations = 999)
gts_gar_WGP_18_18_proc_test = protest(CSWYC_WGP_18_18_cov_mat, 
                                     GTSW_WGP_18_18_cov_mat, 
                                     permutations = 999)


gts_gar_WGP_12_12_residuals = residuals(gts_gar_WGP_12_12_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

gts_gar_WGP_12_18_residuals = residuals(gts_gar_WGP_12_18_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

gts_gar_WGP_18_12_residuals = residuals(gts_gar_WGP_18_12_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

gts_gar_WGP_18_18_residuals = residuals(gts_gar_WGP_18_18_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)


gts_gar_WGP_12_12_resid_plot = ggplot(gts_gar_WGP_12_12_residuals, 
                                     aes(x = factor(trait), 
                                         y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

gts_gar_WGP_12_18_resid_plot = ggplot(gts_gar_WGP_12_18_residuals, 
                                     aes(x = factor(trait), 
                                         y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

gts_gar_WGP_18_12_resid_plot = ggplot(gts_gar_WGP_18_12_residuals, 
                                     aes(x = factor(trait), 
                                         y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

gts_gar_WGP_18_18_resid_plot = ggplot(gts_gar_WGP_18_18_residuals, 
                                     aes(x = factor(trait), 
                                         y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

gts_gar_exp_treat_combo = (gts_gar_WGP_12_12_resid_plot+gts_gar_WGP_12_18_resid_plot)/(gts_gar_WGP_18_12_resid_plot+gts_gar_WGP_18_18_resid_plot)



# Parallelism in plasticity divergence in trait covariance ----------------

## kendals test of concordance 
##Rows are objects and columns are judges (variables)
## Traits in rows and pop pairs in columns

ashn_WGP_12_12_kendal = ashn_WGP_12_12_residuals %>% 
  dplyr::rename(ASHN = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

myv_WGP_12_12_kendal = myv_WGP_12_12_residuals%>% 
  dplyr::rename(MYV = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

skr_WGP_12_12_kendal = skr_WGP_12_12_residuals%>% 
  dplyr::rename(SKR = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

gts_gar_WGP_12_12_kendal = gts_gar_WGP_12_12_residuals%>% 
  dplyr::rename(GTS_GAR = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

WGP_12_12_kendal_df = bind_cols(
  # ashn_WGP_12_12_kendal, 
                               myv_WGP_12_12_kendal, 
                               skr_WGP_12_12_kendal, 
                               gts_gar_WGP_12_12_kendal)

# group = c('ASHN','MYV','SKR','GTS_GAR')
group = c('MYV','SKR','GTS_GAR')

WGP_12_12_kendal = kendall.global(WGP_12_12_kendal_df,
                                 # group = pops,
                                 nperm = 999, 
                                 mult = 'hochberg')
WGP_12_12_kendal_post = kendall.post(WGP_12_12_kendal_df, 
                                    mult = 'hochberg', 
                                    nperm = 999)


ashn_WGP_12_18_kendal = ashn_WGP_12_18_residuals %>% 
  dplyr::rename(ASHN = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

myv_WGP_12_18_kendal = myv_WGP_12_18_residuals%>% 
  dplyr::rename(MYV = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

skr_WGP_12_18_kendal = skr_WGP_12_18_residuals%>% 
  dplyr::rename(SKR = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

gts_gar_WGP_12_18_kendal = gts_gar_WGP_12_18_residuals%>% 
  dplyr::rename(GTS_GAR = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

WGP_12_18_kendal_df = bind_cols(ashn_WGP_12_18_kendal, 
                               myv_WGP_12_18_kendal, 
                               skr_WGP_12_18_kendal, 
                               gts_gar_WGP_12_18_kendal)


WGP_12_18_kendal = kendall.global(WGP_12_18_kendal_df,
                                 # group = pops,
                                 nperm = 999, 
                                 mult = 'hochberg')
WGP_12_18_kendal_post = kendall.post(WGP_12_18_kendal_df, 
                                    mult = 'hochberg', 
                                    nperm = 999)



ashn_WGP_18_12_kendal = ashn_WGP_18_12_residuals %>% 
  dplyr::rename(ASHN = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

myv_WGP_18_12_kendal = myv_WGP_18_12_residuals%>% 
  dplyr::rename(MYV = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

skr_WGP_18_12_kendal = skr_WGP_18_12_residuals%>% 
  dplyr::rename(SKR = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

gts_gar_WGP_18_12_kendal = gts_gar_WGP_18_12_residuals%>% 
  dplyr::rename(GTS_GAR = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

WGP_18_12_kendal_df = bind_cols(ashn_WGP_18_12_kendal, 
                               myv_WGP_18_12_kendal, 
                               skr_WGP_18_12_kendal, 
                               gts_gar_WGP_18_12_kendal)


WGP_18_12_kendal = kendall.global(WGP_18_12_kendal_df,
                                 # group = pops,
                                 nperm = 999, 
                                 mult = 'hochberg')
WGP_18_12_kendal_post = kendall.post(WGP_18_12_kendal_df, 
                                    mult = 'hochberg', 
                                    nperm = 999)




ashn_WGP_18_18_kendal = ashn_WGP_18_18_residuals %>% 
  dplyr::rename(ASHN = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

myv_WGP_18_18_kendal = myv_WGP_18_18_residuals%>% 
  dplyr::rename(MYV = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

skr_WGP_18_18_kendal = skr_WGP_18_18_residuals%>% 
  dplyr::rename(SKR = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

gts_gar_WGP_18_18_kendal = gts_gar_WGP_18_18_residuals%>% 
  dplyr::rename(GTS_GAR = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

WGP_18_18_kendal_df = bind_cols(ashn_WGP_18_18_kendal, 
                               myv_WGP_18_18_kendal, 
                               skr_WGP_18_18_kendal, 
                               gts_gar_WGP_18_18_kendal)

WGP_18_18_kendal = kendall.global(WGP_18_18_kendal_df,
                                 # group = pops,
                                 nperm = 999, 
                                 mult = 'hochberg')
WGP_18_18_kendal_post = kendall.post(WGP_18_18_kendal_df, 
                                    mult = 'hochberg', 
                                    nperm = 999)


WGP_12_12_kendal_post
WGP_12_18_kendal_post
WGP_18_12_kendal_post
WGP_18_18_kendal_post




# Nemeny post hoc test ----------------------------------------------------

library(PMCMRplus)


ashn_WGP_12_12_nemeny_df = ashn_WGP_12_12_residuals %>% 
  mutate(pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10)

myv_WGP_12_12_nemeny_df = myv_WGP_12_12_residuals %>% 
  mutate(pop_pair = 'MYV') %>% 
  dplyr::select(-resid_10)

skr_WGP_12_12_nemeny_df = skr_WGP_12_12_residuals %>% 
  mutate(pop_pair = 'SKR') %>% 
  dplyr::select(-resid_10)

gts_gar_WGP_12_12_nemeny_df = gts_gar_WGP_12_12_residuals %>% 
  mutate(pop_pair = 'GTS/GAR') %>% 
  dplyr::select(-resid_10)

WGP_12_12_nemeny_df = bind_rows(
  # ashn_WGP_12_12_nemeny_df, 
                               myv_WGP_12_12_nemeny_df, 
                               skr_WGP_12_12_nemeny_df, 
                               gts_gar_WGP_12_12_nemeny_df)


WGP_12_12_nemeny_test = frdAllPairsNemenyiTest(y = nemeny_df$residual, 
                                              groups = nemeny_df$trait, 
                                              blocks = nemeny_df$pop_pair, 
                                              p.adjust.method = 'hochberg')

WGP_12_12_nemeny_test$p.value %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname) %>% 
  # filter(across(c(2:33), function(x) x<0.05))
  write_csv('WGP_1212_residual_NemenyTest_results.csv')



ashn_WGP_12_18_nemeny_df = ashn_WGP_12_18_residuals %>% 
  mutate(pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10)

myv_WGP_12_18_nemeny_df = myv_WGP_12_18_residuals %>% 
  mutate(pop_pair = 'MYV') %>% 
  dplyr::select(-resid_10)

skr_WGP_12_18_nemeny_df = skr_WGP_12_18_residuals %>% 
  mutate(pop_pair = 'SKR') %>% 
  dplyr::select(-resid_10)

gts_gar_WGP_12_18_nemeny_df = gts_gar_WGP_12_18_residuals %>% 
  mutate(pop_pair = 'GTS/GAR') %>% 
  dplyr::select(-resid_10)

WGP_12_18_nemeny_df = bind_rows(ashn_WGP_12_18_nemeny_df, 
                               myv_WGP_12_18_nemeny_df, 
                               skr_WGP_12_18_nemeny_df, 
                               gts_gar_WGP_12_18_nemeny_df)


WGP_12_18_nemeny_test = frdAllPairsNemenyiTest(y = WGP_12_18_nemeny_df$residual, 
                                              groups = WGP_12_18_nemeny_df$trait, 
                                              blocks = WGP_12_18_nemeny_df$pop_pair, 
                                              p.adjust.method = 'hochberg')

WGP_12_18_nemeny_test$p.value %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname) %>% 
  # filter(across(c(2:33), function(x) x<0.05))
write_csv('WGP_1218_residual_NemenyTest_results.csv')

ashn_WGP_18_12_nemeny_df = ashn_WGP_18_12_residuals %>% 
  mutate(pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10)

myv_WGP_18_12_nemeny_df = myv_WGP_18_12_residuals %>% 
  mutate(pop_pair = 'MYV') %>% 
  dplyr::select(-resid_10)

skr_WGP_18_12_nemeny_df = skr_WGP_18_12_residuals %>% 
  mutate(pop_pair = 'SKR') %>% 
  dplyr::select(-resid_10)

gts_gar_WGP_18_12_nemeny_df = gts_gar_WGP_18_12_residuals %>% 
  mutate(pop_pair = 'GTS/GAR') %>% 
  dplyr::select(-resid_10)

WGP_18_12_nemeny_df = bind_rows(ashn_WGP_18_12_nemeny_df, 
                               myv_WGP_18_12_nemeny_df, 
                               skr_WGP_18_12_nemeny_df, 
                               gts_gar_WGP_18_12_nemeny_df)


WGP_18_12_nemeny_test = frdAllPairsNemenyiTest(y = WGP_18_12_nemeny_df$residual, 
                                              groups = WGP_18_12_nemeny_df$trait, 
                                              blocks = WGP_18_12_nemeny_df$pop_pair, 
                                              p.adjust.method = 'hochberg')


WGP_18_12_nemeny_test$p.value %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname) %>% 
  # filter(across(c(2:33), function(x) x<0.05))
write_csv('WGP_1812_residual_NemenyTest_results.csv')



ashn_WGP_18_18_nemeny_df = ashn_WGP_18_18_residuals %>% 
  mutate(pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10)

myv_WGP_18_18_nemeny_df = myv_WGP_18_18_residuals %>% 
  mutate(pop_pair = 'MYV') %>% 
  dplyr::select(-resid_10)

skr_WGP_18_18_nemeny_df = skr_WGP_18_18_residuals %>% 
  mutate(pop_pair = 'SKR') %>% 
  dplyr::select(-resid_10)

gts_gar_WGP_18_18_nemeny_df = gts_gar_WGP_18_18_residuals %>% 
  mutate(pop_pair = 'GTS/GAR') %>% 
  dplyr::select(-resid_10)

WGP_18_18_nemeny_df = bind_rows(ashn_WGP_18_18_nemeny_df, 
                               myv_WGP_18_18_nemeny_df, 
                               skr_WGP_18_18_nemeny_df, 
                               gts_gar_WGP_18_18_nemeny_df)


WGP_18_18_nemeny_test = frdAllPairsNemenyiTest(y = WGP_18_18_nemeny_df$residual, 
                                              groups = WGP_18_18_nemeny_df$trait, 
                                              blocks = WGP_18_18_nemeny_df$pop_pair, 
                                              p.adjust.method = 'hochberg')


WGP_18_18_nemeny_test$p.value %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname) %>% 
  # filter(across(c(2:33), function(x) x<0.05))
write_csv('WGP_1818_residual_NemenyTest_results.csv')

# parallelism test - eye ball method --------------------------------------


(ashn_WGP_12_12_resid_plot+myv_WGP_12_12_resid_plot)/(skr_WGP_12_12_resid_plot+gts_gar_WGP_12_12_resid_plot)

(ashn_WGP_12_18_resid_plot+myv_WGP_12_18_resid_plot)/(skr_WGP_12_18_resid_plot+gts_gar_WGP_12_18_resid_plot)
## Opercular KT divergence plastically divergent between geo and amb cov mats 
## across all pops
(ashn_WGP_18_12_resid_plot+myv_WGP_18_12_resid_plot)/(skr_WGP_18_12_resid_plot+gts_gar_WGP_18_12_resid_plot)

(ashn_WGP_18_18_resid_plot+myv_WGP_18_18_resid_plot)/(skr_WGP_18_18_resid_plot+gts_gar_WGP_18_18_resid_plot)


# Combine trait covariance matrices ---------------------------------------


# ASHN covariance mats ----------------------------------------------------


ashn_WGP_cov_matrices <- list(
  # Temp1_Ambient = ashnc_WGP_12_12_cov_mat,
  # Temp1_Geothermal = ashnw_WGP_12_12_cov_mat,
  Temp2_Ambient = ashnc_WGP_12_18_cov_mat,
  Temp2_Geothermal = ashnw_WGP_12_18_cov_mat,
  Temp3_Ambient = ashnc_WGP_18_12_cov_mat,
  Temp3_Geothermal = ashnw_WGP_18_12_cov_mat,
  Temp4_Ambient = ashnc_WGP_18_18_cov_mat,
  Temp4_Geothermal = ashnw_WGP_18_18_cov_mat)

for (i in 2:4) {
  ashn_WGP_cov_matrices[[paste0("Temp", i, "_Difference")]] <-
    ashn_WGP_cov_matrices[[paste0("Temp", i, "_Geothermal")]] -
    ashn_WGP_cov_matrices[[paste0("Temp", i, "_Ambient")]]
}

ashn_WGP_cov_long <- purrr::imap_dfr(ashn_WGP_cov_matrices, function(mat, name) {
  as.data.frame(as.table(mat)) %>%
    setNames(c("Trait1", "Trait2", "Covariance")) %>%
    mutate(
      Temp = sub("_.*", "", name),
      Type = sub(".*_", "", name)
    )
})

ashn_WGP_cov_long$Type <- factor(
  ashn_WGP_cov_long$Type,
  levels = c("Ambient", "Geothermal", "Difference")
)

ashn_WGP_cov_long = ashn_WGP_cov_long %>% 
  mutate(Temp_info = case_when(
    # str_detect(Temp, 'Temp1') ~ '12_12',
    str_detect(Temp, 'Temp2') ~ '12_18',
    str_detect(Temp, 'Temp3') ~ '18_12',
    str_detect(Temp, 'Temp4') ~ '18_18'))

ashn_WGP_cov_plot = ggplot(ashn_WGP_cov_long, 
                       aes(Trait1, 
                           Trait2, 
                           fill = Covariance)) +
  geom_tile() +
  facet_grid(Temp_info ~ Type) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0) +
  coord_equal() +
  labs(x = NULL, 
       y = NULL, 
       fill = "Covariance", 
       title = 'ASHN WGP trait covariance')+
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1),
        panel.grid = element_blank(), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'), 
        axis.text = element_text(size = 6), 
        legend.position = 'none')

# MYV covariance mats ----------------------------------------------------


myv_WGP_cov_matrices <- list(
  Temp1_Ambient = myvc_WGP_12_12_cov_mat,
  Temp1_Geothermal = myvw_WGP_12_12_cov_mat,
  Temp2_Ambient = myvc_WGP_12_18_cov_mat,
  Temp2_Geothermal = myvw_WGP_12_18_cov_mat,
  Temp3_Ambient = myvc_WGP_18_12_cov_mat,
  Temp3_Geothermal = myvw_WGP_18_12_cov_mat,
  Temp4_Ambient = myvc_WGP_18_18_cov_mat,
  Temp4_Geothermal = myvw_WGP_18_18_cov_mat)

for (i in 1:4) {
  myv_WGP_cov_matrices[[paste0("Temp", i, "_Difference")]] <-
    myv_WGP_cov_matrices[[paste0("Temp", i, "_Geothermal")]] -
    myv_WGP_cov_matrices[[paste0("Temp", i, "_Ambient")]]
}

myv_WGP_cov_long <- purrr::imap_dfr(myv_WGP_cov_matrices, function(mat, name) {
  as.data.frame(as.table(mat)) %>%
    setNames(c("Trait1", "Trait2", "Covariance")) %>%
    mutate(
      Temp = sub("_.*", "", name),
      Type = sub(".*_", "", name)
    )
})

myv_WGP_cov_long$Type <- factor(
  myv_WGP_cov_long$Type,
  levels = c("Ambient", "Geothermal", "Difference")
)

myv_WGP_cov_long = myv_WGP_cov_long %>% 
  mutate(Temp_info = case_when(
    str_detect(Temp, 'Temp1') ~ '12_12',
    str_detect(Temp, 'Temp2') ~ '12_18',
    str_detect(Temp, 'Temp3') ~ '18_12',
    str_detect(Temp, 'Temp4') ~ '18_18'))

myv_WGP_cov_plot = ggplot(myv_WGP_cov_long, 
                      aes(Trait1, 
                          Trait2, 
                          fill = Covariance)) +
  geom_tile() +
  facet_grid(Temp_info ~ Type) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0) +
  coord_equal() +
  labs(x = NULL, 
       y = NULL, 
       fill = "Covariance", 
       title = 'MYV WGP trait covariance')+
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1),
        panel.grid = element_blank(), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'), 
        axis.text = element_text(size = 6), 
        legend.position = 'none')

# skr covariance mats ----------------------------------------------------


skr_WGP_cov_matrices <- list(
  Temp1_Ambient = skrc_WGP_12_12_cov_mat,
  Temp1_Geothermal = skrw_WGP_12_12_cov_mat,
  Temp2_Ambient = skrc_WGP_12_18_cov_mat,
  Temp2_Geothermal = skrw_WGP_12_18_cov_mat,
  Temp3_Ambient = skrc_WGP_18_12_cov_mat,
  Temp3_Geothermal = skrw_WGP_18_12_cov_mat,
  Temp4_Ambient = skrc_WGP_18_18_cov_mat,
  Temp4_Geothermal = skrw_WGP_18_18_cov_mat)

for (i in 1:4) {
  skr_WGP_cov_matrices[[paste0("Temp", i, "_Difference")]] <-
    skr_WGP_cov_matrices[[paste0("Temp", i, "_Geothermal")]] -
    skr_WGP_cov_matrices[[paste0("Temp", i, "_Ambient")]]
}

skr_WGP_cov_long <- purrr::imap_dfr(skr_WGP_cov_matrices, function(mat, name) {
  as.data.frame(as.table(mat)) %>%
    setNames(c("Trait1", "Trait2", "Covariance")) %>%
    mutate(
      Temp = sub("_.*", "", name),
      Type = sub(".*_", "", name)
    )
})

skr_WGP_cov_long$Type <- factor(
  skr_WGP_cov_long$Type,
  levels = c("Ambient", "Geothermal", "Difference")
)

skr_WGP_cov_long = skr_WGP_cov_long %>% 
  mutate(Temp_info = case_when(
    str_detect(Temp, 'Temp1') ~ '12_12',
    str_detect(Temp, 'Temp2') ~ '12_18',
    str_detect(Temp, 'Temp3') ~ '18_12',
    str_detect(Temp, 'Temp4') ~ '18_18'))

skr_WGP_cov_plot = ggplot(skr_WGP_cov_long, 
                      aes(Trait1, 
                          Trait2, 
                          fill = Covariance)) +
  geom_tile() +
  facet_grid(Temp_info ~ Type) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0) +
  coord_equal() +
  labs(x = NULL, 
       y = NULL, 
       fill = "Covariance", 
       title = 'SKR WGP trait covariance')+
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1),
        panel.grid = element_blank(), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'), 
        axis.text = element_text(size = 6), 
        legend.position = 'none')

# gts_gar covariance mats ----------------------------------------------------


gts_gar_WGP_cov_matrices <- list(
  Temp1_Ambient = CSWYC_WGP_12_12_cov_mat,
  Temp1_Geothermal = GTSW_WGP_12_12_cov_mat,
  Temp2_Ambient = CSWYC_WGP_12_18_cov_mat,
  Temp2_Geothermal = GTSW_WGP_12_18_cov_mat,
  Temp3_Ambient = CSWYC_WGP_18_12_cov_mat,
  Temp3_Geothermal = GTSW_WGP_18_12_cov_mat,
  Temp4_Ambient = CSWYC_WGP_18_18_cov_mat,
  Temp4_Geothermal = GTSW_WGP_18_18_cov_mat)

for (i in 1:4) {
  gts_gar_WGP_cov_matrices[[paste0("Temp", i, "_Difference")]] <-
    gts_gar_WGP_cov_matrices[[paste0("Temp", i, "_Geothermal")]] -
    gts_gar_WGP_cov_matrices[[paste0("Temp", i, "_Ambient")]]
}

gts_gar_WGP_cov_long <- purrr::imap_dfr(gts_gar_WGP_cov_matrices, function(mat, name) {
  as.data.frame(as.table(mat)) %>%
    setNames(c("Trait1", "Trait2", "Covariance")) %>%
    mutate(
      Temp = sub("_.*", "", name),
      Type = sub(".*_", "", name)
    )
})

gts_gar_WGP_cov_long$Type <- factor(
  gts_gar_WGP_cov_long$Type,
  levels = c("Ambient", "Geothermal", "Difference")
)

gts_gar_WGP_cov_long = gts_gar_WGP_cov_long %>% 
  mutate(Temp_info = case_when(
    str_detect(Temp, 'Temp1') ~ '12_12',
    str_detect(Temp, 'Temp2') ~ '12_18',
    str_detect(Temp, 'Temp3') ~ '18_12',
    str_detect(Temp, 'Temp4') ~ '18_18'))

gts_gar_WGP_cov_plot = ggplot(gts_gar_WGP_cov_long, 
                          aes(Trait1, 
                              Trait2, 
                              fill = Covariance)) +
  geom_tile() +
  facet_grid(Temp_info ~ Type) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0) +
  coord_equal() +
  labs(x = NULL, 
       y = NULL, 
       fill = "Covariance", 
       title = 'GTS/GAR WGP trait covariance')+
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1),
        panel.grid = element_blank(), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'), 
        axis.text = element_text(size = 6), 
        legend.position = 'none')



(ashn_cov_plot + myv_cov_plot)/(skr_cov_plot + gts_gar_cov_plot)

(ashn_cov_plot + myv_cov_plot + skr_cov_plot + gts_gar_cov_plot)
