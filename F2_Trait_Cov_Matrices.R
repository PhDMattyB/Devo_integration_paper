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



f2_orig_data = read_csv("F2_Original_univariate_traits_FIXED_11.02.2026.csv") %>% 
  dplyr::select(Lake_morph, 
                Ecotype_pair, 
                lake_morph_Pair_Full_Temp, 
                Morph, 
                Offspring_temp, 
                Parent_temp, 
                jaw_length:caudal2_15_17) %>% 
  unite(col = Exp_treats, 
        c('Parent_temp', 
          'Offspring_temp'), 
        sep = '_')


# Redoing the analysis due to bad instructions but I hate deleting --------


# Trait cov matrices ------------------------------------------------------


f2_cov_matrices <- list(
  Pair1_Ambient = ashnc_f2_cov_mat,
  Pair1_Geothermal = ashnw_f2_cov_mat,
  Pair2_Ambient = myvc_f2_cov_mat,
  Pair2_Geothermal = myvw_f2_cov_mat,
  Pair3_Ambient = skrc_f2_cov_mat,
  Pair3_Geothermal = skrw_f2_cov_mat,
  Pair4_Ambient = gar_f2_cov_mat,
  Pair4_Geothermal = gts_f2_cov_mat)

for (i in 1:4) {
  f2_cov_matrices[[paste0("Pair", i, "_Difference")]] <-
    f2_cov_matrices[[paste0("Pair", i, "_Geothermal")]] -
    f2_cov_matrices[[paste0("Pair", i, "_Ambient")]]
}

f2_cov_long <- purrr::imap_dfr(f2_cov_matrices, function(mat, name) {
  as.data.frame(as.table(mat)) %>%
    setNames(c("Trait1", "Trait2", "Covariance")) %>%
    mutate(
      Temp = sub("_.*", "", name),
      Type = sub(".*_", "", name)
    )
})

f2_cov_long$Type <- factor(
  f2_cov_long$Type,
  levels = c("Ambient", "Geothermal", "Difference")
)

f2_cov_long = f2_cov_long %>% 
  mutate(pop_pairs = case_when(
    str_detect(Temp, 'Pair1') ~ 'ASHN',
    str_detect(Temp, 'Pair2') ~ 'MYV',
    str_detect(Temp, 'Pair3') ~ 'SKR',
    str_detect(Temp, 'Pair4') ~ 'GTS/GAR'))

f2_cov_plot = ggplot(f2_cov_long, 
                       aes(Trait1, 
                           Trait2, 
                           fill = Covariance)) +
  geom_tile() +
  facet_grid(pop_pairs ~ Type) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0) +
  coord_equal() +
  labs(x = NULL, 
       y = NULL, 
       fill = "Covariance", 
       title = 'F2 trait covariance matrices')+
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1),
        panel.grid = element_blank(), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'), 
        axis.text = element_text(size = 6), 
        legend.position = 'none')


# Procrustes superimposition - cov matrices -------------------------------



# ASHN --------------------------------------------------------------------


ashnc_f2_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'ASHNC') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

ashnw_f2_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'ASHNW') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

ashnc_f2_cov_mat = cov(ashnc_f2_trait_mat, 
                             use = 'complete.obs')
ashnw_f2_cov_mat = cov(ashnw_f2_trait_mat, 
                             use = 'complete.obs')
ashn_f2_krzcov = evolqg::KrzCor(ashnc_f2_cov_mat, 
                                      ashnw_f2_cov_mat)
ashn_f2_proc_res = procrustes(ashnc_f2_cov_mat, 
                                    ashnw_f2_cov_mat)

plot(ashn_f2_proc_res, 
     kind = 1)

ashn_f2_proc_test = protest(ashnc_f2_cov_mat, 
                                  ashnw_f2_cov_mat, 
                                  permutations = 999)


ashn_f2_residuals = residuals(ashn_f2_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                cov_residual = 2)

# %>% 
#   mutate(resid_10 = residual*10)

ashn_f2_resid_plot = ggplot(ashn_f2_residuals, 
                                  aes(x = factor(trait), 
                                      y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  # geom_hline(yintercept = 1.5, 
  #            linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))


# MYV ---------------------------------------------------------------------

myvc_f2_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'MYVC') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

myvw_f2_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'MYVW') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

myvc_f2_cov_mat = cov(myvc_f2_trait_mat, 
                       use = 'complete.obs')
myvw_f2_cov_mat = cov(myvw_f2_trait_mat, 
                       use = 'complete.obs')
myv_f2_krzcov = evolqg::KrzCor(myvc_f2_cov_mat, 
                                myvw_f2_cov_mat)
myv_f2_proc_res = procrustes(myvc_f2_cov_mat, 
                              myvw_f2_cov_mat)

plot(myv_f2_proc_res, 
     kind = 1)

myv_f2_proc_test = protest(myvc_f2_cov_mat, 
                            myvw_f2_cov_mat, 
                            permutations = 999)


myv_f2_residuals = residuals(myv_f2_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                cov_residual = 2)

myv_f2_resid_plot = ggplot(myv_f2_residuals, 
                            aes(x = factor(trait), 
                                y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  # geom_hline(yintercept = 1.5, 
  #            linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))


# SKR ---------------------------------------------------------------------

skrc_f2_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'SKRC') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

skrw_f2_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'SKRW') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

skrc_f2_cov_mat = cov(skrc_f2_trait_mat, 
                       use = 'complete.obs')
skrw_f2_cov_mat = cov(skrw_f2_trait_mat, 
                       use = 'complete.obs')
skr_f2_krzcov = evolqg::KrzCor(skrc_f2_cov_mat, 
                                skrw_f2_cov_mat)
skr_f2_proc_res = procrustes(skrc_f2_cov_mat, 
                              skrw_f2_cov_mat)

plot(skr_f2_proc_res, 
     kind = 1)

skr_f2_proc_test = protest(skrc_f2_cov_mat, 
                            skrw_f2_cov_mat, 
                            permutations = 999)


skr_f2_residuals = residuals(skr_f2_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                cov_residual = 2)

skr_f2_resid_plot = ggplot(skr_f2_residuals, 
                            aes(x = factor(trait), 
                                y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  # geom_hline(yintercept = 1.5, 
  #            linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))


# GTS/GAR -----------------------------------------------------------------

gar_f2_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'CSWYC') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

gts_f2_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'GTSW') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

gar_f2_cov_mat = cov(gar_f2_trait_mat, 
                       use = 'complete.obs')
gts_f2_cov_mat = cov(gts_f2_trait_mat, 
                       use = 'complete.obs')
gts_gar_f2_krzcov = evolqg::KrzCor(gar_f2_cov_mat, 
                                gts_f2_cov_mat)
gts_gar_f2_proc_res = procrustes(gar_f2_cov_mat, 
                              gts_f2_cov_mat)

plot(gts_gar_f2_proc_res, 
     kind = 1)

gts_gar_f2_proc_test = protest(gar_f2_cov_mat, 
                            gts_f2_cov_mat, 
                            permutations = 999)


gts_gar_f2_residuals = residuals(gts_gar_f2_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                cov_residual = 2)

gts_gar_f2_resid_plot = ggplot(gts_gar_f2_residuals, 
                            aes(x = factor(trait), 
                                y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  # geom_hline(yintercept = 1.5, 
  #            linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))








# Procrustes superimposition - eigenvectors -------------------------------


# ASHN ---------------------------------------------------------------------


ashnc_f2_eigen = eigen(ashnc_f2_cov_mat)
ashnw_f2_eigen = eigen(ashnw_f2_cov_mat)

# Extract leading eigenvectors
ashnc_f2_vectors = ashnc_f2_eigen$vectors[, 1:17, drop = FALSE]
ashnw_f2_vectors = ashnw_f2_eigen$vectors[, 1:17, drop = FALSE]

# Preserve trait identities
rownames(ashnc_f2_vectors) = colnames(ashnc_f2_cov_mat)
rownames(ashnw_f2_vectors) = colnames(ashnw_f2_cov_mat)


ashn_f2_proc_res_eigen = procrustes(ashnc_f2_vectors, 
                                    ashnw_f2_vectors, 
                                    symmetric = T)
summary(ashn_f2_proc_res_eigen)
plot(ashn_f2_proc_res_eigen, 
     kind = 1)

ashn_f2_eigen_test = protest(ashnc_f2_vectors, 
                             ashnw_f2_vectors, 
                             permutations = 999, 
                             scale = F, 
                             symmetric = T)

ashn_f2_eigen_residuals = residuals(ashn_f2_eigen_test)

ashn_f2_eigen_residuals = ashn_f2_eigen_residuals%>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                eigen_residual = 2) 
# %>%  
#   write_csv('ASHN_F2_12_12_procrustes_eigen_cov_mat_residuals.csv')


bind_cols(ashn_f2_residuals, 
          ashn_f2_eigen_residuals)

ashn_f2_residuals_combo = inner_join(ashn_f2_residuals, 
           ashn_f2_eigen_residuals, 
           by = 'trait') %>% 
  mutate(pop_pair = 'ASHN')

# %>%
#   write_csv('ASHN_F2_procrustes_eigen_cov_mat_residuals.csv')


# MYV ---------------------------------------------------------------------

myvc_f2_eigen = eigen(myvc_f2_cov_mat)
myvw_f2_eigen = eigen(myvw_f2_cov_mat)

# Extract leading eigenvectors
myvc_f2_vectors = myvc_f2_eigen$vectors[, 1:17, drop = FALSE]
myvw_f2_vectors = myvw_f2_eigen$vectors[, 1:17, drop = FALSE]

# Preserve trait identities
rownames(myvc_f2_vectors) = colnames(myvc_f2_cov_mat)
rownames(myvw_f2_vectors) = colnames(myvw_f2_cov_mat)


myv_f2_proc_res_eigen = procrustes(myvc_f2_vectors, 
                                    myvw_f2_vectors, 
                                    symmetric = T)
summary(myv_f2_proc_res_eigen)
plot(myv_f2_proc_res_eigen, 
     kind = 1)

myv_f2_eigen_test = protest(myvc_f2_vectors, 
                             myvw_f2_vectors, 
                             permutations = 999, 
                             scale = F, 
                             symmetric = T)

myv_f2_eigen_residuals = residuals(myv_f2_eigen_test)

myv_f2_eigen_residuals = myv_f2_eigen_residuals%>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                eigen_residual = 2) 
# %>%  
#   write_csv('myv_F2_12_12_procrustes_eigen_cov_mat_residuals.csv')



myv_f2_residuals_combo = inner_join(myv_f2_residuals, 
                                     myv_f2_eigen_residuals, 
                                     by = 'trait') %>% 
  mutate(pop_pair = 'MYV')


# SKR ---------------------------------------------------------------------

skrc_f2_eigen = eigen(skrc_f2_cov_mat)
skrw_f2_eigen = eigen(skrw_f2_cov_mat)

# Extract leading eigenvectors
skrc_f2_vectors = skrc_f2_eigen$vectors[, 1:17, drop = FALSE]
skrw_f2_vectors = skrw_f2_eigen$vectors[, 1:17, drop = FALSE]

# Preserve trait identities
rownames(skrc_f2_vectors) = colnames(skrc_f2_cov_mat)
rownames(skrw_f2_vectors) = colnames(skrw_f2_cov_mat)


skr_f2_proc_res_eigen = procrustes(skrc_f2_vectors, 
                                   skrw_f2_vectors, 
                                   symmetric = T)
summary(skr_f2_proc_res_eigen)
plot(skr_f2_proc_res_eigen, 
     kind = 1)

skr_f2_eigen_test = protest(skrc_f2_vectors, 
                            skrw_f2_vectors, 
                            permutations = 999, 
                            scale = F, 
                            symmetric = T)

skr_f2_eigen_residuals = residuals(skr_f2_eigen_test)

skr_f2_eigen_residuals = skr_f2_eigen_residuals%>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                eigen_residual = 2) 
# %>%  
#   write_csv('skr_F2_12_12_procrustes_eigen_cov_mat_residuals.csv')



skr_f2_residuals_combo = inner_join(skr_f2_residuals, 
                                    skr_f2_eigen_residuals, 
                                    by = 'trait') %>% 
  mutate(pop_pair = 'SKR')


# GTS/GAR -----------------------------------------------------------------

gar_f2_eigen = eigen(gar_f2_cov_mat)
gts_f2_eigen = eigen(gts_f2_cov_mat)

# Extract leading eigenvectors
gar_f2_vectors = gar_f2_eigen$vectors[, 1:17, drop = FALSE]
gts_f2_vectors = gts_f2_eigen$vectors[, 1:17, drop = FALSE]

# Preserve trait identities
rownames(gar_f2_vectors) = colnames(gar_f2_cov_mat)
rownames(gts_f2_vectors) = colnames(gts_f2_cov_mat)


gts_gar_f2_proc_res_eigen = procrustes(gar_f2_vectors, 
                                   gts_f2_vectors, 
                                   symmetric = T)
summary(gts_gar_f2_proc_res_eigen)
plot(gts_gar_f2_proc_res_eigen, 
     kind = 1)

gts_gar_f2_eigen_test = protest(gar_f2_vectors, 
                            gts_f2_vectors, 
                            permutations = 999, 
                            scale = F, 
                            symmetric = T)

gts_gar_f2_eigen_residuals = residuals(gts_gar_f2_eigen_test)

gts_gar_f2_eigen_residuals = gts_gar_f2_eigen_residuals%>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                eigen_residual = 2) 
# %>%  
#   write_csv('gts_gar_F2_12_12_procrustes_eigen_cov_mat_residuals.csv')



gts_gar_f2_residuals_combo = inner_join(gts_gar_f2_residuals, 
                                    gts_gar_f2_eigen_residuals, 
                                    by = 'trait') %>% 
  mutate(pop_pair = 'GTS/GAR')


Proc_test_residual_combo = bind_rows(ashn_f2_residuals_combo,
                                     myv_f2_residuals_combo, 
                                     skr_f2_residuals_combo, 
                                     gts_gar_f2_residuals_combo)


# ProcTest residual plots -------------------------------------------------

cov_residual_plot = ggplot(Proc_test_residual_combo, 
       aes(x = factor(trait), 
           y = cov_residual))+
  geom_col()+
  ylim(0.0, 0.3)+
  facet_grid(~pop_pair)+
  labs(y = 'Covariance residuals')+
  # geom_hline(yintercept = 1.5, 
  #            linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold', 
                                  size = 14))

eigen_residual_plot = ggplot(Proc_test_residual_combo, 
                           aes(x = factor(trait), 
                               y = eigen_residual))+
  geom_col()+
  ylim(0.0, 0.3)+
  facet_grid(~pop_pair)+
  labs(y = 'Eigenvector residuals')+
  # geom_hline(yintercept = 1.5, 
  #            linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold', 
                                  size = 14))


# f2_orig data - KrzCov and procrustes ---------------------------------------
# ASHN --------------------------------------------------------------------

ashnc_f2_12_12_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'ASHNC', 
         Exp_treats == '12_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
ashnc_f2_12_18_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'ASHNC', 
         Exp_treats == '12_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
ashnc_f2_18_12_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'ASHNC', 
         Exp_treats == '18_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
ashnc_f2_18_18_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'ASHNC', 
         Exp_treats == '18_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)


ashnw_f2_12_12_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'ASHNW', 
         Exp_treats == '12_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
ashnw_f2_12_18_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'ASHNW', 
         Exp_treats == '12_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
ashnw_f2_18_12_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'ASHNW', 
         Exp_treats == '18_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
ashnw_f2_18_18_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'ASHNW', 
         Exp_treats == '18_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

ashnc_f2_12_12_cov_mat = cov(ashnc_f2_12_12_trait_mat, 
                             use = 'complete.obs')
ashnc_f2_12_18_cov_mat = cov(ashnc_f2_12_18_trait_mat, 
                             use = 'complete.obs')
ashnc_f2_18_12_cov_mat = cov(ashnc_f2_18_12_trait_mat, 
                             use = 'complete.obs')
ashnc_f2_18_18_cov_mat = cov(ashnc_f2_18_18_trait_mat, 
                             use = 'complete.obs')

ashnw_f2_12_12_cov_mat = cov(ashnw_f2_12_12_trait_mat, 
                             use = 'complete.obs')
ashnw_f2_12_18_cov_mat = cov(ashnw_f2_12_18_trait_mat, 
                             use = 'complete.obs')
ashnw_f2_18_12_cov_mat = cov(ashnw_f2_18_12_trait_mat, 
                             use = 'complete.obs')
ashnw_f2_18_18_cov_mat = cov(ashnw_f2_18_18_trait_mat, 
                             use = 'complete.obs')

ashn_f2_12_12_krzcov = evolqg::KrzCor(ashnc_f2_12_12_cov_mat, 
                                      ashnw_f2_12_12_cov_mat)

ashn_f2_12_18_krzcov = evolqg::KrzCor(ashnc_f2_12_18_cov_mat, 
                                      ashnw_f2_12_18_cov_mat)

ashn_f2_18_12_krzcov = evolqg::KrzCor(ashnc_f2_18_12_cov_mat, 
                                      ashnw_f2_18_12_cov_mat)

ashn_f2_18_18_krzcov = evolqg::KrzCor(ashnc_f2_18_18_cov_mat, 
                                      ashnw_f2_18_18_cov_mat)


ashn_f2_12_12_proc_res = procrustes(ashnc_f2_12_12_cov_mat, 
                                    ashnw_f2_12_12_cov_mat)
ashn_f2_12_18_proc_res = procrustes(ashnc_f2_12_18_cov_mat, 
                                    ashnw_f2_12_18_cov_mat)
ashn_f2_18_12_proc_res = procrustes(ashnc_f2_18_12_cov_mat, 
                                    ashnw_f2_18_12_cov_mat)
ashn_f2_18_18_proc_res = procrustes(ashnc_f2_18_18_cov_mat, 
                                    ashnw_f2_18_18_cov_mat)

summary(ashn_f2_12_12_proc_res)
summary(ashn_f2_12_18_proc_res)
summary(ashn_f2_18_12_proc_res)
summary(ashn_f2_18_18_proc_res)

ashn_f2_12_12_proc_test = protest(ashnc_f2_12_12_cov_mat, 
                                  ashnw_f2_12_12_cov_mat, 
                                  permutations = 999)
ashn_f2_12_18_proc_test = protest(ashnc_f2_12_18_cov_mat, 
                                  ashnw_f2_12_18_cov_mat, 
                                  permutations = 999)
ashn_f2_18_12_proc_test = protest(ashnc_f2_18_12_cov_mat, 
                                  ashnw_f2_18_12_cov_mat, 
                                  permutations = 999)
ashn_f2_18_18_proc_test = protest(ashnc_f2_18_18_cov_mat, 
                                  ashnw_f2_18_18_cov_mat, 
                                  permutations = 999)


ashn_f2_12_12_residuals = residuals(ashn_f2_12_12_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

ashn_f2_12_18_residuals = residuals(ashn_f2_12_18_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

ashn_f2_18_12_residuals = residuals(ashn_f2_18_12_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

ashn_f2_18_18_residuals = residuals(ashn_f2_18_18_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)


ashn_f2_12_12_resid_plot = ggplot(ashn_f2_12_12_residuals, 
       aes(x = factor(trait), 
           y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

ashn_f2_12_18_resid_plot = ggplot(ashn_f2_12_18_residuals, 
                                  aes(x = factor(trait), 
                                      y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

ashn_f2_18_12_resid_plot = ggplot(ashn_f2_18_12_residuals, 
                                  aes(x = factor(trait), 
                                      y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

ashn_f2_18_18_resid_plot = ggplot(ashn_f2_18_18_residuals, 
                                  aes(x = factor(trait), 
                                      y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

ashn_exp_treat_combo = (ashn_f2_12_12_resid_plot+ashn_f2_12_18_resid_plot)/(ashn_f2_18_12_resid_plot+ashn_f2_18_18_resid_plot)


# MYV ---------------------------------------------------------------------
myvc_f2_12_12_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'MYVC', 
         Exp_treats == '12_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
myvc_f2_12_18_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'MYVC', 
         Exp_treats == '12_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
myvc_f2_18_12_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'MYVC', 
         Exp_treats == '18_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
myvc_f2_18_18_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'MYVC', 
         Exp_treats == '18_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)


myvw_f2_12_12_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'MYVW', 
         Exp_treats == '12_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
myvw_f2_12_18_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'MYVW', 
         Exp_treats == '12_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
myvw_f2_18_12_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'MYVW', 
         Exp_treats == '18_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
myvw_f2_18_18_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'MYVW', 
         Exp_treats == '18_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

myvc_f2_12_12_cov_mat = cov(myvc_f2_12_12_trait_mat, 
                             use = 'complete.obs')
myvc_f2_12_18_cov_mat = cov(myvc_f2_12_18_trait_mat, 
                             use = 'complete.obs')
myvc_f2_18_12_cov_mat = cov(myvc_f2_18_12_trait_mat, 
                             use = 'complete.obs')
myvc_f2_18_18_cov_mat = cov(myvc_f2_18_18_trait_mat, 
                             use = 'complete.obs')

myvw_f2_12_12_cov_mat = cov(myvw_f2_12_12_trait_mat, 
                             use = 'complete.obs')
myvw_f2_12_18_cov_mat = cov(myvw_f2_12_18_trait_mat, 
                             use = 'complete.obs')
myvw_f2_18_12_cov_mat = cov(myvw_f2_18_12_trait_mat, 
                             use = 'complete.obs')
myvw_f2_18_18_cov_mat = cov(myvw_f2_18_18_trait_mat, 
                             use = 'complete.obs')

myv_f2_12_12_krzcov = evolqg::KrzCor(myvc_f2_12_12_cov_mat, 
                                      myvw_f2_12_12_cov_mat)

myv_f2_12_18_krzcov = evolqg::KrzCor(myvc_f2_12_18_cov_mat, 
                                      myvw_f2_12_18_cov_mat)

myv_f2_18_12_krzcov = evolqg::KrzCor(myvc_f2_18_12_cov_mat, 
                                      myvw_f2_18_12_cov_mat)

myv_f2_18_18_krzcov = evolqg::KrzCor(myvc_f2_18_18_cov_mat, 
                                      myvw_f2_18_18_cov_mat)


myv_f2_12_12_proc_res = procrustes(myvc_f2_12_12_cov_mat, 
                                    myvw_f2_12_12_cov_mat)
myv_f2_12_18_proc_res = procrustes(myvc_f2_12_18_cov_mat, 
                                    myvw_f2_12_18_cov_mat)
myv_f2_18_12_proc_res = procrustes(myvc_f2_18_12_cov_mat, 
                                    myvw_f2_18_12_cov_mat)
myv_f2_18_18_proc_res = procrustes(myvc_f2_18_18_cov_mat, 
                                    myvw_f2_18_18_cov_mat)

summary(myv_f2_12_12_proc_res)
summary(myv_f2_12_18_proc_res)
summary(myv_f2_18_12_proc_res)
summary(myv_f2_18_18_proc_res)

myv_f2_12_12_proc_test = protest(myvc_f2_12_12_cov_mat, 
                                  myvw_f2_12_12_cov_mat, 
                                  permutations = 999)
myv_f2_12_18_proc_test = protest(myvc_f2_12_18_cov_mat, 
                                  myvw_f2_12_18_cov_mat, 
                                  permutations = 999)
myv_f2_18_12_proc_test = protest(myvc_f2_18_12_cov_mat, 
                                  myvw_f2_18_12_cov_mat, 
                                  permutations = 999)
myv_f2_18_18_proc_test = protest(myvc_f2_18_18_cov_mat, 
                                  myvw_f2_18_18_cov_mat, 
                                  permutations = 999)


myv_f2_12_12_residuals = residuals(myv_f2_12_12_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

myv_f2_12_18_residuals = residuals(myv_f2_12_18_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

myv_f2_18_12_residuals = residuals(myv_f2_18_12_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

myv_f2_18_18_residuals = residuals(myv_f2_18_18_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)


myv_f2_12_12_resid_plot = ggplot(myv_f2_12_12_residuals, 
                                  aes(x = factor(trait), 
                                      y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

myv_f2_12_18_resid_plot = ggplot(myv_f2_12_18_residuals, 
                                  aes(x = factor(trait), 
                                      y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

myv_f2_18_12_resid_plot = ggplot(myv_f2_18_12_residuals, 
                                  aes(x = factor(trait), 
                                      y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

myv_f2_18_18_resid_plot = ggplot(myv_f2_18_18_residuals, 
                                  aes(x = factor(trait), 
                                      y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

myv_exp_treat_combo = (myv_f2_12_12_resid_plot+myv_f2_12_18_resid_plot)/(myv_f2_18_12_resid_plot+myv_f2_18_18_resid_plot)

# SKR ---------------------------------------------------------------------
skrc_f2_12_12_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'SKRC', 
         Exp_treats == '12_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
skrc_f2_12_18_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'SKRC', 
         Exp_treats == '12_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
skrc_f2_18_12_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'SKRC', 
         Exp_treats == '18_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
skrc_f2_18_18_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'SKRC', 
         Exp_treats == '18_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)


skrw_f2_12_12_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'SKRW', 
         Exp_treats == '12_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
skrw_f2_12_18_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'SKRW', 
         Exp_treats == '12_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
skrw_f2_18_12_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'SKRW', 
         Exp_treats == '18_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
skrw_f2_18_18_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'SKRW', 
         Exp_treats == '18_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

skrc_f2_12_12_cov_mat = cov(skrc_f2_12_12_trait_mat, 
                            use = 'complete.obs')
skrc_f2_12_18_cov_mat = cov(skrc_f2_12_18_trait_mat, 
                            use = 'complete.obs')
skrc_f2_18_12_cov_mat = cov(skrc_f2_18_12_trait_mat, 
                            use = 'complete.obs')
skrc_f2_18_18_cov_mat = cov(skrc_f2_18_18_trait_mat, 
                            use = 'complete.obs')

skrw_f2_12_12_cov_mat = cov(skrw_f2_12_12_trait_mat, 
                            use = 'complete.obs')
skrw_f2_12_18_cov_mat = cov(skrw_f2_12_18_trait_mat, 
                            use = 'complete.obs')
skrw_f2_18_12_cov_mat = cov(skrw_f2_18_12_trait_mat, 
                            use = 'complete.obs')
skrw_f2_18_18_cov_mat = cov(skrw_f2_18_18_trait_mat, 
                            use = 'complete.obs')

skr_f2_12_12_krzcov = evolqg::KrzCor(skrc_f2_12_12_cov_mat, 
                                     skrw_f2_12_12_cov_mat)

skr_f2_12_18_krzcov = evolqg::KrzCor(skrc_f2_12_18_cov_mat, 
                                     skrw_f2_12_18_cov_mat)

skr_f2_18_12_krzcov = evolqg::KrzCor(skrc_f2_18_12_cov_mat, 
                                     skrw_f2_18_12_cov_mat)

skr_f2_18_18_krzcov = evolqg::KrzCor(skrc_f2_18_18_cov_mat, 
                                     skrw_f2_18_18_cov_mat)


skr_f2_12_12_proc_res = procrustes(skrc_f2_12_12_cov_mat, 
                                   skrw_f2_12_12_cov_mat)
skr_f2_12_18_proc_res = procrustes(skrc_f2_12_18_cov_mat, 
                                   skrw_f2_12_18_cov_mat)
skr_f2_18_12_proc_res = procrustes(skrc_f2_18_12_cov_mat, 
                                   skrw_f2_18_12_cov_mat)
skr_f2_18_18_proc_res = procrustes(skrc_f2_18_18_cov_mat, 
                                   skrw_f2_18_18_cov_mat)

summary(skr_f2_12_12_proc_res)
summary(skr_f2_12_18_proc_res)
summary(skr_f2_18_12_proc_res)
summary(skr_f2_18_18_proc_res)

skr_f2_12_12_proc_test = protest(skrc_f2_12_12_cov_mat, 
                                 skrw_f2_12_12_cov_mat, 
                                 permutations = 999)
skr_f2_12_18_proc_test = protest(skrc_f2_12_18_cov_mat, 
                                 skrw_f2_12_18_cov_mat, 
                                 permutations = 999)
skr_f2_18_12_proc_test = protest(skrc_f2_18_12_cov_mat, 
                                 skrw_f2_18_12_cov_mat, 
                                 permutations = 999)
skr_f2_18_18_proc_test = protest(skrc_f2_18_18_cov_mat, 
                                 skrw_f2_18_18_cov_mat, 
                                 permutations = 999)


skr_f2_12_12_residuals = residuals(skr_f2_12_12_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

skr_f2_12_18_residuals = residuals(skr_f2_12_18_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

skr_f2_18_12_residuals = residuals(skr_f2_18_12_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

skr_f2_18_18_residuals = residuals(skr_f2_18_18_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)


skr_f2_12_12_resid_plot = ggplot(skr_f2_12_12_residuals, 
                                 aes(x = factor(trait), 
                                     y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

skr_f2_12_18_resid_plot = ggplot(skr_f2_12_18_residuals, 
                                 aes(x = factor(trait), 
                                     y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

skr_f2_18_12_resid_plot = ggplot(skr_f2_18_12_residuals, 
                                 aes(x = factor(trait), 
                                     y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

skr_f2_18_18_resid_plot = ggplot(skr_f2_18_18_residuals, 
                                 aes(x = factor(trait), 
                                     y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

skr_exp_treat_combo = (skr_f2_12_12_resid_plot+skr_f2_12_18_resid_plot)/(skr_f2_18_12_resid_plot+skr_f2_18_18_resid_plot)



# GAR/GTS -----------------------------------------------------------------
CSWYC_f2_12_12_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'CSWYC', 
         Exp_treats == '12_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
CSWYC_f2_12_18_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'CSWYC', 
         Exp_treats == '12_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
CSWYC_f2_18_12_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'CSWYC', 
         Exp_treats == '18_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
CSWYC_f2_18_18_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'CSWYC', 
         Exp_treats == '18_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)


GTSW_f2_12_12_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'GTSW', 
         Exp_treats == '12_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
GTSW_f2_12_18_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'GTSW', 
         Exp_treats == '12_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
GTSW_f2_18_12_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'GTSW', 
         Exp_treats == '18_12') %>% 
  dplyr::select(jaw_length:caudal2_15_17)
GTSW_f2_18_18_trait_mat = f2_orig_data %>% 
  filter(Lake_morph == 'GTSW', 
         Exp_treats == '18_18') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

CSWYC_f2_12_12_cov_mat = cov(CSWYC_f2_12_12_trait_mat, 
                            use = 'complete.obs')
CSWYC_f2_12_18_cov_mat = cov(CSWYC_f2_12_18_trait_mat, 
                            use = 'complete.obs')
CSWYC_f2_18_12_cov_mat = cov(CSWYC_f2_18_12_trait_mat, 
                            use = 'complete.obs')
CSWYC_f2_18_18_cov_mat = cov(CSWYC_f2_18_18_trait_mat, 
                            use = 'complete.obs')

GTSW_f2_12_12_cov_mat = cov(GTSW_f2_12_12_trait_mat, 
                            use = 'complete.obs')
GTSW_f2_12_18_cov_mat = cov(GTSW_f2_12_18_trait_mat, 
                            use = 'complete.obs')
GTSW_f2_18_12_cov_mat = cov(GTSW_f2_18_12_trait_mat, 
                            use = 'complete.obs')
GTSW_f2_18_18_cov_mat = cov(GTSW_f2_18_18_trait_mat, 
                            use = 'complete.obs')

gts_gar_f2_12_12_krzcov = evolqg::KrzCor(CSWYC_f2_12_12_cov_mat, 
                                     GTSW_f2_12_12_cov_mat)

gts_gar_f2_12_18_krzcov = evolqg::KrzCor(CSWYC_f2_12_18_cov_mat, 
                                     GTSW_f2_12_18_cov_mat)

gts_gar_f2_18_12_krzcov = evolqg::KrzCor(CSWYC_f2_18_12_cov_mat, 
                                     GTSW_f2_18_12_cov_mat)

gts_gar_f2_18_18_krzcov = evolqg::KrzCor(CSWYC_f2_18_18_cov_mat, 
                                     GTSW_f2_18_18_cov_mat)


gts_gar_f2_12_12_proc_res = procrustes(CSWYC_f2_12_12_cov_mat, 
                                   GTSW_f2_12_12_cov_mat)
gts_gar_f2_12_18_proc_res = procrustes(CSWYC_f2_12_18_cov_mat, 
                                   GTSW_f2_12_18_cov_mat)
gts_gar_f2_18_12_proc_res = procrustes(CSWYC_f2_18_12_cov_mat, 
                                   GTSW_f2_18_12_cov_mat)
gts_gar_f2_18_18_proc_res = procrustes(CSWYC_f2_18_18_cov_mat, 
                                   GTSW_f2_18_18_cov_mat)

summary(gts_gar_f2_12_12_proc_res)
summary(gts_gar_f2_12_18_proc_res)
summary(gts_gar_f2_18_12_proc_res)
summary(gts_gar_f2_18_18_proc_res)

gts_gar_f2_12_12_proc_test = protest(CSWYC_f2_12_12_cov_mat, 
                                 GTSW_f2_12_12_cov_mat, 
                                 permutations = 999)
gts_gar_f2_12_18_proc_test = protest(CSWYC_f2_12_18_cov_mat, 
                                 GTSW_f2_12_18_cov_mat, 
                                 permutations = 999)
gts_gar_f2_18_12_proc_test = protest(CSWYC_f2_18_12_cov_mat, 
                                 GTSW_f2_18_12_cov_mat, 
                                 permutations = 999)
gts_gar_f2_18_18_proc_test = protest(CSWYC_f2_18_18_cov_mat, 
                                 GTSW_f2_18_18_cov_mat, 
                                 permutations = 999)


gts_gar_f2_12_12_residuals = residuals(gts_gar_f2_12_12_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

gts_gar_f2_12_18_residuals = residuals(gts_gar_f2_12_18_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

gts_gar_f2_18_12_residuals = residuals(gts_gar_f2_18_12_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)

gts_gar_f2_18_18_residuals = residuals(gts_gar_f2_18_18_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2)%>% 
  mutate(resid_10 = residual*10)


gts_gar_f2_12_12_resid_plot = ggplot(gts_gar_f2_12_12_residuals, 
                                 aes(x = factor(trait), 
                                     y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

gts_gar_f2_12_18_resid_plot = ggplot(gts_gar_f2_12_18_residuals, 
                                 aes(x = factor(trait), 
                                     y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

gts_gar_f2_18_12_resid_plot = ggplot(gts_gar_f2_18_12_residuals, 
                                 aes(x = factor(trait), 
                                     y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

gts_gar_f2_18_18_resid_plot = ggplot(gts_gar_f2_18_18_residuals, 
                                 aes(x = factor(trait), 
                                     y = resid_10))+
  geom_col()+
  ylim(0.0, 3.0)+
  geom_hline(yintercept = 1.5, 
             linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

gts_gar_exp_treat_combo = (gts_gar_f2_12_12_resid_plot+gts_gar_f2_12_18_resid_plot)/(gts_gar_f2_18_12_resid_plot+gts_gar_f2_18_18_resid_plot)



# Parallelism in plasticity divergence in trait covariance ----------------

## kendals test of concordance 
##Rows are objects and columns are judges (variables)
## Traits in rows and pop pairs in columns

ashn_f2_12_12_kendal = ashn_f2_12_12_residuals %>% 
  dplyr::rename(ASHN = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

myv_f2_12_12_kendal = myv_f2_12_12_residuals%>% 
  dplyr::rename(MYV = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

skr_f2_12_12_kendal = skr_f2_12_12_residuals%>% 
  dplyr::rename(SKR = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

gts_gar_f2_12_12_kendal = gts_gar_f2_12_12_residuals%>% 
  dplyr::rename(GTS_GAR = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

f2_12_12_kendal_df = bind_cols(ashn_f2_12_12_kendal, 
                            myv_f2_12_12_kendal, 
                            skr_f2_12_12_kendal, 
                            gts_gar_f2_12_12_kendal)

group = c('ASHN','MYV','SKR','GTS_GAR')


f2_12_12_kendal = kendall.global(f2_12_12_kendal_df,
                                 # group = pops,
                                 nperm = 999, 
                                 mult = 'hochberg')
f2_12_12_kendal_post = kendall.post(f2_12_12_kendal_df, 
                                   mult = 'hochberg', 
                                   nperm = 999)


ashn_f2_12_18_kendal = ashn_f2_12_18_residuals %>% 
  dplyr::rename(ASHN = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

myv_f2_12_18_kendal = myv_f2_12_18_residuals%>% 
  dplyr::rename(MYV = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

skr_f2_12_18_kendal = skr_f2_12_18_residuals%>% 
  dplyr::rename(SKR = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

gts_gar_f2_12_18_kendal = gts_gar_f2_12_18_residuals%>% 
  dplyr::rename(GTS_GAR = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

f2_12_18_kendal_df = bind_cols(ashn_f2_12_18_kendal, 
                               myv_f2_12_18_kendal, 
                               skr_f2_12_18_kendal, 
                               gts_gar_f2_12_18_kendal)


f2_12_18_kendal = kendall.global(f2_12_18_kendal_df,
                                 # group = pops,
                                 nperm = 999, 
                                 mult = 'hochberg')
f2_12_18_kendal_post = kendall.post(f2_12_18_kendal_df, 
                                    mult = 'hochberg', 
                                    nperm = 999)



ashn_f2_18_12_kendal = ashn_f2_18_12_residuals %>% 
  dplyr::rename(ASHN = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

myv_f2_18_12_kendal = myv_f2_18_12_residuals%>% 
  dplyr::rename(MYV = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

skr_f2_18_12_kendal = skr_f2_18_12_residuals%>% 
  dplyr::rename(SKR = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

gts_gar_f2_18_12_kendal = gts_gar_f2_18_12_residuals%>% 
  dplyr::rename(GTS_GAR = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

f2_18_12_kendal_df = bind_cols(ashn_f2_18_12_kendal, 
                               myv_f2_18_12_kendal, 
                               skr_f2_18_12_kendal, 
                               gts_gar_f2_18_12_kendal)


f2_18_12_kendal = kendall.global(f2_18_12_kendal_df,
                                 # group = pops,
                                 nperm = 999, 
                                 mult = 'hochberg')
f2_18_12_kendal_post = kendall.post(f2_18_12_kendal_df, 
                                    mult = 'hochberg', 
                                    nperm = 999)




ashn_f2_18_18_kendal = ashn_f2_18_18_residuals %>% 
  dplyr::rename(ASHN = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

myv_f2_18_18_kendal = myv_f2_18_18_residuals%>% 
  dplyr::rename(MYV = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

skr_f2_18_18_kendal = skr_f2_18_18_residuals%>% 
  dplyr::rename(SKR = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

gts_gar_f2_18_18_kendal = gts_gar_f2_18_18_residuals%>% 
  dplyr::rename(GTS_GAR = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

f2_18_18_kendal_df = bind_cols(ashn_f2_18_18_kendal, 
                               myv_f2_18_18_kendal, 
                               skr_f2_18_18_kendal, 
                               gts_gar_f2_18_18_kendal)

f2_18_18_kendal = kendall.global(f2_18_18_kendal_df,
                                 # group = pops,
                                 nperm = 999, 
                                 mult = 'hochberg')
f2_18_18_kendal_post = kendall.post(f2_18_18_kendal_df, 
                                    mult = 'hochberg', 
                                    nperm = 999)


f2_12_12_kendal_post
f2_12_18_kendal_post
f2_18_12_kendal_post
f2_18_18_kendal_post




# Nemeny post hoc test ----------------------------------------------------

library(PMCMRplus)


ashn_f2_12_12_nemeny_df = ashn_f2_12_12_residuals %>% 
  mutate(pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10)

myv_f2_12_12_nemeny_df = myv_f2_12_12_residuals %>% 
  mutate(pop_pair = 'MYV') %>% 
  dplyr::select(-resid_10)

skr_f2_12_12_nemeny_df = skr_f2_12_12_residuals %>% 
  mutate(pop_pair = 'SKR') %>% 
  dplyr::select(-resid_10)

gts_gar_f2_12_12_nemeny_df = gts_gar_f2_12_12_residuals %>% 
  mutate(pop_pair = 'GTS/GAR') %>% 
  dplyr::select(-resid_10)

f2_12_12_nemeny_df = bind_rows(ashn_f2_12_12_nemeny_df, 
                      myv_f2_12_12_nemeny_df, 
                      skr_f2_12_12_nemeny_df, 
                      gts_gar_f2_12_12_nemeny_df)


f2_12_12_nemeny_test = frdAllPairsNemenyiTest(y = nemeny_df$residual, 
                       groups = nemeny_df$trait, 
                       blocks = nemeny_df$pop_pair, 
                       p.adjust.method = 'hochberg')

f2_12_12_nemeny_test$p.value %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname) %>% 
  # filter(across(c(2:33), function(x) x<0.05))
write_csv('F2_1212_residual_NemenyTest_results.csv')



ashn_f2_12_18_nemeny_df = ashn_f2_12_18_residuals %>% 
  mutate(pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10)

myv_f2_12_18_nemeny_df = myv_f2_12_18_residuals %>% 
  mutate(pop_pair = 'MYV') %>% 
  dplyr::select(-resid_10)

skr_f2_12_18_nemeny_df = skr_f2_12_18_residuals %>% 
  mutate(pop_pair = 'SKR') %>% 
  dplyr::select(-resid_10)

gts_gar_f2_12_18_nemeny_df = gts_gar_f2_12_18_residuals %>% 
  mutate(pop_pair = 'GTS/GAR') %>% 
  dplyr::select(-resid_10)

f2_12_18_nemeny_df = bind_rows(ashn_f2_12_18_nemeny_df, 
                               myv_f2_12_18_nemeny_df, 
                               skr_f2_12_18_nemeny_df, 
                               gts_gar_f2_12_18_nemeny_df)


f2_12_18_nemeny_test = frdAllPairsNemenyiTest(y = f2_12_18_nemeny_df$residual, 
                                              groups = f2_12_18_nemeny_df$trait, 
                                              blocks = f2_12_18_nemeny_df$pop_pair, 
                                              p.adjust.method = 'hochberg')

f2_12_18_nemeny_test$p.value %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname) %>% 
  filter(across(c(2:33), function(x) x<0.05))
  # write_csv('F2_1218_residual_NemenyTest_results.csv')

ashn_f2_18_12_nemeny_df = ashn_f2_18_12_residuals %>% 
  mutate(pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10)

myv_f2_18_12_nemeny_df = myv_f2_18_12_residuals %>% 
  mutate(pop_pair = 'MYV') %>% 
  dplyr::select(-resid_10)

skr_f2_18_12_nemeny_df = skr_f2_18_12_residuals %>% 
  mutate(pop_pair = 'SKR') %>% 
  dplyr::select(-resid_10)

gts_gar_f2_18_12_nemeny_df = gts_gar_f2_18_12_residuals %>% 
  mutate(pop_pair = 'GTS/GAR') %>% 
  dplyr::select(-resid_10)

f2_18_12_nemeny_df = bind_rows(ashn_f2_18_12_nemeny_df, 
                               myv_f2_18_12_einemeny_df, 
                               skr_f2_18_12_nemeny_df, 
                               gts_gar_f2_18_12_nemeny_df)


f2_18_12_nemeny_test = frdAllPairsNemenyiTest(y = f2_18_12_nemeny_df$residual, 
                                              groups = f2_18_12_nemeny_df$trait, 
                                              blocks = f2_18_12_nemeny_df$pop_pair, 
                                              p.adjust.method = 'hochberg')


f2_18_12_nemeny_test$p.value %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname) %>% 
  filter(across(c(2:33), function(x) x<0.05))
  # write_csv('F2_1812_residual_NemenyTest_results.csv')



ashn_f2_18_18_nemeny_df = ashn_f2_18_18_residuals %>% 
  mutate(pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10)

myv_f2_18_18_nemeny_df = myv_f2_18_18_residuals %>% 
  mutate(pop_pair = 'MYV') %>% 
  dplyr::select(-resid_10)

skr_f2_18_18_nemeny_df = skr_f2_18_18_residuals %>% 
  mutate(pop_pair = 'SKR') %>% 
  dplyr::select(-resid_10)

gts_gar_f2_18_18_nemeny_df = gts_gar_f2_18_18_residuals %>% 
  mutate(pop_pair = 'GTS/GAR') %>% 
  dplyr::select(-resid_10)

f2_18_18_nemeny_df = bind_rows(ashn_f2_18_18_nemeny_df, 
                               myv_f2_18_18_nemeny_df, 
                               skr_f2_18_18_nemeny_df, 
                               gts_gar_f2_18_18_nemeny_df)


f2_18_18_nemeny_test = frdAllPairsNemenyiTest(y = f2_18_18_nemeny_df$residual, 
                                              groups = f2_18_18_nemeny_df$trait, 
                                              blocks = f2_18_18_nemeny_df$pop_pair, 
                                              p.adjust.method = 'hochberg')


f2_18_18_nemeny_test$p.value %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname) %>% 
  filter(across(c(2:33), function(x) x<0.05))
  # %>% write_csv('F2_1818_residual_NemenyTest_results.csv')

# parallelism test - eye ball method --------------------------------------


(ashn_f2_12_12_resid_plot+myv_f2_12_12_resid_plot)/(skr_f2_12_12_resid_plot+gts_gar_f2_12_12_resid_plot)

(ashn_f2_12_18_resid_plot+myv_f2_12_18_resid_plot)/(skr_f2_12_18_resid_plot+gts_gar_f2_12_18_resid_plot)
## Opercular KT divergence plastically divergent between geo and amb cov mats 
## across all pops
(ashn_f2_18_12_resid_plot+myv_f2_18_12_resid_plot)/(skr_f2_18_12_resid_plot+gts_gar_f2_18_12_resid_plot)

(ashn_f2_18_18_resid_plot+myv_f2_18_18_resid_plot)/(skr_f2_18_18_resid_plot+gts_gar_f2_18_18_resid_plot)


# Combine trait covariance matrices ---------------------------------------


# ASHN covariance mats ----------------------------------------------------


ashn_cov_matrices <- list(
  Temp1_Ambient = ashnc_f2_12_12_cov_mat,
  Temp1_Geothermal = ashnw_f2_12_12_cov_mat,
  Temp2_Ambient = ashnc_f2_12_18_cov_mat,
  Temp2_Geothermal = ashnw_f2_12_18_cov_mat,
  Temp3_Ambient = ashnc_f2_18_12_cov_mat,
  Temp3_Geothermal = ashnw_f2_18_12_cov_mat,
  Temp4_Ambient = ashnc_f2_18_18_cov_mat,
  Temp4_Geothermal = ashnw_f2_18_18_cov_mat)

for (i in 1:4) {
  ashn_cov_matrices[[paste0("Temp", i, "_Difference")]] <-
    ashn_cov_matrices[[paste0("Temp", i, "_Geothermal")]] -
    ashn_cov_matrices[[paste0("Temp", i, "_Ambient")]]
}

ashn_cov_long <- purrr::imap_dfr(ashn_cov_matrices, function(mat, name) {
  as.data.frame(as.table(mat)) %>%
    setNames(c("Trait1", "Trait2", "Covariance")) %>%
    mutate(
      Temp = sub("_.*", "", name),
      Type = sub(".*_", "", name)
    )
})

ashn_cov_long$Type <- factor(
  ashn_cov_long$Type,
  levels = c("Ambient", "Geothermal", "Difference")
)

ashn_cov_long = ashn_cov_long %>% 
  mutate(Temp_info = case_when(
    str_detect(Temp, 'Temp1') ~ '12_12',
    str_detect(Temp, 'Temp2') ~ '12_18',
    str_detect(Temp, 'Temp3') ~ '18_12',
    str_detect(Temp, 'Temp4') ~ '18_18'))

ashn_cov_plot = ggplot(ashn_cov_long, 
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
       title = 'ASHN F2 trait covariance')+
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1),
        panel.grid = element_blank(), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'), 
        axis.text = element_text(size = 6), 
        legend.position = 'none')

# MYV covariance mats ----------------------------------------------------


myv_cov_matrices <- list(
  Temp1_Ambient = myvc_f2_12_12_cov_mat,
  Temp1_Geothermal = myvw_f2_12_12_cov_mat,
  Temp2_Ambient = myvc_f2_12_18_cov_mat,
  Temp2_Geothermal = myvw_f2_12_18_cov_mat,
  Temp3_Ambient = myvc_f2_18_12_cov_mat,
  Temp3_Geothermal = myvw_f2_18_12_cov_mat,
  Temp4_Ambient = myvc_f2_18_18_cov_mat,
  Temp4_Geothermal = myvw_f2_18_18_cov_mat)

for (i in 1:4) {
  myv_cov_matrices[[paste0("Temp", i, "_Difference")]] <-
    myv_cov_matrices[[paste0("Temp", i, "_Geothermal")]] -
    myv_cov_matrices[[paste0("Temp", i, "_Ambient")]]
}

myv_cov_long <- purrr::imap_dfr(myv_cov_matrices, function(mat, name) {
  as.data.frame(as.table(mat)) %>%
    setNames(c("Trait1", "Trait2", "Covariance")) %>%
    mutate(
      Temp = sub("_.*", "", name),
      Type = sub(".*_", "", name)
    )
})

myv_cov_long$Type <- factor(
  myv_cov_long$Type,
  levels = c("Ambient", "Geothermal", "Difference")
)

myv_cov_long = myv_cov_long %>% 
  mutate(Temp_info = case_when(
    str_detect(Temp, 'Temp1') ~ '12_12',
    str_detect(Temp, 'Temp2') ~ '12_18',
    str_detect(Temp, 'Temp3') ~ '18_12',
    str_detect(Temp, 'Temp4') ~ '18_18'))

myv_cov_plot = ggplot(myv_cov_long, 
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
       title = 'MYV F2 trait covariance')+
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1),
        panel.grid = element_blank(), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'), 
        axis.text = element_text(size = 6), 
        legend.position = 'none')

# skr covariance mats ----------------------------------------------------


skr_cov_matrices <- list(
  Temp1_Ambient = skrc_f2_12_12_cov_mat,
  Temp1_Geothermal = skrw_f2_12_12_cov_mat,
  Temp2_Ambient = skrc_f2_12_18_cov_mat,
  Temp2_Geothermal = skrw_f2_12_18_cov_mat,
  Temp3_Ambient = skrc_f2_18_12_cov_mat,
  Temp3_Geothermal = skrw_f2_18_12_cov_mat,
  Temp4_Ambient = skrc_f2_18_18_cov_mat,
  Temp4_Geothermal = skrw_f2_18_18_cov_mat)

for (i in 1:4) {
  skr_cov_matrices[[paste0("Temp", i, "_Difference")]] <-
    skr_cov_matrices[[paste0("Temp", i, "_Geothermal")]] -
    skr_cov_matrices[[paste0("Temp", i, "_Ambient")]]
}

skr_cov_long <- purrr::imap_dfr(skr_cov_matrices, function(mat, name) {
  as.data.frame(as.table(mat)) %>%
    setNames(c("Trait1", "Trait2", "Covariance")) %>%
    mutate(
      Temp = sub("_.*", "", name),
      Type = sub(".*_", "", name)
    )
})

skr_cov_long$Type <- factor(
  skr_cov_long$Type,
  levels = c("Ambient", "Geothermal", "Difference")
)

skr_cov_long = skr_cov_long %>% 
  mutate(Temp_info = case_when(
    str_detect(Temp, 'Temp1') ~ '12_12',
    str_detect(Temp, 'Temp2') ~ '12_18',
    str_detect(Temp, 'Temp3') ~ '18_12',
    str_detect(Temp, 'Temp4') ~ '18_18'))

skr_cov_plot = ggplot(skr_cov_long, 
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
       title = 'SKR F2 trait covariance')+
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1),
        panel.grid = element_blank(), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'), 
        axis.text = element_text(size = 6), 
        legend.position = 'none')

# gts_gar covariance mats ----------------------------------------------------


gts_gar_cov_matrices <- list(
  Temp1_Ambient = CSWYC_f2_12_12_cov_mat,
  Temp1_Geothermal = GTSW_f2_12_12_cov_mat,
  Temp2_Ambient = CSWYC_f2_12_18_cov_mat,
  Temp2_Geothermal = GTSW_f2_12_18_cov_mat,
  Temp3_Ambient = CSWYC_f2_18_12_cov_mat,
  Temp3_Geothermal = GTSW_f2_18_12_cov_mat,
  Temp4_Ambient = CSWYC_f2_18_18_cov_mat,
  Temp4_Geothermal = GTSW_f2_18_18_cov_mat)

for (i in 1:4) {
  gts_gar_cov_matrices[[paste0("Temp", i, "_Difference")]] <-
    gts_gar_cov_matrices[[paste0("Temp", i, "_Geothermal")]] -
    gts_gar_cov_matrices[[paste0("Temp", i, "_Ambient")]]
}

gts_gar_cov_long <- purrr::imap_dfr(gts_gar_cov_matrices, function(mat, name) {
  as.data.frame(as.table(mat)) %>%
    setNames(c("Trait1", "Trait2", "Covariance")) %>%
    mutate(
      Temp = sub("_.*", "", name),
      Type = sub(".*_", "", name)
    )
})

gts_gar_cov_long$Type <- factor(
  gts_gar_cov_long$Type,
  levels = c("Ambient", "Geothermal", "Difference")
)

gts_gar_cov_long = gts_gar_cov_long %>% 
  mutate(Temp_info = case_when(
    str_detect(Temp, 'Temp1') ~ '12_12',
    str_detect(Temp, 'Temp2') ~ '12_18',
    str_detect(Temp, 'Temp3') ~ '18_12',
    str_detect(Temp, 'Temp4') ~ '18_18'))

gts_gar_cov_plot = ggplot(gts_gar_cov_long, 
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
       title = 'GTS/GAR F2 trait covariance')+
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1),
        panel.grid = element_blank(), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'), 
        axis.text = element_text(size = 6), 
        legend.position = 'none')



(ashn_cov_plot + myv_cov_plot)/(skr_cov_plot + gts_gar_cov_plot)

(ashn_cov_plot + myv_cov_plot + skr_cov_plot + gts_gar_cov_plot)

# KrzCov null distribution resampling -------------------------------------


# ASHN --------------------------------------------------------------------

# Combine individuals from both ecotypes
ashn_combined = bind_rows(ashnc_f2_orig_trait_mat,
                          ashnw_f2_orig_trait_mat)

n_ashn_cold = nrow(ashnc_f2_orig_trait_mat)
n_ashn_warm = nrow(ashnw_f2_orig_trait_mat)

# Observed similarity
ashn_observed = evolqg::KrzCor(cov(ashnc_f2_orig_trait_mat),
                               cov(ashnw_f2_orig_trait_mat))

# Permutation test
set.seed(1738)

n_perm = 9999

ashn_null_krz = replicate(n_perm, {
  
  idx = sample(nrow(ashn_combined))
  
  ashn_cold_perm = ashn_combined[idx[1:n_ashn_cold], , drop = FALSE]
  
  ashn_warm_perm = ashn_combined[
    idx[(n_ashn_cold + 1):(n_ashn_cold + n_ashn_warm)],
    ,
    drop = FALSE]
  
  evolqg::KrzCor(cov(ashn_cold_perm),
                 cov(ashn_warm_perm))
})

# Lower similarity indicates greater divergence
ashn_p_value = (sum(ashn_null_krz <= ashn_observed) + 1)/(n_perm + 1)

ggplot(data.frame(ashn_KrzCor = ashn_null_krz),
       aes(x = ashn_KrzCor)) +
  geom_histogram(bins = 40,
                 fill = "grey75",
                 colour = "white") +
  geom_vline(xintercept = ashn_observed,
             colour = "#C44E52",
             linewidth = 1.2) +
  annotate("text",
           x = ashn_observed,
           y = Inf,
           label = paste0(
             "Observed = ", round(ashn_observed, 3),
             "\np = ", signif(ashn_p_value, 3)),
           vjust = 1.5,
           hjust = 1.1) +
  labs(x = "Krzanowski subspace similarity",
       y = "Permutation frequency",
       title = "ASHN geothermal–ambient comparison")


# MYV ---------------------------------------------------------------------


# Combine individuals from both ecotypes
myv_combined = bind_rows(myvc_f2_orig_trait_mat,
                         myvw_f2_orig_trait_mat)

n_myv_cold = nrow(myvc_f2_orig_trait_mat)
n_myv_warm = nrow(myvw_f2_orig_trait_mat)

# Observed similarity
myv_observed = evolqg::KrzCor(cov(myvc_f2_orig_trait_mat),
                              cov(myvw_f2_orig_trait_mat))

# Permutation test
set.seed(1738)

n_perm = 9999

myv_null_krz = replicate(n_perm, {
  
  idx = sample(nrow(myv_combined))
  
  myv_cold_perm = myv_combined[idx[1:n_myv_cold], , drop = FALSE]
  
  myv_warm_perm = myv_combined[
    idx[(n_myv_cold + 1):(n_myv_cold + n_myv_warm)],
    ,
    drop = FALSE]
  
  evolqg::KrzCor(cov(myv_cold_perm),
                 cov(myv_warm_perm))
})

# Lower similarity indicates greater divergence
myv_p_value = (sum(myv_null_krz <= myv_observed) + 1)/(n_perm + 1)

ggplot(data.frame(myv_KrzCor = myv_null_krz),
       aes(x = myv_KrzCor)) +
  geom_histogram(bins = 40,
                 fill = "grey75",
                 colour = "white") +
  geom_vline(xintercept = myv_observed,
             colour = "#C44E52",
             linewidth = 1.2) +
  annotate("text",
           x = myv_observed,
           y = Inf,
           label = paste0(
             "Observed = ", round(myv_observed, 3),
             "\np = ", signif(myv_p_value, 3)),
           vjust = 1.5,
           hjust = 1.1) +
  labs(x = "Krzanowski subspace similarity",
       y = "Permutation frequency",
       title = "MYV geothermal–ambient comparison")


# SKR ---------------------------------------------------------------------

# Combine individuals from both ecotypes
skr_combined = bind_rows(skrc_f2_orig_trait_mat,
                         skrw_f2_orig_trait_mat)

n_skr_cold = nrow(skrc_f2_orig_trait_mat)
n_skr_warm = nrow(skrw_f2_orig_trait_mat)

# Observed similarity
skr_observed = evolqg::KrzCor(cov(skrc_f2_orig_trait_mat),
                              cov(skrw_f2_orig_trait_mat))

# Permutation test
set.seed(1738)

n_perm = 9999

skr_null_krz = replicate(n_perm, {
  
  idx = sample(nrow(skr_combined))
  
  skr_cold_perm = skr_combined[idx[1:n_skr_cold], , drop = FALSE]
  
  skr_warm_perm = skr_combined[
    idx[(n_skr_cold + 1):(n_skr_cold + n_skr_warm)],
    ,
    drop = FALSE]
  
  evolqg::KrzCor(cov(skr_cold_perm),
                 cov(skr_warm_perm))
})

# Lower similarity indicates greater divergence
skr_p_value = (sum(skr_null_krz <= skr_observed) + 1)/(n_perm + 1)

ggplot(data.frame(skr_KrzCor = skr_null_krz),
       aes(x = skr_KrzCor)) +
  geom_histogram(bins = 40,
                 fill = "grey75",
                 colour = "white") +
  geom_vline(xintercept = skr_observed,
             colour = "#C44E52",
             linewidth = 1.2) +
  annotate("text",
           x = skr_observed,
           y = Inf,
           label = paste0(
             "Observed = ", round(skr_observed, 3),
             "\np = ", signif(skr_p_value, 3)),
           vjust = 1.5,
           hjust = 1.1) +
  labs(x = "Krzanowski subspace similarity",
       y = "Permutation frequency",
       title = "SKR geothermal–ambient comparison")


# GTS/GAR -----------------------------------------------------------------


# Combine individuals from both ecotypes
gts_gar_combined = bind_rows(GAR_f2_orig_trait_mat,
                             GTS_f2_orig_trait_mat)

n_GAR_cold = nrow(GAR_f2_orig_trait_mat)
n_GTS_warm = nrow(GTS_f2_orig_trait_mat)

# Observed similarity
gts_gar_observed = evolqg::KrzCor(cov(GAR_f2_orig_trait_mat),
                                  cov(GTS_f2_orig_trait_mat))

# Permutation test
set.seed(1738)

n_perm = 9999

gts_gar_null_krz = replicate(n_perm, {
  
  idx = sample(nrow(gts_gar_combined))
  
  GAR_cold_perm = gts_gar_combined[idx[1:n_GAR_cold], , drop = FALSE]
  
  GTS_warm_perm = gts_gar_combined[
    idx[(n_GAR_cold + 1):(n_GAR_cold + n_GTS_warm)],
    ,
    drop = FALSE]
  
  evolqg::KrzCor(cov(GAR_cold_perm),
                 cov(GTS_warm_perm))
})

# Lower similarity indicates greater divergence
gts_gar_p_value = (sum(gts_gar_null_krz <= gts_gar_observed) + 1)/(n_perm + 1)

ggplot(data.frame(gts_gar_KrzCor = gts_gar_null_krz),
       aes(x = gts_gar_KrzCor)) +
  geom_histogram(bins = 40,
                 fill = "grey75",
                 colour = "white") +
  geom_vline(xintercept = gts_gar_observed,
             colour = "#C44E52",
             linewidth = 1.2) +
  annotate("text",
           x = gts_gar_observed,
           y = Inf,
           label = paste0(
             "Observed = ", round(gts_gar_observed, 3),
             "\np = ", signif(gts_gar_p_value, 3)),
           vjust = 1.5,
           hjust = 1.1) +
  labs(x = "Krzanowski subspace similarity",
       y = "Permutation frequency",
       title = "GTS/GAR geothermal–ambient comparison")


# procrustes matrix superimposition ---------------------------------------
## tests differences in trait configuration derived from eigenvectors


# ASHN --------------------------------------------------------------------

ashnc_f2_12_12_eigen = eigen(ashnc_f2_12_12_cov_mat)
ashnw_f2_12_12_eigen = eigen(ashnw_f2_12_12_cov_mat)

ashnc_f2_12_18_eigen = eigen(ashnc_f2_12_18_cov_mat)
ashnw_f2_12_18_eigen = eigen(ashnw_f2_12_18_cov_mat)

ashnc_f2_18_12_eigen = eigen(ashnc_f2_18_12_cov_mat)
ashnw_f2_18_12_eigen = eigen(ashnw_f2_18_12_cov_mat)

ashnc_f2_18_18_eigen = eigen(ashnc_f2_18_18_cov_mat)
ashnw_f2_18_18_eigen = eigen(ashnw_f2_18_18_cov_mat)

# Extract leading eigenvectors
ashnc_f2_12_12_vectors = ashnc_f2_12_12_eigen$vectors[, 1:17, drop = FALSE]
ashnw_f2_12_12_vectors = ashnw_f2_12_12_eigen$vectors[, 1:17, drop = FALSE]

ashnc_f2_12_18_vectors = ashnc_f2_12_18_eigen$vectors[, 1:17, drop = FALSE]
ashnw_f2_12_18_vectors = ashnw_f2_12_18_eigen$vectors[, 1:17, drop = FALSE]

ashnc_f2_18_12_vectors = ashnc_f2_18_12_eigen$vectors[, 1:17, drop = FALSE]
ashnw_f2_18_12_vectors = ashnw_f2_18_12_eigen$vectors[, 1:17, drop = FALSE]

ashnc_f2_18_18_vectors = ashnc_f2_18_18_eigen$vectors[, 1:17, drop = FALSE]
ashnw_f2_18_18_vectors = ashnw_f2_18_18_eigen$vectors[, 1:17, drop = FALSE]


# Preserve trait identities
rownames(ashnc_f2_12_12_vectors) = colnames(ashnc_f2_12_12_cov_mat)
rownames(ashnw_f2_12_12_vectors) = colnames(ashnw_f2_12_12_cov_mat)

rownames(ashnc_f2_12_18_vectors) = colnames(ashnc_f2_12_18_cov_mat)
rownames(ashnw_f2_12_18_vectors) = colnames(ashnw_f2_12_18_cov_mat)

rownames(ashnc_f2_18_12_vectors) = colnames(ashnc_f2_18_12_cov_mat)
rownames(ashnw_f2_18_12_vectors) = colnames(ashnw_f2_18_12_cov_mat)

rownames(ashnc_f2_18_18_vectors) = colnames(ashnc_f2_18_18_cov_mat)
rownames(ashnw_f2_18_18_vectors) = colnames(ashnw_f2_18_18_cov_mat)


ashn_f2_12_12_proc_res = procrustes(ashnc_f2_12_12_vectors, 
                                    ashnw_f2_12_12_vectors, 
                                    symmetric = T)
summary(ashn_f2_12_12_proc_res)
plot(ashn_f2_12_12_proc_res, 
     kind = 2)
ashn_f2_12_12_test = protest(ashnc_f2_12_12_vectors, 
                        ashnw_f2_12_12_vectors, 
                        permutations = 999, 
                        scale = F, 
                        symmetric = T)

ashn_f2_12_12_residuals = residuals(ashn_f2_12_12_test)

ashn_f2_12_12_residuals%>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2) %>%  
  write_csv('ASHN_F2_12_12_procrustes_eigen_cov_mat_residuals.csv')


ashn_f2_12_18_proc_res = procrustes(ashnc_f2_12_18_vectors, 
                                    ashnw_f2_12_18_vectors, 
                                    symmetric = T)
summary(ashn_f2_12_18_proc_res)
plot(ashn_f2_12_18_proc_res, 
     kind = 2)
ashn_f2_12_18_test = protest(ashnc_f2_12_18_vectors, 
                             ashnw_f2_12_18_vectors, 
                             permutations = 999, 
                             scale = F, 
                             symmetric = T)
ashn_f2_12_18_residuals = residuals(ashn_f2_12_18_test)

ashn_f2_12_18_residuals%>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2) %>%  
  write_csv('ASHN_F2_12_18_procrustes_eigen_cov_mat_residuals.csv')


ashn_f2_18_12_proc_res = procrustes(ashnc_f2_18_12_vectors, 
                                    ashnw_f2_18_12_vectors, 
                                    symmetric = T)
summary(ashn_f2_18_12_proc_res)
plot(ashn_f2_18_12_proc_res, 
     kind = 2)

ashn_f2_18_12_test = protest(ashnc_f2_18_12_vectors, 
                             ashnw_f2_18_12_vectors, 
                             permutations = 999, 
                             scale = F, 
                             symmetric = T)
ashn_f2_18_12_residuals = residuals(ashn_f2_18_12_test)

ashn_f2_18_12_residuals%>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2) %>%  
  write_csv('ASHN_F2_18_12_procrustes_eigen_cov_mat_residuals.csv')


ashn_f2_18_18_proc_res = procrustes(ashnc_f2_18_18_vectors, 
                                    ashnw_f2_18_18_vectors, 
                                    symmetric = T)
summary(ashn_f2_18_18_proc_res)
plot(ashn_f2_18_18_proc_res, 
     kind = 2)
ashn_f2_18_18_test = protest(ashnc_f2_18_18_vectors, 
                             ashnw_f2_18_18_vectors, 
                             permutations = 999, 
                             scale = F, 
                             symmetric = T)

ashn_f2_18_18_residuals = residuals(ashn_f2_18_18_test)

ashn_f2_18_18_residuals%>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2) %>%  
  write_csv('ASHN_F2_18_18_procrustes_eigen_cov_mat_residuals.csv')


# MYV ---------------------------------------------------------------------


myvc_eigen = eigen(myvc_trait_cov_mat)
myvw_eigen = eigen(myvw_trait_cov_mat)

# Number of principal components to retain
k = floor(ncol(myvc_trait_cov_mat) / 2)

# Extract leading eigenvectors
myvc_vectors = myvc_eigen$vectors[, 1:3, drop = FALSE]
myvw_vectors = myvw_eigen$vectors[, 1:3, drop = FALSE]

# Preserve trait identities
rownames(myvc_vectors) = colnames(myvc_trait_cov_mat)
rownames(myvw_vectors) = colnames(myvw_trait_cov_mat)

myv_proc_res = procrustes(myvc_vectors, 
                          myvw_vectors, 
                          symmetric = T)

summary(myv_proc_res)

plot(myv_proc_res)
plot(myv_proc_res, 
     kind = 2)

myv_test_res = protest(myvc_vectors, 
                       myvw_vectors, 
                       permutations = 999)



# SKR ---------------------------------------------------------------------

skrc_eigen = eigen(skrc_trait_cov_mat)
skrw_eigen = eigen(skrw_trait_cov_mat)

# Number of principal components to retain
k = floor(ncol(skrc_trait_cov_mat) / 2)

# Extract leading eigenvectors
skrc_vectors = skrc_eigen$vectors[, 1:3, drop = FALSE]
skrw_vectors = skrw_eigen$vectors[, 1:3, drop = FALSE]

# Preserve trait identities
rownames(skrc_vectors) = colnames(skrc_trait_cov_mat)
rownames(skrw_vectors) = colnames(skrw_trait_cov_mat)

skr_proc_res = procrustes(skrc_vectors, 
                          skrw_vectors, 
                          symmetric = T)

summary(skr_proc_res)

plot(skr_proc_res)
plot(skr_proc_res, 
     kind = 2)

skr_test_res = protest(skrc_vectors, 
                       skrw_vectors, 
                       permutations = 999)



# GTS/GAR -----------------------------------------------------------------

gar_eigen = eigen(GAR_trait_cov_mat)
gts_eigen = eigen(GTS_trait_cov_mat)


# Extract leading eigenvectors
gar_vectors = gar_eigen$vectors[, 1:3, drop = FALSE]
gts_vectors = gts_eigen$vectors[, 1:3, drop = FALSE]

# Preserve trait identities
rownames(gar_vectors) = colnames(GAR_trait_cov_mat)
rownames(gts_vectors) = colnames(GTS_trait_cov_mat)

gts_gar_proc_res = procrustes(gar_vectors, 
                              gts_vectors, 
                              symmetric = T)

summary(gts_gar_proc_res)

plot(gts_gar_proc_res)
plot(gts_gar_proc_res, 
     kind = 2)

gts_gar_test_res = protest(gar_vectors, 
                           gts_vectors, 
                           permutations = 999)

# Principal axes of divergence?  ------------------------------------------

ashn_eig_amb <- eigen(ashnc_trait_cov_mat)
ashn_eig_geo <- eigen(ashnw_trait_cov_mat)

ashn_loadings <- tibble(
  Trait = rownames(ashnc_trait_cov_mat),
  Ambient_PC1 = ashn_eig_amb$vectors[, 1],
  Ambient_PC2 = ashn_eig_amb$vectors[, 2],
  Ambient_PC3 = ashn_eig_amb$vectors[, 3],
  Geothermal_PC1 = ashn_eig_geo$vectors[, 1], 
  Geothermal_PC2 = ashn_eig_geo$vectors[, 2], 
  Geothermal_PC3 = ashn_eig_geo$vectors[, 3])

# Align arbitrary eigenvector signs
if (sum(ashn_loadings$Ambient_PC1 * ashn_loadings$Geothermal_PC1) < 0) {
  ashn_loadings$Geothermal_PC1 <- -ashn_loadings$Geothermal_PC1
}

if (sum(ashn_loadings$Ambient_PC2 * ashn_loadings$Geothermal_PC2) < 0) {
  ashn_loadings$Geothermal_PC2 <- -ashn_loadings$Geothermal_PC2
}
if (sum(ashn_loadings$Ambient_PC3 * ashn_loadings$Geothermal_PC3) < 0) {
  ashn_loadings$Geothermal_PC3 <- -ashn_loadings$Geothermal_PC3
}


ashn_pc_loadings = ggplot(ashn_loadings, 
                          aes(Ambient_PC1, 
                              Geothermal_PC1)) +
  # geom_abline(slope = 1, 
  #             intercept = 0, 
  #             linetype = 2) +
  geom_hline(yintercept = 0.0, 
             linetype = 2)+
  geom_vline(xintercept = 0.0, 
             linetype = 2)+
  geom_point(size = 3) +
  geom_text(aes(label = Trait), 
            vjust = -0.7) +
  coord_equal() +
  theme_classic() +
  labs(
    x = "Ambient PC1 loading",
    y = "Geothermal PC1 loading", 
    title = 'A) ASHN')+
  
  ggplot(ashn_loadings, 
         aes(Ambient_PC2, 
             Geothermal_PC2)) +
  # geom_abline(slope = 1, 
  #             intercept = 0, 
  #             linetype = 2) +
  geom_hline(yintercept = 0.0, 
             linetype = 2)+
  geom_vline(xintercept = 0.0, 
             linetype = 2)+
  geom_point(size = 3) +
  geom_text(aes(label = Trait), 
            vjust = -0.7) +
  coord_equal() +
  theme_classic() +
  labs(
    x = "Ambient PC2 loading",
    y = "Geothermal PC2 loading"
  )+
  
  ggplot(ashn_loadings, 
         aes(Ambient_PC3, 
             Geothermal_PC3)) +
  # geom_abline(slope = 1, 
  #             intercept = 0, 
  #             linetype = 2) +
  geom_hline(yintercept = 0.0, 
             linetype = 2)+
  geom_vline(xintercept = 0.0, 
             linetype = 2)+
  geom_point(size = 3) +
  geom_text(aes(label = Trait), 
            vjust = -0.7) +
  coord_equal() +
  theme_classic() +
  labs(
    x = "Ambient PC3 loading",
    y = "Geothermal PC3 loading"
  )


ashnc_loading_pcs = ashn_loadings %>%
  dplyr::select(Ambient_PC1, 
                Ambient_PC2, 
                Ambient_PC3)
ashnw_loading_pcs = ashn_loadings %>%
  dplyr::select(Geothermal_PC1, 
                Geothermal_PC2, 
                Geothermal_PC3)

ashn_loading_proc_res = protest(ashnc_loading_pcs, 
                                ashnw_loading_pcs, 
                                permutations = 999)


# MYV pc loadings ---------------------------------------------------------


myv_eig_amb <- eigen(myvc_trait_cov_mat)
myv_eig_geo <- eigen(myvw_trait_cov_mat)

myv_loadings <- tibble(
  Trait = rownames(myvc_trait_cov_mat),
  Ambient_PC1 = myv_eig_amb$vectors[, 1],
  Ambient_PC2 = myv_eig_amb$vectors[, 2],
  Ambient_PC3 = myv_eig_amb$vectors[, 3],
  Geothermal_PC1 = myv_eig_geo$vectors[, 1], 
  Geothermal_PC2 = myv_eig_geo$vectors[, 2], 
  Geothermal_PC3 = myv_eig_geo$vectors[, 3])

# Align arbitrary eigenvector signs
if (sum(myv_loadings$Ambient_PC1 * myv_loadings$Geothermal_PC1) < 0) {
  myv_loadings$Geothermal_PC1 <- -myv_loadings$Geothermal_PC1
}

if (sum(myv_loadings$Ambient_PC2 * myv_loadings$Geothermal_PC2) < 0) {
  myv_loadings$Geothermal_PC2 <- -myv_loadings$Geothermal_PC2
}
if (sum(myv_loadings$Ambient_PC3 * myv_loadings$Geothermal_PC3) < 0) {
  myv_loadings$Geothermal_PC3 <- -myv_loadings$Geothermal_PC3
}


myv_pc_loadings = ggplot(myv_loadings, 
                         aes(Ambient_PC1, 
                             Geothermal_PC1)) +
  # geom_abline(slope = 1, 
  #             intercept = 0, 
  #             linetype = 2) +
  geom_hline(yintercept = 0.0, 
             linetype = 2)+
  geom_vline(xintercept = 0.0, 
             linetype = 2)+
  geom_point(size = 3) +
  geom_text(aes(label = Trait), 
            vjust = -0.7) +
  coord_equal() +
  theme_classic() +
  labs(
    x = "Ambient PC1 loading",
    y = "Geothermal PC1 loading", 
    title = 'B) MYV'
  )+
  
  ggplot(myv_loadings, 
         aes(Ambient_PC2, 
             Geothermal_PC2)) +
  # geom_abline(slope = 1, 
  #             intercept = 0, 
  #             linetype = 2) +
  geom_hline(yintercept = 0.0, 
             linetype = 2)+
  geom_vline(xintercept = 0.0, 
             linetype = 2)+
  geom_point(size = 3) +
  geom_text(aes(label = Trait), 
            vjust = -0.7) +
  coord_equal() +
  theme_classic() +
  labs(
    x = "Ambient PC2 loading",
    y = "Geothermal PC2 loading"
  )+
  
  ggplot(myv_loadings, 
         aes(Ambient_PC3, 
             Geothermal_PC3)) +
  # geom_abline(slope = 1, 
  #             intercept = 0, 
  #             linetype = 2) +
  geom_hline(yintercept = 0.0, 
             linetype = 2)+
  geom_vline(xintercept = 0.0, 
             linetype = 2)+
  geom_point(size = 3) +
  geom_text(aes(label = Trait), 
            vjust = -0.7) +
  coord_equal() +
  theme_classic() +
  labs(
    x = "Ambient PC3 loading",
    y = "Geothermal PC3 loading"
  )



# SKR pc loadings ---------------------------------------------------------


skr_eig_amb <- eigen(skrc_trait_cov_mat)
skr_eig_geo <- eigen(skrw_trait_cov_mat)

skr_loadings <- tibble(
  Trait = rownames(skrc_trait_cov_mat),
  Ambient_PC1 = skr_eig_amb$vectors[, 1],
  Ambient_PC2 = skr_eig_amb$vectors[, 2],
  Ambient_PC3 = skr_eig_amb$vectors[, 3],
  Geothermal_PC1 = skr_eig_geo$vectors[, 1], 
  Geothermal_PC2 = skr_eig_geo$vectors[, 2], 
  Geothermal_PC3 = skr_eig_geo$vectors[, 3])

# Align arbitrary eigenvector signs
if (sum(skr_loadings$Ambient_PC1 * skr_loadings$Geothermal_PC1) < 0) {
  skr_loadings$Geothermal_PC1 <- -skr_loadings$Geothermal_PC1
}

if (sum(skr_loadings$Ambient_PC2 * skr_loadings$Geothermal_PC2) < 0) {
  skr_loadings$Geothermal_PC2 <- -skr_loadings$Geothermal_PC2
}
if (sum(skr_loadings$Ambient_PC3 * skr_loadings$Geothermal_PC3) < 0) {
  skr_loadings$Geothermal_PC3 <- -skr_loadings$Geothermal_PC3
}


skr_pc_loadings = ggplot(skr_loadings, 
                         aes(Ambient_PC1, 
                             Geothermal_PC1)) +
  # geom_abline(slope = 1, 
  #             intercept = 0, 
  #             linetype = 2) +
  geom_hline(yintercept = 0.0, 
             linetype = 2)+
  geom_vline(xintercept = 0.0, 
             linetype = 2)+
  geom_point(size = 3) +
  geom_text(aes(label = Trait), 
            vjust = -0.7) +
  coord_equal() +
  theme_classic() +
  labs(
    x = "Ambient PC1 loading",
    y = "Geothermal PC1 loading", 
    title = 'C) SKR'
  )+
  
  ggplot(skr_loadings, 
         aes(Ambient_PC2, 
             Geothermal_PC2)) +
  # geom_abline(slope = 1, 
  #             intercept = 0, 
  #             linetype = 2) +
  geom_hline(yintercept = 0.0, 
             linetype = 2)+
  geom_vline(xintercept = 0.0, 
             linetype = 2)+
  geom_point(size = 3) +
  geom_text(aes(label = Trait), 
            vjust = -0.7) +
  coord_equal() +
  theme_classic() +
  labs(
    x = "Ambient PC2 loading",
    y = "Geothermal PC2 loading"
  )+
  
  ggplot(skr_loadings, 
         aes(Ambient_PC3, 
             Geothermal_PC3)) +
  # geom_abline(slope = 1, 
  #             intercept = 0, 
  #             linetype = 2) +
  geom_hline(yintercept = 0.0, 
             linetype = 2)+
  geom_vline(xintercept = 0.0, 
             linetype = 2)+
  geom_point(size = 3) +
  geom_text(aes(label = Trait), 
            vjust = -0.7) +
  coord_equal() +
  theme_classic() +
  labs(
    x = "Ambient PC3 loading",
    y = "Geothermal PC3 loading"
  )


# GTS/GAR -----------------------------------------------------------------


GAR_eig_amb <- eigen(GAR_trait_cov_mat)
GTS_eig_geo <- eigen(GTS_trait_cov_mat)

gts_gar_loadings <- tibble(
  Trait = rownames(GAR_trait_cov_mat),
  Ambient_PC1 = GAR_eig_amb$vectors[, 1],
  Ambient_PC2 = GAR_eig_amb$vectors[, 2],
  Ambient_PC3 = GAR_eig_amb$vectors[, 3],
  Geothermal_PC1 = GTS_eig_geo$vectors[, 1], 
  Geothermal_PC2 = GTS_eig_geo$vectors[, 2], 
  Geothermal_PC3 = GTS_eig_geo$vectors[, 3])

# Align arbitrary eigenvector signs
if (sum(gts_gar_loadings$Ambient_PC1 * gts_gar_loadings$Geothermal_PC1) < 0) {
  gts_gar_loadings$Geothermal_PC1 <- -gts_gar_loadings$Geothermal_PC1
}

if (sum(gts_gar_loadings$Ambient_PC2 * gts_gar_loadings$Geothermal_PC2) < 0) {
  gts_gar_loadings$Geothermal_PC2 <- -gts_gar_loadings$Geothermal_PC2
}
if (sum(gts_gar_loadings$Ambient_PC3 * gts_gar_loadings$Geothermal_PC3) < 0) {
  gts_gar_loadings$Geothermal_PC3 <- -gts_gar_loadings$Geothermal_PC3
}


gts_gar_pc_loadings = ggplot(gts_gar_loadings, 
                             aes(Ambient_PC1, 
                                 Geothermal_PC1)) +
  # geom_abline(slope = 1, 
  #             intercept = 0, 
  #             linetype = 2) +
  geom_hline(yintercept = 0.0, 
             linetype = 2)+
  geom_vline(xintercept = 0.0, 
             linetype = 2)+
  geom_point(size = 3) +
  geom_text(aes(label = Trait), 
            vjust = -0.7) +
  coord_equal() +
  theme_classic() +
  labs(
    x = "Ambient PC1 loading",
    y = "Geothermal PC1 loading", 
    title = 'D) GTS/GAR'
  )+
  
  ggplot(gts_gar_loadings, 
         aes(Ambient_PC2, 
             Geothermal_PC2)) +
  # geom_abline(slope = 1, 
  #             intercept = 0, 
  #             linetype = 2) +
  geom_hline(yintercept = 0.0, 
             linetype = 2)+
  geom_vline(xintercept = 0.0, 
             linetype = 2)+
  geom_point(size = 3) +
  geom_text(aes(label = Trait), 
            vjust = -0.7) +
  coord_equal() +
  theme_classic() +
  labs(
    x = "Ambient PC2 loading",
    y = "Geothermal PC2 loading"
  )+
  
  ggplot(gts_gar_loadings, 
         aes(Ambient_PC3, 
             Geothermal_PC3)) +
  # geom_abline(slope = 1, 
  #             intercept = 0, 
  #             linetype = 2) +
  geom_hline(yintercept = 0.0, 
             linetype = 2)+
  geom_vline(xintercept = 0.0, 
             linetype = 2)+
  geom_point(size = 3) +
  geom_text(aes(label = Trait), 
            vjust = -0.7) +
  coord_equal() +
  theme_classic() +
  labs(
    x = "Ambient PC3 loading",
    y = "Geothermal PC3 loading"
  )


ashn_pc_loadings/myv_pc_loadings/skr_pc_loadings/gts_gar_pc_loadings



# Frobenius distance between covariance matrices --------------------------

# The two statistics address slightly different aspects of 
# morphological covariance divergence:

# Krzanowski similarity: Do geothermal and ambient ecotypes share 
# the same major axes of morphological variation?

# Frobenius distance: How different are their overall covariance matrices, 
# accounting for both the magnitude and structure of trait covariation?

# cov_distance <- function(A, B) {
#   sqrt(sum((A - B)^2))
# }
# 
# distances <- tibble(
#   Pair = paste0("Pair", 1:4),
#   Distance = map_dbl(1:4, function(i) {
#     cov_distance(
#       cov_matrices[[paste0("Pair", i, "_Ambient")]],
#       cov_matrices[[paste0("Pair", i, "_Geothermal")]]
#     )
#   })
# ) %>% 
#   mutate(Ecotype_pair = case_when(
#     str_detect(Pair, 'Pair1') ~ 'ASHN',
#     str_detect(Pair, 'Pair2') ~ 'MYV',
#     str_detect(Pair, 'Pair3') ~ 'SKR',
#     str_detect(Pair, 'Pair4') ~ 'GTS/GAR'))
# 
# ggplot(distances, 
#        aes(Ecotype_pair, 
#            Distance)) +
#   geom_col(fill = "steelblue") +
#   labs(
#     x = "Population pair",
#     y = "Covariance matrix distance")


# ASHN --------------------------------------------------------------------

# Function to calculate Frobenius distance
frob_distance = function(A, B) {
  sqrt(sum((A - B)^2))
}

# Observed Frobenius distance
ashn_observed_frob = frob_distance(ashnc_trait_cov_mat,
                                   ashnw_trait_cov_mat)


# Number of permutations
n_perm <- 9999

set.seed(8008)

# Generate null distribution
ashn_null_frob <- replicate(n_perm, {
  
  # Randomly shuffle individuals
  idx <- sample(nrow(ashn_combined))
  
  # Reassign individuals to groups
  ashn_cold_perm <- ashn_combined[
    idx[1:n_ashn_cold],
    ,
    drop = FALSE
  ]
  
  ashn_warm_perm <- ashn_combined[
    idx[(n_ashn_cold + 1):(n_ashn_cold + n_ashn_warm)],
    ,
    drop = FALSE
  ]
  
  # Calculate covariance matrices
  cov_ashn_cold <- cov(ashn_cold_perm)
  cov_ashn_warm <- cov(ashn_warm_perm)
  
  # Frobenius distance
  frob_distance(cov_ashn_cold, cov_ashn_warm)
})

# One-sided permutation test
# Greater distance = greater covariance divergence
ashn_p_frob = (sum(ashn_null_frob >= ashn_observed_frob) + 1) /
  (n_perm + 1)

ggplot(
  data.frame(ashn_Frobenius = ashn_null_frob),
  aes(x = ashn_Frobenius)) +
  geom_histogram(bins = 40,
                 fill = "grey75",
                 colour = "white") +
  geom_vline(xintercept = ashn_observed_frob,
             colour = "#C44E52",
             linewidth = 1.2) +
  annotate("text",
           x = ashn_observed_frob,
           y = Inf,
           label = paste0(
             "Observed = ", round(ashn_observed_frob, 3),
             "\np = ", signif(ashn_p_frob, 3)),
           vjust = 1.5,
           hjust = 1.1) +
  labs(x = "Frobenius distance",
       y = "Permutation frequency",
       title = "ASHN geothermal–ambient covariance divergence")

# MYV ---------------------------------------------------------------------

# Observed Frobenius distance
myv_observed_frob = frob_distance(myvc_trait_cov_mat,
                                  myvw_trait_cov_mat)


# Number of permutations
n_perm <- 9999

set.seed(8008)

# Generate null distribution
myv_null_frob <- replicate(n_perm, {
  
  # Randomly shuffle individuals
  idx <- sample(nrow(myv_combined))
  
  # Reassign individuals to groups
  myv_cold_perm <- myv_combined[
    idx[1:n_myv_cold],
    ,
    drop = FALSE
  ]
  
  myv_warm_perm <- myv_combined[
    idx[(n_myv_cold + 1):(n_myv_cold + n_myv_warm)],
    ,
    drop = FALSE
  ]
  
  # Calculate covariance matrices
  cov_myv_cold <- cov(myv_cold_perm)
  cov_myv_warm <- cov(myv_warm_perm)
  
  # Frobenius distance
  frob_distance(cov_myv_cold, cov_myv_warm)
})

# One-sided permutation test
# Greater distance = greater covariance divergence
myv_p_frob = (sum(myv_null_frob >= myv_observed_frob) + 1) /
  (n_perm + 1)

ggplot(
  data.frame(myv_Frobenius = myv_null_frob),
  aes(x = myv_Frobenius)) +
  geom_histogram(bins = 40,
                 fill = "grey75",
                 colour = "white") +
  geom_vline(xintercept = myv_observed_frob,
             colour = "#C44E52",
             linewidth = 1.2) +
  annotate("text",
           x = myv_observed_frob,
           y = Inf,
           label = paste0(
             "Observed = ", round(myv_observed_frob, 3),
             "\np = ", signif(myv_p_frob, 3)),
           vjust = 1.5,
           hjust = 1.1) +
  labs(x = "Frobenius distance",
       y = "Permutation frequency",
       title = "MYV geothermal–ambient covariance divergence")

# SKR ---------------------------------------------------------------------

# Observed Frobenius distance
skr_observed_frob = frob_distance(skrc_trait_cov_mat,
                                  skrw_trait_cov_mat)


# Number of permutations
n_perm <- 9999

set.seed(8008)

# Generate null distribution
skr_null_frob <- replicate(n_perm, {
  
  # Randomly shuffle individuals
  idx <- sample(nrow(skr_combined))
  
  # Reassign individuals to groups
  skr_cold_perm <- skr_combined[
    idx[1:n_skr_cold],
    ,
    drop = FALSE
  ]
  
  skr_warm_perm <- skr_combined[
    idx[(n_skr_cold + 1):(n_skr_cold + n_skr_warm)],
    ,
    drop = FALSE
  ]
  
  # Calculate covariance matrices
  cov_skr_cold <- cov(skr_cold_perm)
  cov_skr_warm <- cov(skr_warm_perm)
  
  # Frobenius distance
  frob_distance(cov_skr_cold, cov_skr_warm)
})

# One-sided permutation test
# Greater distance = greater covariance divergence
skr_p_frob = (sum(skr_null_frob >= skr_observed_frob) + 1) /
  (n_perm + 1)

ggplot(
  data.frame(skr_Frobenius = skr_null_frob),
  aes(x = skr_Frobenius)) +
  geom_histogram(bins = 40,
                 fill = "grey75",
                 colour = "white") +
  geom_vline(xintercept = skr_observed_frob,
             colour = "#C44E52",
             linewidth = 1.2) +
  annotate("text",
           x = skr_observed_frob,
           y = Inf,
           label = paste0(
             "Observed = ", round(skr_observed_frob, 3),
             "\np = ", signif(skr_p_frob, 3)),
           vjust = 1.5,
           hjust = 1.1) +
  labs(x = "Frobenius distance",
       y = "Permutation frequency",
       title = "SKR geothermal–ambient covariance divergence")

# GTS/GAR -----------------------------------------------------------------

# Observed Frobenius distance
gts_gar_observed_frob = frob_distance(GAR_trait_cov_mat,
                                      GTS_trait_cov_mat)


# Number of permutations
n_perm <- 9999

set.seed(8008)

# Generate null distribution
gts_gar_null_frob <- replicate(n_perm, {
  
  # Randomly shuffle individuals
  idx <- sample(nrow(gts_gar_combined))
  
  # Reassign individuals to groups
  GAR_cold_perm <- gts_gar_combined[
    idx[1:n_GAR_cold],
    ,
    drop = FALSE
  ]
  
  GTS_warm_perm <- gts_gar_combined[
    idx[(n_GAR_cold + 1):(n_GAR_cold + n_GTS_warm)],
    ,
    drop = FALSE
  ]
  
  # Calculate covariance matrices
  cov_GAR_cold <- cov(GAR_cold_perm)
  cov_GTS_warm <- cov(GTS_warm_perm)
  
  # Frobenius distance
  frob_distance(cov_GAR_cold, cov_GTS_warm)
})

# One-sided permutation test
# Greater distance = greater covariance divergence
gts_gar_p_frob = (sum(gts_gar_null_frob >= gts_gar_observed_frob) + 1) /
  (n_perm + 1)

ggplot(
  data.frame(gts_gar_Frobenius = gts_gar_null_frob),
  aes(x = gts_gar_Frobenius)) +
  geom_histogram(bins = 40,
                 fill = "grey75",
                 colour = "white") +
  geom_vline(xintercept = gts_gar_observed_frob,
             colour = "#C44E52",
             linewidth = 1.2) +
  annotate("text",
           x = gts_gar_observed_frob,
           y = Inf,
           label = paste0(
             "Observed = ", round(gts_gar_observed_frob, 3),
             "\np = ", signif(gts_gar_p_frob, 3)),
           vjust = 1.5,
           hjust = 1.1) +
  labs(x = "Frobenius distance",
       y = "Permutation frequency",
       title = "GTS/GAR geothermal–ambient covariance divergence")

