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


wgp_orig_data = read_csv("WGP_TRAITS_SCALED_FIXED_11.02.2026.csv") %>% 
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



# Procrustes superimposition - cov matrices -------------------------------



# ASHN --------------------------------------------------------------------


ashnc_wgp_trait_mat = wgp_orig_data %>% 
  filter(Lake_morph == 'ASHNC') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

ashnw_wgp_trait_mat = wgp_orig_data %>% 
  filter(Lake_morph == 'ASHNW') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

ashnc_wgp_cov_mat = cov(ashnc_wgp_trait_mat, 
                       use = 'complete.obs')
ashnw_wgp_cov_mat = cov(ashnw_wgp_trait_mat, 
                       use = 'complete.obs')
ashn_wgp_krzcov = evolqg::KrzCor(ashnc_wgp_cov_mat, 
                                ashnw_wgp_cov_mat)
ashn_wgp_proc_res = procrustes(ashnc_wgp_cov_mat, 
                              ashnw_wgp_cov_mat)

# plot(ashn_wgp_proc_res, 
#      kind = 1)

ashn_wgp_proc_test = protest(ashnc_wgp_cov_mat, 
                            ashnw_wgp_cov_mat, 
                            permutations = 999)


ashn_wgp_residuals = residuals(ashn_wgp_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                cov_residual = 2)

# %>% 
#   mutate(resid_10 = residual*10)

ashn_wgp_resid_plot = ggplot(ashn_wgp_residuals, 
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

myvc_wgp_trait_mat = wgp_orig_data %>% 
  filter(Lake_morph == 'MYVC') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

myvw_wgp_trait_mat = wgp_orig_data %>% 
  filter(Lake_morph == 'MYVW') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

myvc_wgp_cov_mat = cov(myvc_wgp_trait_mat, 
                      use = 'complete.obs')
myvw_wgp_cov_mat = cov(myvw_wgp_trait_mat, 
                      use = 'complete.obs')
myv_wgp_krzcov = evolqg::KrzCor(myvc_wgp_cov_mat, 
                               myvw_wgp_cov_mat)
myv_wgp_proc_res = procrustes(myvc_wgp_cov_mat, 
                             myvw_wgp_cov_mat)

plot(myv_wgp_proc_res, 
     kind = 1)

myv_wgp_proc_test = protest(myvc_wgp_cov_mat, 
                           myvw_wgp_cov_mat, 
                           permutations = 999)


myv_wgp_residuals = residuals(myv_wgp_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                cov_residual = 2)

myv_wgp_resid_plot = ggplot(myv_wgp_residuals, 
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

skrc_wgp_trait_mat = wgp_orig_data %>% 
  filter(Lake_morph == 'SKRC') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

skrw_wgp_trait_mat = wgp_orig_data %>% 
  filter(Lake_morph == 'SKRW') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

skrc_wgp_cov_mat = cov(skrc_wgp_trait_mat, 
                      use = 'complete.obs')
skrw_wgp_cov_mat = cov(skrw_wgp_trait_mat, 
                      use = 'complete.obs')
skr_wgp_krzcov = evolqg::KrzCor(skrc_wgp_cov_mat, 
                               skrw_wgp_cov_mat)
skr_wgp_proc_res = procrustes(skrc_wgp_cov_mat, 
                             skrw_wgp_cov_mat)

plot(skr_wgp_proc_res, 
     kind = 1)

skr_wgp_proc_test = protest(skrc_wgp_cov_mat, 
                           skrw_wgp_cov_mat, 
                           permutations = 999)


skr_wgp_residuals = residuals(skr_wgp_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                cov_residual = 2)

skr_wgp_resid_plot = ggplot(skr_wgp_residuals, 
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

gar_wgp_trait_mat = wgp_orig_data %>% 
  filter(Lake_morph == 'CSWYC') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

gts_wgp_trait_mat = wgp_orig_data %>% 
  filter(Lake_morph == 'GTSW') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

gar_wgp_cov_mat = cov(gar_wgp_trait_mat, 
                     use = 'complete.obs')
gts_wgp_cov_mat = cov(gts_wgp_trait_mat, 
                     use = 'complete.obs')
gts_gar_wgp_krzcov = evolqg::KrzCor(gar_wgp_cov_mat, 
                                   gts_wgp_cov_mat)
gts_gar_wgp_proc_res = procrustes(gar_wgp_cov_mat, 
                                 gts_wgp_cov_mat)

plot(gts_gar_wgp_proc_res, 
     kind = 1)

gts_gar_wgp_proc_test = protest(gar_wgp_cov_mat, 
                               gts_wgp_cov_mat, 
                               permutations = 999)


gts_gar_wgp_residuals = residuals(gts_gar_wgp_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                cov_residual = 2)

gts_gar_wgp_resid_plot = ggplot(gts_gar_wgp_residuals, 
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


# PCA find significant eigenvectors ---------------------------------------

ashnc_wgp_pca = prcomp(ashnc_wgp_cov_mat, 
                      center = F, 
                      scale. = F)
paran(ashnc_wgp_cov_mat, 
      iterations = 1000, 
      graph = T)
## 4 PC axes retained
ashnw_wgp_pca = prcomp(ashnw_wgp_cov_mat, 
                      center = F, 
                      scale. = F)
paran(ashnw_wgp_cov_mat, 
      iterations = 1000, 
      graph = T)
## 3 axes

myvc_wgp_pca = prcomp(myvc_wgp_cov_mat, 
                     center = F, 
                     scale. = F)
paran(myvc_wgp_cov_mat, 
      iterations = 1000, 
      graph = T)
## 3 PC axes retained
myvw_wgp_pca = prcomp(myvw_wgp_cov_mat, 
                     center = F, 
                     scale. = F)
paran(myvw_wgp_cov_mat, 
      iterations = 1000, 
      graph = T)
## 3 axes

skrc_wgp_pca = prcomp(skrc_wgp_cov_mat, 
                     center = F, 
                     scale. = F)
paran(skrc_wgp_cov_mat, 
      iterations = 1000, 
      graph = T)
## 3 PC axes retained
skrw_wgp_pca = prcomp(skrw_wgp_cov_mat, 
                     center = F, 
                     scale. = F)
paran(skrw_wgp_cov_mat, 
      iterations = 1000, 
      graph = T)
## 3 axes

gar_wgp_pca = prcomp(gar_wgp_cov_mat, 
                    center = F, 
                    scale. = F)
paran(gar_wgp_cov_mat, 
      iterations = 1000, 
      graph = T)
## 3 PC axes retained
gts_wgp_pca = prcomp(gts_wgp_cov_mat, 
                    center = F, 
                    scale. = F)
paran(gts_wgp_cov_mat, 
      iterations = 1000, 
      graph = T)
## 4 axes

# ASHN ---------------------------------------------------------------------


ashnc_wgp_eigen = eigen(ashnc_wgp_cov_mat)
ashnw_wgp_eigen = eigen(ashnw_wgp_cov_mat)

# Extract leading eigenvectors
ashnc_wgp_vectors = ashnc_wgp_eigen$vectors[, 1:4, drop = FALSE]
ashnw_wgp_vectors = ashnw_wgp_eigen$vectors[, 1:4, drop = FALSE]

# Preserve trait identities
rownames(ashnc_wgp_vectors) = colnames(ashnc_wgp_cov_mat)
rownames(ashnw_wgp_vectors) = colnames(ashnw_wgp_cov_mat)


ashn_wgp_proc_res_eigen = procrustes(ashnc_wgp_vectors, 
                                    ashnw_wgp_vectors, 
                                    symmetric = T)
summary(ashn_wgp_proc_res_eigen)
# plot(ashn_wgp_proc_res_eigen, 
#      kind = 1)

ashn_wgp_eigen_test = protest(ashnc_wgp_vectors, 
                             ashnw_wgp_vectors, 
                             permutations = 999, 
                             scale = F, 
                             symmetric = T)

ashn_wgp_eigen_residuals = residuals(ashn_wgp_eigen_test)

ashn_wgp_eigen_residuals = ashn_wgp_eigen_residuals%>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                eigen_residual = 2) 
# %>%  
#   write_csv('ASHN_wgp_12_12_procrustes_eigen_cov_mat_residuals.csv')


# bind_cols(ashn_wgp_residuals, 
#           ashn_wgp_eigen_residuals)
# 
ashn_wgp_residuals_combo = inner_join(ashn_wgp_residuals,
                                     ashn_wgp_eigen_residuals,
                                     by = 'trait') %>%
  mutate(pop_pair = 'ASHN')

# %>%
#   write_csv('ASHN_wgp_procrustes_eigen_cov_mat_residuals.csv')

# MYV ---------------------------------------------------------------------

myvc_wgp_eigen = eigen(myvc_wgp_cov_mat)
myvw_wgp_eigen = eigen(myvw_wgp_cov_mat)

# Extract leading eigenvectors
myvc_wgp_vectors = myvc_wgp_eigen$vectors[, 1:4, drop = FALSE]
myvw_wgp_vectors = myvw_wgp_eigen$vectors[, 1:4, drop = FALSE]

# Preserve trait identities
rownames(myvc_wgp_vectors) = colnames(myvc_wgp_cov_mat)
rownames(myvw_wgp_vectors) = colnames(myvw_wgp_cov_mat)


myv_wgp_proc_res_eigen = procrustes(myvc_wgp_vectors, 
                                   myvw_wgp_vectors, 
                                   symmetric = T)
summary(myv_wgp_proc_res_eigen)
# plot(myv_wgp_proc_res_eigen, 
#      kind = 1)

myv_wgp_eigen_test = protest(myvc_wgp_vectors, 
                            myvw_wgp_vectors, 
                            permutations = 999, 
                            scale = F, 
                            symmetric = T)

myv_wgp_eigen_residuals = residuals(myv_wgp_eigen_test)

myv_wgp_eigen_residuals = myv_wgp_eigen_residuals%>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                eigen_residual = 2) 
# %>%  
#   write_csv('myv_wgp_12_12_procrustes_eigen_cov_mat_residuals.csv')



myv_wgp_residuals_combo = inner_join(myv_wgp_residuals, 
                                    myv_wgp_eigen_residuals, 
                                    by = 'trait') %>% 
  mutate(pop_pair = 'MYV')


# SKR ---------------------------------------------------------------------

skrc_wgp_eigen = eigen(skrc_wgp_cov_mat)
skrw_wgp_eigen = eigen(skrw_wgp_cov_mat)

# Extract leading eigenvectors
skrc_wgp_vectors = skrc_wgp_eigen$vectors[, 1:4, drop = FALSE]
skrw_wgp_vectors = skrw_wgp_eigen$vectors[, 1:4, drop = FALSE]

# Preserve trait identities
rownames(skrc_wgp_vectors) = colnames(skrc_wgp_cov_mat)
rownames(skrw_wgp_vectors) = colnames(skrw_wgp_cov_mat)


skr_wgp_proc_res_eigen = procrustes(skrc_wgp_vectors, 
                                   skrw_wgp_vectors, 
                                   symmetric = T)
summary(skr_wgp_proc_res_eigen)
# plot(skr_wgp_proc_res_eigen, 
#      kind = 1)

skr_wgp_eigen_test = protest(skrc_wgp_vectors, 
                            skrw_wgp_vectors, 
                            permutations = 999, 
                            scale = F, 
                            symmetric = T)

skr_wgp_eigen_residuals = residuals(skr_wgp_eigen_test)

skr_wgp_eigen_residuals = skr_wgp_eigen_residuals%>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                eigen_residual = 2) 
# %>%  
#   write_csv('skr_wgp_12_12_procrustes_eigen_cov_mat_residuals.csv')



skr_wgp_residuals_combo = inner_join(skr_wgp_residuals, 
                                    skr_wgp_eigen_residuals, 
                                    by = 'trait') %>% 
  mutate(pop_pair = 'SKR')

# GTS/GAR -----------------------------------------------------------------

gar_wgp_eigen = eigen(gar_wgp_cov_mat)
gts_wgp_eigen = eigen(gts_wgp_cov_mat)

# Extract leading eigenvectors
gar_wgp_vectors = gar_wgp_eigen$vectors[, 1:4, drop = FALSE]
gts_wgp_vectors = gts_wgp_eigen$vectors[, 1:4, drop = FALSE]

# Preserve trait identities
rownames(gar_wgp_vectors) = colnames(gar_wgp_cov_mat)
rownames(gts_wgp_vectors) = colnames(gts_wgp_cov_mat)

gts_gar_wgp_proc_res_eigen = procrustes(gar_wgp_vectors, 
                                        gts_wgp_vectors, 
                                        symmetric = T)

# gts_gar_wgp_proc_res_eigen = procrustes(gts_wgp_vectors, 
#                                        gar_wgp_vectors, 
#                                        symmetric = T)
summary(gts_gar_wgp_proc_res_eigen)
# plot(gts_gar_wgp_proc_res_eigen, 
#      kind = 1)


gts_gar_wgp_eigen_test = protest(gar_wgp_vectors, 
                                 gts_wgp_vectors, 
                                 permutations = 999, 
                                 scale = F, 
                                 symmetric = T)

# gts_gar_wgp_eigen_test = protest(gts_wgp_vectors, 
#                                 gar_wgp_vectors, 
#                                 permutations = 999, 
#                                 scale = F, 
#                                 symmetric = T)

gts_gar_wgp_eigen_residuals = residuals(gts_gar_wgp_eigen_test)

gts_gar_wgp_eigen_residuals = gts_gar_wgp_eigen_residuals%>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                eigen_residual = 2) 
# %>%  
#   write_csv('gts_gar_wgp_12_12_procrustes_eigen_cov_mat_residuals.csv')



gts_gar_wgp_residuals_combo = inner_join(gts_gar_wgp_residuals, 
                                        gts_gar_wgp_eigen_residuals, 
                                        by = 'trait') %>% 
  mutate(pop_pair = 'GTS/GAR')


wgp_proc_test_residual_combo = bind_rows(ashn_wgp_residuals_combo,
                                     myv_wgp_residuals_combo, 
                                     skr_wgp_residuals_combo, 
                                     gts_gar_wgp_residuals_combo)

##
# Trait cov matrices ------------------------------------------------------


wgp_cov_matrices <- list(
  Pair1_Ambient = ashnc_wgp_cov_mat,
  Pair1_Geothermal = ashnw_wgp_cov_mat,
  Pair2_Ambient = myvc_wgp_cov_mat,
  Pair2_Geothermal = myvw_wgp_cov_mat,
  Pair3_Ambient = skrc_wgp_cov_mat,
  Pair3_Geothermal = skrw_wgp_cov_mat,
  Pair4_Ambient = gar_wgp_cov_mat,
  Pair4_Geothermal = gts_wgp_cov_mat)

for (i in 1:4) {
  wgp_cov_matrices[[paste0("Pair", i, "_Difference")]] <-
    wgp_cov_matrices[[paste0("Pair", i, "_Geothermal")]] -
    wgp_cov_matrices[[paste0("Pair", i, "_Ambient")]]
}

wgp_cov_long <- purrr::imap_dfr(wgp_cov_matrices, function(mat, name) {
  as.data.frame(as.table(mat)) %>%
    setNames(c("Trait1", "Trait2", "Covariance")) %>%
    mutate(
      Temp = sub("_.*", "", name),
      Type = sub(".*_", "", name)
    )
})

wgp_cov_long$Type <- factor(
  wgp_cov_long$Type,
  levels = c("Ambient", "Geothermal", "Difference")
)

wgp_cov_long = wgp_cov_long %>% 
  mutate(pop_pairs = case_when(
    str_detect(Temp, 'Pair1') ~ 'ASHN',
    str_detect(Temp, 'Pair2') ~ 'MYV',
    str_detect(Temp, 'Pair3') ~ 'SKR',
    str_detect(Temp, 'Pair4') ~ 'GTS/GAR'))

wgp_cov_plot = ggplot(wgp_cov_long, 
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
       title = 'WGP trait covariance matrices')+
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1),
        panel.grid = element_blank(), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'), 
        axis.text = element_text(size = 6), 
        legend.position = 'none')

ggsave('WGP_Geo_Amb_Covariance_matrix.svg', 
       plot = wgp_cov_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 20)


# ProcTest residual plots -------------------------------------------------

resid_col_pal = c('#f72585', 
                  '#7209b7', 
                  '#480ca8', 
                  '#3f37c9')

wgp_cov_residual_plot = ggplot(wgp_proc_test_residual_combo, 
                           aes(x = factor(trait), 
                               y = cov_residual))+
  geom_col(aes(fill = pop_pair), 
           col = 'black')+
  ylim(0.0, 0.3)+
  facet_grid(~pop_pair)+
  labs(y = 'Covariance residuals')+
  scale_fill_manual(values = resid_col_pal)+
  # geom_hline(yintercept = 1.5, 
  #            linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold', 
                                  size = 14), 
        legend.position = 'none')

wgp_eigen_residual_plot = ggplot(wgp_proc_test_residual_combo, 
                             aes(x = factor(trait), 
                                 y = eigen_residual))+
  geom_col(aes(fill = pop_pair), 
           col = 'black')+
  ylim(0.0, 0.3)+
  facet_grid(~pop_pair)+
  labs(y = 'Eigenvector residuals')+
  scale_fill_manual(values = resid_col_pal)+
  # geom_hline(yintercept = 1.5, 
  #            linetype = 2)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold', 
                                  size = 14), 
        legend.position = 'none')


WGP_Trait_residual_plots = wgp_cov_residual_plot/wgp_eigen_residual_plot

ggsave('WGP_First_4_Eigenvectors_trait_residual_plots.svg', 
       plot = WGP_Trait_residual_plots, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 15)


##
# Kendalls concordance plot - Covariance matrices -------------------------

ashn_wgp_kendal = ashn_wgp_residuals %>% 
  dplyr::rename(ASHN = cov_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

myv_wgp_kendal = myv_wgp_residuals%>% 
  dplyr::rename(MYV = cov_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

skr_wgp_kendal = skr_wgp_residuals%>% 
  dplyr::rename(SKR = cov_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

gts_gar_wgp_kendal = gts_gar_wgp_residuals%>% 
  dplyr::rename(GTS_GAR = cov_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

wgp_kendal_df = bind_cols(ashn_wgp_kendal, 
                         myv_wgp_kendal, 
                         skr_wgp_kendal, 
                         gts_gar_wgp_kendal)

# group = c('ASHN','MYV','SKR','GTS_GAR')


wgp_kendal = kendall.global(wgp_kendal_df,
                           # group = pops,
                           nperm = 999, 
                           mult = 'hochberg')
wgp_kendal_post = kendall.post(wgp_kendal_df, 
                              mult = 'hochberg', 
                              nperm = 999)



# Kendall concordance - eigen vector traits -------------------------------
ashn_wgp_eigen_kendal = ashn_wgp_eigen_residuals %>% 
  dplyr::rename(ASHN = eigen_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

myv_wgp_eigen_kendal = myv_wgp_eigen_residuals%>% 
  dplyr::rename(MYV = eigen_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

skr_wgp_eigen_kendal = skr_wgp_eigen_residuals%>% 
  dplyr::rename(SKR = eigen_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

gts_gar_wgp_eigen_kendal = gts_gar_wgp_eigen_residuals%>% 
  dplyr::rename(GTS_GAR = eigen_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

wgp_eigen_kendal_df = bind_cols(ashn_wgp_eigen_kendal, 
                               myv_wgp_eigen_kendal, 
                               skr_wgp_eigen_kendal, 
                               gts_gar_wgp_eigen_kendal)

# group = c('ASHN','MYV','SKR','GTS_GAR')


wgp_eigen_kendal = kendall.global(wgp_eigen_kendal_df,
                                 # group = pops,
                                 nperm = 999, 
                                 mult = 'hochberg')
wgp_eigen_kendal_post = kendall.post(wgp_eigen_kendal_df, 
                                    mult = 'hochberg', 
                                    nperm = 999)

##
# Nemeny post hoc - trait covarianec --------------------------------------
library(PMCMRplus)


ashn_wgp_nemeny_df = ashn_wgp_residuals %>% 
  mutate(pop_pair = 'ASHN') 

myv_wgp_nemeny_df = myv_wgp_residuals %>% 
  mutate(pop_pair = 'MYV')

skr_wgp_nemeny_df = skr_wgp_residuals %>% 
  mutate(pop_pair = 'SKR')

gts_gar_wgp_nemeny_df = gts_gar_wgp_residuals %>% 
  mutate(pop_pair = 'GTS/GAR') 

wgp_nemeny_df = bind_rows(ashn_wgp_nemeny_df, 
                         myv_wgp_nemeny_df, 
                         skr_wgp_nemeny_df, 
                         gts_gar_wgp_nemeny_df)


wgp_nemeny_test = frdAllPairsNemenyiTest(y = wgp_nemeny_df$cov_residual, 
                                        groups = wgp_nemeny_df$trait, 
                                        blocks = wgp_nemeny_df$pop_pair, 
                                        p.adjust.method = 'hochberg')

wgp_nemeny_test$p.value %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname) %>% 
  # View() 
filter(if_any(c(2:33), function(x) x<0.05))
write_csv('wgp_1212_residual_NemenyTest_results.csv')



# nemeny post hoc - eigen vectors -----------------------------------------

library(PMCMRplus)


ashn_wgp_eigen_nemeny_df = ashn_wgp_eigen_residuals %>% 
  mutate(pop_pair = 'ASHN') 

myv_wgp_eigen_nemeny_df = myv_wgp_eigen_residuals %>% 
  mutate(pop_pair = 'MYV') 

skr_wgp_eigen_nemeny_df = skr_wgp_eigen_residuals %>% 
  mutate(pop_pair = 'SKR') 

gts_gar_wgp_eigen_nemeny_df = gts_gar_wgp_eigen_residuals %>% 
  mutate(pop_pair = 'GTS/GAR') 

wgp_eigen_nemeny_df = bind_rows(ashn_wgp_eigen_nemeny_df, 
                               myv_wgp_eigen_nemeny_df, 
                               skr_wgp_eigen_nemeny_df, 
                               gts_gar_wgp_eigen_nemeny_df)


wgp_eigen_nemeny_test = frdAllPairsNemenyiTest(y = wgp_eigen_nemeny_df$eigen_residual, 
                                              groups = wgp_eigen_nemeny_df$trait, 
                                              blocks = wgp_eigen_nemeny_df$pop_pair, 
                                              p.adjust.method = 'hochberg')



wgp_eigen_nemeny_test$p.value %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname) %>% 
  # View()
filter(if_any(c(2:33), function(x) x<0.05))
write_csv('wgp_1212_residual_NemenyTest_results.csv')

