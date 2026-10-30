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

tgp_orig_data = read_csv("TGP_TRAITS_SCALED_FIXED_11.02.2026.csv") %>% 
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


ashnc_tgp_trait_mat = tgp_orig_data %>% 
  filter(Lake_morph == 'ASHNC') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

ashnw_tgp_trait_mat = tgp_orig_data %>% 
  filter(Lake_morph == 'ASHNW') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

ashnc_tgp_cov_mat = cov(ashnc_tgp_trait_mat, 
                        use = 'complete.obs')
ashnw_tgp_cov_mat = cov(ashnw_tgp_trait_mat, 
                        use = 'complete.obs')
ashn_tgp_krzcov = evolqg::KrzCor(ashnc_tgp_cov_mat, 
                                 ashnw_tgp_cov_mat)
ashn_tgp_proc_res = procrustes(ashnc_tgp_cov_mat, 
                               ashnw_tgp_cov_mat)

# plot(ashn_tgp_proc_res, 
#      kind = 1)

ashn_tgp_proc_test = protest(ashnc_tgp_cov_mat, 
                             ashnw_tgp_cov_mat, 
                             permutations = 999)


ashn_tgp_residuals = residuals(ashn_tgp_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                cov_residual = 2)

# %>% 
#   mutate(resid_10 = residual*10)

ashn_tgp_resid_plot = ggplot(ashn_tgp_residuals, 
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

myvc_tgp_trait_mat = tgp_orig_data %>% 
  filter(Lake_morph == 'MYVC') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

myvw_tgp_trait_mat = tgp_orig_data %>% 
  filter(Lake_morph == 'MYVW') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

myvc_tgp_cov_mat = cov(myvc_tgp_trait_mat, 
                       use = 'complete.obs')
myvw_tgp_cov_mat = cov(myvw_tgp_trait_mat, 
                       use = 'complete.obs')
myv_tgp_krzcov = evolqg::KrzCor(myvc_tgp_cov_mat, 
                                myvw_tgp_cov_mat)
myv_tgp_proc_res = procrustes(myvc_tgp_cov_mat, 
                              myvw_tgp_cov_mat)

# plot(myv_tgp_proc_res, 
#      kind = 1)

myv_tgp_proc_test = protest(myvc_tgp_cov_mat, 
                            myvw_tgp_cov_mat, 
                            permutations = 999)


myv_tgp_residuals = residuals(myv_tgp_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                cov_residual = 2)

myv_tgp_resid_plot = ggplot(myv_tgp_residuals, 
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

skrc_tgp_trait_mat = tgp_orig_data %>% 
  filter(Lake_morph == 'SKRC') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

skrw_tgp_trait_mat = tgp_orig_data %>% 
  filter(Lake_morph == 'SKRW') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

skrc_tgp_cov_mat = cov(skrc_tgp_trait_mat, 
                       use = 'complete.obs')
skrw_tgp_cov_mat = cov(skrw_tgp_trait_mat, 
                       use = 'complete.obs')
skr_tgp_krzcov = evolqg::KrzCor(skrc_tgp_cov_mat, 
                                skrw_tgp_cov_mat)
skr_tgp_proc_res = procrustes(skrc_tgp_cov_mat, 
                              skrw_tgp_cov_mat)

# plot(skr_tgp_proc_res, 
#      kind = 1)

skr_tgp_proc_test = protest(skrc_tgp_cov_mat, 
                            skrw_tgp_cov_mat, 
                            permutations = 999)


skr_tgp_residuals = residuals(skr_tgp_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                cov_residual = 2)

skr_tgp_resid_plot = ggplot(skr_tgp_residuals, 
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

gar_tgp_trait_mat = tgp_orig_data %>% 
  filter(Lake_morph == 'CSWYC') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

gts_tgp_trait_mat = tgp_orig_data %>% 
  filter(Lake_morph == 'GTSW') %>% 
  dplyr::select(jaw_length:caudal2_15_17)

gar_tgp_cov_mat = cov(gar_tgp_trait_mat, 
                      use = 'complete.obs')
gts_tgp_cov_mat = cov(gts_tgp_trait_mat, 
                      use = 'complete.obs')
gts_gar_tgp_krzcov = evolqg::KrzCor(gar_tgp_cov_mat, 
                                    gts_tgp_cov_mat)
gts_gar_tgp_proc_res = procrustes(gar_tgp_cov_mat, 
                                  gts_tgp_cov_mat)

# plot(gts_gar_tgp_proc_res, 
#      kind = 1)

gts_gar_tgp_proc_test = protest(gar_tgp_cov_mat, 
                                gts_tgp_cov_mat, 
                                permutations = 999)


gts_gar_tgp_residuals = residuals(gts_gar_tgp_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                cov_residual = 2)

gts_gar_tgp_resid_plot = ggplot(gts_gar_tgp_residuals, 
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

ashnc_tgp_pca = prcomp(ashnc_tgp_cov_mat, 
                       center = F, 
                       scale. = F)
paran(ashnc_tgp_cov_mat, 
      iterations = 1000, 
      graph = T)
## 4 PC axes retained
ashnw_tgp_pca = prcomp(ashnw_tgp_cov_mat, 
                       center = F, 
                       scale. = F)
paran(ashnw_tgp_cov_mat, 
      iterations = 1000, 
      graph = T)
## 3 axes

myvc_tgp_pca = prcomp(myvc_tgp_cov_mat, 
                      center = F, 
                      scale. = F)
paran(myvc_tgp_cov_mat, 
      iterations = 1000, 
      graph = T)
## 3 PC axes retained
myvw_tgp_pca = prcomp(myvw_tgp_cov_mat, 
                      center = F, 
                      scale. = F)
paran(myvw_tgp_cov_mat, 
      iterations = 1000, 
      graph = T)
## 3 axes

skrc_tgp_pca = prcomp(skrc_tgp_cov_mat, 
                      center = F, 
                      scale. = F)
paran(skrc_tgp_cov_mat, 
      iterations = 1000, 
      graph = T)
## 3 PC axes retained
skrw_tgp_pca = prcomp(skrw_tgp_cov_mat, 
                      center = F, 
                      scale. = F)
paran(skrw_tgp_cov_mat, 
      iterations = 1000, 
      graph = T)
## 2 axes

gar_tgp_pca = prcomp(gar_tgp_cov_mat, 
                     center = F, 
                     scale. = F)
paran(gar_tgp_cov_mat, 
      iterations = 1000, 
      graph = T)
## 3 PC axes retained
gts_tgp_pca = prcomp(gts_tgp_cov_mat, 
                     center = F, 
                     scale. = F)
paran(gts_tgp_cov_mat, 
      iterations = 1000, 
      graph = T)
## 3 axes


# ASHN ---------------------------------------------------------------------


ashnc_tgp_eigen = eigen(ashnc_tgp_cov_mat)
ashnw_tgp_eigen = eigen(ashnw_tgp_cov_mat)

# Extract leading eigenvectors
ashnc_tgp_vectors = ashnc_tgp_eigen$vectors[, 1:4, drop = FALSE]
ashnw_tgp_vectors = ashnw_tgp_eigen$vectors[, 1:4, drop = FALSE]

# Preserve trait identities
rownames(ashnc_tgp_vectors) = colnames(ashnc_tgp_cov_mat)
rownames(ashnw_tgp_vectors) = colnames(ashnw_tgp_cov_mat)


ashn_tgp_proc_res_eigen = procrustes(ashnc_tgp_vectors, 
                                     ashnw_tgp_vectors, 
                                     symmetric = T)
summary(ashn_tgp_proc_res_eigen)
# plot(ashn_tgp_proc_res_eigen, 
#      kind = 1)

ashn_tgp_eigen_test = protest(ashnc_tgp_vectors, 
                              ashnw_tgp_vectors, 
                              permutations = 999, 
                              scale = F, 
                              symmetric = T)

ashn_tgp_eigen_residuals = residuals(ashn_tgp_eigen_test)

ashn_tgp_eigen_residuals = ashn_tgp_eigen_residuals%>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                eigen_residual = 2) 
# %>%  
#   write_csv('ASHN_tgp_12_12_procrustes_eigen_cov_mat_residuals.csv')


ashn_tgp_residuals_combo = inner_join(ashn_tgp_residuals,
                                      ashn_tgp_eigen_residuals,
                                      by = 'trait') %>%
  mutate(pop_pair = 'ASHN')

# %>%
#   write_csv('ASHN_tgp_procrustes_eigen_cov_mat_residuals.csv')


# MYV ---------------------------------------------------------------------

myvc_tgp_eigen = eigen(myvc_tgp_cov_mat)
myvw_tgp_eigen = eigen(myvw_tgp_cov_mat)

# Extract leading eigenvectors
myvc_tgp_vectors = myvc_tgp_eigen$vectors[, 1:4, drop = FALSE]
myvw_tgp_vectors = myvw_tgp_eigen$vectors[, 1:4, drop = FALSE]

# Preserve trait identities
rownames(myvc_tgp_vectors) = colnames(myvc_tgp_cov_mat)
rownames(myvw_tgp_vectors) = colnames(myvw_tgp_cov_mat)


myv_tgp_proc_res_eigen = procrustes(myvc_tgp_vectors, 
                                    myvw_tgp_vectors, 
                                    symmetric = T)
summary(myv_tgp_proc_res_eigen)
# plot(myv_tgp_proc_res_eigen, 
#      kind = 1)

myv_tgp_eigen_test = protest(myvc_tgp_vectors, 
                             myvw_tgp_vectors, 
                             permutations = 999, 
                             scale = F, 
                             symmetric = T)

myv_tgp_eigen_residuals = residuals(myv_tgp_eigen_test)

myv_tgp_eigen_residuals = myv_tgp_eigen_residuals%>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                eigen_residual = 2) 
# %>%  
#   write_csv('myv_tgp_12_12_procrustes_eigen_cov_mat_residuals.csv')



myv_tgp_residuals_combo = inner_join(myv_tgp_residuals, 
                                     myv_tgp_eigen_residuals, 
                                     by = 'trait') %>% 
  mutate(pop_pair = 'MYV')


# SKR ---------------------------------------------------------------------

skrc_tgp_eigen = eigen(skrc_tgp_cov_mat)
skrw_tgp_eigen = eigen(skrw_tgp_cov_mat)

# Extract leading eigenvectors
skrc_tgp_vectors = skrc_tgp_eigen$vectors[, 1:4, drop = FALSE]
skrw_tgp_vectors = skrw_tgp_eigen$vectors[, 1:4, drop = FALSE]

# Preserve trait identities
rownames(skrc_tgp_vectors) = colnames(skrc_tgp_cov_mat)
rownames(skrw_tgp_vectors) = colnames(skrw_tgp_cov_mat)


skr_tgp_proc_res_eigen = procrustes(skrc_tgp_vectors, 
                                    skrw_tgp_vectors, 
                                    symmetric = T)
summary(skr_tgp_proc_res_eigen)
# plot(skr_tgp_proc_res_eigen, 
#      kind = 1)

skr_tgp_eigen_test = protest(skrc_tgp_vectors, 
                             skrw_tgp_vectors, 
                             permutations = 999, 
                             scale = F, 
                             symmetric = T)

skr_tgp_eigen_residuals = residuals(skr_tgp_eigen_test)

skr_tgp_eigen_residuals = skr_tgp_eigen_residuals%>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                eigen_residual = 2) 
# %>%  
#   write_csv('skr_tgp_12_12_procrustes_eigen_cov_mat_residuals.csv')



skr_tgp_residuals_combo = inner_join(skr_tgp_residuals, 
                                     skr_tgp_eigen_residuals, 
                                     by = 'trait') %>% 
  mutate(pop_pair = 'SKR')

# GTS/GAR -----------------------------------------------------------------

gar_tgp_eigen = eigen(gar_tgp_cov_mat)
gts_tgp_eigen = eigen(gts_tgp_cov_mat)

# Extract leading eigenvectors
gar_tgp_vectors = gar_tgp_eigen$vectors[, 1:4, drop = FALSE]
gts_tgp_vectors = gts_tgp_eigen$vectors[, 1:4, drop = FALSE]

# Preserve trait identities
rownames(gar_tgp_vectors) = colnames(gar_tgp_cov_mat)
rownames(gts_tgp_vectors) = colnames(gts_tgp_cov_mat)


gts_gar_tgp_proc_res_eigen = procrustes(gar_tgp_vectors, 
                                        gts_tgp_vectors, 
                                        symmetric = T)
summary(gts_gar_tgp_proc_res_eigen)
# plot(gts_gar_tgp_proc_res_eigen, 
#      kind = 1)

gts_gar_tgp_eigen_test = protest(gar_tgp_vectors, 
                                 gts_tgp_vectors, 
                                 permutations = 999, 
                                 scale = F, 
                                 symmetric = T)

gts_gar_tgp_eigen_residuals = residuals(gts_gar_tgp_eigen_test)

gts_gar_tgp_eigen_residuals = gts_gar_tgp_eigen_residuals%>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                eigen_residual = 2) 
# %>%  
#   write_csv('gts_gar_tgp_12_12_procrustes_eigen_cov_mat_residuals.csv')



gts_gar_tgp_residuals_combo = inner_join(gts_gar_tgp_residuals, 
                                         gts_gar_tgp_eigen_residuals, 
                                         by = 'trait') %>% 
  mutate(pop_pair = 'GTS/GAR')


tgp_proc_test_residual_combo = bind_rows(ashn_tgp_residuals_combo,
                                         myv_tgp_residuals_combo, 
                                         skr_tgp_residuals_combo, 
                                         gts_gar_tgp_residuals_combo)

##
# Trait cov matrices ------------------------------------------------------


tgp_cov_matrices <- list(
  Pair1_Ambient = ashnc_tgp_cov_mat,
  Pair1_Geothermal = ashnw_tgp_cov_mat,
  Pair2_Ambient = myvc_tgp_cov_mat,
  Pair2_Geothermal = myvw_tgp_cov_mat,
  Pair3_Ambient = skrc_tgp_cov_mat,
  Pair3_Geothermal = skrw_tgp_cov_mat,
  Pair4_Ambient = gar_tgp_cov_mat,
  Pair4_Geothermal = gts_tgp_cov_mat)

for (i in 1:4) {
  tgp_cov_matrices[[paste0("Pair", i, "_Difference")]] <-
    tgp_cov_matrices[[paste0("Pair", i, "_Geothermal")]] -
    tgp_cov_matrices[[paste0("Pair", i, "_Ambient")]]
}

tgp_cov_long <- purrr::imap_dfr(tgp_cov_matrices, function(mat, name) {
  as.data.frame(as.table(mat)) %>%
    setNames(c("Trait1", "Trait2", "Covariance")) %>%
    mutate(
      Temp = sub("_.*", "", name),
      Type = sub(".*_", "", name)
    )
})

tgp_cov_long$Type <- factor(
  tgp_cov_long$Type,
  levels = c("Ambient", "Geothermal", "Difference")
)

tgp_cov_long = tgp_cov_long %>% 
  mutate(pop_pairs = case_when(
    str_detect(Temp, 'Pair1') ~ 'ASHN',
    str_detect(Temp, 'Pair2') ~ 'MYV',
    str_detect(Temp, 'Pair3') ~ 'SKR',
    str_detect(Temp, 'Pair4') ~ 'GTS/GAR'))

tgp_cov_plot = ggplot(tgp_cov_long, 
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
       title = 'tgp trait covariance matrices')+
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1),
        panel.grid = element_blank(), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'), 
        axis.text = element_text(size = 6), 
        legend.position = 'none')

ggsave('tgp_Geo_Amb_Covariance_matrix.svg', 
       plot = tgp_cov_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 20)


# ProcTest residual plots -------------------------------------------------

resid_col_pal = c('#f72585', 
                  '#7209b7', 
                  '#480ca8', 
                  '#3f37c9')

tgp_cov_residual_plot = ggplot(tgp_proc_test_residual_combo, 
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

tgp_eigen_residual_plot = ggplot(tgp_proc_test_residual_combo, 
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


tgp_Trait_residual_plots = tgp_cov_residual_plot/tgp_eigen_residual_plot

ggsave('tgp_First_4_Eigenvectors_trait_residual_plots.svg', 
       plot = tgp_Trait_residual_plots, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 15)


##
# Kendalls concordance plot - Covariance matrices -------------------------

ashn_tgp_kendal = ashn_tgp_residuals %>% 
  dplyr::rename(ASHN = cov_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

myv_tgp_kendal = myv_tgp_residuals%>% 
  dplyr::rename(MYV = cov_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

skr_tgp_kendal = skr_tgp_residuals%>% 
  dplyr::rename(SKR = cov_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

gts_gar_tgp_kendal = gts_gar_tgp_residuals%>% 
  dplyr::rename(GTS_GAR = cov_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

tgp_kendal_df = bind_cols(ashn_tgp_kendal, 
                          myv_tgp_kendal, 
                          skr_tgp_kendal, 
                          gts_gar_tgp_kendal)

# group = c('ASHN','MYV','SKR','GTS_GAR')


tgp_kendal = kendall.global(tgp_kendal_df,
                            # group = pops,
                            nperm = 999, 
                            mult = 'hochberg')
tgp_kendal_post = kendall.post(tgp_kendal_df, 
                               mult = 'hochberg', 
                               nperm = 999)



# Kendall concordance - eigen vector traits -------------------------------
ashn_tgp_eigen_kendal = ashn_tgp_eigen_residuals %>% 
  dplyr::rename(ASHN = eigen_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

myv_tgp_eigen_kendal = myv_tgp_eigen_residuals%>% 
  dplyr::rename(MYV = eigen_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

skr_tgp_eigen_kendal = skr_tgp_eigen_residuals%>% 
  dplyr::rename(SKR = eigen_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

gts_gar_tgp_eigen_kendal = gts_gar_tgp_eigen_residuals%>% 
  dplyr::rename(GTS_GAR = eigen_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

tgp_eigen_kendal_df = bind_cols(ashn_tgp_eigen_kendal, 
                                myv_tgp_eigen_kendal, 
                                skr_tgp_eigen_kendal, 
                                gts_gar_tgp_eigen_kendal)

# group = c('ASHN','MYV','SKR','GTS_GAR')


tgp_eigen_kendal = kendall.global(tgp_eigen_kendal_df,
                                  # group = pops,
                                  nperm = 999, 
                                  mult = 'hochberg')
tgp_eigen_kendal_post = kendall.post(tgp_eigen_kendal_df, 
                                     mult = 'hochberg', 
                                     nperm = 999)

##
# Nemeny post hoc - trait covarianec --------------------------------------
library(PMCMRplus)


ashn_tgp_nemeny_df = ashn_tgp_residuals %>% 
  mutate(pop_pair = 'ASHN') 

myv_tgp_nemeny_df = myv_tgp_residuals %>% 
  mutate(pop_pair = 'MYV')

skr_tgp_nemeny_df = skr_tgp_residuals %>% 
  mutate(pop_pair = 'SKR')

gts_gar_tgp_nemeny_df = gts_gar_tgp_residuals %>% 
  mutate(pop_pair = 'GTS/GAR') 

tgp_nemeny_df = bind_rows(ashn_tgp_nemeny_df, 
                          myv_tgp_nemeny_df, 
                          skr_tgp_nemeny_df, 
                          gts_gar_tgp_nemeny_df)


tgp_nemeny_test = frdAllPairsNemenyiTest(y = tgp_nemeny_df$cov_residual, 
                                         groups = tgp_nemeny_df$trait, 
                                         blocks = tgp_nemeny_df$pop_pair, 
                                         p.adjust.method = 'hochberg')

tgp_nemeny_test$p.value %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname) %>% 
  # View()
  filter(if_any(c(2:33), function(x) x<0.05))
write_csv('tgp_1212_residual_NemenyTest_results.csv')



# nemeny post hoc - eigen vectors -----------------------------------------

library(PMCMRplus)


ashn_tgp_eigen_nemeny_df = ashn_tgp_eigen_residuals %>% 
  mutate(pop_pair = 'ASHN') 

myv_tgp_eigen_nemeny_df = myv_tgp_eigen_residuals %>% 
  mutate(pop_pair = 'MYV') 

skr_tgp_eigen_nemeny_df = skr_tgp_eigen_residuals %>% 
  mutate(pop_pair = 'SKR') 

gts_gar_tgp_eigen_nemeny_df = gts_gar_tgp_eigen_residuals %>% 
  mutate(pop_pair = 'GTS/GAR') 

tgp_eigen_nemeny_df = bind_rows(ashn_tgp_eigen_nemeny_df, 
                                myv_tgp_eigen_nemeny_df, 
                                skr_tgp_eigen_nemeny_df, 
                                gts_gar_tgp_eigen_nemeny_df)


tgp_eigen_nemeny_test = frdAllPairsNemenyiTest(y = tgp_eigen_nemeny_df$eigen_residual, 
                                               groups = tgp_eigen_nemeny_df$trait, 
                                               blocks = tgp_eigen_nemeny_df$pop_pair, 
                                               p.adjust.method = 'hochberg')



tgp_eigen_nemeny_test$p.value %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname) %>% 
  # View()
  filter(if_any(c(2:33), function(x) x<0.05))  %>% View()
write_csv('tgp_1212_residual_NemenyTest_results.csv')

tgp_significant_traits =  c("lm_19_18",
                           "fbar_23_24", 
                           'fbar_8_27', 
                           'head_depth', 
                           'lm_6_12') # Replace with actual Nemenyi-identified traits


tgp_eigen_nemeny_df %>%
  filter(trait %in% tgp_significant_traits) %>%
  ggplot(aes(x = trait,
             y = eigen_residual,
             group = pop_pair,
             colour = pop_pair)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 3) +
  labs(x = "Morphological trait",
       y = "Procrustes residual distance",
       colour = "Population pair")


## Zoom in on trait covariance structure. 
tgp_ashn_cov_differences = ashnw_tgp_cov_mat - ashnc_tgp_cov_mat 
tgp_ashn_cov_diff = tgp_ashn_cov_differences%>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname) %>% 
  dplyr::filter(trait %in% c(tgp_significant_traits)) %>% 
  dplyr::select(trait, 
                lm_19_18, 
                fbar_23_24, 
                fbar_8_27, 
                head_depth, 
                lm_6_12) %>% 
  mutate(pop_pair = 'ASHN')

tgp_myv_cov_differences = myvw_tgp_cov_mat - myvc_tgp_cov_mat 
tgp_myv_cov_diff = tgp_myv_cov_differences%>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname)%>% 
  dplyr::filter(trait %in% c(tgp_significant_traits)) %>% 
  dplyr::select(trait, 
                lm_19_18, 
                fbar_23_24, 
                fbar_8_27, 
                head_depth, 
                lm_6_12) %>% 
  mutate(pop_pair = 'MYV')

tgp_skr_cov_differences = skrw_tgp_cov_mat - skrc_tgp_cov_mat 
tgp_skr_cov_diff = tgp_skr_cov_differences%>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname)%>% 
  dplyr::filter(trait %in% c(tgp_significant_traits)) %>% 
  dplyr::select(trait, 
                lm_19_18, 
                fbar_23_24, 
                fbar_8_27, 
                head_depth, 
                lm_6_12) %>% 
  mutate(pop_pair = 'SKR')

tgp_gts_gar_cov_differences = gts_tgp_cov_mat - gar_tgp_cov_mat 
tgp_gts_gar_cov_diff = tgp_gts_gar_cov_differences%>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname)%>% 
  dplyr::filter(trait %in% c(tgp_significant_traits)) %>% 
  dplyr::select(trait, 
                lm_19_18, 
                fbar_23_24, 
                fbar_8_27, 
                head_depth, 
                lm_6_12) %>% 
  mutate(pop_pair = 'GTS/GAR')


tgp_sig_cov_mat = bind_rows(tgp_ashn_cov_diff, 
                            tgp_myv_cov_diff, 
                            tgp_skr_cov_diff, 
                            tgp_gts_gar_cov_diff)

tgp_sig_cov_long <- tgp_sig_cov_mat %>%
  pivot_longer(
    cols = all_of(tgp_significant_traits),
    names_to = "Trait2",
    values_to = "Covariance") %>%
  dplyr::rename(Trait1 = trait) %>%
  mutate(Trait1 = factor(Trait1, levels = tgp_significant_traits),
         Trait2 = factor(Trait2, levels = tgp_significant_traits),
         pop_pair = factor(
           pop_pair,
           levels = c("ASHN", "MYV", "SKR", "GTS/GAR")))


tgp_sig_cov_triangle = tgp_sig_cov_long %>%
  mutate(Trait1_num = as.integer(Trait1),
         Trait2_num = as.integer(Trait2)) %>%
  filter(Trait1_num > Trait2_num)

tgp_sig_covariance_traits = ggplot(tgp_sig_cov_triangle,
                                  aes(x = Trait2, y = Trait1, fill = Covariance)) +
  geom_tile(colour = "black", 
            linewidth = 0.5) +
  facet_wrap(~ pop_pair, nrow = 2) +
  scale_fill_gradient2(
    low = "#0f4c5c",
    mid = "white",
    high = "#9a031e",
    midpoint = 0,
    limits = c(-1.05, 1.05)) +
  coord_equal() +
  labs(x = NULL,
       y = NULL,
       fill = expression(Delta~'Covariance'), 
       title = 'C) TGP')+
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold"), 
        strip.background = element_rect(fill = 'white'), 
        strip.text.x = element_text(face = 'bold', 
                                    size = 14), 
        plot.title = element_text(size = 16),
        legend.position = 'none')



sig_cov_plots = f2_sig_covariance_traits+wgp_sig_covariance_traits+tgp_sig_covariance_traits

ggsave('Covariance_plots_Significant_Nemenyi_test.svg', 
       plot = sig_cov_plots, 
       dpi = 'retina', 
       units = 'cm', 
       width = 25, 
       height = 15)
