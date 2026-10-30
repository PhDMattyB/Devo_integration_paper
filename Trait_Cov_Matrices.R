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



wild_orig_data = read_csv("WILD_SCALED_FIXED_27.05.2026.csv") %>% 
  dplyr::select(Lake_morph, 
                rowname, 
                Lake, 
                Morph, 
                5:39)


# Streamlining the analysis to match F2 generation ------------------------

# Procrustes superimposition - cov matrices -------------------------------
# ASHN --------------------------------------------------------------------


ashnc_wild_trait_mat = wild_orig_data %>% 
  filter(Lake_morph == 'ASHNC') %>% 
  dplyr::select(jaw_length:39)

ashnw_wild_trait_mat = wild_orig_data %>% 
  filter(Lake_morph == 'ASHNW') %>% 
  dplyr::select(jaw_length:39)

ashnc_wild_cov_mat = cov(ashnc_wild_trait_mat, 
                       use = 'complete.obs')
ashnw_wild_cov_mat = cov(ashnw_wild_trait_mat, 
                       use = 'complete.obs')
ashn_wild_krzcov = evolqg::KrzCor(ashnc_wild_cov_mat, 
                                ashnw_wild_cov_mat)
ashn_wild_proc_res = procrustes(ashnc_wild_cov_mat, 
                              ashnw_wild_cov_mat)

plot(ashn_wild_proc_res, 
     kind = 1)

ashn_wild_proc_test = protest(ashnc_wild_cov_mat, 
                            ashnw_wild_cov_mat, 
                            permutations = 999)


ashn_wild_residuals = residuals(ashn_wild_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                cov_residual = 2)

# %>% 
#   mutate(resid_10 = residual*10)

ashn_wild_resid_plot = ggplot(ashn_wild_residuals, 
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

myvc_wild_trait_mat = wild_orig_data %>% 
  filter(Lake_morph == 'MYVC') %>% 
  dplyr::select(jaw_length:39)

myvw_wild_trait_mat = wild_orig_data %>% 
  filter(Lake_morph == 'MYVW') %>% 
  dplyr::select(jaw_length:39)

myvc_wild_cov_mat = cov(myvc_wild_trait_mat, 
                      use = 'complete.obs')
myvw_wild_cov_mat = cov(myvw_wild_trait_mat, 
                      use = 'complete.obs')
myv_wild_krzcov = evolqg::KrzCor(myvc_wild_cov_mat, 
                               myvw_wild_cov_mat)
myv_wild_proc_res = procrustes(myvc_wild_cov_mat, 
                             myvw_wild_cov_mat)

plot(myv_wild_proc_res, 
     kind = 1)

myv_wild_proc_test = protest(myvc_wild_cov_mat, 
                           myvw_wild_cov_mat, 
                           permutations = 999)


myv_wild_residuals = residuals(myv_wild_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                cov_residual = 2)

myv_wild_resid_plot = ggplot(myv_wild_residuals, 
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

skrc_wild_trait_mat = wild_orig_data %>% 
  filter(Lake_morph == 'SKRC') %>% 
  dplyr::select(jaw_length:39)

skrw_wild_trait_mat = wild_orig_data %>% 
  filter(Lake_morph == 'SKRW') %>% 
  dplyr::select(jaw_length:39)

skrc_wild_cov_mat = cov(skrc_wild_trait_mat, 
                      use = 'complete.obs')
skrw_wild_cov_mat = cov(skrw_wild_trait_mat, 
                      use = 'complete.obs')
skr_wild_krzcov = evolqg::KrzCor(skrc_wild_cov_mat, 
                               skrw_wild_cov_mat)
skr_wild_proc_res = procrustes(skrc_wild_cov_mat, 
                             skrw_wild_cov_mat)

plot(skr_wild_proc_res, 
     kind = 1)

skr_wild_proc_test = protest(skrc_wild_cov_mat, 
                           skrw_wild_cov_mat, 
                           permutations = 999)


skr_wild_residuals = residuals(skr_wild_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                cov_residual = 2)

skr_wild_resid_plot = ggplot(skr_wild_residuals, 
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

gar_wild_trait_mat = wild_orig_data %>% 
  filter(Lake_morph == 'CSWY') %>% 
  dplyr::select(jaw_length:39)

gts_wild_trait_mat = wild_orig_data %>% 
  filter(Lake_morph == 'GTS') %>% 
  dplyr::select(jaw_length:39)

gar_wild_cov_mat = cov(gar_wild_trait_mat, 
                     use = 'complete.obs')
gts_wild_cov_mat = cov(gts_wild_trait_mat, 
                     use = 'complete.obs')
gts_gar_wild_krzcov = evolqg::KrzCor(gar_wild_cov_mat, 
                                   gts_wild_cov_mat)
gts_gar_wild_proc_res = procrustes(gar_wild_cov_mat, 
                                 gts_wild_cov_mat)

plot(gts_gar_wild_proc_res, 
     kind = 1)

gts_gar_wild_proc_test = protest(gar_wild_cov_mat, 
                               gts_wild_cov_mat, 
                               permutations = 999)


gts_gar_wild_residuals = residuals(gts_gar_wild_proc_test) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                cov_residual = 2)

gts_gar_wild_resid_plot = ggplot(gts_gar_wild_residuals, 
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

ashnc_wild_pca = prcomp(ashnc_wild_cov_mat, 
                      center = F, 
                      scale. = F)
paran(ashnc_wild_cov_mat, 
      iterations = 1000, 
      graph = T)
## 3 PC axes retained
ashnw_wild_pca = prcomp(ashnw_wild_cov_mat, 
                      center = F, 
                      scale. = F)
paran(ashnw_wild_cov_mat, 
      iterations = 1000, 
      graph = T)
## 2 axes

myvc_wild_pca = prcomp(myvc_wild_cov_mat, 
                     center = F, 
                     scale. = F)
paran(myvc_wild_cov_mat, 
      iterations = 1000, 
      graph = T)
## 4 PC axes retained
myvw_wild_pca = prcomp(myvw_wild_cov_mat, 
                     center = F, 
                     scale. = F)
paran(myvw_wild_cov_mat, 
      iterations = 1000, 
      graph = T)
## 3 axes

skrc_wild_pca = prcomp(skrc_wild_cov_mat, 
                     center = F, 
                     scale. = F)
paran(skrc_wild_cov_mat, 
      iterations = 1000, 
      graph = T)
## 3 PC axes retained
skrw_wild_pca = prcomp(skrw_wild_cov_mat, 
                     center = F, 
                     scale. = F)
paran(skrw_wild_cov_mat, 
      iterations = 1000, 
      graph = T)
## 2 axes

gar_wild_pca = prcomp(gar_wild_cov_mat, 
                    center = F, 
                    scale. = F)
paran(gar_wild_cov_mat, 
      iterations = 1000, 
      graph = T)
## 3 PC axes retained
gts_wild_pca = prcomp(gts_wild_cov_mat, 
                    center = F, 
                    scale. = F)
paran(gts_wild_cov_mat, 
      iterations = 1000, 
      graph = T)
## 4 axes

variance_explained <- function(cov_mat) {
  
  eig <- eigen(cov_mat)$values
  
  cumulative <- cumsum(eig) / sum(eig)
  
  data.frame(
    PC = seq_along(eig),
    Variance = eig / sum(eig),
    Cumulative = cumulative
  )
}

variance_explained(ashnc_wild_cov_mat)[1:6, ]
variance_explained(ashnw_wild_cov_mat)[1:6, ]

variance_explained(myvc_wild_cov_mat)[1:6, ]
variance_explained(myvw_wild_cov_mat)[1:6, ]

variance_explained(skrc_wild_cov_mat)[1:6, ]
variance_explained(skrw_wild_cov_mat)[1:6, ]

variance_explained(gar_wild_cov_mat)[1:6, ]
variance_explained(gts_wild_cov_mat)[1:6, ]

# ASHN ---------------------------------------------------------------------


ashnc_wild_eigen = eigen(ashnc_wild_cov_mat)
ashnw_wild_eigen = eigen(ashnw_wild_cov_mat)

# Extract leading eigenvectors
ashnc_wild_vectors = ashnc_wild_eigen$vectors[, 1:4, drop = FALSE]
ashnw_wild_vectors = ashnw_wild_eigen$vectors[, 1:4, drop = FALSE]

# Preserve trait identities
rownames(ashnc_wild_vectors) = colnames(ashnc_wild_cov_mat)
rownames(ashnw_wild_vectors) = colnames(ashnw_wild_cov_mat)


ashn_wild_proc_res_eigen = procrustes(ashnc_wild_vectors, 
                                    ashnw_wild_vectors, 
                                    symmetric = T)
summary(ashn_wild_proc_res_eigen)
# plot(ashn_wild_proc_res_eigen, 
#      kind = 1)

ashn_wild_eigen_test = protest(ashnc_wild_vectors, 
                             ashnw_wild_vectors, 
                             permutations = 999, 
                             scale = F, 
                             symmetric = T)

ashn_wild_eigen_residuals = residuals(ashn_wild_eigen_test)

ashn_wild_eigen_residuals = ashn_wild_eigen_residuals%>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                eigen_residual = 2) 
# %>%  
#   write_csv('ASHN_wild_12_12_procrustes_eigen_cov_mat_residuals.csv')


ashn_wild_residuals_combo = inner_join(ashn_wild_residuals, 
                                     ashn_wild_eigen_residuals, 
                                     by = 'trait') %>% 
  mutate(pop_pair = 'ASHN')

# %>%
#   write_csv('ASHN_wild_procrustes_eigen_cov_mat_residuals.csv')


# MYV ---------------------------------------------------------------------

myvc_wild_eigen = eigen(myvc_wild_cov_mat)
myvw_wild_eigen = eigen(myvw_wild_cov_mat)

# Extract leading eigenvectors
myvc_wild_vectors = myvc_wild_eigen$vectors[, 1:4, drop = FALSE]
myvw_wild_vectors = myvw_wild_eigen$vectors[, 1:4, drop = FALSE]

# Preserve trait identities
rownames(myvc_wild_vectors) = colnames(myvc_wild_cov_mat)
rownames(myvw_wild_vectors) = colnames(myvw_wild_cov_mat)


myv_wild_proc_res_eigen = procrustes(myvc_wild_vectors, 
                                   myvw_wild_vectors, 
                                   symmetric = T)
summary(myv_wild_proc_res_eigen)
# plot(myv_wild_proc_res_eigen, 
#      kind = 1)

myv_wild_eigen_test = protest(myvc_wild_vectors, 
                            myvw_wild_vectors, 
                            permutations = 999, 
                            scale = F, 
                            symmetric = T)

myv_wild_eigen_residuals = residuals(myv_wild_eigen_test)

myv_wild_eigen_residuals = myv_wild_eigen_residuals%>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                eigen_residual = 2) 
# %>%  
#   write_csv('myv_wild_12_12_procrustes_eigen_cov_mat_residuals.csv')



myv_wild_residuals_combo = inner_join(myv_wild_residuals, 
                                    myv_wild_eigen_residuals, 
                                    by = 'trait') %>% 
  mutate(pop_pair = 'MYV')


# SKR ---------------------------------------------------------------------

skrc_wild_eigen = eigen(skrc_wild_cov_mat)
skrw_wild_eigen = eigen(skrw_wild_cov_mat)

# Extract leading eigenvectors
skrc_wild_vectors = skrc_wild_eigen$vectors[, 1:4, drop = FALSE]
skrw_wild_vectors = skrw_wild_eigen$vectors[, 1:4, drop = FALSE]

# Preserve trait identities
rownames(skrc_wild_vectors) = colnames(skrc_wild_cov_mat)
rownames(skrw_wild_vectors) = colnames(skrw_wild_cov_mat)


skr_wild_proc_res_eigen = procrustes(skrc_wild_vectors, 
                                   skrw_wild_vectors, 
                                   symmetric = T)
summary(skr_wild_proc_res_eigen)
# plot(skr_wild_proc_res_eigen, 
#      kind = 1)

skr_wild_eigen_test = protest(skrc_wild_vectors, 
                            skrw_wild_vectors, 
                            permutations = 999, 
                            scale = F, 
                            symmetric = T)

skr_wild_eigen_residuals = residuals(skr_wild_eigen_test)

skr_wild_eigen_residuals = skr_wild_eigen_residuals%>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                eigen_residual = 2) 
# %>%  
#   write_csv('skr_wild_12_12_procrustes_eigen_cov_mat_residuals.csv')



skr_wild_residuals_combo = inner_join(skr_wild_residuals, 
                                    skr_wild_eigen_residuals, 
                                    by = 'trait') %>% 
  mutate(pop_pair = 'SKR')


# GTS/GAR -----------------------------------------------------------------

gar_wild_eigen = eigen(gar_wild_cov_mat)
gts_wild_eigen = eigen(gts_wild_cov_mat)

# Extract leading eigenvectors
gar_wild_vectors = gar_wild_eigen$vectors[, 1:4, drop = FALSE]
gts_wild_vectors = gts_wild_eigen$vectors[, 1:4, drop = FALSE]

# Preserve trait identities
rownames(gar_wild_vectors) = colnames(gar_wild_cov_mat)
rownames(gts_wild_vectors) = colnames(gts_wild_cov_mat)


gts_gar_wild_proc_res_eigen = procrustes(gar_wild_vectors, 
                                       gts_wild_vectors, 
                                       symmetric = T)
summary(gts_gar_wild_proc_res_eigen)
# plot(gts_gar_wild_proc_res_eigen, 
#      kind = 1)

gts_gar_wild_eigen_test = protest(gar_wild_vectors, 
                                gts_wild_vectors, 
                                permutations = 999, 
                                scale = F, 
                                symmetric = T)

gts_gar_wild_eigen_residuals = residuals(gts_gar_wild_eigen_test)

gts_gar_wild_eigen_residuals = gts_gar_wild_eigen_residuals%>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                eigen_residual = 2) 
# %>%  
#   write_csv('gts_gar_wild_12_12_procrustes_eigen_cov_mat_residuals.csv')



gts_gar_wild_residuals_combo = inner_join(gts_gar_wild_residuals, 
                                        gts_gar_wild_eigen_residuals, 
                                        by = 'trait') %>% 
  mutate(pop_pair = 'GTS/GAR')


wild_Proc_test_residual_combo = bind_rows(ashn_wild_residuals_combo,
                                     myv_wild_residuals_combo, 
                                     skr_wild_residuals_combo, 
                                     gts_gar_wild_residuals_combo)


# ProcTest residual plots -------------------------------------------------

resid_col_pal = c('#f72585', 
                  '#7209b7', 
                  '#480ca8', 
                  '#3f37c9')

wild_cov_residual_plot = ggplot(wild_Proc_test_residual_combo, 
                           aes(x = factor(trait), 
                               y = cov_residual))+
  geom_col(aes(fill = pop_pair), 
           col = 'black')+
  ylim(0.0, 0.5)+
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

wild_eigen_residual_plot = ggplot(wild_Proc_test_residual_combo, 
                             aes(x = factor(trait), 
                                 y = eigen_residual))+
  geom_col(aes(fill = pop_pair), 
           col = 'black')+
  ylim(0.0, 0.5)+
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


wild_Trait_residual_plots = wild_cov_residual_plot/wild_eigen_residual_plot

ggsave('wild_First_4Eigenvectors_trait_residual_plots.svg', 
       plot = wild_Trait_residual_plots, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 15)

# Trait cov matrices ------------------------------------------------------


wild_cov_matrices <- list(
  Pair1_Ambient = ashnc_wild_cov_mat,
  Pair1_Geothermal = ashnw_wild_cov_mat,
  Pair2_Ambient = myvc_wild_cov_mat,
  Pair2_Geothermal = myvw_wild_cov_mat,
  Pair3_Ambient = skrc_wild_cov_mat,
  Pair3_Geothermal = skrw_wild_cov_mat,
  Pair4_Ambient = gar_wild_cov_mat,
  Pair4_Geothermal = gts_wild_cov_mat)

for (i in 1:4) {
  wild_cov_matrices[[paste0("Pair", i, "_Difference")]] <-
    wild_cov_matrices[[paste0("Pair", i, "_Geothermal")]] -
    wild_cov_matrices[[paste0("Pair", i, "_Ambient")]]
}

wild_cov_long <- purrr::imap_dfr(wild_cov_matrices, function(mat, name) {
  as.data.frame(as.table(mat)) %>%
    setNames(c("Trait1", "Trait2", "Covariance")) %>%
    mutate(
      Temp = sub("_.*", "", name),
      Type = sub(".*_", "", name)
    )
})

wild_cov_long$Type <- factor(
  wild_cov_long$Type,
  levels = c("Ambient", "Geothermal", "Difference")
)

wild_cov_long = wild_cov_long %>% 
  mutate(pop_pairs = case_when(
    str_detect(Temp, 'Pair1') ~ 'ASHN',
    str_detect(Temp, 'Pair2') ~ 'MYV',
    str_detect(Temp, 'Pair3') ~ 'SKR',
    str_detect(Temp, 'Pair4') ~ 'GTS/GAR'))

wild_cov_plot = ggplot(wild_cov_long, 
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
       title = 'Wild trait covariance matrices')+
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1),
        panel.grid = element_blank(), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'), 
        axis.text = element_text(size = 6), 
        legend.position = 'none')

ggsave('wild_Geo_Amb_Covariance_matrix.svg', 
       plot = wild_cov_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 20)


 ##
# Kendalls concordance plot - Covariance matrices -------------------------

ashn_wild_kendal = ashn_wild_residuals %>% 
  dplyr::rename(ASHN = cov_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

myv_wild_kendal = myv_wild_residuals%>% 
  dplyr::rename(MYV = cov_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

skr_wild_kendal = skr_wild_residuals%>% 
  dplyr::rename(SKR = cov_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

gts_gar_wild_kendal = gts_gar_wild_residuals%>% 
  dplyr::rename(GTS_GAR = cov_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

wild_kendal_df = bind_cols(ashn_wild_kendal, 
                         myv_wild_kendal, 
                         skr_wild_kendal, 
                         gts_gar_wild_kendal)

# group = c('ASHN','MYV','SKR','GTS_GAR')


wild_kendal = kendall.global(wild_kendal_df,
                           # group = pops,
                           nperm = 999, 
                           mult = 'hochberg')
wild_kendal_post = kendall.post(wild_kendal_df, 
                              mult = 'hochberg', 
                              nperm = 999)



# Kendall concordance - eigen vector traits -------------------------------
ashn_wild_eigen_kendal = ashn_wild_eigen_residuals %>% 
  dplyr::rename(ASHN = eigen_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

myv_wild_eigen_kendal = myv_wild_eigen_residuals%>% 
  dplyr::rename(MYV = eigen_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

skr_wild_eigen_kendal = skr_wild_eigen_residuals%>% 
  dplyr::rename(SKR = eigen_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

gts_gar_wild_eigen_kendal = gts_gar_wild_eigen_residuals%>% 
  dplyr::rename(GTS_GAR = eigen_residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  # dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

wild_eigen_kendal_df = bind_cols(ashn_wild_eigen_kendal, 
                               myv_wild_eigen_kendal, 
                               skr_wild_eigen_kendal, 
                               gts_gar_wild_eigen_kendal)

# group = c('ASHN','MYV','SKR','GTS_GAR')


wild_eigen_kendal = kendall.global(wild_eigen_kendal_df,
                                 # group = pops,
                                 nperm = 999, 
                                 mult = 'hochberg')
wild_eigen_kendal_post = kendall.post(wild_eigen_kendal_df, 
                                    mult = 'hochberg', 
                                    nperm = 999)
###

# Nemeny post hoc - trait covarianec --------------------------------------
library(PMCMRplus)


ashn_wild_nemeny_df = ashn_wild_residuals %>% 
  mutate(pop_pair = 'ASHN') 

myv_wild_nemeny_df = myv_wild_residuals %>% 
  mutate(pop_pair = 'MYV')

skr_wild_nemeny_df = skr_wild_residuals %>% 
  mutate(pop_pair = 'SKR')

gts_gar_wild_nemeny_df = gts_gar_wild_residuals %>% 
  mutate(pop_pair = 'GTS/GAR') 

wild_nemeny_df = bind_rows(ashn_wild_nemeny_df, 
                         myv_wild_nemeny_df, 
                         skr_wild_nemeny_df, 
                         gts_gar_wild_nemeny_df)


wild_nemeny_test = frdAllPairsNemenyiTest(y = wild_nemeny_df$cov_residual, 
                                        groups = wild_nemeny_df$trait, 
                                        blocks = wild_nemeny_df$pop_pair, 
                                        p.adjust.method = 'hochberg')

wild_nemeny_test$p.value %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname) %>% View() 
filter(if_any(c(2:33), function(x) x<0.05))
write_csv('wild_1212_residual_NemenyTest_results.csv')



# nemeny post hoc - eigen vectors -----------------------------------------

library(PMCMRplus)


ashn_wild_eigen_nemeny_df = ashn_wild_eigen_residuals %>% 
  mutate(pop_pair = 'ASHN') 

myv_wild_eigen_nemeny_df = myv_wild_eigen_residuals %>% 
  mutate(pop_pair = 'MYV') 

skr_wild_eigen_nemeny_df = skr_wild_eigen_residuals %>% 
  mutate(pop_pair = 'SKR') 

gts_gar_wild_eigen_nemeny_df = gts_gar_wild_eigen_residuals %>% 
  mutate(pop_pair = 'GTS/GAR') 

wild_eigen_nemeny_df = bind_rows(ashn_wild_eigen_nemeny_df, 
                               myv_wild_eigen_nemeny_df, 
                               skr_wild_eigen_nemeny_df, 
                               gts_gar_wild_eigen_nemeny_df)


wild_eigen_nemeny_test = frdAllPairsNemenyiTest(y = wild_eigen_nemeny_df$eigen_residual, 
                                              groups = wild_eigen_nemeny_df$trait, 
                                              blocks = wild_eigen_nemeny_df$pop_pair, 
                                              p.adjust.method = 'hochberg')


wild_eigen_nemeny_test$p.value %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname) %>% 
  # View()
  filter(if_any(c(2:33), function(x) x<0.05)) %>% View()
# %>%
# write_csv('wild_4eigen_ALL_residual_NemenyTest_results.csv')

Proc_test_residual_combo

wild_significant_traits =  c("lm_19_18",
                           "fbar_23_24", 
                           'fbar_8_27', 
                           'head_depth', 
                           'lm_6_12') # Replace with actual Nemenyi-identified traits


wild_eigen_nemeny_df %>%
  filter(trait %in% wild_significant_traits) %>%
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
wild_ashn_cov_differences = ashnw_wild_cov_mat - ashnc_wild_cov_mat 
wild_ashn_cov_diff = wild_ashn_cov_differences%>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname) %>% 
  dplyr::filter(trait %in% c(wild_significant_traits)) %>% 
  dplyr::select(trait, 
                lm_19_18, 
                fbar_23_24, 
                fbar_8_27, 
                head_depth, 
                lm_6_12) %>% 
  mutate(pop_pair = 'ASHN')

wild_myv_cov_differences = myvw_wild_cov_mat - myvc_wild_cov_mat 
wild_myv_cov_diff = wild_myv_cov_differences%>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname)%>% 
  dplyr::filter(trait %in% c(wild_significant_traits)) %>% 
  dplyr::select(trait, 
                lm_19_18, 
                fbar_23_24, 
                fbar_8_27, 
                head_depth, 
                lm_6_12) %>% 
  mutate(pop_pair = 'MYV')

wild_skr_cov_differences = skrw_wild_cov_mat - skrc_wild_cov_mat 
wild_skr_cov_diff = wild_skr_cov_differences%>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname)%>% 
  dplyr::filter(trait %in% c(wild_significant_traits)) %>% 
  dplyr::select(trait, 
                lm_19_18, 
                fbar_23_24, 
                fbar_8_27, 
                head_depth, 
                lm_6_12) %>% 
  mutate(pop_pair = 'SKR')

wild_gts_gar_cov_differences = gts_wild_cov_mat - gar_wild_cov_mat 
wild_gts_gar_cov_diff = wild_gts_gar_cov_differences%>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname)%>% 
  dplyr::filter(trait %in% c(wild_significant_traits)) %>% 
  dplyr::select(trait, 
                lm_19_18, 
                fbar_23_24, 
                fbar_8_27, 
                head_depth, 
                lm_6_12) %>% 
  mutate(pop_pair = 'GTS/GAR')


wild_sig_cov_mat = bind_rows(wild_ashn_cov_diff, 
                           wild_myv_cov_diff, 
                           wild_skr_cov_diff, 
                           wild_gts_gar_cov_diff)

wild_sig_cov_long <- wild_sig_cov_mat %>%
  pivot_longer(
    cols = all_of(wild_significant_traits),
    names_to = "Trait2",
    values_to = "Covariance") %>%
  dplyr::rename(Trait1 = trait) %>%
  mutate(Trait1 = factor(Trait1, levels = wild_significant_traits),
         Trait2 = factor(Trait2, levels = wild_significant_traits),
         pop_pair = factor(
           pop_pair,
           levels = c("ASHN", "MYV", "SKR", "GTS/GAR")))


wild_sig_cov_triangle = wild_sig_cov_long %>%
  mutate(Trait1_num = as.integer(Trait1),
         Trait2_num = as.integer(Trait2)) %>%
  filter(Trait1_num > Trait2_num)

wild_sig_covariance_traits = ggplot(wild_sig_cov_triangle,
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
       fill = expression(Delta~'Covariance'),, 
       title = 'A) F2 generation')+
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold"), 
        strip.background = element_rect(fill = 'white'), 
        strip.text.x = element_text(face = 'bold', 
                                    size = 14),
        plot.title = element_text(size = 16),
        legend.position = 'none')








# Code I dont need anymore ------------------------------------------------


# wild data - KrzCov and procrustes ---------------------------------------
# ASHN --------------------------------------------------------------------


ashnc_wild_trait_mat = wild_data %>% 
  filter(Lake_morph %in% c('ASHNC')) %>% 
  dplyr::select(5:39) 
# %>% 
#   as.matrix()

ashnw_wild_trait_mat = wild_data %>% 
  filter(Lake_morph %in% c('ASHNW')) %>% 
  dplyr::select(5:39) 

ashnc_trait_cov_mat = cov(ashnc_wild_trait_mat, 
                          use = 'complete.obs')

ashnw_trait_cov_mat = cov(ashnw_wild_trait_mat, 
                          use = 'complete.obs')

ashn_krzcov = evolqg::KrzCor(ashnc_trait_cov_mat, 
                     ashnw_trait_cov_mat)


ashn_proc_res = procrustes(ashnc_trait_cov_mat, 
                          ashnw_trait_cov_mat)

summary(ashn_proc_res)

plot(ashn_proc_res)
plot(ashn_proc_res, 
     kind = 2)

ashn_test_res <- protest(ashnc_trait_cov_mat, 
                        ashnw_trait_cov_mat, 
                        permutations = 999)

ashn_wild_residuals = residuals(ashn_test_res)

ashn_wild_residuals %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
         residual = 2) %>% 
  write_csv('ASHN_WILD_procrustes_cov_mat_residuals.csv')

ashn_wild_residuals %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2) 

# %>% View() 
#   filter(residual > 0.14)
# MYV ---------------------------------------------------------------------
myvc_wild_trait_mat = wild_data %>% 
  filter(Lake_morph %in% c('MYVC')) %>% 
  dplyr::select(5:39) 
# %>% 
#   as.matrix()

myvw_wild_trait_mat = wild_data %>% 
  filter(Lake_morph %in% c('MYVW')) %>% 
  dplyr::select(5:39) 

myvc_trait_cov_mat = cov(myvc_wild_trait_mat, 
                         use = 'complete.obs')

myvw_trait_cov_mat = cov(myvw_wild_trait_mat, 
                         use = 'complete.obs')

myv_krzcov = KrzCor(myvc_trait_cov_mat, 
                    myvw_trait_cov_mat)

myv_proc_res = procrustes(myvc_trait_cov_mat, 
                          myvw_trait_cov_mat)

summary(myv_proc_res)

plot(myv_proc_res)
plot(myv_proc_res, 
     kind = 2)

myv_test_res <- protest(myvc_trait_cov_mat, 
                        myvw_trait_cov_mat, 
                        permutations = 999)

myv_wild_residuals = residuals(myv_test_res)

myv_wild_residuals %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2) %>% 
  write_csv('MYV_WILD_procrustes_cov_mat_residuals.csv')

# myv_wild_residuals %>% 
#   as.data.frame() %>%
#   rownames_to_column() %>% 
#   as_tibble() %>% 
#   dplyr::rename(trait = rowname, 
#                 residual = 2) %>%View()
# SKR ---------------------------------------------------------------------

skrc_wild_trait_mat = wild_data %>% 
  filter(Lake_morph %in% c('SKRC')) %>% 
  dplyr::select(5:39) 
# %>% 
#   as.matrix()

skrw_wild_trait_mat = wild_data %>% 
  filter(Lake_morph %in% c('SKRW')) %>% 
  dplyr::select(5:39) 

skrc_trait_cov_mat = cov(skrc_wild_trait_mat, 
                         use = 'complete.obs')

skrw_trait_cov_mat = cov(skrw_wild_trait_mat, 
                         use = 'complete.obs')

skr_krzcov = KrzCor(skrc_trait_cov_mat, 
                    skrw_trait_cov_mat)

skr_proc_res = procrustes(skrc_trait_cov_mat, 
                          skrw_trait_cov_mat)

summary(skr_proc_res)

plot(skr_proc_res)
plot(skr_proc_res, 
     kind = 2)

skr_test_res <- protest(skrc_trait_cov_mat, 
                        skrw_trait_cov_mat, 
                        permutations = 999)

skr_wild_residuals = residuals(skr_test_res)

skr_wild_residuals %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2) %>% View()
  write_csv('SKR_WILD_procrustes_cov_mat_residuals.csv')


# GAR/GTS -----------------------------------------------------------------

GAR_wild_trait_mat = wild_data %>% 
  filter(Lake_morph %in% c('CSWY')) %>% 
  dplyr::select(5:39) 
# %>% 
#   as.matrix()

GTS_wild_trait_mat = wild_data %>% 
  filter(Lake_morph %in% c('GTS')) %>% 
  dplyr::select(5:39) 

GAR_trait_cov_mat = cov(GAR_wild_trait_mat, 
                        use = 'complete.obs')

GTS_trait_cov_mat = cov(GTS_wild_trait_mat, 
                        use = 'complete.obs')

gts_gar_krzcov = KrzCor(GAR_trait_cov_mat, 
                        GTS_trait_cov_mat)

gts_gar_proc_res = procrustes(GAR_trait_cov_mat, 
                              GTS_trait_cov_mat)

summary(gts_gar_proc_res)

plot(gts_gar_proc_res)
plot(gts_gar_proc_res, 
     kind = 2)

gts_gar_test_res <- protest(GAR_trait_cov_mat, 
                            GTS_trait_cov_mat, 
                            permutations = 999)

gts_gar_wild_residuals = residuals(gts_gar_test_res)

gts_gar_wild_residuals %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname, 
                residual = 2) %>%  
  write_csv('GTS_GAR_WILD_procrustes_cov_mat_residuals.csv')


# protest the residuals ---------------------------------------------------

ashn_wild_residuals = read_csv('ASHN_WILD_procrustes_cov_mat_residuals.csv') %>% 
  arrange(desc(residual)) %>% 
  mutate(resid_10 = residual*10)
myv_wild_residuals = read_csv('MYV_WILD_procrustes_cov_mat_residuals.csv')%>% 
  arrange(desc(residual)) %>% 
  mutate(resid_10 = residual*10)
skr_wild_residuals = read_csv('SKR_WILD_procrustes_cov_mat_residuals.csv')%>% 
  arrange(desc(residual)) %>% 
  mutate(resid_10 = residual*10)
gts_gar_wild_residuals = read_csv('GTS_GAR_WILD_procrustes_cov_mat_residuals.csv')%>% 
  arrange(desc(residual)) %>% 
  mutate(resid_10 = residual*10)


ggplot(ashn_wild_residuals, 
       aes(x = factor(trait), 
           y = resid_10))+
  geom_col()+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

ggplot(myv_wild_residuals, 
       aes(x = factor(trait), 
           y = resid_10))+
  geom_col()+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

ggplot(skr_wild_residuals, 
       aes(x = factor(trait), 
           y = resid_10))+
  geom_col()+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

ggplot(gts_gar_wild_residuals, 
       aes(x = factor(trait), 
           y = resid_10))+
  geom_col()+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 6, 
                                   angle = 90))

# Combine trait covariance matrices ---------------------------------------


cov_matrices <- list(
  Pair1_Ambient = ashnc_trait_cov_mat,
  Pair1_Geothermal = ashnw_trait_cov_mat,
  Pair2_Ambient = myvc_trait_cov_mat,
  Pair2_Geothermal = myvw_trait_cov_mat,
  Pair3_Ambient = skrc_trait_cov_mat,
  Pair3_Geothermal = skrw_trait_cov_mat,
  Pair4_Ambient = GAR_trait_cov_mat,
  Pair4_Geothermal = GTS_trait_cov_mat
)

for (i in 1:4) {
  cov_matrices[[paste0("Pair", i, "_Difference")]] <-
    cov_matrices[[paste0("Pair", i, "_Geothermal")]] -
    cov_matrices[[paste0("Pair", i, "_Ambient")]]
}

cov_long <- purrr::imap_dfr(cov_matrices, function(mat, name) {
  as.data.frame(as.table(mat)) %>%
    setNames(c("Trait1", "Trait2", "Covariance")) %>%
    mutate(
      Pair = sub("_.*", "", name),
      Type = sub(".*_", "", name)
    )
})

cov_long$Type <- factor(
  cov_long$Type,
  levels = c("Ambient", "Geothermal", "Difference")
)

cov_long = cov_long %>% 
  mutate(Ecotype_pair = case_when(
    str_detect(Pair, 'Pair1') ~ 'ASHN',
    str_detect(Pair, 'Pair2') ~ 'MYV',
    str_detect(Pair, 'Pair3') ~ 'SKR',
    str_detect(Pair, 'Pair4') ~ 'GTS/GAR'))

ggplot(cov_long, 
       aes(Trait1, 
           Trait2, 
           fill = Covariance)) +
  geom_tile() +
  facet_grid(Ecotype_pair ~ Type) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0) +
  coord_equal() +
  labs(x = NULL, 
       y = NULL, 
       fill = "Covariance")+
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1),
        panel.grid = element_blank(), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'), 
        axis.text = element_text(size = 6))




# KrzCov null distribution resampling -------------------------------------


# ASHN --------------------------------------------------------------------


# Combine individuals from both ecotypes
ashn_combined = bind_rows(ashnc_wild_trait_mat,
  ashnw_wild_trait_mat)

n_ashn_cold = nrow(ashnc_wild_trait_mat)
n_ashn_warm = nrow(ashnw_wild_trait_mat)

# Observed similarity
ashn_observed = evolqg::KrzCor(cov(ashnc_wild_trait_mat),
  cov(ashnw_wild_trait_mat))

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
myv_combined = bind_rows(myvc_wild_trait_mat,
                          myvw_wild_trait_mat)

n_myv_cold = nrow(myvc_wild_trait_mat)
n_myv_warm = nrow(myvw_wild_trait_mat)

# Observed similarity
myv_observed = evolqg::KrzCor(cov(myvc_wild_trait_mat),
                                cov(myvw_wild_trait_mat))

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
skr_combined = bind_rows(skrc_wild_trait_mat,
                         skrw_wild_trait_mat)

n_skr_cold = nrow(skrc_wild_trait_mat)
n_skr_warm = nrow(skrw_wild_trait_mat)

# Observed similarity
skr_observed = evolqg::KrzCor(cov(skrc_wild_trait_mat),
                              cov(skrw_wild_trait_mat))

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
gts_gar_combined = bind_rows(GAR_wild_trait_mat,
                         GTS_wild_trait_mat)

n_GAR_cold = nrow(GAR_wild_trait_mat)
n_GTS_warm = nrow(GTS_wild_trait_mat)

# Observed similarity
gts_gar_observed = evolqg::KrzCor(cov(GAR_wild_trait_mat),
                              cov(GTS_wild_trait_mat))

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


# Eigendecomposition of covariance matrices
ashnc_eigen = eigen(ashnc_trait_cov_mat)
ashnw_eigen = eigen(ashnw_trait_cov_mat)

# Number of principal components to retain
# k = floor(ncol(ashnc_trait_cov_mat) / 2)

# Extract leading eigenvectors
ashnc_vectors = ashnc_eigen$vectors[, 1:3, drop = FALSE]
ashnw_vectors = ashnw_eigen$vectors[, 1:3, drop = FALSE]

# Preserve trait identities
rownames(ashnc_vectors) = colnames(ashnc_trait_cov_mat)
rownames(ashnw_vectors) = colnames(ashnw_trait_cov_mat)

ashn_proc_res = procrustes(ashnc_vectors, 
                           ashnw_vectors, 
                           symmetric = T)

summary(ashn_proc_res)

plot(ashn_proc_res)
plot(ashn_proc_res, 
     kind = 2)

ashn_test_res = protest(ashnc_vectors, 
                         ashnw_vectors, 
                         permutations = 999, 
                        scale = F, 
                        symmetric = T)


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

