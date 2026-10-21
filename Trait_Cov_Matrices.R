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



wild_data = read_csv("WILD_SCALED_FIXED_27.05.2026.csv") %>% 
  dplyr::select(Lake_morph, 
                rowname, 
                Lake, 
                Morph, 
                5:39)


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

