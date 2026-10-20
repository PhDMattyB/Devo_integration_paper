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

ashn_krzcov = KrzCor(ashnc_trait_cov_mat, 
                     ashnw_trait_cov_mat)

ashn_proc_res = procrustes(ashnc_trait_cov_mat, 
                           ashnw_trait_cov_mat)

summary(ashn_proc_res)

plot(ashn_proc_res, 
     kind = 1)

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

# plot(myv_proc_res, 
#      kind = 1)

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

# plot(skr_proc_res, 
#      kind = 1)

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
# 
# plot(ashn_proc_res, 
#      kind = 1)

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
  
