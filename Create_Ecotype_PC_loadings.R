##############################
## PC loading matrices per ecotype pair
##
## Matt Brachmann (PhDMattyB)
##
## 27.05.2026
##
##############################

setwd("~/Parsons_postdoc/Stickleback_Morphometric_data/Updated Landmarks/")

library(tidyverse)
library(lme4)
library(lmerTest)
library(brms)
library(paran)
library(broom.mixed)
library(MatrixCorrelation)


# functions ---------------------------------------------------------------

cor_fun <- function(data, pc_cols, trait_cols) {
  
  expand_grid(
    pc = pc_cols,
    trait = trait_cols) %>%
    mutate(test = map2(pc,trait,
                       \(x, y)
                       cor.test(
                         x = data[[x]],
                         y = data[[y]],
                         method = "pearson") %>%
                         tidy())) %>%
    unnest(test) %>%
    transmute(
      pc,
      trait,
      r = estimate,
      p = p.value)
}

# wild data ----------------------------------------------------------------
wild_data = read_csv("WILD_SCALED_FIXED_27.05.2026.csv") %>% 
  dplyr::select(Lake_morph, 
                rowname, 
                Lake, 
                Morph, 
                5:39)
ashn_wild_data = wild_data %>% 
  filter(Lake_morph %in% c('ASHNC', 
                           'ASHNW'))
ashn_wild_trait_mat = ashn_wild_data %>% 
  dplyr::select(5:39) %>% 
  as.matrix()

myv_wild_data = wild_data %>% 
  filter(Lake_morph %in% c('MYVC', 
                           'MYVW'))
myv_wild_trait_mat = myv_wild_data %>% 
  dplyr::select(5:39) %>% 
  as.matrix()

skr_wild_data = wild_data %>% 
  filter(Lake_morph %in% c('SKRC', 
                           'SKRW'))
skr_wild_trait_mat = skr_wild_data %>% 
  dplyr::select(5:39) %>% 
  as.matrix()

gts_gar_wild_data = wild_data %>% 
  filter(Lake_morph %in% c('CSWY', 
                           'GTS'))
gts_gar_wild_trait_mat = gts_gar_wild_data %>% 
  dplyr::select(5:39) %>% 
  as.matrix()


# PCA analysis on each wild ecotype pair ----------------------------------

ashn_wild_pca = prcomp(ashn_wild_trait_mat, 
                       center = F, 
                       scale. = F)
paran(ashn_wild_trait_mat, 
      iterations = 1000, 
      graph = T)
## 5 pc axes retained

ashn_wild_scores = as.data.frame(ashn_wild_pca$x[,1:5])
ashn_wild_PCA = bind_cols(ashn_wild_data, 
                          ashn_wild_scores)


myv_wild_pca = prcomp(myv_wild_trait_mat, 
                       center = F, 
                       scale. = F)
paran(myv_wild_trait_mat, 
      iterations = 1000, 
      graph = T)
## 5 pc axes retained

myv_wild_scores = as.data.frame(myv_wild_pca$x[,1:5])
myv_wild_PCA = bind_cols(myv_wild_data, 
                          myv_wild_scores)


skr_wild_pca = prcomp(skr_wild_trait_mat, 
                       center = F, 
                       scale. = F)
paran(skr_wild_trait_mat, 
      iterations = 1000, 
      graph = T)
## 5 pc axes retained

skr_wild_scores = as.data.frame(skr_wild_pca$x[,1:5])
skr_wild_PCA = bind_cols(skr_wild_data, 
                          skr_wild_scores)


gts_gar_wild_pca = prcomp(gts_gar_wild_trait_mat, 
                       center = F, 
                       scale. = F)
paran(gts_gar_wild_trait_mat, 
      iterations = 1000, 
      graph = T)
## 6 pc axes retained

gts_gar_wild_scores = as.data.frame(gts_gar_wild_pca$x[,1:6])
gts_gar_wild_PCA = bind_cols(gts_gar_wild_data, 
                          gts_gar_wild_scores)


