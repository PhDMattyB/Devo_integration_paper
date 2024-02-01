##############################
## Variance-Covariance wild stickleback
##
## Matt Brachmann (PhDMattyB)
##
## 25.01.2024
##
##############################

setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/')

# install.packages('geomorph')
# install.packages('vcvComp')

library(geomorph)
library(vcvComp)
library(tidyverse)

theme_set(theme_bw())
## Example data from cichlids
data("Tropheus")

landmarks = read_csv('allometry minimised data (XY) with ID (6 population pairs).csv')


identifiers = landmarks %>% 
  select(ID, 
         POP, 
         Morph, 
         CS)

LM_data = landmarks %>% 
  select(-starts_with('LMS'))

phenotypes = as.matrix(LM_data[which(names(LM_data) == 'LM1X'):
                    which(names(LM_data) == 'LM22Y')])

rownames(phenotypes) = LM_data$ID
# dim(phenotypes)

phenotype_array = arrayspecs(phenotypes, p = 22, k = 2)

phenotype_gpa = gpagen(phenotype_array, 
                       print.progress = F)

proc_coord = two.d.array(phenotype_gpa$coords)
colnames(proc_coord) = colnames(phenotypes)

phenotype_pca = prcomp(proc_coord, 
                       rank. = 5, 
                       tol = sqrt(.Machine$double.eps))
pca_scores = phenotype_pca$x


## This is for within population analyses. 
## Testing for differences between phenotype variance-covariance 
## relationships between morphs within each population
## Need to test for overall differences in cold vs warm morphs

phenotypes_pooled_var = cov.group(pca_scores, 
                                  groups = LM_data$POP)


phenotype_eigen_vals = mat.sq.dist(phenotypes_pooled_var, 
            dist. = 'Riemannian')

prcoa = pr.coord(phenotype_eigen_vals)
prcoa$Variance

## change to ggplot format for better graphs
## don't do this now and my brain isn't working for graphs
barplot(prcoa$Variance$exVar, las = 1, col = "darkblue",
        names.arg = 1:nrow(prcoa$Variance), cex.axis = 0.8, cex  = 0.8,
        xlab = "Dimensions", ylab = "Variance explained")

# plot(prcoa$PCoords[, 1], prcoa$PCoords[, 2])
# abline(h = 0) ; abline(v = 0)
# text(prcoa$PCoords[, 1], 
#      prcoa$PCoords[, 1], 
#      labels = rownames(prcoa$PCoords))
# 

## Need to plot point estimates for each population
## we will plot the pooled phenotypic variances 
## Plot using $PCoords


# Pooled principal coordinate analysis ------------------------------------


pooled_pc_coords = prcoa$PCoords %>% 
  as.data.frame() %>% 
  rownames_to_column('POP') %>% 
  as.tibble()


## Need to colour code the populations a bit better
## ASHN = Dark blue
## MYV = Red
## GTS = BLACK
## CSWY = Dark grey
## SKR = Light blue
## RKL = Hot pink
## STN = Green - blue


pooled_pc_coords = mutate(.data = pooled_pc_coords, 
                      Morph = as.factor(case_when(
                        POP == 'ASHNC' ~ 'Ambient',
                        POP == 'ASHNW' ~ 'Geothermal',
                        POP == 'CSWY' ~ 'Ambient',
                        POP == 'GTS' ~ 'Geothermal',
                        POP == 'MYVC' ~ 'Ambient',
                        POP == 'MYVW' ~ 'Geothermal',
                        POP == 'SKRC' ~ 'Ambient',
                        POP == 'SKRW' ~ 'Geothermal',
                        POP == 'RKLTC' ~ 'Ambient', 
                        POP == 'RKLTW' ~ 'Geothermal', 
                        POP == 'STNC' ~ 'Ambient', 
                        POP == 'STNW' ~ 'Geothermal'
                        
                      )))

pooled_pc_coords = mutate(.data = pooled_pc_coords, 
                          POP_only = as.factor(case_when(
                            POP == 'ASHNC' ~ 'ASHN',
                            POP == 'ASHNW' ~ 'ASHN',
                            POP == 'CSWY' ~ 'CSWY',
                            POP == 'GTS' ~ 'GTS',
                            POP == 'MYVC' ~ 'MYV',
                            POP == 'MYVW' ~ 'MYV',
                            POP == 'SKRC' ~ 'SKR',
                            POP == 'SKRW' ~ 'SKR',
                            POP == 'RKLTC' ~ 'RKLT', 
                            POP == 'RKLTW' ~ 'RKLT', 
                            POP == 'STNC' ~ 'STN', 
                            POP == 'STNW' ~ 'STN'
                            
                          )))

pops_col_pal = c('#277da1',
                 '#023e8a',
                 '#0d1b2a', 
                 '#415a77',
                 '#780000', 
                 '#c1121f', 
                 '#fb6f92',
                 '#edafb8', 
                 '#00b4d8', 
                 '#90e0ef', 
                 '#4d908e', 
                 '#43aa8b')

pop_only_pal = c('#277da1',
                 '#415a77',
                 '#0d1b2a',
                 '#c1121f',
                 '#fb6f92',
                 '#43aa8b',
                 '#90e0ef')

## Good enough for now
principal_coord_analysis = ggplot(data = pooled_pc_coords, 
       aes(x = PCo1, 
           y = PCo2, 
           col = POP_only, 
           group = POP_only))+
  geom_line(col = 'black')+
  geom_point(size = 3, 
             aes(shape = Morph))+
  scale_color_manual(values = pop_only_pal)+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

# Per population relative variance analyses -------------------------------


## Need to test for differences in variance-covariance relationships
table(LM_data$POP)
## Maximum likelihood test


# ASHN relative eigenvectors ----------------------------------------------

## ASHNC cold vs warm
prop.vcv.test(n = c(30,30), 
              phenotypes_pooled_var[,,'ASHNC'], 
              phenotypes_pooled_var[,,'ASHNW'])
## 0.34 covariance matrices not different
relGV.multi(phenotypes_pooled_var[,,c('ASHNC', 'ASHNW')], 
            logGV = F)
ashn_rel_eign = relative.eigen(phenotypes_pooled_var[,,'ASHNC'], 
                               phenotypes_pooled_var[,,'ASHNW'])

plot(ashn_rel_eign$relValues[1:ashn_rel_eign$q], 
     log = 'y', 
     las = 1, 
     col = 'blue', 
     type = 'b', 
     main = 'ASHNC relative to ASHNW', 
     cex = 0.8, 
     cex.main = 1, 
     cex.axis = 0.8, 
     cex.sub = 0.7, 
     sub = paste('Relative generalized variance =', ashn_rel_eign$relGV), 
     xlab = NA, 
     ylab = 'Relative eigenvalues')
abline(h = 1)

## deformation grids to show relative differences in covariation
## between the cold vs warm morphs

ashn = c(which(LM_data$POP %in% 'ASHNC'), 
         which(LM_data$POP %in% 'ASHNW'))
## calculate avg shape for the population
ref_ashn = mshape(phenotype_gpa$coords[,,ashn])

ashn_shape = arrayspecs(t(phenotype_pca$rotation %*% ashn_rel_eign$relVectors), 
                        p = 22, k = 2)
## Need to think about the linkages between the landmarks below
## ask kevin about this and see if someone has used this
# WF <- cbind(c(1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 1, 12, 14, 14), 
#             c(19, 18, 18, 3, 4, 5, 6, 7, 8, 9, 10, 10, 11, 11, 15))

gp3 = gridPar(grid.col = "grey", 
               tar.link.col = "blue", 
               tar.pt.size = 0.7, 
               tar.pt.bg = "blue")

# Visualization of the first dimension
par(new = FALSE, 
    mfrow = c(1, 2), 
    mar = c(0.5, 0.5, 0.5, 0.5))

plotRefToTarget(ref_ashn, 
                (ref_ashn - 0.01 * ashn_shape[,,1]), 
                mag = 4, 
                method = 'TPS', 
                gridPar = gp3)
title(main = "-", line = -1)
plotRefToTarget(ref_ashn, 
                (ref_ashn + 0.01 * ashn_shape[,,1]), 
                mag = 4, 
                method = 'TPS', 
                gridPar = gp3)
title(main = "+", line = -1)
title("First relative eigenvector - ASHN", outer = TRUE, line = - 1)
plotRefToTarget(ref_ashn, 
                (ref_ashn - 0.01 * ashn_shape[,,5]), 
                mag = 4, 
                method = 'TPS', 
                gridPar = gp3)
title(main = "-", line = -1)
plotRefToTarget(ref_ashn, 
                (ref_ashn + 0.01 * ashn_shape[,,5]), 
                mag = 4, 
                method = 'TPS', 
                gridPar = gp3)
title(main = "+", line = -1)
title("Last relative eigenvector - ASHN", outer = TRUE, line = - 1)

# Visualization of the last dimension (fig. 8: bottom)

## These are the features with max excess of variance in 
## ASHNC relative to ASHNW
## OR these are the shape features maximally canalized in 
## the warm populations. Selection is acting on head shape
## Canalization has occured in the warm morphs??

# MYV relative eigenvectors -----------------------------------------------

## MYV cold vs warm
prop.vcv.test(n = c(30,30), 
              phenotypes_pooled_var[,,'MYVC'], 
              phenotypes_pooled_var[,,'MYVW'])
## 0.020 covariance matrices different

relGV.multi(phenotypes_pooled_var[,,c('MYVC', 'MYVW')], 
            logGV = F)
myv_rel_eign = relative.eigen(phenotypes_pooled_var[,,'MYVC'], 
                               phenotypes_pooled_var[,,'MYVW'])

plot(myv_rel_eign$relValues[1:myv_rel_eign$q], 
     log = 'y', 
     las = 1, 
     col = 'blue', 
     type = 'b', 
     main = 'MYVC relative to MYVW', 
     cex = 0.8, 
     cex.main = 1, 
     cex.axis = 0.8, 
     cex.sub = 0.7, 
     sub = paste('Relative generalized variance =', myv_rel_eign$relGV), 
     xlab = NA, 
     ylab = 'Relative eigenvalues')
abline(h = 1)

## deformation grids to show relative differences in covariation
## between the cold vs warm morphs

myv = c(which(LM_data$POP %in% 'MYVC'), 
         which(LM_data$POP %in% 'MYVW'))
## calculate avg shape for the population
ref_myv = mshape(phenotype_gpa$coords[,,myv])

myv_shape = arrayspecs(t(phenotype_pca$rotation %*% myv_rel_eign$relVectors), 
                        p = 22, k = 2)
## Need to think about the linkages between the landmarks below
## ask kevin about this and see if someone has used this
# WF <- cbind(c(1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 1, 12, 14, 14), 
#             c(19, 18, 18, 3, 4, 5, 6, 7, 8, 9, 10, 10, 11, 11, 15))

gp3 = gridPar(grid.col = "grey", 
              tar.link.col = "blue", 
              tar.pt.size = 0.7, 
              tar.pt.bg = "blue")

# Visualization of the first dimension
par(new = FALSE, 
    mfrow = c(1, 2), 
    mar = c(0.5, 0.5, 0.5, 0.5))

plotRefToTarget(ref_myv, 
                (ref_myv - 0.01 * myv_shape[,,1]), 
                mag = 4, 
                method = 'TPS', 
                gridPar = gp3)
title(main = "-", line = -1)
plotRefToTarget(ref_myv, 
                (ref_myv + 0.01 * myv_shape[,,1]), 
                mag = 4, 
                method = 'TPS', 
                gridPar = gp3)
title(main = "+", line = -1)
title("First relative eigenvector - MYV", outer = TRUE, line = - 1)


# SKR relative eigenvectors -----------------------------------------------

## SKR cold vs warm
prop.vcv.test(n = c(31,29), 
              phenotypes_pooled_var[,,'SKRC'], 
              phenotypes_pooled_var[,,'SKRW'])
## 0.00029 covariance matrices different
relGV.multi(phenotypes_pooled_var[,,c('SKRC', 'SKRW')], 
            logGV = F)
skr_rel_eign = relative.eigen(phenotypes_pooled_var[,,'SKRC'], 
                              phenotypes_pooled_var[,,'SKRW'])

plot(skr_rel_eign$relValues[1:skr_rel_eign$q], 
     log = 'y', 
     las = 1, 
     col = 'blue', 
     type = 'b', 
     main = 'SKRC relative to SKRW', 
     cex = 0.8, 
     cex.main = 1, 
     cex.axis = 0.8, 
     cex.sub = 0.7, 
     sub = paste('Relative generalized variance =', skr_rel_eign$relGV), 
     xlab = NA, 
     ylab = 'Relative eigenvalues')
abline(h = 1)

## deformation grids to show relative differences in covariation
## between the cold vs warm morphs

skr = c(which(LM_data$POP %in% 'SKRC'), 
         which(LM_data$POP %in% 'SKRW'))
## calculate avg shape for the population
ref_skr = mshape(phenotype_gpa$coords[,,skr])

skr_shape = arrayspecs(t(phenotype_pca$rotation %*% skr_rel_eign$relVectors), 
                        p = 22, k = 2)
## Need to think about the linkages between the landmarks below
## ask kevin about this and see if someone has used this
# WF <- cbind(c(1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 1, 12, 14, 14), 
#             c(19, 18, 18, 3, 4, 5, 6, 7, 8, 9, 10, 10, 11, 11, 15))

gp3 = gridPar(grid.col = "grey", 
              tar.link.col = "blue", 
              tar.pt.size = 0.7, 
              tar.pt.bg = "blue")

# Visualization of the first dimension
par(new = FALSE, 
    mfrow = c(1, 2), 
    mar = c(0.5, 0.5, 0.5, 0.5))

plotRefToTarget(ref_skr, 
                (ref_skr - 0.01 * skr_shape[,,1]), 
                mag = 4, 
                method = 'TPS', 
                gridPar = gp3)
title(main = "-", line = -1)
plotRefToTarget(ref_skr, 
                (ref_skr + 0.01 * skr_shape[,,1]), 
                mag = 4, 
                method = 'TPS', 
                gridPar = gp3)
title(main = "+", line = -1)
title("First relative eigenvector - SKR", outer = TRUE, line = - 1)

# RKLT relative eigenvectors ----------------------------------------------

## RKLT cold vs warm
prop.vcv.test(n = c(18,14), 
              phenotypes_pooled_var[,,'RKLTC'], 
              phenotypes_pooled_var[,,'RKLTW'])
## 0.0037 covariance matrices different

relGV.multi(phenotypes_pooled_var[,,c('RKLTC', 'RKLTW')], 
            logGV = F)
rklt_rel_eign = relative.eigen(phenotypes_pooled_var[,,'RKLTC'], 
                              phenotypes_pooled_var[,,'RKLTW'])

plot(rklt_rel_eign$relValues[1:rklt_rel_eign$q], 
     log = 'y', 
     las = 1, 
     col = 'blue', 
     type = 'b', 
     main = 'RKLTC relative to RKLTW', 
     cex = 0.8, 
     cex.main = 1, 
     cex.axis = 0.8, 
     cex.sub = 0.7, 
     sub = paste('Relative generalized variance =', rklt_rel_eign$relGV), 
     xlab = NA, 
     ylab = 'Relative eigenvalues')
abline(h = 1)

## deformation grids to show relative differences in covariation
## between the cold vs warm morphs

RKLT = c(which(LM_data$POP %in% 'RKLTC'), 
         which(LM_data$POP %in% 'RKLTW'))
## calculate avg shape for the population
ref_RKLT = mshape(phenotype_gpa$coords[,,RKLT])

RKLT_shape = arrayspecs(t(phenotype_pca$rotation %*% rklt_rel_eign$relVectors), 
                        p = 22, k = 2)
## Need to think about the linkages between the landmarks below
## ask kevin about this and see if someone has used this
# WF <- cbind(c(1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 1, 12, 14, 14), 
#             c(19, 18, 18, 3, 4, 5, 6, 7, 8, 9, 10, 10, 11, 11, 15))

gp3 = gridPar(grid.col = "grey", 
              tar.link.col = "blue", 
              tar.pt.size = 0.7, 
              tar.pt.bg = "blue")

# Visualization of the first dimension
par(new = FALSE, 
    mfrow = c(1, 2), 
    mar = c(0.5, 0.5, 0.5, 0.5))

plotRefToTarget(ref_RKLT, 
                (ref_RKLT - 0.01 * RKLT_shape[,,1]), 
                mag = 4, 
                method = 'TPS', 
                gridPar = gp3)
title(main = "-", line = -1)
plotRefToTarget(ref_RKLT, 
                (ref_RKLT + 0.01 * RKLT_shape[,,1]), 
                mag = 4, 
                method = 'TPS', 
                gridPar = gp3)
title(main = "+", line = -1)
title("First relative eigenvector - RKLT", outer = TRUE, line = - 1)

# STN relative eigenvectors -----------------------------------------------

## STN cold vs warm
prop.vcv.test(n = c(32,28), 
              phenotypes_pooled_var[,,'STNC'], 
              phenotypes_pooled_var[,,'STNW'])
## 0.117 covariance matrices not different 
relGV.multi(phenotypes_pooled_var[,,c('STNC', 'STNW')], 
            logGV = F)
stn_rel_eign = relative.eigen(phenotypes_pooled_var[,,'STNC'], 
                              phenotypes_pooled_var[,,'STNW'])

plot(stn_rel_eign$relValues[1:stn_rel_eign$q], 
     log = 'y', 
     las = 1, 
     col = 'blue', 
     type = 'b', 
     main = 'STNC relative to STNW', 
     cex = 0.8, 
     cex.main = 1, 
     cex.axis = 0.8, 
     cex.sub = 0.7, 
     sub = paste('Relative generalized variance =', stn_rel_eign$relGV), 
     xlab = NA, 
     ylab = 'Relative eigenvalues')
abline(h = 1)

## deformation grids to show relative differences in covariation
## between the cold vs warm morphs

STN = c(which(LM_data$POP %in% 'STNC'), 
         which(LM_data$POP %in% 'STNW'))
## calculate avg shape for the population
ref_STN = mshape(phenotype_gpa$coords[,,STN])

STN_shape = arrayspecs(t(phenotype_pca$rotation %*% stn_rel_eign$relVectors), 
                        p = 22, k = 2)
## Need to think about the linkages between the landmarks below
## ask kevin about this and see if someone has used this
# WF <- cbind(c(1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 1, 12, 14, 14), 
#             c(19, 18, 18, 3, 4, 5, 6, 7, 8, 9, 10, 10, 11, 11, 15))

gp3 = gridPar(grid.col = "grey", 
              tar.link.col = "blue", 
              tar.pt.size = 0.7, 
              tar.pt.bg = "blue")

# Visualization of the first dimension
par(new = FALSE, 
    mfrow = c(1, 2), 
    mar = c(0.5, 0.5, 0.5, 0.5))

plotRefToTarget(ref_STN, 
                (ref_STN - 0.01 * STN_shape[,,1]), 
                mag = 4, 
                method = 'TPS', 
                gridPar = gp3)
title(main = "-", line = -1)
plotRefToTarget(ref_STN, 
                (ref_STN + 0.01 * STN_shape[,,1]), 
                mag = 4, 
                method = 'TPS', 
                gridPar = gp3)
title(main = "+", line = -1)
title("First relative eigenvector - STN", outer = TRUE, line = - 1)

# GTS-CSWY relative eigenvectors ------------------------------------------

## GTS vs CSWY
prop.vcv.test(n = c(30,29), 
              phenotypes_pooled_var[,,'CSWY'], 
              phenotypes_pooled_var[,,'GTS'])

## 0.033 covariance matrices different

relGV.multi(phenotypes_pooled_var[,,c('CSWY', 'GTS')], 
            logGV = F)
gts_cswy_rel_eign = relative.eigen(phenotypes_pooled_var[,,'CSWY'], 
                              phenotypes_pooled_var[,,'GTS'])

plot(gts_cswy_rel_eign$relValues[1:gts_cswy_rel_eign$q], 
     log = 'y', 
     las = 1, 
     col = 'blue', 
     type = 'b', 
     main = 'CSWY relative to GTS', 
     cex = 0.8, 
     cex.main = 1, 
     cex.axis = 0.8, 
     cex.sub = 0.7, 
     sub = paste('Relative generalized variance =', gts_cswy_rel_eign$relGV), 
     xlab = NA, 
     ylab = 'Relative eigenvalues')
abline(h = 1)


## deformation grids to show relative differences in covariation
## between the cold vs warm morphs

gts_cswy = c(which(LM_data$POP %in% 'CSWY'), 
         which(LM_data$POP %in% 'GTS'))
## calculate avg shape for the population
ref_gts_cswy = mshape(phenotype_gpa$coords[,,gts_cswy])

gts_cswy_shape = arrayspecs(t(phenotype_pca$rotation %*% gts_cswy_rel_eign$relVectors), 
                        p = 22, k = 2)
## Need to think about the linkages between the landmarks below
## ask kevin about this and see if someone has used this
# WF <- cbind(c(1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 1, 12, 14, 14), 
#             c(19, 18, 18, 3, 4, 5, 6, 7, 8, 9, 10, 10, 11, 11, 15))

gp3 = gridPar(grid.col = "grey", 
              tar.link.col = "blue", 
              tar.pt.size = 0.7, 
              tar.pt.bg = "blue")

# Visualization of the first dimension
par(new = FALSE, 
    mfrow = c(1, 2), 
    mar = c(0.5, 0.5, 0.5, 0.5))

plotRefToTarget(ref_gts_cswy, 
                (ref_gts_cswy - 0.01 * gts_cswy_shape[,,1]), 
                mag = 4, 
                method = 'TPS', 
                gridPar = gp3)
title(main = "-", line = -1)
plotRefToTarget(ref_gts_cswy, 
                (ref_gts_cswy + 0.01 * gts_cswy_shape[,,1]), 
                mag = 4, 
                method = 'TPS', 
                gridPar = gp3)
title(main = "+", line = -1)
title("First relative eigenvector - CSWY vs GTS", outer = TRUE, line = - 1)



# General warm-cold comparison --------------------------------------------

phenotypes_WC = cov.group(pca_scores, 
                                  groups = LM_data$Morph)
WC_eigen_vals = mat.sq.dist(phenotypes_WC, 
                                   dist. = 'Riemannian')
WC_prcoa = pr.coord(WC_eigen_vals)
WC_prcoa$Variance

## This won't work as there aren't enough dimensions



# between vs within population covariance matrices ------------------------

between_cov = cov.B(pca_scores, 
      groups = LM_data$POP)
dim(LM_data)
within_cov = cov.W(pca_scores, 
                   groups = LM_data$POP)

prop.vcv.test(n = c(12, 331), 
              between_cov, 
              within_cov)


# two block PLS -----------------------------------------------------------
landmarks = read_csv('allometry minimised data (XY) with ID (6 population pairs).csv')


identifiers = landmarks %>% 
  select(ID, 
         POP, 
         Morph, 
         CS)
landmarks = mutate(.data = landmarks, 
                          POP_only = as.factor(case_when(
                            POP == 'ASHNC' ~ 'ASHN',
                            POP == 'ASHNW' ~ 'ASHN',
                            POP == 'CSWY' ~ 'CSWY',
                            POP == 'GTS' ~ 'GTS',
                            POP == 'MYVC' ~ 'MYV',
                            POP == 'MYVW' ~ 'MYV',
                            POP == 'SKRC' ~ 'SKR',
                            POP == 'SKRW' ~ 'SKR',
                            POP == 'RKLTC' ~ 'RKLT', 
                            POP == 'RKLTW' ~ 'RKLT', 
                            POP == 'STNC' ~ 'STN', 
                            POP == 'STNW' ~ 'STN'
                            
                          )))

ashnc_lms = landmarks %>%
  filter(POP == 'ASHNC') %>% 
  select(-starts_with('LMS'))

ashnw_lms = landmarks %>%
  filter(POP == 'ASHNW') %>% 
  select(-starts_with('LMS'))

ashn_pheno = as.matrix(ashnc_lms[which(names(ashnc_lms) == 'LM1X'):
                                  which(names(ashnc_lms) == 'LM22Y')])   

ashnw_pheno = as.matrix(ashnw_lms[which(names(ashnw_lms) == 'LM1X'):
                                  which(names(ashnw_lms) == 'LM22Y')])   


rownames(ashnc_pheno) = ashnc_lms$ID
rownames(ashnw_pheno) = ashnw_lms$ID


# dim(phenotypes)
ashnc_array = arrayspecs(ashnc_pheno, 
                         p = 22, 
                         k = 2)
ashnw_array = arrayspecs(ashnw_pheno, 
                         p = 22, 
                         k = 2)

ashnc_gpa = gpagen(ashnc_array)
ashnw_gpa = gpagen(ashnw_array)


## Hmm this doesn't seem to work
## Need to find another function to compare two sets of shape data
## We can compare integration among landmarks from the same
## Individual though!
PLS = two.b.pls(ashnc_gpa$coords,
                ashnw_gpa$coords,
                iter=999)
summary(PLS)
plot(PLS)


# morphological disparity -------------------------------------------------


# ashn --------------------------------------------------------------------


## Need to use morphol.disparity to look for patterns of variation
ashn = landmarks %>%
  filter(POP_only == 'ASHN')


ashn_pheno = as.matrix(ashn[which(names(ashn) == 'LM1X'):
                                    which(names(ashn) == 'LM22Y')])   

rownames(ashn_pheno) = ashn$ID

# dim(phenotypes)
ashn_array = arrayspecs(ashn_pheno, 
                         p = 22, 
                         k = 2)

ashn_gpa = gpagen(ashn_array)

ashn_df = geomorph.data.frame(ashn_gpa, 
                           # species = plethodon$species, 
                           species = ashn$POP)

morphol.disparity(coords ~ 1,
                  groups = ~species,
                  data = ashn_df, 
                  iter = 999)


# myv ---------------------------------------------------------------------
myv = landmarks %>%
  filter(POP_only == 'MYV')


myv_pheno = as.matrix(myv[which(names(myv) == 'LM1X'):
                              which(names(myv) == 'LM22Y')])   

rownames(myv_pheno) = myv$ID

# dim(phenotypes)
myv_array = arrayspecs(myv_pheno, 
                        p = 22, 
                        k = 2)

myv_gpa = gpagen(myv_array)

myv_df = geomorph.data.frame(myv_gpa, 
                              # species = plethodon$species, 
                              species = myv$POP)

morphol.disparity(coords ~ 1,
                  groups = ~species,
                  data = myv_df,
                  iter = 999)

# skr ---------------------------------------------------------------------
skr = landmarks %>%
  filter(POP_only == 'SKR')


skr_pheno = as.matrix(skr[which(names(skr) == 'LM1X'):
                            which(names(skr) == 'LM22Y')])   

rownames(skr_pheno) = skr$ID

# dim(phenotypes)
skr_array = arrayspecs(skr_pheno, 
                       p = 22, 
                       k = 2)

skr_gpa = gpagen(skr_array)

skr_df = geomorph.data.frame(skr_gpa, 
                             # species = plethodon$species, 
                             species = skr$POP)

morphol.disparity(coords ~ 1,
                  groups = ~species,
                  data = skr_df,
                  iter = 999)

# rklt --------------------------------------------------------------------
rklt = landmarks %>%
  filter(POP_only == 'RKLT')


rklt_pheno = as.matrix(rklt[which(names(rklt) == 'LM1X'):
                            which(names(rklt) == 'LM22Y')])   

rownames(rklt_pheno) = rklt$ID

# dim(phenotypes)
rklt_array = arrayspecs(rklt_pheno, 
                       p = 22, 
                       k = 2)

rklt_gpa = gpagen(rklt_array)

rklt_df = geomorph.data.frame(rklt_gpa, 
                             # species = plethodon$species, 
                             species = rklt$POP)

morphol.disparity(coords ~ 1,
                  groups = ~species,
                  data = rklt_df,
                  iter = 999)

# stn ---------------------------------------------------------------------

stn = landmarks %>%
  filter(POP_only == 'STN')


stn_pheno = as.matrix(stn[which(names(stn) == 'LM1X'):
                            which(names(stn) == 'LM22Y')])   

rownames(stn_pheno) = stn$ID

# dim(phenotypes)
stn_array = arrayspecs(stn_pheno, 
                       p = 22, 
                       k = 2)

stn_gpa = gpagen(stn_array)

stn_df = geomorph.data.frame(stn_gpa, 
                             # species = plethodon$species, 
                             species = stn$POP)

morphol.disparity(coords ~ 1,
                  groups = ~species,
                  data = stn_df,
                  iter = 999)


# gts vs cswy -------------------------------------------------------------

gts_cswy = landmarks %>%
  filter(POP_only %in% c('GTS', 
                         'CSWY'))


gts_cswy_pheno = as.matrix(gts_cswy[which(names(gts_cswy) == 'LM1X'):
                            which(names(gts_cswy) == 'LM22Y')])   

rownames(gts_cswy_pheno) = gts_cswy$ID

# dim(phenotypes)
gts_cswy_array = arrayspecs(gts_cswy_pheno, 
                       p = 22, 
                       k = 2)

gts_cswy_gpa = gpagen(gts_cswy_array)

gts_cswy_df = geomorph.data.frame(gts_cswy_gpa, 
                             # species = plethodon$species, 
                             species = gts_cswy$POP)

morphol.disparity(coords ~ 1,
                  groups = ~species,
                  data = gts_cswy_df,
                  iter = 999)



# Phenotypic modularity ---------------------------------------------------


## Modularity.test might work to find patterns of modularity in warm vs cold

# ASNW modularity ---------------------------------------------------------

ashnw = landmarks %>%
  filter(POP == 'ASHNW')

ashnw_pheno = as.matrix(ashnw[which(names(ashnw) == 'LM1X'):
                              which(names(ashnw) == 'LM22Y')])   

rownames(ashnw_pheno) = ashnw$ID

# dim(phenotypes)
ashnw_array = arrayspecs(ashnw_pheno, 
                        p = 22, 
                        k = 2)
ashnw_gpa = gpagen(ashnw_array)
land.gps = rep('a',22); land.gps[1:6]<-'b'
ashnw_mt = modularity.test(ashnw_gpa$coords,
                           land.gps,
                           CI=T,
                           iter=999)
summary(ashnw_mt) # Test summary
plot(ashnw_mt) # Histogram of CR sampling distribution


# ASHNC modularity --------------------------------------------------------

ashnc = landmarks %>%
  filter(POP == 'ASHNC')

ashnc_pheno = as.matrix(ashnc[which(names(ashnc) == 'LM1X'):
                                which(names(ashnc) == 'LM22Y')])   

rownames(ashnc_pheno) = ashnc$ID

# dim(phenotypes)
ashnc_array = arrayspecs(ashnc_pheno, 
                         p = 22, 
                         k = 2)
ashnc_gpa = gpagen(ashnc_array)
land.gps = rep('a',22); land.gps[1:6]<-'b'
ashnc_mt = modularity.test(ashnc_gpa$coords,
                           land.gps,
                           CI=T,
                           iter=999)
summary(ashnc_mt) # Test summary
plot(ashnc_mt) 


# MYVW modularity ---------------------------------------------------------

myvw = landmarks %>%
  filter(POP == 'MYVW')

myvw_pheno = as.matrix(myvw[which(names(myvw) == 'LM1X'):
                                which(names(myvw) == 'LM22Y')])   

rownames(myvw_pheno) = myvw$ID

# dim(phenotypes)
myvw_array = arrayspecs(myvw_pheno, 
                         p = 22, 
                         k = 2)
myvw_gpa = gpagen(myvw_array)
land.gps = rep('a',22); land.gps[1:6]<-'b'
myvw_mt = modularity.test(myvw_gpa$coords,
                           land.gps,
                           CI=T,
                           iter=999)
summary(myvw_mt) # Test summary
plot(myvw_mt) # Histogram of CR sampling distribution

# MYVC modularity ---------------------------------------------------------

MYVC = landmarks %>%
  filter(POP == 'MYVC')

MYVC_pheno = as.matrix(MYVC[which(names(MYVC) == 'LM1X'):
                              which(names(MYVC) == 'LM22Y')])   

rownames(MYVC_pheno) = MYVC$ID

# dim(phenotypes)
MYVC_array = arrayspecs(MYVC_pheno, 
                        p = 22, 
                        k = 2)
MYVC_gpa = gpagen(MYVC_array)
land.gps = rep('a',22); land.gps[1:6]<-'b'
MYVC_mt = modularity.test(MYVC_gpa$coords,
                          land.gps,
                          CI=T,
                          iter=999)
summary(MYVC_mt) # Test summary
plot(MYVC_mt) # Histogram of CR sampling distribution

# SRKW modularity ---------------------------------------------------------

SKRW = landmarks %>%
  filter(POP == 'SKRW')

SKRW_pheno = as.matrix(SKRW[which(names(SKRW) == 'LM1X'):
                              which(names(SKRW) == 'LM22Y')])   

rownames(SKRW_pheno) = SKRW$ID

# dim(phenotypes)
SKRW_array = arrayspecs(SKRW_pheno, 
                        p = 22, 
                        k = 2)
SKRW_gpa = gpagen(SKRW_array)
land.gps = rep('a',22); land.gps[1:6]<-'b'
SKRW_mt = modularity.test(SKRW_gpa$coords,
                          land.gps,
                          CI=T,
                          iter=999)
summary(SKRW_mt) # Test summary
plot(SKRW_mt) # Histogram of CR sampling distribution

# SRKC modularity ---------------------------------------------------------

SKRC = landmarks %>%
  filter(POP == 'SKRC')

SKRC_pheno = as.matrix(SKRC[which(names(SKRC) == 'LM1X'):
                              which(names(SKRC) == 'LM22Y')])   

rownames(SKRC_pheno) = SKRC$ID

# dim(phenotypes)
SKRC_array = arrayspecs(SKRC_pheno, 
                        p = 22, 
                        k = 2)
SKRC_gpa = gpagen(SKRC_array)
land.gps = rep('a',22); land.gps[1:6]<-'b'
SKRC_mt = modularity.test(SKRC_gpa$coords,
                          land.gps,
                          CI=T,
                          iter=999)
summary(SKRC_mt) # Test summary
plot(SKRC_mt) #

# RKLTW modularity ---------------------------------------------------------

RKLTW = landmarks %>%
  filter(POP == 'RKLTW')

RKLTW_pheno = as.matrix(RKLTW[which(names(RKLTW) == 'LM1X'):
                              which(names(RKLTW) == 'LM22Y')])   

rownames(RKLTW_pheno) = RKLTW$ID

# dim(phenotypes)
RKLTW_array = arrayspecs(RKLTW_pheno, 
                        p = 22, 
                        k = 2)
RKLTW_gpa = gpagen(RKLTW_array)
land.gps = rep('a',22); land.gps[1:6]<-'b'
RKLTW_mt = modularity.test(RKLTW_gpa$coords,
                          land.gps,
                          CI=T,
                          iter=999)
summary(RKLTW_mt) # Test summary
plot(RKLTW_mt) #

# RKLTC modularity ---------------------------------------------------------

RKLTC = landmarks %>%
  filter(POP == 'RKLTC')

RKLTC_pheno = as.matrix(RKLTC[which(names(RKLTC) == 'LM1X'):
                                which(names(RKLTC) == 'LM22Y')])   

rownames(RKLTC_pheno) = RKLTC$ID

# dim(phenotypes)
RKLTC_array = arrayspecs(RKLTC_pheno, 
                         p = 22, 
                         k = 2)
RKLTC_gpa = gpagen(RKLTC_array)
land.gps = rep('a',22); land.gps[1:6]<-'b'
RKLTC_mt = modularity.test(RKLTC_gpa$coords,
                           land.gps,
                           CI=T,
                           iter=999)
summary(RKLTC_mt) # Test summary
plot(RKLTC_mt) #

# STNW modularity ---------------------------------------------------------

STNW = landmarks %>%
  filter(POP == 'STNW')

STNW_pheno = as.matrix(STNW[which(names(STNW) == 'LM1X'):
                                which(names(STNW) == 'LM22Y')])   

rownames(STNW_pheno) = STNW$ID

# dim(phenotypes)
STNW_array = arrayspecs(STNW_pheno, 
                         p = 22, 
                         k = 2)
STNW_gpa = gpagen(STNW_array)
land.gps = rep('a',22); land.gps[1:6]<-'b'
STNW_mt = modularity.test(STNW_gpa$coords,
                           land.gps,
                           CI=T,
                           iter=999)
summary(STNW_mt) # Test summary
plot(STNW_mt) #

# STNC modularity ---------------------------------------------------------

STNC = landmarks %>%
  filter(POP == 'STNC')

STNC_pheno = as.matrix(STNC[which(names(STNC) == 'LM1X'):
                                which(names(STNC) == 'LM22Y')])   

rownames(STNC_pheno) = STNC$ID

# dim(phenotypes)
STNC_array = arrayspecs(STNC_pheno, 
                         p = 22, 
                         k = 2)
STNC_gpa = gpagen(STNC_array)
land.gps = rep('a',22); land.gps[1:6]<-'b'
STNC_mt = modularity.test(STNC_gpa$coords,
                           land.gps,
                           CI=T,
                           iter=999)
summary(STNC_mt) # Test summary
plot(STNC_mt) #

# GTS modularity ---------------------------------------------------------

GTS = landmarks %>%
  filter(POP == 'GTS')

GTS_pheno = as.matrix(GTS[which(names(GTS) == 'LM1X'):
                              which(names(GTS) == 'LM22Y')])   

rownames(GTS_pheno) = GTS$ID

# dim(phenotypes)
GTS_array = arrayspecs(GTS_pheno, 
                        p = 22, 
                        k = 2)
GTS_gpa = gpagen(GTS_array)
land.gps = rep('a',22); land.gps[1:6]<-'b'
GTS_mt = modularity.test(GTS_gpa$coords,
                          land.gps,
                          CI=T,
                          iter=999)
summary(GTS_mt) # Test summary
plot(GTS_mt) #

# CSWY modularity ---------------------------------------------------------

CSWY = landmarks %>%
  filter(POP == 'CSWY')

CSWY_pheno = as.matrix(CSWY[which(names(CSWY) == 'LM1X'):
                              which(names(CSWY) == 'LM22Y')])   

rownames(CSWY_pheno) = CSWY$ID

# dim(phenotypes)
CSWY_array = arrayspecs(CSWY_pheno, 
                        p = 22, 
                        k = 2)
CSWY_gpa = gpagen(CSWY_array)
land.gps = rep('a',22); land.gps[1:6]<-'b'
CSWY_mt = modularity.test(CSWY_gpa$coords,
                          land.gps,
                          CI=T,
                          iter=999)
summary(CSWY_mt) # Test summary
plot(CSWY_mt) #

# ASHN modularity ---------------------------------------------------------

ASHN = landmarks %>%
  filter(POP_only == 'ASHN')

ASHN_pheno = as.matrix(ASHN[which(names(ASHN) == 'LM1X'):
                                which(names(ASHN) == 'LM22Y')])   

rownames(ASHN_pheno) = ASHN$ID

# dim(phenotypes)
ASHN_array = arrayspecs(ASHN_pheno, 
                         p = 22, 
                         k = 2)
ASHN_gpa = gpagen(ASHN_array)
land.gps = rep('a',22); land.gps[1:6]<-'b'
ASHN_mt = modularity.test(ASHN_gpa$coords,
                           land.gps,
                           CI=T,
                           iter=999)
summary(ASHN_mt) # Test summary
plot(ASHN_mt) 

# MYV modularity ---------------------------------------------------------

MYV = landmarks %>%
  filter(POP_only == 'MYV')

MYV_pheno = as.matrix(MYV[which(names(MYV) == 'LM1X'):
                              which(names(MYV) == 'LM22Y')])   

rownames(MYV_pheno) = MYV$ID

# dim(phenotypes)
MYV_array = arrayspecs(MYV_pheno, 
                        p = 22, 
                        k = 2)
MYV_gpa = gpagen(MYV_array)
land.gps = rep('a',22); land.gps[1:6]<-'b'
MYV_mt = modularity.test(MYV_gpa$coords,
                          land.gps,
                          CI=T,
                          iter=999)
summary(MYV_mt) # Test summary
plot(MYV_mt) 

# SKR modularity ---------------------------------------------------------

SKR = landmarks %>%
  filter(POP_only == 'SKR')

SKR_pheno = as.matrix(SKR[which(names(SKR) == 'LM1X'):
                              which(names(SKR) == 'LM22Y')])   

rownames(SKR_pheno) = SKR$ID

# dim(phenotypes)
SKR_array = arrayspecs(SKR_pheno, 
                        p = 22, 
                        k = 2)
SKR_gpa = gpagen(SKR_array)
land.gps = rep('a',22); land.gps[1:6]<-'b'
SKR_mt = modularity.test(SKR_gpa$coords,
                          land.gps,
                          CI=T,
                          iter=999)
summary(SKR_mt) # Test summary
plot(SKR_mt) 

# RKLT modularity ---------------------------------------------------------

RKLT = landmarks %>%
  filter(POP_only == 'RKLT')

RKLT_pheno = as.matrix(RKLT[which(names(RKLT) == 'LM1X'):
                              which(names(RKLT) == 'LM22Y')])   

rownames(RKLT_pheno) = RKLT$ID

# dim(phenotypes)
RKLT_array = arrayspecs(RKLT_pheno, 
                        p = 22, 
                        k = 2)
RKLT_gpa = gpagen(RKLT_array)
land.gps = rep('a',22); land.gps[1:6]<-'b'
RKLT_mt = modularity.test(RKLT_gpa$coords,
                          land.gps,
                          CI=T,
                          iter=999)
summary(RKLT_mt) # Test summary
plot(RKLT_mt) 

# STN modularity ---------------------------------------------------------

STN = landmarks %>%
  filter(POP_only == 'STN')

STN_pheno = as.matrix(STN[which(names(STN) == 'LM1X'):
                              which(names(STN) == 'LM22Y')])   

rownames(STN_pheno) = STN$ID

# dim(phenotypes)
STN_array = arrayspecs(STN_pheno, 
                        p = 22, 
                        k = 2)
STN_gpa = gpagen(STN_array)
land.gps = rep('a',22); land.gps[1:6]<-'b'
STN_mt = modularity.test(STN_gpa$coords,
                          land.gps,
                          CI=T,
                          iter=999)
summary(STN_mt) # Test summary
plot(STN_mt) 
# Phenotypic integration between modules ----------------------------------




