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




