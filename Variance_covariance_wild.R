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

## Morph works but not pop, why? 
phenotypes_pooled_var = cov.group(pca_scores, 
                                  groups = LM_data$POP)

pca_scores = as.matrix(pca_scores)
mat.sq.dist(pca_scores, 
            dist. = 'Riemannian')

pr.coord(pca_scores)


eigen.phen <- mat.sq.dist(S.phen.pooled, dist. = "Riemannian")  # Riemannian distances
prcoa <- pr.coord(eigen.phen)  # ordination
prcoa$Variance  # variance explained