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

## Good enough for now
ggplot(data = pooled_pc_coords, 
       aes(x = PCo1, 
           y = PCo2, 
           col = POP))+
  geom_point(size = 3)+
  scale_color_manual(values = pops_col_pal)

ggplot(data = pooled_pc_coords, 
       aes(x = PCo1, 
           y = PCo3, 
           col = POP))+
  geom_point(size = 3)+
  scale_color_manual(values = pops_col_pal)



## Need to test for differences in variance-covariance relationships
table(LM_data$POP)
## Maximum likelihood test
## ASHNC cold vs warm
prop.vcv.test(n = c(30,30), 
              phenotypes_pooled_var[,,'ASHNC'], 
              phenotypes_pooled_var[,,'ASHNW'])
## 0.34

## MYV cold vs warm
prop.vcv.test(n = c(30,30), 
              phenotypes_pooled_var[,,'MYVC'], 
              phenotypes_pooled_var[,,'MYVW'])
## 0.020

## SKR cold vs warm
prop.vcv.test(n = c(31,29), 
              phenotypes_pooled_var[,,'SKRC'], 
              phenotypes_pooled_var[,,'SKRW'])
## 0.00029

## RKLT cold vs warm
prop.vcv.test(n = c(18,14), 
              phenotypes_pooled_var[,,'RKLTC'], 
              phenotypes_pooled_var[,,'RKLTW'])
## 0.0037

## STN cold vs warm
prop.vcv.test(n = c(32,28), 
              phenotypes_pooled_var[,,'STNC'], 
              phenotypes_pooled_var[,,'STNW'])
## 0.117

## GTS vs CSWY
prop.vcv.test(n = c(30,29), 
              phenotypes_pooled_var[,,'CSWY'], 
              phenotypes_pooled_var[,,'GTS'])

## 0.033





