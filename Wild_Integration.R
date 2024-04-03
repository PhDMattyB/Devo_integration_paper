##############################
## Wild fish integration analysis
##
## Matt Brachmann (PhDMattyB)
##
## 04.03.2024
##
##############################


setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/Updated Landmarks/')

library(geomorph)
library(vcvComp)
library(factoextra)
library(tidyverse)

wild_tps = readland.tps('Wild_Final.TPS', 
                        specID = 'imageID')

identifiers = read_csv('TPS_Wild_metadata.csv') 

## superimposition on the entire dataset
wild_gpa = gpagen(wild_tps, 
                       print.progress = F)

wild_gpa$coords

dim(wild_gpa$coords)



# Allometry ---------------------------------------------------------------

allometry_model1 = procD.lm(wild_gpa$coords ~ log(wild_gpa$Csize), 
         iter = 999, 
         RRPP = T)
summary(allometry_model1)


shape_resid = arrayspecs(allometry_model1$residuals, 
                         p = dim(wild_gpa$coords)[1], 
                         k = dim(wild_gpa$coords)[2])
allometry_adj_shape = shape_resid + array(wild_gpa$consensus, 
                                dim(shape_resid))


# allometry_model2 = procD.lm(wild_gpa$coords ~ wild_gpa$Csize * identifiers$Lake, 
#                             iter = 999, 
#                             RRPP = T)
# summary(allometry_model2)
# 
# 
# allometry_model3 = procD.lm(wild_gpa$coords ~ wild_gpa$Csize * identifiers$Morph, 
#                             iter = 999, 
#                             RRPP = T)
# summary(allometry_model3)
# 

allometry_model4 = procD.lm(wild_gpa$coords ~ log(wild_gpa$Csize) * identifiers$Lake_morph, 
                            iter = 999, 
                            RRPP = T)
summary(allometry_model4)



# PCA ---------------------------------------------------------------------


## Need to figure out a way to make this ggplot able
PCA_allometry = gm.prcomp(wild_gpa$coords)
plot(PCA_allometry,
     pch=21,
     # bg=identifiers$Lake_morph,
     cex=1.5)


PCA_allometry_adj = gm.prcomp(allometry_adj_shape)
plot(PCA_allometry_adj)
