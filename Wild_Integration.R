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

allometry_model1 = procD.lm(wild_gpa$coords ~ wild_gpa$Csize, 
         iter = 999, 
         RRPP = T)
summary(allometry_model1)


allometry_model2 = procD.lm(wild_gpa$coords ~ wild_gpa$Csize * identifiers$Lake, 
                            iter = 999, 
                            RRPP = T)
summary(allometry_model2)


allometry_model3 = procD.lm(wild_gpa$coords ~ wild_gpa$Csize * identifiers$Morph, 
                            iter = 999, 
                            RRPP = T)
summary(allometry_model3)


allometry_model4 = procD.lm(wild_gpa$coords ~ wild_gpa$Csize * identifiers$Lake_morph, 
                            iter = 999, 
                            RRPP = T)
summary(allometry_model4)
