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

wild_identifiers = read_csv('TPS_Wild_metadata.csv') 

wild_identifiers %>% 
  group_by(Lake_morph) %>% 
  summarize(n = n())

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


## Need to figure out a way to make this ggplotable
PCA_allometry = gm.prcomp(wild_gpa$coords)
plot(PCA_allometry,
     pch=21,
     # bg=identifiers$Lake_morph,
     cex=1.5)


PCA_allometry_adj = gm.prcomp(allometry_adj_shape)
plot(PCA_allometry_adj)

# Integration analyses ----------------------------------------------------

# wild integration analysis -----------------------------------------------

## integration with allometric variation

subset_wild_coords = coords.subset(wild_gpa$coords, 
                                   wild_identifiers$Lake_morph)

vrel_wild_coords = Map(function(x) integration.Vrel(x), 
                       subset_wild_coords)

ASHN_compare = compare.ZVrel(vrel_wild_coords$ASHNC, 
                             vrel_wild_coords$ASHNW)

MYV_compare = compare.ZVrel(vrel_wild_coords$MYVC, 
                            vrel_wild_coords$MYVW)

SKR_compare = compare.ZVrel(vrel_wild_coords$SKRC, 
                            vrel_wild_coords$SKRW)

# STN_compare = compare.ZVrel(vrel_wild_coords$STNC, 
#                             vrel_wild_coords$STNW)
# 
# RKLT_compare = compare.ZVrel(vrel_wild_coords$RKLTC, 
#                             vrel_wild_coords$RKLTW)

GTS_CSWY_compare = compare.ZVrel(vrel_wild_coords$CSWY, 
                                 vrel_wild_coords$GTS)



## integration with allometric variation removed
## standardized for allometry

subset_wild_coords_allo_adj = coords.subset(allometry_adj_shape, 
                                   wild_identifiers$Lake_morph)

vrel_wild_coords_allo_adj = Map(function(x) integration.Vrel(x), 
                       subset_wild_coords_allo_adj)

ASHN_noallo_compare = compare.ZVrel(vrel_wild_coords_allo_adj$ASHNC, 
                             vrel_wild_coords_allo_adj$ASHNW)

MYV_noallo_compare = compare.ZVrel(vrel_wild_coords_allo_adj$MYVC, 
                            vrel_wild_coords_allo_adj$MYVW)

SKR_noallo_compare = compare.ZVrel(vrel_wild_coords_allo_adj$SKRC, 
                            vrel_wild_coords_allo_adj$SKRW)

# STN_noallo_compare = compare.ZVrel(vrel_wild_coords_allo_adj$STNC, 
#                             vrel_wild_coords_allo_adj$STNW)
# 
# RKLT_noallo_compare = compare.ZVrel(vrel_wild_coords_allo_adj$RKLTC, 
#                              vrel_wild_coords_allo_adj$RKLTW)

GTS_CSWY_noallo_compare = compare.ZVrel(vrel_wild_coords_allo_adj$CSWY, 
                                 vrel_wild_coords_allo_adj$GTS)



# wild craniofacial integration -------------------------------------------

## integration of craniofacial patterns with allometric variation included

wild_craniofacial = readland.tps('Wild_Craniofacial_LM.TPS',
                        specID = 'imageID')

identifiers = read_csv('TPS_Wild_metadata.csv')

## superimposition on the entire dataset
wild_craniofacial_gpa = gpagen(wild_craniofacial,
                  print.progress = F)


subset_wild_craniofacial_coords = coords.subset(wild_craniofacial_gpa$coords,
                                   identifiers$Lake_morph)

vrel_wild_craniofacial = Map(function(x) integration.Vrel(x),
                       subset_wild_craniofacial_coords)

ASHN_craniofacial = compare.ZVrel(vrel_wild_craniofacial$ASHNC,
                             vrel_wild_craniofacial$ASHNW)

MYV_craniofacial = compare.ZVrel(vrel_wild_craniofacial$MYVC,
                            vrel_wild_craniofacial$MYVW)

SKR_craniofacial = compare.ZVrel(vrel_wild_craniofacial$SKRC,
                            vrel_wild_craniofacial$SKRW)

STN_craniofacial = compare.ZVrel(vrel_wild_craniofacial$STNC,
                            vrel_wild_craniofacial$STNW)

RKLT_craniofacial = compare.ZVrel(vrel_wild_craniofacial$RKLTC,
                             vrel_wild_craniofacial$RKLTW)

GTS_CSWY_craniofacial = compare.ZVrel(vrel_wild_craniofacial$CSWY,
                                 vrel_wild_craniofacial$GTS)


# wild craniofacial integration no allometry -------------------------------------------


wild_craniofacial = readland.tps('Wild_Craniofacial_LM.TPS', 
                                 specID = 'imageID')

identifiers = read_csv('TPS_Wild_metadata.csv') 

## superimposition on the entire dataset
wild_craniofacial_gpa = gpagen(wild_craniofacial, 
                               print.progress = F)


craniofacial_allometry = procD.lm(wild_craniofacial_gpa$coords ~ log(wild_craniofacial_gpa$Csize), 
                            iter = 999, 
                            RRPP = T)
summary(craniofacial_allometry)


craniofacial_resid = arrayspecs(craniofacial_allometry$residuals, 
                         p = dim(wild_craniofacial_gpa$coords)[1], 
                         k = dim(wild_craniofacial_gpa$coords)[2])
craniofacial_allo_adj = craniofacial_resid + array(wild_craniofacial_gpa$consensus, 
                                          dim(craniofacial_resid))



wild_craniofacial_coords_noallo = coords.subset(craniofacial_allo_adj, 
                                                identifiers$Lake_morph)

vrel_wild_craniofacial_noallo = Map(function(x) integration.Vrel(x), 
                             wild_craniofacial_coords_noallo)

ASHN_noallo_craniofacial = compare.ZVrel(vrel_wild_craniofacial_noallo$ASHNC, 
                                  vrel_wild_craniofacial_noallo$ASHNW)

MYV_noallo_craniofacial = compare.ZVrel(vrel_wild_craniofacial_noallo$MYVC, 
                                 vrel_wild_craniofacial_noallo$MYVW)

## This is the only model that shows that there is significantly different
## integration in craniofacial shape between the two ecotypes
SKR_noallo_craniofacial = compare.ZVrel(vrel_wild_craniofacial_noallo$SKRC, 
                                 vrel_wild_craniofacial_noallo$SKRW)

STN_noallo_craniofacial = compare.ZVrel(vrel_wild_craniofacial_noallo$STNC, 
                                 vrel_wild_craniofacial_noallo$STNW)

RKLT_noallo_craniofacial = compare.ZVrel(vrel_wild_craniofacial_noallo$RKLTC, 
                                  vrel_wild_craniofacial_noallo$RKLTW)

GTS_CSWY_noallo_craniofacial = compare.ZVrel(vrel_wild_craniofacial_noallo$CSWY, 
                                      vrel_wild_craniofacial_noallo$GTS)




# Integration between ecotypes ---------------------------------------------

wild_craniofacial = readland.tps('Wild_Craniofacial_LM.TPS',
                                 specID = 'imageID')

identifiers = read_csv('TPS_Wild_metadata.csv')

## superimposition on the entire dataset
wild_craniofacial_gpa = gpagen(wild_craniofacial,
                               print.progress = F)


subset_wild_craniofacial_coords = coords.subset(wild_craniofacial_gpa$coords,
                                                identifiers$Morph)

vrel_wild_craniofacial = Map(function(x) integration.Vrel(x),
                             subset_wild_craniofacial_coords)

Ecotype_craniofacial = compare.ZVrel(vrel_wild_craniofacial$Cold,
                                  vrel_wild_craniofacial$Warm)



## Analysis with no allometric variation

wild_craniofacial = readland.tps('Wild_Craniofacial_LM.TPS', 
                                 specID = 'imageID')

identifiers = read_csv('TPS_Wild_metadata.csv') 

## superimposition on the entire dataset
wild_craniofacial_gpa = gpagen(wild_craniofacial, 
                               print.progress = F)


craniofacial_allometry = procD.lm(wild_craniofacial_gpa$coords ~ log(wild_craniofacial_gpa$Csize), 
                                  iter = 999, 
                                  RRPP = T)
summary(craniofacial_allometry)


craniofacial_resid = arrayspecs(craniofacial_allometry$residuals, 
                                p = dim(wild_craniofacial_gpa$coords)[1], 
                                k = dim(wild_craniofacial_gpa$coords)[2])
craniofacial_allo_adj = craniofacial_resid + array(wild_craniofacial_gpa$consensus, 
                                                   dim(craniofacial_resid))



wild_craniofacial_coords_noallo = coords.subset(craniofacial_allo_adj, 
                                                identifiers$Morph)

vrel_wild_craniofacial_noallo = Map(function(x) integration.Vrel(x), 
                                    wild_craniofacial_coords_noallo)

Ecotype_noallo_craniofacial = compare.ZVrel(vrel_wild_craniofacial_noallo$Cold, 
                                         vrel_wild_craniofacial_noallo$Warm)





