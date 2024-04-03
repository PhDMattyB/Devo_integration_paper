
##############################
## F2 fish integration analysis
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

F2_tps = readland.tps('F2_No_GT.TPS', 
                        specID = 'imageID')

identifiers = read_csv('F2_Metadata.CSV', 
                       col_names = T) 

## superimposition on the entire dataset
F2_gpa = gpagen(F2_tps, 
                  print.progress = F)

allometry_model1 = procD.lm(F2_gpa$coords ~ log(F2_gpa$Csize), 
                            iter = 999, 
                            RRPP = T)
summary(allometry_model1)


F2_shape_resid = arrayspecs(allometry_model1$residuals, 
                         p = dim(F2_gpa$coords)[1], 
                         k = dim(F2_gpa$coords)[2])
F2_allometry_adj_shape = F2_shape_resid + array(F2_gpa$consensus, 
                                          dim(F2_shape_resid))


# F2 integration analysis -----------------------------------------------

## integration with allometric variation

subset_F2_coords = coords.subset(F2_gpa$coords, 
                                   identifiers$Lake_morph)

vrel_F2_coords = Map(function(x) integration.Vrel(x), 
                       subset_F2_coords)

ASHN_compare = compare.ZVrel(vrel_F2_coords$ASHNC, 
                             vrel_F2_coords$ASHNW)

MYV_compare = compare.ZVrel(vrel_F2_coords$MYVC, 
                            vrel_F2_coords$MYVW)

SKR_compare = compare.ZVrel(vrel_F2_coords$SKRC, 
                            vrel_F2_coords$SKRW)

GTS_CSWY_compare = compare.ZVrel(vrel_F2_coords$CSWY, 
                                 vrel_F2_coords$GTS)



## integration with allometric variation removed
## standardized for allometry

subset_F2_coords_allo_adj = coords.subset(F2_allometry_adj_shape, 
                                            identifiers$Lake_morph)

vrel_F2_coords_allo_adj = Map(function(x) integration.Vrel(x), 
                                subset_F2_coords_allo_adj)

ASHN_noallo_compare = compare.ZVrel(vrel_F2_coords_allo_adj$ASHNC, 
                                    vrel_F2_coords_allo_adj$ASHNW)

MYV_noallo_compare = compare.ZVrel(vrel_F2_coords_allo_adj$MYVC, 
                                   vrel_F2_coords_allo_adj$MYVW)

SKR_noallo_compare = compare.ZVrel(vrel_F2_coords_allo_adj$SKRC, 
                                   vrel_F2_coords_allo_adj$SKRW)

GTS_CSWY_noallo_compare = compare.ZVrel(vrel_F2_coords_allo_adj$CSWY, 
                                        vrel_F2_coords_allo_adj$GTS)



# F2 craniofacial integration -------------------------------------------

## integration of craniofacial patterns with allometric variation included

F2_craniofacial = readland.tps('F2_Craniofacial_LM.TPS',
                                 specID = 'imageID')

identifiers = read_csv('F2_Metadata.CSV')

## superimposition on the entire dataset
F2_craniofacial_gpa = gpagen(F2_craniofacial,
                               print.progress = F)


subset_F2_craniofacial_coords = coords.subset(F2_craniofacial_gpa$coords,
                                                identifiers$Lake_morph)

vrel_F2_craniofacial = Map(function(x) integration.Vrel(x),
                             subset_F2_craniofacial_coords)

ASHN_craniofacial = compare.ZVrel(vrel_F2_craniofacial$ASHNC,
                                  vrel_F2_craniofacial$ASHNW)

MYV_craniofacial = compare.ZVrel(vrel_F2_craniofacial$MYVC,
                                 vrel_F2_craniofacial$MYVW)

SKR_craniofacial = compare.ZVrel(vrel_F2_craniofacial$SKRC,
                                 vrel_F2_craniofacial$SKRW)

GTS_CSWY_craniofacial = compare.ZVrel(vrel_F2_craniofacial$CSWY,
                                      vrel_F2_craniofacial$GTS)


# F2 craniofacial integration no allometry -------------------------------------------


F2_craniofacial = readland.tps('F2_Craniofacial_LM.TPS', 
                                 specID = 'imageID')

identifiers = read_csv('F2_metadata.csv') 

## superimposition on the entire dataset
F2_craniofacial_gpa = gpagen(F2_craniofacial, 
                               print.progress = F)


craniofacial_allometry = procD.lm(F2_craniofacial_gpa$coords ~ log(F2_craniofacial_gpa$Csize), 
                                  iter = 999, 
                                  RRPP = T)
summary(craniofacial_allometry)


craniofacial_resid = arrayspecs(craniofacial_allometry$residuals, 
                                p = dim(F2_craniofacial_gpa$coords)[1], 
                                k = dim(F2_craniofacial_gpa$coords)[2])
craniofacial_allo_adj = craniofacial_resid + array(F2_craniofacial_gpa$consensus, 
                                                   dim(craniofacial_resid))



F2_craniofacial_coords_noallo = coords.subset(craniofacial_allo_adj, 
                                                identifiers$Lake_morph)

vrel_F2_craniofacial_noallo = Map(function(x) integration.Vrel(x), 
                                    F2_craniofacial_coords_noallo)

ASHN_noallo_craniofacial = compare.ZVrel(vrel_F2_craniofacial_noallo$ASHNC, 
                                         vrel_F2_craniofacial_noallo$ASHNW)

MYV_noallo_craniofacial = compare.ZVrel(vrel_F2_craniofacial_noallo$MYVC, 
                                        vrel_F2_craniofacial_noallo$MYVW)

SKR_noallo_craniofacial = compare.ZVrel(vrel_F2_craniofacial_noallo$SKRC, 
                                        vrel_F2_craniofacial_noallo$SKRW)

GTS_CSWY_noallo_craniofacial = compare.ZVrel(vrel_F2_craniofacial_noallo$CSWY, 
                                             vrel_F2_craniofacial_noallo$GTS)




# Integration plasticity  ---------------------------------------------

F2_craniofacial = readland.tps('F2_Craniofacial_LM.TPS',
                                 specID = 'imageID')

identifiers = read_csv('F2_metadata.csv')

## superimposition on the entire dataset
F2_craniofacial_gpa = gpagen(F2_craniofacial,
                               print.progress = F)


subset_F2_craniofacial_coords = coords.subset(F2_craniofacial_gpa$coords,
                                                identifiers$Full_temp)

vrel_F2_craniofacial = Map(function(x) integration.Vrel(x),
                             subset_F2_craniofacial_coords)

Cold_Plasticity_craniofacial = compare.ZVrel(vrel_F2_craniofacial$`12@12`,
                                     vrel_F2_craniofacial$`12@18`)

Warm_Plasticity_craniofacial = compare.ZVrel(vrel_F2_craniofacial$`18@12`,
                                        vrel_F2_craniofacial$`18@18`)

Cold_off_Plasticity_craniofacial = compare.ZVrel(vrel_F2_craniofacial$`12@12`,
                                             vrel_F2_craniofacial$`18@12`)

Warm_off_Plasticity_craniofacial = compare.ZVrel(vrel_F2_craniofacial$`12@18`,
                                             vrel_F2_craniofacial$`18@18`)


## Analysis with no allometric variation

F2_craniofacial = readland.tps('F2_Craniofacial_LM.TPS', 
                                 specID = 'imageID')

identifiers = read_csv('F2_metadata.csv') 

## superimposition on the entire dataset
F2_craniofacial_gpa = gpagen(F2_craniofacial, 
                               print.progress = F)


craniofacial_allometry = procD.lm(F2_craniofacial_gpa$coords ~ log(F2_craniofacial_gpa$Csize), 
                                  iter = 999, 
                                  RRPP = T)
summary(craniofacial_allometry)


craniofacial_resid = arrayspecs(craniofacial_allometry$residuals, 
                                p = dim(F2_craniofacial_gpa$coords)[1], 
                                k = dim(F2_craniofacial_gpa$coords)[2])
craniofacial_allo_adj = craniofacial_resid + array(F2_craniofacial_gpa$consensus, 
                                                   dim(craniofacial_resid))



F2_craniofacial_coords_noallo = coords.subset(craniofacial_allo_adj, 
                                                identifiers$Full_temp)

vrel_F2_craniofacial_noallo = Map(function(x) integration.Vrel(x), 
                                    F2_craniofacial_coords_noallo)

NoAllo_Cold_Plasticity_craniofacial = compare.ZVrel(vrel_F2_craniofacial_noallo$`12@12`,
                                             vrel_F2_craniofacial_noallo$`12@18`)

NoAllo_Warm_Plasticity_craniofacial = compare.ZVrel(vrel_F2_craniofacial_noallo$`18@12`,
                                             vrel_F2_craniofacial_noallo$`18@18`)

NoAllo_Cold_off_Plasticity_craniofacial = compare.ZVrel(vrel_F2_craniofacial_noallo$`12@12`,
                                                 vrel_F2_craniofacial_noallo$`18@12`)

NoAllo_Warm_off_Plasticity_craniofacial = compare.ZVrel(vrel_F2_craniofacial_noallo$`12@18`,
                                                 vrel_F2_craniofacial_noallo$`18@18`)



# Going deeper into lake differences in plasticity -----------------------------------------------------

F2_craniofacial_gpa = gpagen(F2_craniofacial,
                             print.progress = F)

## This doesn't work for some reason, can't subset on multiple categories?
subset_F2_craniofacial_coords = coords.subset(F2_craniofacial_gpa$coords,
                                              identifiers$Full_temp)

subset_F2_craniofacial_coords = coords.subset(subset_F2_craniofacial_coords$coords,
                                              identifiers$Lake_morph)

vrel_F2_craniofacial = Map(function(x) integration.Vrel(x),
                           subset_F2_craniofacial_coords)

Cold_Plasticity_craniofacial = compare.ZVrel(vrel_F2_craniofacial$`12@12`,
                                             vrel_F2_craniofacial$`12@18`)

Warm_Plasticity_craniofacial = compare.ZVrel(vrel_F2_craniofacial$`18@12`,
                                             vrel_F2_craniofacial$`18@18`)

Cold_off_Plasticity_craniofacial = compare.ZVrel(vrel_F2_craniofacial$`12@12`,
                                                 vrel_F2_craniofacial$`18@12`)

Warm_off_Plasticity_craniofacial = compare.ZVrel(vrel_F2_craniofacial$`12@18`,
                                                 vrel_F2_craniofacial$`18@18`)


