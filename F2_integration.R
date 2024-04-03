
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

