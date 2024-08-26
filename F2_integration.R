
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
identifiers %>% 
  group_by(Lake_morph, 
           Parent_temp, 
           Offspring_temp) %>% 
  summarize(n = n()) 

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

## Integration between whole body shape between the ecotypes 
## in the F2 generation does not appear to be different

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

## Difference in craniofacial integration between ASHN morphs
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

## Difference in craniofacial integration between ASHN ecotypes
ASHN_noallo_craniofacial = compare.ZVrel(vrel_F2_craniofacial_noallo$ASHNC, 
                                         vrel_F2_craniofacial_noallo$ASHNW)

MYV_noallo_craniofacial = compare.ZVrel(vrel_F2_craniofacial_noallo$MYVC, 
                                        vrel_F2_craniofacial_noallo$MYVW)

SKR_noallo_craniofacial = compare.ZVrel(vrel_F2_craniofacial_noallo$SKRC, 
                                        vrel_F2_craniofacial_noallo$SKRW)

GTS_CSWY_noallo_craniofacial = compare.ZVrel(vrel_F2_craniofacial_noallo$CSWY, 
                                             vrel_F2_craniofacial_noallo$GTS)

## patterns of integration between F2 ecotypes hold up regardless of 
## allometry. Craniofacial integration only differs in ASHN ecotypes


# 4bar linkage integration with ALLOMETRY ---------------------------------
F2_4bar = readland.tps('F2_4bar_linkage.TPS', 
                       specID = 'imageID')

identifiers = read_csv('F2_metadata.csv') 

## superimposition on the entire dataset
F2_4bar_gpa = gpagen(F2_4bar, 
                     print.progress = F)


F2_4bar_coords = coords.subset(F2_4bar_gpa$coords, 
                                      identifiers$Lake_morph)

vrel_F2_4bar = Map(function(x) integration.Vrel(x), 
                          F2_4bar_coords)

ASHN_4bar = compare.ZVrel(vrel_F2_4bar$ASHNC, 
                                 vrel_F2_4bar$ASHNW)

MYV_4bar = compare.ZVrel(vrel_F2_4bar$MYVC, 
                                vrel_F2_4bar$MYVW)

SKR_4bar = compare.ZVrel(vrel_F2_4bar$SKRC, 
                                vrel_F2_4bar$SKRW)

GTS_CSWY_4bar = compare.ZVrel(vrel_F2_4bar$CSWY, 
                                     vrel_F2_4bar$GTS)

## 4bar linkage integration on different between GTS and CSWY

# 4-bar linkage integration NO ALLOMETRY -----------------------------------------------
F2_4bar = readland.tps('F2_4bar_linkage.TPS', 
                               specID = 'imageID')

identifiers = read_csv('F2_metadata.csv') 

## superimposition on the entire dataset
F2_4bar_gpa = gpagen(F2_4bar, 
                             print.progress = F)


F2_4bar_allometry = procD.lm(F2_4bar_gpa$coords ~ log(F2_4bar_gpa$Csize), 
                                  iter = 999, 
                                  RRPP = T)
summary(F2_4bar_allometry)


F2_4bar_resid = arrayspecs(F2_4bar_allometry$residuals, 
                                p = dim(F2_4bar_gpa$coords)[1], 
                                k = dim(F2_4bar_gpa$coords)[2])
F2_4bar_allo_adj = F2_4bar_resid + array(F2_4bar_gpa$consensus, 
                                                   dim(F2_4bar_resid))



F2_4bar_coords_noallo = coords.subset(F2_4bar_allo_adj, 
                                              identifiers$Lake_morph)

vrel_F2_4bar_noallo = Map(function(x) integration.Vrel(x), 
                                  F2_4bar_coords_noallo)

ASHN_noallo_4bar = compare.ZVrel(vrel_F2_4bar_noallo$ASHNC, 
                                         vrel_F2_4bar_noallo$ASHNW)

MYV_noallo_4bar = compare.ZVrel(vrel_F2_4bar_noallo$MYVC, 
                                        vrel_F2_4bar_noallo$MYVW)

SKR_noallo_4bar = compare.ZVrel(vrel_F2_4bar_noallo$SKRC, 
                                        vrel_F2_4bar_noallo$SKRW)

GTS_CSWY_noallo_4bar = compare.ZVrel(vrel_F2_4bar_noallo$CSWY, 
                                             vrel_F2_4bar_noallo$GTS)

## Getting rid of allometry got rid of the differences in integration
## between GTS and CSWY


# 4bar integration plasticity ---------------------------------------------

F2_4bar = readland.tps('F2_4bar_linkage.TPS', 
                       specID = 'imageID')

identifiers = read_csv('F2_metadata.csv') 

## superimposition on the entire dataset
F2_4bar_gpa = gpagen(F2_4bar, 
                     print.progress = F)

F2_4bar_allometry = procD.lm(F2_4bar_gpa$coords ~ log(F2_4bar_gpa$Csize), 
                             iter = 999, 
                             RRPP = T)
summary(F2_4bar_allometry)


F2_4bar_resid = arrayspecs(F2_4bar_allometry$residuals, 
                           p = dim(F2_4bar_gpa$coords)[1], 
                           k = dim(F2_4bar_gpa$coords)[2])
F2_4bar_allo_adj = F2_4bar_resid + array(F2_4bar_gpa$consensus, 
                                         dim(F2_4bar_resid))

F2_4bar_coords = coords.subset(F2_4bar_allo_adj, 
                               identifiers$Full_temp)

# F2_4bar_coords = coords.subset(F2_4bar_gpa$coords, 
#                                identifiers$Full_temp)


vrel_F2_4bar = Map(function(x) integration.Vrel(x),
                           F2_4bar_coords)

Cold_Plasticity_4bar = compare.ZVrel(vrel_F2_4bar$`12@12`,
                                             vrel_F2_4bar$`12@18`)

Warm_Plasticity_4bar = compare.ZVrel(vrel_F2_4bar$`18@12`,
                                             vrel_F2_4bar$`18@18`)

Cold_off_Plasticity_4bar = compare.ZVrel(vrel_F2_4bar$`12@12`,
                                                 vrel_F2_4bar$`18@12`)

Warm_off_Plasticity_4bar = compare.ZVrel(vrel_F2_4bar$`12@18`,
                                                 vrel_F2_4bar$`18@18`)



# 4bar integration interlandmark dist -------------------------------------

# F2_4bar = readland.tps('F2_4bar_linkage.TPS', 
#                        specID = 'imageID')

F2_4bar = readland.tps('F2_No_GT.TPS', 
                      specID = 'imageID')

identifiers = read_csv('F2_metadata.csv') 

## superimposition on the entire dataset
F2_4bar_gpa = gpagen(F2_4bar,
                     print.progress = F)

# lmks = matrix(c(8,24,8,27,23,27,23,24), 
#               ncol = 2, 
#               byrow = T, 
#               dimnames = list(c('824_dist', '827_dist', '2327_dist', '2324_dist'), c('start', 'end')))

lmks = data.frame(dist1_824 = c(8, 24), 
                  dist2_827 = c(8, 27), 
                  dist3_2327 = c(23, 27), 
                  dist4_2324 = c(23, 24), 
                  row.names = c('start', 
                                'end'))

A = F2_4bar_gpa$coords
F2_4bar_lineardist = interlmkdist(A, 
                                  lmks)

arrayspecs(F2_4bar_lineardist, 
           4, 
           3)

## linear distances not in 3d array
F2_4bar_dist = coords.subset(F2_4bar_lineardist,
                             identifiers$Full_temp)

vrel_F2_craniofacial = Map(function(x) integration.Vrel(x),
                           subset_F2_craniofacial_coords)

# Craniofacial Integration plasticity  ---------------------------------------------

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

## integration differs between F2 12 and 18 while parents from 12 degrees
Cold_Plasticity_craniofacial = compare.ZVrel(vrel_F2_craniofacial$`12@12`,
                                     vrel_F2_craniofacial$`12@18`)

## Integration differs between F2 12 and 18 while parents from 18
Warm_Plasticity_craniofacial = compare.ZVrel(vrel_F2_craniofacial$`18@12`,
                                        vrel_F2_craniofacial$`18@18`)

## Effect size of plastic response is greater when parents are from
## 12 degrees and not 18 degrees


## When parents are from different temps and F2 at 12, 
## there is no difference in integration
Cold_off_Plasticity_craniofacial = compare.ZVrel(vrel_F2_craniofacial$`12@12`,
                                             vrel_F2_craniofacial$`18@12`)

## When parents are from different temps and F2 at 12, 
## there is no difference in integration
Warm_off_Plasticity_craniofacial = compare.ZVrel(vrel_F2_craniofacial$`12@18`,
                                             vrel_F2_craniofacial$`18@18`)

compare.ZVrel(vrel_F2_craniofacial$`12@18`,
              vrel_F2_craniofacial$`18@12`)

compare.ZVrel(vrel_F2_craniofacial$`12@12`,
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

compare.ZVrel(vrel_F2_craniofacial_noallo$`12@18`,
              vrel_F2_craniofacial_noallo$`18@12`)

compare.ZVrel(vrel_F2_craniofacial_noallo$`12@12`,
              vrel_F2_craniofacial_noallo$`18@18`)

## When allometry is removed, the pattern of integration plasticity stays
## When parental temps are from 12 or 18 and the offspring are raised
## under different temps then the patterns of integration are different

## If parents are from different temperatures and offspring are raised 
## under the same temperature, then there is no difference in integration
## Essentially the parental temperature doesn't seem to be doing to much

# Ecotype plasticity and transgenerational integration -----------------------------------------------------

F2_craniofacial = readland.tps('F2_Craniofacial_LM.TPS',
                               specID = 'imageID')

identifiers = read_csv('F2_metadata.csv') %>% 
  unite('Lake_Morph_Full_Temp', 
        Lake_morph, 
        Full_temp, 
        sep = '_', 
        remove = F)


F2_craniofacial_gpa = gpagen(F2_craniofacial,
                             print.progress = F)

subset_F2_craniofacial_coords = coords.subset(F2_craniofacial_gpa$coords,
                                              identifiers$Lake_Morph_Full_Temp)

vrel_F2_craniofacial = Map(function(x) integration.Vrel(x),
                           subset_F2_craniofacial_coords)

ASHNC_12_plasticity = compare.ZVrel(vrel_F2_craniofacial$`ASHNC_12@12`,
                                             vrel_F2_craniofacial$`ASHNC_12@18`)

ASHNC_18_plasticity = compare.ZVrel(vrel_F2_craniofacial$`ASHNC_18@12`,
                                    vrel_F2_craniofacial$`ASHNC_18@18`)

ASHNC_12_Trans = compare.ZVrel(vrel_F2_craniofacial$`ASHNC_12@12`,
                                    vrel_F2_craniofacial$`ASHNC_18@12`)

ASHNC_18_trans = compare.ZVrel(vrel_F2_craniofacial$`ASHNC_12@18`,
                                    vrel_F2_craniofacial$`ASHNC_18@18`)


ASHNW_12_plasticity = compare.ZVrel(vrel_F2_craniofacial$`ASHNW_12@12`,
                                    vrel_F2_craniofacial$`ASHNW_12@18`)

ASHNW_18_plasticity = compare.ZVrel(vrel_F2_craniofacial$`ASHNW_18@12`,
                                    vrel_F2_craniofacial$`ASHNW_18@18`)

ASHNW_12_Trans = compare.ZVrel(vrel_F2_craniofacial$`ASHNW_12@12`,
                               vrel_F2_craniofacial$`ASHNW_18@12`)

ASHNW_18_trans = compare.ZVrel(vrel_F2_craniofacial$`ASHNW_12@18`,
                               vrel_F2_craniofacial$`ASHNW_18@18`)


SKRC_12_plasticity = compare.ZVrel(vrel_F2_craniofacial$`SKRC_12@12`,
                                    vrel_F2_craniofacial$`SKRC_12@18`)

SKRC_18_plasticity = compare.ZVrel(vrel_F2_craniofacial$`SKRC_18@12`,
                                    vrel_F2_craniofacial$`SKRC_18@18`)

SKRC_12_Trans = compare.ZVrel(vrel_F2_craniofacial$`SKRC_12@12`,
                               vrel_F2_craniofacial$`SKRC_18@12`)

SKRC_18_trans = compare.ZVrel(vrel_F2_craniofacial$`SKRC_12@18`,
                               vrel_F2_craniofacial$`SKRC_18@18`)


SKRW_12_plasticity = compare.ZVrel(vrel_F2_craniofacial$`SKRW_12@12`,
                                   vrel_F2_craniofacial$`SKRW_12@18`)

SKRW_18_plasticity = compare.ZVrel(vrel_F2_craniofacial$`SKRW_18@12`,
                                   vrel_F2_craniofacial$`SKRW_18@18`)

SKRW_12_Trans = compare.ZVrel(vrel_F2_craniofacial$`SKRW_12@12`,
                              vrel_F2_craniofacial$`SKRW_18@12`)

SKRW_18_trans = compare.ZVrel(vrel_F2_craniofacial$`SKRW_12@18`,
                              vrel_F2_craniofacial$`SKRW_18@18`)


MYVC_12_plasticity = compare.ZVrel(vrel_F2_craniofacial$`MYVC_12@12`,
                                   vrel_F2_craniofacial$`MYVC_12@18`)

MYVC_18_plasticity = compare.ZVrel(vrel_F2_craniofacial$`MYVC_18@12`,
                                   vrel_F2_craniofacial$`MYVC_18@18`)

MYVC_12_Trans = compare.ZVrel(vrel_F2_craniofacial$`MYVC_12@12`,
                              vrel_F2_craniofacial$`MYVC_18@12`)

MYVC_18_trans = compare.ZVrel(vrel_F2_craniofacial$`MYVC_12@18`,
                              vrel_F2_craniofacial$`MYVC_18@18`)

MYVW_12_plasticity = compare.ZVrel(vrel_F2_craniofacial$`MYVW_12@12`,
                                   vrel_F2_craniofacial$`MYVW_12@18`)

MYVW_18_plasticity = compare.ZVrel(vrel_F2_craniofacial$`MYVW_18@12`,
                                   vrel_F2_craniofacial$`MYVW_18@18`)

MYVW_12_Trans = compare.ZVrel(vrel_F2_craniofacial$`MYVW_12@12`,
                              vrel_F2_craniofacial$`MYVW_18@12`)

MYVW_18_trans = compare.ZVrel(vrel_F2_craniofacial$`MYVW_12@18`,
                              vrel_F2_craniofacial$`MYVW_18@18`)


GTS_12_plasticity = compare.ZVrel(vrel_F2_craniofacial$`GTSW_12@12`,
                                   vrel_F2_craniofacial$`GTSW_12@18`)

GTS_18_plasticity = compare.ZVrel(vrel_F2_craniofacial$`GTSW_18@12`,
                                   vrel_F2_craniofacial$`GTSW_18@18`)

GTS_12_Trans = compare.ZVrel(vrel_F2_craniofacial$`GTSW_12@12`,
                              vrel_F2_craniofacial$`GTSW_18@12`)

GTS_18_trans = compare.ZVrel(vrel_F2_craniofacial$`GTSW_12@18`,
                              vrel_F2_craniofacial$`GTSW_18@18`)

CSWYC_12_plasticity = compare.ZVrel(vrel_F2_craniofacial$`CSWYC_12@12`,
                                   vrel_F2_craniofacial$`CSWYC_12@18`)

CSWYC_18_plasticity = compare.ZVrel(vrel_F2_craniofacial$`CSWYC_18@12`,
                                   vrel_F2_craniofacial$`CSWYC_18@18`)

CSWYC_12_Trans = compare.ZVrel(vrel_F2_craniofacial$`CSWYC_12@12`,
                              vrel_F2_craniofacial$`CSWYC_18@12`)

CSWYC_18_trans = compare.ZVrel(vrel_F2_craniofacial$`CSWYC_12@18`,
                              vrel_F2_craniofacial$`CSWYC_18@18`)

## There appears to be no differences in integration between 
## offspring raised at different temps or if their parents were
## raised at different temps. 


# F2 Integration ----------------------------------------------------------

F2_tps = readland.tps('F2_No_GT.TPS', 
                      specID = 'imageID')

identifiers = read_csv('F2_Metadata.CSV', 
                       col_names = T) 

## superimposition on the entire dataset
F2_gpa = gpagen(F2_tps, 
                print.progress = F)

subset_F2_coords = coords.subset(F2_gpa$coords,
                                              identifiers$lake_morph_Pair_Full_Temp)

vrel_F2_ = Map(function(x) integration.Vrel(x),
                           subset_F2_coords)

compare.ZVrel(vrel_F2$`ASHNC_12@12`,
              vrel_F2$`ASHNW_12@12`)

compare.ZVrel(vrel_F2_$`ASHNC_12@18`,
              vrel_F2_$`ASHNW_12@18`)

compare.ZVrel(vrel_F2_$`ASHNC_18@18`,
              vrel_F2_$`ASHNW_18@18`)

compare.ZVrel(vrel_F2_$`ASHNC_18@12`,
              vrel_F2_$`ASHNW_18@12`)


compare.ZVrel(vrel_F2_$`SKRC_12@12`,
              vrel_F2_$`SKRW_12@12`)

compare.ZVrel(vrel_F2_$`SKRC_12@18`,
              vrel_F2_$`SKRW_12@18`)

compare.ZVrel(vrel_F2_$`SKRC_18@18`,
              vrel_F2_$`SKRW_18@18`)

compare.ZVrel(vrel_F2_$`SKRC_18@12`,
              vrel_F2_$`SKRW_18@12`)


compare.ZVrel(vrel_F2_$`MYVC_12@12`,
              vrel_F2_$`MYVW_12@12`)

compare.ZVrel(vrel_F2_$`MYVC_12@18`,
              vrel_F2_$`MYVW_12@18`)

compare.ZVrel(vrel_F2_$`MYVC_18@18`,
              vrel_F2_$`MYVW_18@18`)

compare.ZVrel(vrel_F2_$`MYVC_18@12`,
              vrel_F2_$`MYVW_18@12`)


compare.ZVrel(vrel_F2_$`GTSW_12@12`,
              vrel_F2_$`CSWYC_12@12`)

compare.ZVrel(vrel_F2_$`GTSW_12@18`,
              vrel_F2_$`CSWYC_12@18`)

compare.ZVrel(vrel_F2_$`GTSW_18@18`,
              vrel_F2_$`CSWYC_18@18`)

compare.ZVrel(vrel_F2_$`GTSW_18@12`,
              vrel_F2_$`CSWYC_18@12`)

# Lake plasticity and transgenerational integration -----------------------


F2_craniofacial = readland.tps('F2_Craniofacial_LM.TPS',
                               specID = 'imageID')

identifiers = read_csv('F2_metadata.csv') %>%
  unite('Ecotype_Pair_Full_Temp',
        Ecotype_pair,
        Full_temp,
        sep = '_',
        remove = F)

identifiers = read_csv('F2_metadata.csv') %>%
  unite('lake_morph_Pair_Full_Temp',
        Lake_morph,
        Full_temp,
        sep = '_',
        remove = F)



subset_F2_craniofacial_coords = coords.subset(F2_craniofacial_gpa$coords,
                                              identifiers$lake_morph_Pair_Full_Temp)

vrel_F2_craniofacial = Map(function(x) integration.Vrel(x),
                           subset_F2_craniofacial_coords)

compare.ZVrel(vrel_F2_craniofacial$`ASHNC_12@12`,
              vrel_F2_craniofacial$`ASHNW_12@12`)

compare.ZVrel(vrel_F2_craniofacial$`ASHNC_12@18`,
              vrel_F2_craniofacial$`ASHNW_12@18`)

compare.ZVrel(vrel_F2_craniofacial$`ASHNC_18@18`,
              vrel_F2_craniofacial$`ASHNW_18@18`)

compare.ZVrel(vrel_F2_craniofacial$`ASHNC_18@12`,
              vrel_F2_craniofacial$`ASHNW_18@12`)


compare.ZVrel(vrel_F2_craniofacial$`SKRC_12@12`,
              vrel_F2_craniofacial$`SKRW_12@12`)

compare.ZVrel(vrel_F2_craniofacial$`SKRC_12@18`,
              vrel_F2_craniofacial$`SKRW_12@18`)

compare.ZVrel(vrel_F2_craniofacial$`SKRC_18@18`,
              vrel_F2_craniofacial$`SKRW_18@18`)

compare.ZVrel(vrel_F2_craniofacial$`SKRC_18@12`,
              vrel_F2_craniofacial$`SKRW_18@12`)


compare.ZVrel(vrel_F2_craniofacial$`MYVC_12@12`,
              vrel_F2_craniofacial$`MYVW_12@12`)

compare.ZVrel(vrel_F2_craniofacial$`MYVC_12@18`,
              vrel_F2_craniofacial$`MYVW_12@18`)

compare.ZVrel(vrel_F2_craniofacial$`MYVC_18@18`,
              vrel_F2_craniofacial$`MYVW_18@18`)

compare.ZVrel(vrel_F2_craniofacial$`MYVC_18@12`,
              vrel_F2_craniofacial$`MYVW_18@12`)


compare.ZVrel(vrel_F2_craniofacial$`GTSW_12@12`,
              vrel_F2_craniofacial$`CSWYC_12@12`)

compare.ZVrel(vrel_F2_craniofacial$`GTSW_12@18`,
              vrel_F2_craniofacial$`CSWYC_12@18`)

compare.ZVrel(vrel_F2_craniofacial$`GTSW_18@18`,
              vrel_F2_craniofacial$`CSWYC_18@18`)

compare.ZVrel(vrel_F2_craniofacial$`GTSW_18@12`,
              vrel_F2_craniofacial$`CSWYC_18@12`)


ASHN_18_plasticity = compare.ZVrel(vrel_F2_craniofacial$`ASHN_18@12`,
                                   vrel_F2_craniofacial$`ASHN_18@18`)

ASHN_12_Trans = compare.ZVrel(vrel_F2_craniofacial$`ASHN_12@12`,
                              vrel_F2_craniofacial$`ASHN_18@12`)

ASHN_18_trans = compare.ZVrel(vrel_F2_craniofacial$`ASHN_12@18`,
                              vrel_F2_craniofacial$`ASHN_18@18`)


compare.ZVrel(vrel_F2_craniofacial_noallo$`12@18`,
              vrel_F2_craniofacial_noallo$`18@12`)

compare.ZVrel(vrel_F2_craniofacial$`12@12`,
              vrel_F2_craniofacial$`18@18`)



subset_F2_craniofacial_coords = coords.subset(F2_craniofacial_gpa$coords,
                                              identifiers$Ecotype_Pair_Full_Temp)

vrel_F2_craniofacial = Map(function(x) integration.Vrel(x),
                           subset_F2_craniofacial_coords)

## ASHN shows a high degree of plasticity in the F2 generation
## and no effect of the partental F1 generation on integration
ASHN_12_plasticity = compare.ZVrel(vrel_F2_craniofacial$`ASHN_12@12`,
                                    vrel_F2_craniofacial$`ASHN_12@18`)

ASHN_18_plasticity = compare.ZVrel(vrel_F2_craniofacial$`ASHN_18@12`,
                                    vrel_F2_craniofacial$`ASHN_18@18`)

ASHN_12_Trans = compare.ZVrel(vrel_F2_craniofacial$`ASHN_12@12`,
                               vrel_F2_craniofacial$`ASHN_18@12`)

ASHN_18_trans = compare.ZVrel(vrel_F2_craniofacial$`ASHN_12@18`,
                               vrel_F2_craniofacial$`ASHN_18@18`)


compare.ZVrel(vrel_F2_craniofacial_noallo$`12@18`,
              vrel_F2_craniofacial_noallo$`18@12`)

compare.ZVrel(vrel_F2_craniofacial$`12@12`,
              vrel_F2_craniofacial$`18@18`)

## ASHN shows the pattern of integration that we've seen in the F2
## with all populations combined. Patterns of craniofacial integration
## changes when the ecotypes are raised under different temperatures 
## but when the parental temps differ there is no difference in integration


SKR_12_plasticity = compare.ZVrel(vrel_F2_craniofacial$`SKR_12@12`,
                                   vrel_F2_craniofacial$`SKR_12@18`)

SKR_18_plasticity = compare.ZVrel(vrel_F2_craniofacial$`SKR_18@12`,
                                   vrel_F2_craniofacial$`SKR_18@18`)

SKR_12_Trans = compare.ZVrel(vrel_F2_craniofacial$`SKR_12@12`,
                              vrel_F2_craniofacial$`SKR_18@12`)

SKR_18_trans = compare.ZVrel(vrel_F2_craniofacial$`SKR_12@18`,
                              vrel_F2_craniofacial$`SKR_18@18`)



MYV_12_plasticity = compare.ZVrel(vrel_F2_craniofacial$`MYV_12@12`,
                                   vrel_F2_craniofacial$`MYV_12@18`)

MYV_18_plasticity = compare.ZVrel(vrel_F2_craniofacial$`MYV_18@12`,
                                   vrel_F2_craniofacial$`MYV_18@18`)

MYV_12_Trans = compare.ZVrel(vrel_F2_craniofacial$`MYV_12@12`,
                              vrel_F2_craniofacial$`MYV_18@12`)

MYV_18_trans = compare.ZVrel(vrel_F2_craniofacial$`MYV_12@18`,
                              vrel_F2_craniofacial$`MYV_18@18`)


GTS_CSWY_12_plasticity = compare.ZVrel(vrel_F2_craniofacial$`GTS_CSWY_12@12`,
                                  vrel_F2_craniofacial$`GTS_CSWY_12@18`)

GTS_CSWY_18_plasticity = compare.ZVrel(vrel_F2_craniofacial$`GTS_CSWY_18@12`,
                                  vrel_F2_craniofacial$`GTS_CSWY_18@18`)

GTS_CSWY_12_Trans = compare.ZVrel(vrel_F2_craniofacial$`GTS_CSWY_12@12`,
                             vrel_F2_craniofacial$`GTS_CSWY_18@12`)

GTS_CSWY_18_trans = compare.ZVrel(vrel_F2_craniofacial$`GTS_CSWY_12@18`,
                             vrel_F2_craniofacial$`GTS_CSWY_18@18`)


## ASHN appears to be driving the differences in craniofacial integration
## plasticity that we've been seeing in the F2 generation
## There is differences in craniofacial integration when the current 
## generation is raised at different temperatures 
## But there doesn't appear to be an influence of parental generation
## on patterns of craniofacial integration
## Or even Grandparent effect of integration when you go into the specific
## ecotype differences. 

## Integration appears to be a plastic trait in ASHN only and 
## only really due to the current temperatures they're exposed to
## transgenerational effects appear to have minor involvement
## However this needs to be tested in more detail
## This is what the meeting is about, how do we actually test this? 