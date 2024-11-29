##############################
## Jaw kinetics (Corins code)
##
## Matt Brachmann (PhDMattyB)
##
## 28.11.2024
##
##############################

setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/Updated Landmarks/')

library(geomorph)
library(tidyverse)
library(linkR)


## landmarks read during the wild_pattern script
## TGP and WPG effects isolated in plasticity_analysis script
## All LMs have been put through a common gpa
## they should be good to go straight into this script


# Wild kinetics -----------------------------------------------------------

###Calculate transmission coefficients for Anterior/Pre-Maxillary linkage (23,27,37,3)
wild_PM_LMs <- c(23,27,28,3)
wild_coords = wild_lmk_coords[PM_LMs,,]


#Create empty table for data
wild_KT_PreMax <- numeric(length=331)
wild_PreMax_Rotation <- numeric(length=331)
for (i in 1:331) {
  #Define linkages, joints etc.
  PM_ind <- wild_coords[,,i]
  lms <- PM_ind
  joint.types <- c('R','S','S','R')
  joint.conn <- rbind(c(0,1),c(1,2),c(2,3),c(3,0))
  joint.cons <- list(NA,NA,NA,NA)
  
  linkage <- defineLinkage(joint.coor = lms, joint.types = joint.types,
                           joint.cons = joint.cons, joint.conn=joint.conn)
  #Overcome error message
  linkage2 <- linkage
  linkage2$joint.cons <- list(c(0,0,1), NA, NA, c(0,0,1))  
  
  #Animate linkage, calculate transmission coefficients, add to table
  drawLinkage(linkage2, animate=FALSE)
  animate <- animateLinkage(linkage2, input.param=c(0,pi/12, length = 15), input.joint=1)
  kine <- linkageKinematics(animate)
  KT <- abs(kine$links.rdis.d['Link2',] / kine$links.rdis.d['Link1',])
  KT_PreMax[i] <- KT[2]
  KT_PM_df <- as.data.frame(KT_PreMax)
  
  #Calculate Output rotation
  PM_Output <- print(kine$links.rdis.d['Link2',])
  PreMax_Rotation[i] <- PM_Output[2]
  PM_Rotation_DF <- as.data.frame(PreMax_Rotation)
}

Pheno_4Bar = cbind(KT_PM_df, PM_Rotation_DF)



###Calculate transmission coefficients for Opercular Linkage (27,23,24,25)

#Define landmarks
OP_LMs <- c(27,23,24,25)
OP_coords <- wild_lmk_coords[OP_LMs,,]

##Calculate transmission coefficients, create dataframe

#Create empty table for data
wild_KT_Opercular <- numeric(length=331)
wild_KT_Output_Rotation <- numeric(length=331)
for (i in 1:331) {
  #Define linkages, joints etc.
  OP_ind <- OP_coords[,,i]
  lms <- OP_ind
  joint.types <- c('R','S','S','R')
  joint.conn <- rbind(c(0,1),c(1,2),c(2,3),c(3,0))
  joint.cons <- list(NA,NA,NA,NA)
  
  linkage <- defineLinkage(joint.coor = lms, joint.types = joint.types,
                           joint.cons = joint.cons, joint.conn=joint.conn)
  
  #Overcome error message
  linkage2 <- linkage
  linkage2$joint.cons <- list(c(0,0,1), NA, NA, c(0,0,1))
  
  #Animate linkage, calculate transmission coefficients, add to table
  drawLinkage(linkage2)
  animate <- animateLinkage(linkage2, input.param=c(0,pi/12, length = 15), input.joint=1)
  kine <- linkageKinematics(animate)
  KT <- abs(kine$links.rdis.d['Link1',] / kine$links.rdis.d['Link3',])
  KT_Opercular[i] <- KT[2]
  KT_df <- as.data.frame(KT_Opercular)
  
  #Calculate Output rotation
  OP_Output <- print(kine$links.rdis.d['Link1',])
  KT_Output_Rotation[i] <- OP_Output[2]
  OP_Rotation_DF <- as.data.frame(KT_Output_Rotation)
}

operculum_df = cbind(KT_df, OP_Rotation_DF)


wild_Pheno_4Bar <- cbind(Pheno_4Bar, operculum_df)
colnames(wild_Pheno_4Bar) <- c("PreMax_KT", "PreMax_Rotation",
                          "Opercular_KT", "Opercular_Rotation")

wild_Pheno_4Bar %>% 
  as_tibble() %>% 
  write_csv('Wild_Jaw_kinetic_traits.csv')



# F2 uncorrected traits ---------------------------------------------------

###Calculate transmission coefficients for Anterior/Pre-Maxillary linkage (23,27,37,3)
PM_LMs <- c(23,27,28,3)
F2_PM_coords = F2_lmk_coords[PM_LMs,,]


#Create empty table for data
KT_PreMax <- numeric(length=931)
PreMax_Rotation <- numeric(length=931)
for (i in 1:931) {
  #Define linkages, joints etc.
  PM_ind <- F2_PM_coords[,,i]
  lms <- PM_ind
  joint.types <- c('R','S','S','R')
  joint.conn <- rbind(c(0,1),c(1,2),c(2,3),c(3,0))
  joint.cons <- list(NA,NA,NA,NA)
  
  linkage <- defineLinkage(joint.coor = lms, joint.types = joint.types,
                           joint.cons = joint.cons, joint.conn=joint.conn)
  #Overcome error message
  linkage2 <- linkage
  linkage2$joint.cons <- list(c(0,0,1), NA, NA, c(0,0,1))  
  
  #Animate linkage, calculate transmission coefficients, add to table
  drawLinkage(linkage2, animate=FALSE)
  animate <- animateLinkage(linkage2, input.param=c(0,pi/12, length = 15), input.joint=1)
  kine <- linkageKinematics(animate)
  KT <- abs(kine$links.rdis.d['Link2',] / kine$links.rdis.d['Link1',])
  KT_PreMax[i] <- KT[2]
  KT_PM_df <- as.data.frame(KT_PreMax)
  
  #Calculate Output rotation
  PM_Output <- print(kine$links.rdis.d['Link2',])
  PreMax_Rotation[i] <- PM_Output[2]
  PM_Rotation_DF <- as.data.frame(PreMax_Rotation)
}

F2_Pheno_4Bar = cbind(KT_PM_df, PM_Rotation_DF)



###Calculate transmission coefficients for Opercular Linkage (27,23,24,25)

#Define landmarks
OP_LMs <- c(27,23,24,25)
OP_coords <- F2_lmk_coords[OP_LMs,,]

##Calculate transmission coefficients, create dataframe

#Create empty table for data
KT_Opercular <- numeric(length=931)
KT_Output_Rotation <- numeric(length=931)
for (i in 1:931) {
  #Define linkages, joints etc.
  OP_ind <- OP_coords[,,i]
  lms <- OP_ind
  joint.types <- c('R','S','S','R')
  joint.conn <- rbind(c(0,1),c(1,2),c(2,3),c(3,0))
  joint.cons <- list(NA,NA,NA,NA)
  
  linkage <- defineLinkage(joint.coor = lms, joint.types = joint.types,
                           joint.cons = joint.cons, joint.conn=joint.conn)
  
  #Overcome error message
  linkage2 <- linkage
  linkage2$joint.cons <- list(c(0,0,1), NA, NA, c(0,0,1))
  
  #Animate linkage, calculate transmission coefficients, add to table
  drawLinkage(linkage2)
  animate <- animateLinkage(linkage2, input.param=c(0,pi/12, length = 15), input.joint=1)
  kine <- linkageKinematics(animate)
  KT <- abs(kine$links.rdis.d['Link1',] / kine$links.rdis.d['Link3',])
  KT_Opercular[i] <- KT[2]
  KT_df <- as.data.frame(KT_Opercular)
  
  #Calculate Output rotation
  OP_Output <- print(kine$links.rdis.d['Link1',])
  KT_Output_Rotation[i] <- OP_Output[2]
  OP_Rotation_DF <- as.data.frame(KT_Output_Rotation)
}

F2_operculum_4bar = cbind(KT_df, OP_Rotation_DF)


F2_Pheno_4Bar_full <- cbind(F2_Pheno_4Bar, F2_operculum_4bar)
colnames(F2_Pheno_4Bar_full) <- c("PreMax_KT", "PreMax_Rotation",
                          "Opercular_KT", "Opercular_Rotation")

F2_Pheno_4Bar_full %>% 
  as_tibble() %>% 
  write_csv('F2_uncorrected_Jaw_kinetic_traits.csv')


# TGP effects jaw kinetics ------------------------------------------------

###Calculate transmission coefficients for Anterior/Pre-Maxillary linkage (23,27,37,3)
PM_LMs <- c(23,27,28,3)
TGP_PM_coords = F2_parent_array_consensus[PM_LMs,,]


#Create empty table for data
KT_PreMax <- numeric(length=931)
PreMax_Rotation <- numeric(length=931)
for (i in 1:931) {
  #Define linkages, joints etc.
  PM_ind <- TGP_PM_coords[,,i]
  lms <- PM_ind
  joint.types <- c('R','S','S','R')
  joint.conn <- rbind(c(0,1),c(1,2),c(2,3),c(3,0))
  joint.cons <- list(NA,NA,NA,NA)
  
  linkage <- defineLinkage(joint.coor = lms, joint.types = joint.types,
                           joint.cons = joint.cons, joint.conn=joint.conn)
  #Overcome error message
  linkage2 <- linkage
  linkage2$joint.cons <- list(c(0,0,1), NA, NA, c(0,0,1))  
  
  #Animate linkage, calculate transmission coefficients, add to table
  drawLinkage(linkage2, animate=FALSE)
  animate <- animateLinkage(linkage2, input.param=c(0,pi/12, length = 15), input.joint=1)
  kine <- linkageKinematics(animate)
  KT <- abs(kine$links.rdis.d['Link2',] / kine$links.rdis.d['Link1',])
  KT_PreMax[i] <- KT[2]
  KT_PM_df <- as.data.frame(KT_PreMax)
  
  #Calculate Output rotation
  PM_Output <- print(kine$links.rdis.d['Link2',])
  PreMax_Rotation[i] <- PM_Output[2]
  PM_Rotation_DF <- as.data.frame(PreMax_Rotation)
}

TGP_Pheno_4Bar = cbind(KT_PM_df, PM_Rotation_DF)



###Calculate transmission coefficients for Opercular Linkage (27,23,24,25)

#Define landmarks
OP_LMs <- c(27,23,24,25)
TGP_OP_coords <- F2_parent_array_consensus[OP_LMs,,]

##Calculate transmission coefficients, create dataframe

#Create empty table for data
KT_Opercular <- numeric(length=931)
KT_Output_Rotation <- numeric(length=931)
for (i in 1:931) {
  #Define linkages, joints etc.
  OP_ind <- TGP_OP_coords[,,i]
  lms <- OP_ind
  joint.types <- c('R','S','S','R')
  joint.conn <- rbind(c(0,1),c(1,2),c(2,3),c(3,0))
  joint.cons <- list(NA,NA,NA,NA)
  
  linkage <- defineLinkage(joint.coor = lms, joint.types = joint.types,
                           joint.cons = joint.cons, joint.conn=joint.conn)
  
  #Overcome error message
  linkage2 <- linkage
  linkage2$joint.cons <- list(c(0,0,1), NA, NA, c(0,0,1))
  
  #Animate linkage, calculate transmission coefficients, add to table
  drawLinkage(linkage2)
  animate <- animateLinkage(linkage2, input.param=c(0,pi/12,length = 15), input.joint=1)
  kine <- linkageKinematics(animate)
  KT <- abs(kine$links.rdis.d['Link1',] / kine$links.rdis.d['Link3',])
  KT_Opercular[i] <- KT[2]
  KT_df <- as.data.frame(KT_Opercular)
  
  #Calculate Output rotation
  OP_Output <- print(kine$links.rdis.d['Link1',])
  KT_Output_Rotation[i] <- OP_Output[2]
  OP_Rotation_DF <- as.data.frame(KT_Output_Rotation)
}

TGP_operculum = cbind(KT_df, OP_Rotation_DF)

TGP_Pheno_4Bar_full <- cbind(TGP_Pheno_4Bar,TGP_operculum)
colnames(TGP_Pheno_4Bar_full) <- c("PreMax_KT", "PreMax_Rotation",
                             "Opercular_KT", "Opercular_Rotation")

TGP_Pheno_4Bar_full %>% 
  as_tibble() %>% 
  write_csv('TGP_Jaw_kinetic_traits.csv')


# WGP effects jaw kinetics ------------------------------------------------

###Calculate transmission coefficients for Anterior/Pre-Maxillary linkage (23,27,37,3)
PM_LMs <- c(23,27,28,3)
WGP_PM_coords = F2_array_consensus[PM_LMs,,]


#Create empty table for data
KT_PreMax <- numeric(length=931)
PreMax_Rotation <- numeric(length=931)
for (i in 1:931) {
  #Define linkages, joints etc.
  PM_ind <- WGP_PM_coords[,,i]
  lms <- PM_ind
  joint.types <- c('R','S','S','R')
  joint.conn <- rbind(c(0,1),c(1,2),c(2,3),c(3,0))
  joint.cons <- list(NA,NA,NA,NA)
  
  linkage <- defineLinkage(joint.coor = lms, joint.types = joint.types,
                           joint.cons = joint.cons, joint.conn=joint.conn)
  #Overcome error message
  linkage2 <- linkage
  linkage2$joint.cons <- list(c(0,0,1), NA, NA, c(0,0,1))  
  
  #Animate linkage, calculate transmission coefficients, add to table
  drawLinkage(linkage2, animate=FALSE)
  animate <- animateLinkage(linkage2, input.param=c(0,pi/12,length=15), input.joint=1)
  kine <- linkageKinematics(animate)
  KT <- abs(kine$links.rdis.d['Link2',] / kine$links.rdis.d['Link1',])
  KT_PreMax[i] <- KT[2]
  KT_PM_df <- as.data.frame(KT_PreMax)
  
  #Calculate Output rotation
  PM_Output <- print(kine$links.rdis.d['Link2',])
  PreMax_Rotation[i] <- PM_Output[2]
  PM_Rotation_DF <- as.data.frame(PreMax_Rotation)
}

WGP_Pheno_4Bar = cbind(KT_PM_df, PM_Rotation_DF)



###Calculate transmission coefficients for Opercular Linkage (27,23,24,25)

#Define landmarks
OP_LMs <- c(27,23,24,25)
WGP_OP_coords <- F2_array_consensus[OP_LMs,,]

##Calculate transmission coefficients, create dataframe

#Create empty table for data
KT_Opercular <- numeric(length=931)
KT_Output_Rotation <- numeric(length=931)
for (i in 1:931) {
  #Define linkages, joints etc.
  OP_ind <- WGP_OP_coords[,,i]
  lms <- OP_ind
  joint.types <- c('R','S','S','R')
  joint.conn <- rbind(c(0,1),c(1,2),c(2,3),c(3,0))
  joint.cons <- list(NA,NA,NA,NA)
  
  linkage <- defineLinkage(joint.coor = lms, joint.types = joint.types,
                           joint.cons = joint.cons, joint.conn=joint.conn)
  
  #Overcome error message
  linkage2 <- linkage
  linkage2$joint.cons <- list(c(0,0,1), NA, NA, c(0,0,1))
  
  #Animate linkage, calculate transmission coefficients, add to table
  drawLinkage(linkage2)
  animate <- animateLinkage(linkage2, input.param=c(0,pi/12,length=15), input.joint=1)
  kine <- linkageKinematics(animate)
  KT <- abs(kine$links.rdis.d['Link1',] / kine$links.rdis.d['Link3',])
  KT_Opercular[i] <- KT[2]
  KT_df <- as.data.frame(KT_Opercular)
  
  #Calculate Output rotation
  OP_Output <- print(kine$links.rdis.d['Link1',])
  KT_Output_Rotation[i] <- OP_Output[2]
  OP_Rotation_DF <- as.data.frame(KT_Output_Rotation)
}

WGP_operculum = cbind(KT_df, OP_Rotation_DF)

WGP_Pheno_4Bar <- cbind(WGP_Pheno_4Bar, WGP_operculum)
colnames(WGP_Pheno_4Bar) <- c("PreMax_KT", "PreMax_Rotation",
                              "Opercular_KT", "Opercular_Rotation")

WGP_Pheno_4Bar %>% 
  as_tibble() %>% 
  write_csv('WGP_Jaw_kinetic_traits.csv')

