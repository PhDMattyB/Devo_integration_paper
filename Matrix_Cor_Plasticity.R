##############################
## matrix correlation - plasticity
##
## Matt Brachmann (PhDMattyB)
##
## 26.08.2024
##
##############################


setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/Updated Landmarks/')

# install.packages('MatrixCorrelation')
library(MatrixCorrelation)


# Mat cor F2 offspring effect -------------------------------------------------------
## Run the plasticity_analysis.R script to get the 
## corrected trait correlations

## if this next object doesn't run, re run the script
## Plasticity due to offspring temperature 
off_plasticity_trait_cor

off_plast_ASHNW = off_plasticity_trait_cor$ASHNW
off_plast_ASHNC = off_plasticity_trait_cor$ASHNC

allCorrelations(X1 = off_plast_ASHNW, 
                X2 = off_plast_ASHNC, 
                ncomp1 = 28,
                ncomp2 = 28)

off_plast_MYVW = off_plasticity_trait_cor$MYVW
off_plast_MYVC = off_plasticity_trait_cor$MYVC


allCorrelations(X1 = off_plast_MYVW, 
                X2 = off_plast_MYVC, 
                ncomp1 = 28,
                ncomp2 = 28)

off_plast_SKRW = off_plasticity_trait_cor$SKRW
off_plast_SKRC = off_plasticity_trait_cor$SKRC

allCorrelations(X1 = off_plast_SKRW, 
                X2 = off_plast_SKRC, 
                ncomp1 = 28,
                ncomp2 = 28)

off_plast_GTSW = off_plasticity_trait_cor$GTSW
off_plast_CSWYC = off_plasticity_trait_cor$CSWYC

allCorrelations(X1 = off_plast_GTSW, 
                X2 = off_plast_CSWYC, 
                ncomp1 = 28,
                ncomp2 = 28)


# Mat cor F1 parental effect -------------------------------------------------

parent_plasticity_trait_cor

parent_plast_ASHNW = parent_plasticity_trait_cor$ASHNW
parent_plast_ASHNC = parent_plasticity_trait_cor$ASHNC

allCorrelations(X1 = parent_plast_ASHNW, 
                X2 = parent_plast_ASHNC, 
                ncomp1 = 28,
                ncomp2 = 28)

parent_plast_MYVW = parent_plasticity_trait_cor$MYVW
parent_plast_MYVC = parent_plasticity_trait_cor$MYVC


allCorrelations(X1 = parent_plast_MYVW, 
                X2 = parent_plast_MYVC, 
                ncomp1 = 28,
                ncomp2 = 28)

parent_plast_SKRW = parent_plasticity_trait_cor$SKRW
parent_plast_SKRC = parent_plasticity_trait_cor$SKRC

allCorrelations(X1 = parent_plast_SKRW, 
                X2 = parent_plast_SKRC, 
                ncomp1 = 28,
                ncomp2 = 28)

parent_plast_GTSW = parent_plasticity_trait_cor$GTSW
parent_plast_CSWYC = parent_plasticity_trait_cor$CSWYC

allCorrelations(X1 = parent_plast_GTSW, 
                X2 = parent_plast_CSWYC, 
                ncomp1 = 28,
                ncomp2 = 28)
