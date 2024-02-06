##############################
##  variance - covariance example data
##
## Matt Brachmann (PhDMattyB)
##
## 06.02.2024
##
##############################

library(geomorph)
library(vcvComp)
library(tidyverse)

theme_set(theme_bw())
## Example data from cichlids
data("Tropheus")

## example code all from the vcvComp package vignette

outliers <- c(18, 56, 155, 351, 624)
Tropheus.IK <- Tropheus[- outliers, ]
# Sample reduced to six populations
Tropheus.IK <- subset(Tropheus.IK, subset = POP.ID %in% levels(POP.ID)[1:6])
Tropheus.IK$POP.ID <- factor(Tropheus.IK$POP.ID)
# New variable combining population and sex
Tropheus.IK$SexPop <- paste(Tropheus.IK$POP.ID, Tropheus.IK$Sex, sep = "_")
Tropheus.IK$SexPop <- as.factor(Tropheus.IK$SexPop)

PHEN <- as.matrix(Tropheus.IK[which(names(Tropheus.IK) == "X1"):
                                which(names(Tropheus.IK) == "Y19")])
rownames(PHEN) <- Tropheus.IK$List_TropheusData_ID


PHEN_array <- arrayspecs(PHEN, p = 19, k = 2)
# Procrustes superimposition
phen.gpa <- gpagen(PHEN_array, print.progress = FALSE)
# conversion array -> matrix of Procrustes coordinates
proc.coord <- two.d.array(phen.gpa$coords)
colnames(proc.coord) <- colnames(PHEN)

phen.pca <- prcomp(proc.coord, rank. = 5, tol = sqrt(.Machine$double.eps))
pc.scores <- phen.pca$x

S.phen.pooled <- cov.group(pc.scores, groups = Tropheus.IK$POP.ID, sex = Tropheus.IK$Sex)

eigen.phen <- mat.sq.dist(S.phen.pooled, dist. = "Riemannian")  # Riemannian distances
prcoa <- pr.coord(eigen.phen)  # ordination
prcoa$Variance  # variance explained

table(Tropheus.IK$POP.ID)  # sample sizes
prop.vcv.test(n = c(69,75), S.phen.pooled[,,"IKA1"], S.phen.pooled[,,"IKS5"])  # ML test

relGV.multi(S.phen.pooled[, , c("IKA1", "IKS5")], logGV = FALSE)
# Relative PCA = relative eigenanalysis
relEigen.a1s5 <- relative.eigen(S.phen.pooled[, , "IKA1"], S.phen.pooled[, , "IKS5"])
relEigen.a1s5$relValues  # relative eigenvalues
# Visualization of the relative eigenvalues (fig. 3)
plot(relEigen.a1s5$relValues[1:relEigen.a1s5$q], 
     log = "y",  las = 1, col = "blue", type = "b", 
     main = "IKA1 relative to IKS5", cex = 0.8, 
     cex.main = 1, cex.axis = 0.8, cex.sub = 0.7, 
     sub = paste("Relative generalized variance =", relEigen.a1s5$relGV), 
     xlab = NA, ylab = "Relative eigenvalues")
abline(h = 1)
