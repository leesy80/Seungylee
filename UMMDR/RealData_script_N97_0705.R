#############
# Clear All #
#############

rm(list=ls())

################
# Load Dataset #
################

tmp <- read.csv('AML_all_dataset.csv', header=TRUE)
str(tmp)
time <- as.numeric(tmp$OS)
censoring <- tmp$isAlive
N <- length(time)
phe <- cbind(time, censoring) # Survival time and Censoring indicators
XX <- cbind(tmp$age,tmp$sex) # covariates
SNPa <- tmp[,-(1:6)] # SNPa
SNPa <- data.matrix(SNPa)
P <- ncol(SNPa) # Number of available SNPs

#################
# Load Packages #
#################

library(survival)
library(compiler)
library(MASS)
library(splines)
library(VennDiagram)

##################
# Load functions #
##################

source('SurvDif.R')
source('SurvMed.R')
# source('Data_Generate.R')
source('Power_Function.R')
source('MDR_KMc_Median_4.R')
source('MDR_Cox_4.R')
source('MDR_Aft.R')
source('MDR_Surv_Median_4.R')
source('UMMDR.R')

## Reduce the SNP pairs that are not significant in coxph model ##
p.value <- betas <- numeric(P)
sur0 <- Surv(time, censoring == 1)
for (i in 1:P){
  sur1 <- coxph(sur0 ~ SNPa[,i] + XX)
  p.value[i] <- summary(sur1)$coef[1,5]
  betas[i] <- summary(sur1)$coef[1,1]
}
varlist <- which(p.value <= 0.05)
length(varlist)

p.value[varlist]
betas[varlist]
colnames(SNPa[,varlist])

# SNPb <- SNPa[,-varlist] ; P2 <- dim(SNPb)[2]
# Ca <- combn(P,2) ; na <- dim(Ca)
# Cb <- combn(P2,2) ; nb <- dim(Cb)

## Real Data Analysis - SNPa ##
snp.all <- SNPa[,-varlist]
pca.mat <- SNPa[,varlist]
nP <- dim(snp.all)[2]

## Obtain PC1 and PC2 ##
t <- 0
for(col in 1:length(varlist)){
  if ( sum( is.na(pca.mat[,col]) ) != 0){
    tm <- which( is.na(pca.mat[,col]) )
    t <- c(t,tm)
  }
}
t <- t[-1] ; ind <- unique(t) 
pca.mat <- pca.mat[-ind,]
pca <- princomp(pca.mat)
XX.pca <- cbind(XX, rep(0,N), rep(0,N))
XX.pca[-ind,3:4] <- pca$scores[,1:2]

## Case I : model = 'umcox', cls.method = 'mg' ## 

cox.ummdr.NoPCA.NoAdj <- UM.MDR(snp.all, phe, XX, cls.method = 'mg', nway = 2, model = 'umcox', adj.main = FALSE, nperm = 20)

cox.ummdr.PCA.NoAdj <- UM.MDR(snp.all, phe, XX.pca, cls.method = 'mg', nway = 2, model = 'umcox', adj.main = FALSE, nperm = 20)

cox.ummdr.NoPCA.Adj <- UM.MDR(snp.all, phe, XX, cls.method = 'mg', nway = 2, model = 'umcox', adj.main = TRUE, nperm = 20)

cox.ummdr.PCA.Adj <- UM.MDR(snp.all, phe, XX.pca, cls.method = 'mg', nway = 2, model = 'umcox', adj.main = TRUE, nperm = 20)

save.image("C:\\Users\\jwl19002\\OneDrive - University of Connecticut\\Desktop\\UMMDR\\N97_0705.RData")

## Case II : model = 'umcox', cls.method = 'km' ## 

km.ummdr.NoPCA.NoAdj <- UM.MDR(snp.all, phe, XX, cls.method = 'km', nway = 2, model = 'umcox', adj.main = FALSE, nperm = 20)

km.ummdr.PCA.NoAdj <- UM.MDR(snp.all, phe, XX.pca, cls.method = 'km', nway = 2, model = 'umcox', adj.main = FALSE, nperm = 20)

km.ummdr.NoPCA.Adj <- UM.MDR(snp.all, phe, XX, cls.method = 'km', nway = 2, model = 'umcox', adj.main = TRUE, nperm = 20)

km.ummdr.PCA.Adj <- UM.MDR(snp.all, phe, XX.pca, cls.method = 'km', nway = 2, model = 'umcox', adj.main = TRUE, nperm = 20)

save.image("C:\\Users\\jwl19002\\OneDrive - University of Connecticut\\Desktop\\UMMDR\\N97_0705.RData")

## Case III : model = 'umcox', cls.method = 'mg2' ## 

cox2.ummdr.NoPCA.NoAdj <- UM.MDR(snp.all, phe, XX, cls.method = 'mg2', nway = 2, model = 'umcox', adj.main = FALSE, nperm = 20)

cox2.ummdr.PCA.NoAdj <- UM.MDR(snp.all, phe, XX.pca, cls.method = 'mg2', nway = 2, model = 'umcox', adj.main = FALSE, nperm = 20)

cox2.ummdr.NoPCA.Adj <- UM.MDR(snp.all, phe, XX, cls.method = 'mg2', nway = 2, model = 'umcox', adj.main = TRUE, nperm = 20)

cox2.ummdr.PCA.Adj <- UM.MDR(snp.all, phe, XX.pca, cls.method = 'mg2', nway = 2, model = 'umcox', adj.main = TRUE, nperm = 20)

save.image("C:\\Users\\jwl19002\\OneDrive - University of Connecticut\\Desktop\\UMMDR\\N97_0705.RData")

	## END of the programm ##

#########################
## Obtain Venn Diagram ##
#########################

## cox.um.mdr ##
cox.ummdr.NoPCA.NoAdj.ind <- which(cox.ummdr.NoPCA.NoAdj$p.value <= 0.05)
cox.ummdr.PCA.Adj.ind <- which(cox.ummdr.PCA.Adj$p.value <= 0.05)
cox.ummdr.NoPCA.Adj.ind <- which(cox.ummdr.NoPCA.Adj$p.value <= 0.05)
cox.ummdr.PCA.NoAdj.ind <- which(cox.ummdr.PCA.NoAdj$p.value <= 0.05)

cox.ummdr.V <- venn.diagram(list(nopca_nomain=cox.ummdr.NoPCA.NoAdj.ind, pca_nomain=cox.ummdr.PCA.NoAdj.ind, nopca_main=cox.ummdr.NoPCA.Adj.ind, pca_main=cox.ummdr.PCA.Adj.ind), filename=NULL, fill=rainbow(4))
grid.newpage()
grid.draw(cox.ummdr.V)

## kmummdr ##
km.ummdr.NoPCA.NoAdj.ind <- which(km.ummdr.NoPCA.NoAdj$p.value <= 0.05)
km.ummdr.PCA.Adj.ind <- which(km.ummdr.PCA.Adj$p.value <= 0.05)
km.ummdr.NoPCA.Adj.ind <- which(km.ummdr.NoPCA.Adj$p.value <= 0.05)
km.ummdr.PCA.NoAdj.ind <- which(km.ummdr.PCA.NoAdj$p.value <= 0.05)

km.ummdr.V <- venn.diagram(list(nopca_nomain=km.ummdr.NoPCA.NoAdj.ind, pca_nomain=km.ummdr.PCA.NoAdj.ind, nopca_main=km.ummdr.NoPCA.Adj.ind, pca_main=km.ummdr.PCA.Adj.ind), filename=NULL, fill=rainbow(4))
grid.newpage()
grid.draw(km.ummdr.V)

## coxummdr2 ##
cox2.ummdr.NoPCA.NoAdj.ind <- which(cox2.ummdr.NoPCA.NoAdj$p.value <= 0.05)
cox2.ummdr.PCA.Adj.ind <- which(cox2.ummdr.PCA.Adj$p.value <= 0.05)
cox2.ummdr.NoPCA.Adj.ind <- which(cox2.ummdr.NoPCA.Adj$p.value <= 0.05)
cox2.ummdr.PCA.NoAdj.ind <- which(cox2.ummdr.PCA.NoAdj$p.value <= 0.05)

cox2.ummdr.V <- venn.diagram(list(nopca_nomain=cox2.ummdr.NoPCA.NoAdj.ind, pca_nomain=cox2.ummdr.PCA.NoAdj.ind, nopca_main=cox2.ummdr.NoPCA.Adj.ind, pca_main=cox2.ummdr.PCA.Adj.ind), filename=NULL, fill=rainbow(4))
grid.newpage()
grid.draw(cox2.ummdr.V)

	## End ##
	## End ##

