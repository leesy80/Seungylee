	#################
	# Load Packages #
	#################

rm(list=ls())
library(survival)
library(compiler)
library(MASS)
library(splines)

	##################
	# Call Functions #
	##################

source('SurvDif.R')
source('SurvMed.R')
source('Power_Function.R')
source('MDR_Surv_Median_4.R')
source('MDR_KMc_Median_4.R')
source('MDR_Cox_4.R')
source('MDR_Aft.R')
source('UMMDR.R')

	###################
	# Calculate Power #
	###################

  ######################################################
  # niter : Number of iteration for calculating power
  # N : Sample Size
  # time.type : Survival time generation option. 'LogN' or 'PH'
  # gamma : gamma value when generating data with covariate
  #       : 1 if you want to include covariate when generating simulation data 
  #       : 0 if you don't want to include covariate when generating simulation data 
  # beta : beta value when generating data with covariate
  # cp : value that controls censoring rate
  # method : AFT, COX, KM, KMC, KMC2, SURV
  # option : 'BA' for Balanced Accuracy, 'LR' for Log Rank Stat
  ######################################################

## Senario... N = 400, beta = -0.8, gamma = 0 ##

## Survival time : 'PH', 'LogN', 'Weibull' ##
time.type <- 'PH'  

## 1. Censoring Rate cutpoint ##

# cp.rate <- 0.1
cp.rate <- 0.3

## 2. adj.main : TRUE or FALSE ##

# gamma <- 0.8
gamma <- 0.0

## 3. model : 'umcox' or 'coxmdr' ## 
model <- 'umcox' ; cls.method <- 'mg2'
# model <- 'umcox' ; cls.method <- 'mg'
# model <- 'umcox' ; cls.method <- 'km'
# model <- 'coxmdr' ; cls.method <- 'mg2'
# model <- 'kmmdr' ; cls.method <- 'km'

niter <- 100
nway <- 2
nsnp <- 10
N <- 400 
beta <- 0.8 
omega <- 0.0
adj.main <- FALSE
delta <- 0.0 # Add non-causal SNP to hazard function
if ( cls.method == 'aft' ){
  beta <- -beta
}

if (cp.rate == 0){
  cp <- 0.0
} else {
  if ( time.type == 'PH' ){
    if( cp.rate == 0.1 ){
      cp <- 7 
    } else if( cp.rate == 0.2 ){
      cp <- 2.45 
    } else if( cp.rate == 0.3 ){
      cp <- 2.20 
    } else if( cp.rate == 0.5 ){
      cp <- 1.45 
    }
  } else if ( time.type == 'LogN') {
    if( cp.rate == 0.1 ){
      cp <- 39.7 
    } else if( cp.rate == 0.2 ){
      cp <- 11.1
    } else if( cp.rate == 0.3 ){
      cp <- 8.7
    } else if( cp.rate == 0.4 ){
      cp <- 5.7
    } else if( cp.rate == 0.5 ){
      cp <- 4.7
    }
  }
}

if (nway == 3){
  Balanced <- TRUE 
  resulttemp <- NULL
  result0 <- NULL
  u <- 0.2
  source('Data_Generate.R')
  for (i in 1){
    result00 <- CalPower2(N = N, nsnp, nway, trait.type = 'time', time.type = time.type, cls.method, model = model,
                          adj.main = adj.main, niter = niter, alpha = 0.05, cand.pmodels = i:i, cp = cp, mu = 1, gamma = gamma, beta = beta, delta = delta, sigma = 1, Balanced)
    resulttemp <- result00$tables[1,]
    save.image(paste0(model,nway,"_N_",N,"_cp_",cp,"_gamma_",gamma,"_beta_",beta,"_delta_",delta,"_cls.method_",cls.method,"_adj.main_",adj.main,"_type_",time.type,"_",".RData"))
    result0 <- rbind(result0,resulttemp)
  } 
  write.csv(result0, paste0(model,nway,"_N_",N,"_cp_",cp,"_gamma_",gamma,"_beta_",beta,"_delta_",delta,"_cls.method_",cls.method,"_adj.main_",adj.main,"_type_",time.type,"_",".csv"))
} else {
  source('Data_Generate.R')
  resulttemp <- NULL
  result0 <- NULL
  Balanced <- FALSE 
  for (i in 1:70){
    result00 <- CalPower2(N = N, nsnp, nway, trait.type = 'time', time.type = time.type, cls.method, model = model,
                          adj.main = adj.main, niter = niter, alpha = 0.05, cand.pmodels = i:i, cp = cp, mu = 1, gamma = gamma, beta = beta, delta = delta, sigma = 1, Balanced)
    resulttemp <- result00$tables[1,]
    save.image(paste0(model,nway,"_N_",N,"_cp_",cp,"_gamma_",gamma,"_beta_",beta,"_delta_",delta,"_cls.method_",cls.method,"_adj.main_",adj.main,"_type_",time.type,"_",".RData"))
    result0 <- rbind(result0,resulttemp)
  } 
  write.csv(result0, paste0(model,nway,"_N_",N,"_cp_",cp,"_gamma_",gamma,"_beta_",beta,"_delta_",delta,"_cls.method_",cls.method,"_adj.main_",adj.main,"_type_",time.type,"_",".csv"))
}
