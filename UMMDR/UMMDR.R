## Classification step, assign high/low for each cell and test ##
class.HL <- function(snp.pair, SS, cls.method, cova){
  # SS[,1] : Survival time, SS[,2] : censoring indicator
  snp.pair <- as.matrix(snp.pair)
  time <- SS[,1]
  censoring <- SS[,2]
  N <- length(time)
  
  # Define a new random variable S : which indicates high/low
  S <- weights <- rep(0, N)

  # Remove missing values
  fids <- which(complete.cases(snp.pair))
  snp.pair <- as.matrix(snp.pair[fids, ])
  k <- dim(snp.pair)[2] 
  
  # Split completed data into cells
  tlist <- vector('list', k)
  for(i in 1:k){
    tlist[[i]] <- snp.pair[, i]
  }
  cells <- split(data.frame(cbind(fids, snp.pair)), tlist)
  
  # delete NULL cells
  obs.cell <- sapply(cells, function(x) nrow(x))
  cell.null <- which(obs.cell == 0)
  if(length(cell.null) > 0 ){
    cells <- cells[-cell.null]
  }
  ncells <- length(cells) # number of cells

  # Classify each cell into H/L by various methods
  if( cls.method == 'mg' ){
    # Use martingale residuals to define H/L, corresponding to time trait
    Ycen <- Surv(time, censoring == 1) # if SS[,1] is 1, then the obs is not censored
    cox <- coxph(Ycen ~ 1) # Fit cox model WITHOUT covariates
    mt.resid <- residuals(cox, type="martingale") # Extract martingale residuals

    # Define H/L values based on mt.residuals
    for(i in 1:ncells){
      temp.ids <- cells[[i]][, 1]
      # Pass the cell whose number of obs is zero
      if(length(temp.ids) == 0){
        next 
      }
      mSS <- sum(na.omit(mt.resid[temp.ids]))
      weights[temp.ids] <- (length(temp.ids))
      # if mss >=0, then the cell is defined to be high risk
      if (mSS >= 0){ 
        S[temp.ids] <- 1
      }
    }
  }
  if( cls.method == 'mg2' ){
    # Use martingale residuals to define H/L, corresponding to time trait
    Ycen <- Surv(time, censoring == 1) # if SS[,1] is 1, then the obs is not censored
    if (is.null(cova) == TRUE){
      cox <- coxph(Ycen ~ 1) 
    } else {
      cox <- coxph(Ycen ~ cova) # Fit cox model USING covariates
    }
    mt.resid <- residuals(cox, type="martingale") # Extract martingale residuals

    # Define H/L values based on mt.residuals
    for(i in 1:ncells){
      temp.ids <- cells[[i]][, 1]
      # Pass the cell whose number of obs is zero
      if(length(temp.ids) == 0){
        next 
      }
      mSS <- sum(na.omit(mt.resid[temp.ids]))
      weights[temp.ids] <- (length(temp.ids))

      # if mss >=0, then the cell is defined to be high risk
      if (mSS >= 0){ 
        S[temp.ids] <- 1
      }
    }
  }
  if( cls.method == 'aft' ){  
    # Use martingale residuals to define H/L, corresponding to time trait
    Ycen <- Surv(time, censoring == 1) # if SS[,1] is 1, then the obs is not censored
 
    # Fit AFT regression WITHOUT covariates
    aft <- survreg(Ycen ~ 1, dist="lognormal") 
    AFT.resi <- (log(time) - aft$coefficients[1]) / aft$scale

    # Standardize the AFT residuals
    AFT.res.sd <- AFT.resi - mean(AFT.resi)

    # Define H/L values based on mt.residuals
    for(i in 1:ncells){
      temp.ids <- cells[[i]][,1]
      if(length(temp.ids) == 0){
        next # Pass the cell whose number of obs is zero
      }
      AFT.ss <- sum(na.omit(AFT.res.sd[temp.ids]))
      weights[temp.ids] <- (length(temp.ids))
      # if AFT.ss <=0, then the cell is defined to be high risk
      if (AFT.ss <= 0){ 
        S[temp.ids] <- 1
      }
    }
  }
  if( cls.method == 'aft2' ){  
    # Use martingale residuals to define H/L, corresponding to time trait
    Ycen <- Surv(time, censoring == 1) # if SS[,1] is 1, then the obs is not censored
    
    # Fit AFT regression WITH covariates
    aft <- survreg(Ycen ~ cova, dist="lognormal") 
    AFT.resi <- (log(time) - aft$coefficients[1]) / aft$scale

    # Standardize the AFT residuals
    AFT.res.sd <- AFT.resi - mean(AFT.resi)

    # Define H/L values based on mt.residuals
    for(i in 1:ncells){
      temp.ids <- cells[[i]][,1]
      if(length(temp.ids) == 0){
        next # Pass the cell whose number of obs is zero
      }
      AFT.ss <- sum(na.omit(AFT.res.sd[temp.ids]))
      weights[temp.ids] <- (length(temp.ids))
      # if AFT.ss <=0, then the cell is defined to be high risk
      if (AFT.ss <= 0){ 
        S[temp.ids] <- 1
      }
    }
  }
  if( cls.method == 'km' ){  
    # Find Median Survival time
    total.med <- SurvMed(time, censoring)[1] 
    # Find survival times within each cells based on KM method. 
    for(i in 1:ncells){
      temp.ids <- cells[[i]][,1]
      if(length(temp.ids) == 0){
        next # Pass the cell whose number of obs is zero
      }
      # Evaluate the median survival time for complement cells
      comp.med <- SurvMed(time[-temp.ids], censoring[-temp.ids])[1]
      weights[temp.ids] <- length(temp.ids)
      if(total.med < comp.med){        
        # This case implies that complement cells are low risk, thus cel-th cells are high.
        S[temp.ids] <- 1
      }
    }
  }
   
  # Missing values in SNPs should be remained as NA in S
  delist <- 0
  for (col in 1:k){
    delist <- c(delist, which(is.na(snp.pair[,k]) == TRUE))
  }   
  delist <- delist[-1]
  delist <- unique(delist) 
  S[delist] <- NA

  output <- list()
  output$S <- S
  output$weights <- weights
  return(output) 
}
class.HL <- cmpfun(class.HL)

## Estimate non-center parameter of the null distribution due to classification ##

nullcenter <- function(SS, cova, snp.pair, model, cls.method, nperm = 20){
  n <- dim(SS)[1]
  nocov <- FALSE
  if( is.null(cova) == TRUE ){
    nocov <- TRUE
  } else if( is.null(dim(cova)) == TRUE ){
    k <- 1
  } else {
    k <- dim(cova)[2]
  }
  center <- rep(0, nperm)
  cls.method <- cls.method

  if (nocov == FALSE){
    for (i in 1:nperm){
      perm.id <- sample(1:n, n, replace = FALSE) # id.permutation in terms of observation 
      perm.SS <- SS[perm.id,] # Permute Traits
      if( k == 1 ){
        perm.cova <- cova[perm.id]
      } else { 
        perm.cova <- cova[perm.id,]
      }
      ## Implement class_HL based on permuted data ##
      cls <- class.HL(snp.pair, perm.SS, cls.method, cova) 
      cls.S <- cls$S # H/L indicator values based on mt.resi   

      ## Run preliminary analysis ##
      Ycen <- Surv(perm.SS[,1], perm.SS[,2] == 1)
      # Combined covariate whose 1st column is H/L indicator
      X <- cbind(cls.S, perm.cova) ; X <- data.matrix(X)
      if (model == 'umcox'){
        sumstat <- summary(coxph(Ycen ~ X))$coef # Extract test statistics regarding H/L indicator
        center[i] <- sumstat[1,4]^2 # Wald statistics related with H/L indicator
      } 
      if (model == 'umaft'){ 
        sumstat <- summary(survreg(Ycen ~ X, dist="lognormal"))$table # Fit AFT Regression with object from Surv
        center[i] <- sumstat[2,3]^2 # Wald statistics related with H/L indicator
      } 
    }
  } else { 
    for (i in 1:nperm){
      perm.id <- sample(1:n, n, replace = FALSE) # id.permutation in terms of observation 
      perm.SS <- SS[perm.id,] # Permute Traits
      ## Implement class_HL based on permuted data ##
      cls <- class.HL(snp.pair, perm.SS, cls.method) 
      cls.S <- cls$S # H/L indicator values based on mt.resi   

      ## Run preliminary analysis ##
      Ycen <- Surv(perm.SS[,1], perm.SS[,2] == 1)
      # Combined covariate whose 1st column is H/L indicator
      X <- data.matrix(cls.S)
      if (model == 'umcox'){
        sumstat <- summary(coxph(Ycen ~ X))$coef # Extract test statistics regarding H/L indicator
        center[i] <- sumstat[1,4]^2 # Wald statistics related with H/L indicator
      } 
      if (model == 'umaft'){ 
        sumstat <- summary(survreg(Ycen ~ X, dist="lognormal"))$table # Fit AFT Regression with object from Surv
        center[i] <- sumstat[2,3]^2 # Wald statistics related with H/L indicator
      } 
    }
  }
  return(mean(na.omit(center)))
}
nullcenter <- cmpfun(nullcenter)

############################################################################################################
## the non-zero center parameter is estimated by a few permutation 
## snp.all ----snp matrix, n by p
## snp.combs ---all snp pairs 2 by (p choose 2)
## phe --- phenotype... must be n by 2 matrix with lst column : Survival time, 2nd column : Censoring Indicator 
## cova -- covariate
## cls.method -- 'mg' / 'mg2' / 'aft' / 'aft2' / 'km' : classification rule for H/L
## model -- 'umcox' / 'umaft' : Determine the model that we use in UM.MDR
## adj.main -- TRUE / FALSE : adjust marginal effect of SNP pairs or not
## nperm -- if um.mdr : small number of permutation time to estimate non-central parameters ... set to be 10
##       -- if cv.mdr : number of cross validations ... set to be 10
## output: --- a list including the coefficient (beta), pvalue (pv),
## ----------- p-values (p.values) raw pvalue (rpv) and statistics(stat) for all snp pairs
############################################################################################################

UM.MDR <- function(snp.all, SS, cova, cls.method, nway, model, adj.main, nperm){
  nSNP <- dim(snp.all)[2] 
  k <- nway
  snp.combs <- combn(nSNP, k) # Create the SNP pairs 
  ns <- dim(snp.combs)[2] # Number of 2-combinations available  
  n <- dim(SS)[1]
  d <- dim(SS)[2]
  pv <- NULL
  cls.method <- cls.method
  p.values <- coefs <- stat <- rpvs<- rep(0, ns) 
  if (is.null(cova) == FALSE){
    cova <- as.matrix(cova)
  }

  ## Select best model(i.e. snp combination)
  Ycen <- Surv(SS[,1], SS[,2] == 1)
  if ( model == 'umcox' ){
    ## Use Coxph to fit a model
    for (j in 1:ns){
      cls <- class.HL(snp.all[, snp.combs[, j]], SS, cls.method, cova)
      cls.S <- cls$S
      if ( adj.main == TRUE ){
        ## Use coxph including cls.S, cova, and SNP pair
        snp1 <- snp.all[, snp.combs[,j]]
        if ( cls.method == 'mg' ){
          if (is.null(cova) == FALSE){
            X <- cbind(cls.S, cova, snp1) ; rm(snp1)
          } else {
            X <- cbind(cls.S, snp1) ; rm(snp1)
          }
          sumstat <- summary(coxph(Ycen ~ X))$coef
        } else {
          X <- cbind(cls.S, snp1) ; rm(snp1)
          sumstat <- summary(coxph(Ycen ~ X))$coef
        }
        null.center <- nullcenter(SS, X, snp.all[, snp.combs[, j]], model, cls.method, nperm = 20)
        lambda <- max(null.center - 1, 0)
      } else {
        ## Use coxph including cls.S, cova
        if ( cls.method == 'mg' ){
          if (is.null(cova) == FALSE){
            X <- cbind(cls.S, cova)
          } else {
            X <- cls.S
          }
        } else {
          X <- cls.S
        }
        sumstat <- summary(coxph(Ycen ~ X))$coef
        null.center <- nullcenter(SS, X, snp.all[, snp.combs[, j]], model, cls.method, nperm = 20)
        lambda <- max(null.center - 1, 0)
      }
      stat[j] <- sumstat[1,4]^2 
      p.values[j] <- pchisq(stat[j], df=1, ncp = lambda, lower.tail = FALSE )
      coefs[j] <- sumstat[1,1]
      rpvs[j] <- pchisq(stat[j], df=1, ncp = 0, lower.tail = FALSE ) # Z^2 statistics or Chisq with df = 1
      cat(j,'th combination is done with p-value : ',p.values[j],'\n') 
    }
  } 
  if ( model == 'umaft' ){
    ## Use lognormal-Aft to fit a model
    for (j in 1:ns){
      cls <- class.HL(snp.all[, snp.combs[, j]], SS, cls.method, cova)
      cls.S <- cls$S
      if ( adj.main == 'TRUE' ){
        ## Use survreg including cls.S, cova, and SNP pair
        snp1 <- snp.all[, snp.combs[,j]]
        if ( cls.method == 'aft' ){
          if (is.null(cova) == FALSE){
            X <- cbind(cls.S, cova, snp1) ; rm(snp1)
          } else {
            X <- cbind(cls.S, snp1) ; rm(snp1)
          }
          tmp <- survreg(Ycen ~ X, dist='lognormal') # Fit AFT Regression with object from Surv
        } else {
          X <- cbind(cls.S, snp1) ; rm(snp1)
          tmp <- survreg(Ycen ~ X, dist='lognormal') # Fit AFT Regression with object from Surv
        }
        null.center <- nullcenter(SS, X, snp.all[, snp.combs[, j]], model, cls.method, nperm = 20) 
        lambda <- max(null.center-1, 0)
      } else {
        # Fit AFT Regression with object from Surv
        if ( cls.method == 'aft' ){
          if (is.null(cova) == FALSE){
            X <- cbind(cls.S, cova) ; rm(snp1)
          } else {
            X <- cls.S ; rm(snp1)
          }
        } else {
          X <- cls.S
        }
        tmp <- survreg(Ycen ~ X, dist='lognormal') 
        null.center <- nullcenter(SS, X, snp.all[, snp.combs[, j]], model, cls.method, nperm = 20) 
        lambda <- max(null.center-1, 0)
      }
      stat[j] <- (summary(tmp)$table[2,3])^2 
      p.values[j] <- pchisq(stat[j], df = 1, ncp = lambda, lower.tail = FALSE )
      coefs[j] <- summary(tmp)$table[2,1]
      rpvs[j] <- pchisq(stat[j], df = 1, ncp = 0, lower.tail = FALSE ) # Z^2 statistics or Chisq with df = 1
      cat(j,'th combination is done with p-value : ',p.values[j],'\n') 
    }
  } 
  output <- list()
  output$cls.method <- cls.method
  output$model <- model
  output$coefficients <- coefs
  output$stat <- stat
  output$rpv <- rpvs
  output$p.values <- p.values
  return(output)
}
UM.MDR <- cmpfun(UM.MDR)

FDR <- function(data, alpha){
  index <- rank(data)
  N <- length(data)
  fdr <- index/N*alpha
  return(fdr)
}
FDR <- cmpfun(FDR)

