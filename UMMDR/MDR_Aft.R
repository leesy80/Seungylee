######################################################################################
# For given snp.dat.part, Implement Cox_MDR approach and return Log_rank_statistics. #
######################################################################################

Aft.MDR <- function(snp.dat.part, SS, X, method = 'BA'){
  death <- SS[,2] 
  time <- SS[,1]
  X <- X
  snp.dat.part <- as.matrix(snp.dat.part)
  if ( is.null(snp.dat.part) == TRUE ){
    k <- 1
  } else {
    k <- dim(snp.dat.part)[2] # Number of SNP factors... likely to be 2 in this case.
  }
  if (method != 'LR'){
    method <- 'BA' # Default method in classification
  }
  result <- list()

  ## Check missing values
  time.fids <- which(is.na(time) == FALSE)
  death.fids <- which(is.na(death) == FALSE)
  snp.fids <- which(complete.cases(snp.dat.part))
  fids <- intersect(intersect(time.fids, death.fids), snp.fids)

  ## remove missing values
  snp.dat.part <- as.matrix(snp.dat.part[fids, ])
  time <- time[fids]
  death <- death[fids]
 
  ## Number of observation without missing
  n <- length(time)

  ## G is a binary variable which has 2 if obs is in high risk, 1 otherwise.
  G <- rep(1, n)

  ## remove missing values from X
  snp.dat.part <- as.matrix(snp.dat.part[fids,])
  time <- time[fids]
  death <- death[fids] 
  if ( is.null(X) == FALSE ){ 
    if ( is.null(dim(X)) == TRUE){
      X <- X[fids]
    } else {
      X <- X[fids,] 
    }
  }

  ## split data into 3^k cells
  tmp.list <- vector('list', k) # make k-dimensional list : tmp.list[[i]]
  for(i in 1:k){
    tmp.list[[i]] <- snp.dat.part[,i] # For each tmp.list[[i]], the i-th column of snp.dat is plugged in
  }
  cells <- split(data.frame(cbind(fids, snp.dat.part)), tmp.list) 
  # Cell becomes list, whose i-th list contains observations whose patterns of tmp.list are same
  # Since tmp.list has 3*3 = 9 patterns, cells should have 9 lists
  # For each list in cells, the observations with same patterns to tmp.list are contained.

  ## delete NULL cells
  obs.cell <- sapply(cells, function(x) nrow(x)) # apply nrow() function for each list in cells.
  cell.null <- which(obs.cell == 0) # Null cells should be deleted.
  if (length(cell.null) > 0 ) {
    cells <- cells[-cell.null]
  }
  cells.num <- length(cells) # Number of cells created

  ## Find the sign of the Standardized residuals for each cells
  sur <- Surv(time, death) # Fit Surv model
  if ( is.null(X) == TRUE ){
    aft <- survreg(sur ~ 1, dist="lognormal") # Fit AFT Regression with object from Surv
    AFT.resi <- ((log(time)-(aft$coefficients[1]))/(aft$scale))
  } else {
    aft <- survreg(sur ~ X, dist="lognormal") # Fit AFT Regression with object from Surv
    AFT.resi <- ((log(time)-(aft$coefficients[1] + X * aft$coefficients[2]))/(aft$scale))
  }

  ## Standardize the AFT residuals
  AFT.res.sd <- AFT.resi - mean(AFT.resi)
  
  ## Find survival times within each cells based on KM method. 
  Cell.membership <- rep(NA, n)
  for(cel in 1:cells.num){
    # Select cell IDs in cel-th cell
    temp.ids <- cells[[cel]][,1] 
    Cell.membership[temp.ids] <- cel 
    # Sum the martingale residuals within the cel-th cells
    Cell.sum <- sum(na.omit(AFT.res.sd[temp.ids]))
    if( Cell.sum < 0 ){        
      # Obs in cel-th cell are listed as High.risk members if their martingale cell.sum < 0
      G[temp.ids] <- 2 # 2 implies High risk  
    }
  }
  if (method == 'BA'){
    ## Calculate balanced accuracy based on the Groups
    TP <- length(intersect(which(AFT.res.sd > 0), which(G == 1))) 
    FN <- length(intersect(which(AFT.res.sd > 0), which(G == 2)))
    FP <- length(intersect(which(AFT.res.sd <= 0), which(G == 1)))
    TN <- length(intersect(which(AFT.res.sd <= 0), which(G == 2)))
    result$stat <- 0.5 * (TP/(TP + FN) + TN/(TN + FP))
  }
  if (method == 'LR'){
    ## Calculate log.rank.statistics based on the Groups
    if (length(levels(as.factor(G))) == 1){ 
      # If observations in training sets are all in 'High' or all in 'Low'
      log.rank.stat <- NA
    } else {	
      # Calculate log.rank statistics  
      log.rank.stat <- SurvDif.LogRank(time, death, G)
    }
    result$stat <- log.rank.stat
  }
  result$HL <- G # High Risk : 2
  return(result)
}

#########################################################
# Find best combination of k=2 model among all SNP data #
#########################################################

CV.Aft.MDR <- function(snp.dat, SS, X, nway, method, nfold){
  snp.dat <- snp.dat
  SS <- SS # SS[,1] : time, SS[,2] : death
  method <- method 
  N <- dim(SS)[1]
  obs.id <- 1:N
  ind <- which(SS[,2] == 1)
  com.id <- obs.id[ind]
  cen.id <- obs.id[-ind]
  N1 <- sum(SS[,2] == 1)
  N2 <- N - N1
  if ( is.null(X) == TRUE ){
    X <- rep(1,N)
  }
  k <- nway
  nfold <- nfold 
  numK <- ncol(snp.dat)
  comb.mat <- combn(numK,k)
  ncombi <- dim(comb.mat)[2]
  CV.train <- matrix(0, ncombi, nfold)
  CV.test <- numeric(ncombi)

  ## Split the whole data into folds .. adjust the cf rates for each fold ##
  cv.len1 <- floor(N1/nfold + 1)
  cv.len2 <- floor(N2/nfold + 1)

  # Investigate all combinations #
  for (combi in 1:ncombi){ 
    tmp.col <- comb.mat[ ,combi]
    tmp.pair <- snp.dat[ ,tmp.col]
    CV.test.HL <- numeric(N)

    # Start CV procedure #
    for (cv in 1:nfold){
      if (cv != nfold){
        test.id1 <- ((cv - 1) * cv.len1 + 1) : (cv * cv.len1) 
        test.id2 <- ((cv - 1) * cv.len2 + 1) : (cv * cv.len2) 
      } else {
        test.id1 <- ((cv - 1) * cv.len1 + 1) : N1
        test.id2 <- ((cv - 1) * cv.len2 + 1) : N2
      }
      test.id <- na.omit(c(com.id[test.id1], cen.id[test.id2]))
  
      # Obtain result from training fold #
      train.SS <- SS[-test.id,]
      if ( is.null(dim(X)) == TRUE ){
        train.X <- X[-test.id] 
      } else {
        train.X <- X[-test.id,] 
      }
      if (k == 1){
        train.snp.pair <- tmp.pair[-test.id]
      } else {
        train.snp.pair <- tmp.pair[-test.id,]  
      }

      # MDR for training #
      tmp <- Aft.MDR(train.snp.pair, train.SS, train.X, method)
      CV.train[combi,cv] <- tmp$stat
 
      # HL membership for testing #
      CV.test.HL[test.id] <- Aft.MDR(tmp.pair, SS, X, method)$HL[test.id]
    }
    CV.test[combi] <- SurvDif.LogRank(SS[,1], SS[,2], CV.test.HL)
    cat(combi,'th combination is done','\n')
  }

  # Obtain CVC # 
  cvc <- apply(CV.train, 2, which.max)

  result <- list()
  result$train.score <- apply(CV.train, 1, mean) # average nfold traning results from CV
  result$test.score <- CV.test    
  result$cvc <- cvc  
  result$best.combi.No <- which.max( na.omit(result$train.score) )[1]
  return(result)
}