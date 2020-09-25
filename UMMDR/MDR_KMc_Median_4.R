###########################################
# Median MDR for single SNP or a SNP pair #
###########################################

Median.MDR <- function(snp.dat.part, SS){
  death <- SS[,2] 
  time <- SS[,1] 
  snp.dat.part <- as.matrix(snp.dat.part)
  if ( is.null(snp.dat.part) == TRUE ){
    k <- 1
  } else {
    k <- dim(snp.dat.part)[2] # Number of SNP factors... likely to be 2 in this case.
  }
  result <- list()

  ## Check missing values
  time.fids <- which(is.na(time) == FALSE)
  death.fids <- which(is.na(death) == FALSE)
  snp.fids <- which(complete.cases(snp.dat.part))
  fids <- intersect(intersect(time.fids, death.fids), snp.fids)

  # Number of observation without missing
  n <- length(time)

  # G is a binary variable which has 2 if obs is in high risk, 1 otherwise.
  G <- rep(1, n)

  ## remove missing values
  snp.dat.part <- as.matrix(snp.dat.part[fids,])
  time <- time[fids]
  death <- death[fids]

  ## Find Median Survival time
  # total.med <- unname(summary(survfit(Surv(time, death) ~ 1))$table["median"])
  total.med <- SurvMed(time, death)[1] # Causes errors too often

  ## If median exist, then implement MDR. If not, return PASS
  if ( is.na(total.med) == FALSE ){
    ## split data into 3^k cells
    tmp.list <- vector('list', k) # make k-dimensional list : tmp.list[[i]]
    for(i in 1:k){
      tmp.list[[i]] <- snp.dat.part[,i] # For each tmp.list[[i]], the i-th column of snp.dat is plugged in
    }
    cells <- split(data.frame(cbind(fids, snp.dat.part)), tmp.list) 
    # Cell becomes list, whose i-th list contains observations whose patterns of tmp.list are same
    # Since tmp.list has 3*3 = 9 patterns, cells should have 9 lists
    # For each list in cells, the observations with same patterns to tmp.list are contained.
  
    ## delete NULL cells ##
    obs.cell <- sapply(cells, function(x) nrow(x)) # apply nrow() function for each list in cells.
    cell.null <- which(obs.cell == 0) # Null cells should be deleted.
    if (length(cell.null) > 0 ) {
      cells <- cells[-cell.null]
    }
    cells.num <- length(cells) # Number of cells created
  
    ## Find survival times within each cells based on KM method ##
    Cell.membership <- rep(NA, n)
    for(cel in 1:cells.num){
      temp.ids <- cells[[cel]][,1] # cell ids in cel-th cell
      Cell.membership[temp.ids] <- cel 
      if(length(temp.ids) != 0){
        # Evaluate the median survival time for complement cells
        # comp.med <- unname(summary(survfit(Surv(time[-temp.ids], death[-temp.ids]) ~ 1))$table["median"])
        comp.med <- SurvMed(time[-temp.ids], death[-temp.ids])[1]
        if (is.na(comp.med) == FALSE){
          if(total.med < comp.med){        
            # This case implies that complement cells are low risk, thus cel-th cell is high.
            G[temp.ids] <- 2  
          }
        } else {
          # If complementary median does not exist, then the cel-th is in High Risk
          G[temp.ids] <- 2
        }
      }
    }

    ## Evaluate each cells based on the overall survival time and cell-survival time
    if (length(levels(as.factor(G))) == 1){ 
      # If observations in training sets are all in 'High' or all in 'Low'
      log.rank.stat <- 0
    } else {	
      # Calculate log.rank statistics  
      log.rank.stat <- SurvDif.LogRank(time, death, G)
    }
    result$stat <- log.rank.stat
    result$HL <- G 
    result$Cell.id <- Cell.membership
  } else {
    # If total median does not exist, then drop the corresponding pair
    result$stat <- 0
  }
  return(result)
}

############################################
# Choose the best SNP via Cross Validation #
############################################

CV.Median.MDR <- function(snp.dat, SS, nway, nfold){
  snp.dat <- snp.dat
  time <- time <- SS[,1]
  death <- death <- SS[,2]
  N <- length(time)
  obs.id <- 1:N
  ind <- which(death == 1)
  com.id <- obs.id[ind]
  cen.id <- obs.id[-ind]
  N1 <- sum(death == 1)
  N2 <- N - N1

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
      test.id <- na.omit( c(com.id[test.id1], cen.id[test.id2]) )
  
      # Obtain result from training fold #
      train.SS <- SS[-test.id, ]
      if (k == 1){
        train.snp.pair <- tmp.pair[-test.id]
      } else {
        train.snp.pair <- tmp.pair[-test.id,]  
      }
     
      # MDR for training #
      tmp <- Median.MDR(train.snp.pair, train.SS)
      CV.train[combi,cv] <- tmp$stat
 
      # HL membership for testing #
      CV.test.HL[test.id] <- Median.MDR(tmp.pair, SS)$HL[test.id]
    }
    CV.test[combi] <- SurvDif.LogRank(time, death, CV.test.HL)
    cat(combi,'th combination is done','\n')
  }

  # Obtain CVC # 
  cvc <- apply(CV.train, 2, which.max)

  result <- list()
  result$train.score <- apply(CV.train, 1, mean) # average nfold traning results from CV
  result$test.score <- CV.test    
  result$cvc <- cvc  
  result$best.combi.No <- as.numeric(most.freq(cvc, n=1, freq = FALSE))
  return(result)
}

#########################
# Impute the HL results #
#########################

Med.MDR.CV.pred <- function(snp.dat.pair, time, death, nfold){
  snp.dat.pair <- snp.dat.pair # Enter selected pairs
  time <- time
  death <- death 
  N <- length(time)
  if (is.null(snp.dat.pair) == TRUE){
    is.vector <- TRUE 
  } else {
    is.vector <- FALSE
  }
  obs.id <- 1:N
  ind <- which(death == 1)
  com.id <- obs.id[ind]
  cen.id <- obs.id[-ind]
  N1 <- sum(death == 1)
  N2 <- N - N1

  nfold <- nfold
  cv.len1 <- floor(N1/nfold + 1)
  cv.len2 <- floor(N2/nfold + 1)

  ## Implement 10-fold model fitting ##
  for (cv in 1:nfold){
    if (cv != nfold){
      test.id1 <- ((cv - 1) * cv.len1 + 1) : (cv * cv.len1) 
      test.id2 <- ((cv - 1) * cv.len2 + 1) : (cv * cv.len2) 
    } else {
      test.id1 <- ((cv - 1) * cv.len1 + 1) : N1
      test.id2 <- ((cv - 1) * cv.len2 + 1) : N2
    }
    test.id <- c(com.id[test.id1], cen.id[test.id2])

    # HL membership for testing #
    CV.test.HL[test.id] <- Median.MDR(tmp.pair, time, death)$HL[test.id]
  }
  return(CV.test.HL)  
}

################################################################### 
# For given SNP pair, calculate permutation p-value without 10-CV #
###################################################################

Med.MDR.Perm.pval.a <- function(snp.dat.part, time, death, nperm){
  snp.dat.pair <- snp.dat.part
  time <- time
  death <- death 
  nperm <- nperm
  if (is.null(dim(snp.dat.pair)) == TRUE){
    N <- length(snp.dat.pair)
  } else {
    N <- dim(snp.dat.pair)[1]
  }

  ## Obtain the obs score ##
  tmp <- Median.MDR(snp.dat.part, time, death)
  obs.score <- tmp$stat 

  ## Calculate the permutational distribution of score ##
  Perm.score <- numeric(nperm)
  for (perm in 1:nperm){
    shuf.id <- sample(1:N, replace=FALSE)
    shuf.time <- time[shuf.id]
    shuf.death <- death[shuf.id] 
    tmp <- Median.MDR(snp.dat.part, shuf.time, shuf.death)
    Perm.score[perm] <- tmp$stat
    cat(perm,'th permutation is done','\n')
  }

  ## Obtain right tail p-value ##
  pvalue <- sum(Perm.score >= obs.score) / nperm

  output <- list()
  output$Perm.score <- Perm.score
  output$pvalue <- pvalue 
  output$obs.score <- obs.score
  return(output)
}

################################################################ 
# For given SNP pair, calculate permutation p-value with 10-CV #
################################################################

Med.MDR.Perm.pval.b <- function(snp.dat.part, time, death, nfold, nperm){
  snp.dat.pair <- snp.dat.part
  time <- time
  death <- death 
  nperm <- nperm
  nfold <- nfold
  if (is.null(snp.dat.pair) == TRUE){
    N <- length(snp.dat.pair)
  } else {
    N <- dim(snp.dat.pair)[1]
  }
  cv.len <- floor(N/nfold + 1)

  ## Obtain the obs score in Cross Validation ##
  obs.cv.score <- numeric(nfold)
  for (cv in 1:nfold){ 
    if (cv != nfold){
      test.id <- ((cv - 1) * cv.len + 1) : (cv * cv.len) 
    } else {
      test.id <- ((cv - 1) * cv.len + 1) : N
    }
    test.time <- time[test.id]   
    test.death <- death[test.id]        
    test.snp.part <- snp.dat.part[test.id,]
    test.tmp <- Median.MDR(test.snp.part, test.time, test.death)
    obs.cv.score[cv] <- test.tmp$stat 
  }
  obs.score <- mean(obs.cv.score)

  ## Calculate the permutational distribution of score ##
  Perm.score <- Perm.train.score <- numeric(nperm)
  for (perm in 1:nperm){
    shuf.id <- sample(1:N, replace=FALSE)
    shuf.time <- time[shuf.id]
    shuf.death <- death[shuf.id]
    train.cv.score <- test.cv.score <- numeric(nfold)
    for (cv in 1:nfold){ 
      if (cv != nfold){
        test.id <- ((cv - 1) * cv.len + 1) : (cv * cv.len) 
      } else {
        test.id <- ((cv - 1) * cv.len + 1) : N
      }
      test.time <- shuf.time[test.id]   
      test.death <- shuf.death[test.id]        
      test.snp.part <- snp.dat.part[test.id,]

      # Calculate testing score #
      test.tmp <- Median.MDR(test.snp.part, test.time, test.death)
      test.cv.score[cv] <- test.tmp$stat 

      # Calculate training score #
      train.tmp <- Median.MDR(snp.dat.part[-test.id,], shuf.time[-test.id], shuf.death[-test.id])
      train.cv.score[cv] <- train.tmp$stat 
      cat('...',cv,'th Cross Validation is done','\n')
    }
    Perm.score[perm] <- mean(test.cv.score)
    Perm.train.score[perm] <- mean(train.cv.score)
    cat(perm,'th Permutation is done','\n')
  }
  
  ## Obtain right tail p-value ##
  pvalue <- sum(Perm.score >= obs.score) / nperm

  output <- list()
  output$Perm.train.score <- Perm.train.score
  output$Perm.score <- Perm.score
  output$pvalue <- pvalue 
  output$obs.score <- obs.score
  return(output)
}

