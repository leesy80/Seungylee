###################################
# Calculate median test statistcs #
###################################

CZ.Median.test <- function(time, death, G){
  time <- time  
  death <- death
  G <- as.factor(G)    
  G.levels <- as.numeric(levels(G)) 
  ngroup <- length(G.levels)

  # Calculate Survivial time
  total.med <- SurvMed(time, death)[1]
  
  ## If median exist, then implement MDR. If not, PASS
  if ( is.na(total.med) == FALSE ){
    # Estimate F_i : survival function by group.
    eta.hat <- theta.med <- Fi.se <- numeric(ngroup)

    ## Calculate median statistics C
    for (g in G.levels){
      # Take only observations in g-th group
      tmp.death <- death[which(G == g)] 
      tmp.time <- time[which(G == g)] 
    
      # Fit Kaplan Meier model within g-th group
      tmp.surv <- survfit(Surv(tmp.time, tmp.death) ~ 1, error='greenwood')

      # Calculate survival probability based on Kaplan Meier
      F.i <- summary(tmp.surv)$surv 

      # Standard error of each survival prob based on Greenwood
      F.se <- summary(tmp.surv)$std.err 

      # Survival times of observations in g-th group
      theta.g <- summary(tmp.surv)$time 
 
      # Median Survival time in g-th group 
      tmp.med <- SurvMed(tmp.time, tmp.death)[1]

      ## If cell median does not exist, then 
      if (is.na(tmp.med) == TRUE){
        cat('      ','Median Survival time in',g,'th cell does not exist !','\n')
        Fi.se[g] <- NA
      } else {    
        # Median survival time in gth group
        theta.med[g] <- tmp.med
  
        # Find the index of survival time which is the closest to the total.med
        index.1 <- which(abs(theta.g - total.med) == min(abs(theta.g - total.med)))[1]  
    
        # Evaluate F_i(total.med)
        eta.hat[g] <- F.i[index.1]
      
        # Find the standard error of F_i(overall.med)   
        Fi.se.1 <- F.se[index.1] # Greenwood error
    
        # Find the index that is closest to the ith group.med among uncensored obs  
        censored <- which(tmp.death == 0) # obs that are censored

        if (sum(theta.med[g] == theta.g) != 0){
          ## CASE 1 : There is any survival time that is same with cell median 

          # Exclude observations who are same with theta.med[g] or censored
          del.index <- union(which(theta.g == theta.med[g]), censored) 
          tmp <- (theta.g - theta.med[g])[-del.index] 
          if (length(tmp) == 0){
            # If there is no observation near cell median
            Fi.se[g] <- Fi.se.1^2 
          } else {
            min.dff <- min(abs(tmp))
            time.dff <- abs(theta.g - theta.med[g])

            # Index of the nearlist theta within gth group  
            index.2 <- which(time.dff == min.dff)[1] 

            # Index of the group median
            index.3 <- which(theta.g == theta.med[g])[1] 
            Fi.se.2 <- (F.i[index.3]-F.i[index.2])^2 / 2  
            Fi.se[g] <- Fi.se.1^2 + Fi.se.2           
          }
        } else if (sum(censored) != 0){
          ## CASE 2 : There is no observation with same cell median, but censored.

          if ( length(theta.g - theta.med[g]) <= length(censored) ){
            ## CASE 2-1 : There are only censored observations left

            # No second component for variance exist.
            Fi.se.2 <- 0

          } else {
            ## CASE 2-1 : There are complete observations as well as censored observations  

            # Exclude observations who are censored
            min.dff <- min(abs((theta.g - theta.med[g])[-censored]))
            time.dff <- abs(theta.g - theta.med[g])  

            # Calculate index of the nearlist theta within gth group
            index.2 <- which(time.dff == min.dff)[1]   

            # Calculate index of the group median
            index.3 <- which(theta.g == theta.med[g])[1] 
            Fi.se.2 <- (F.i[index.3]-F.i[index.2])^2 / 2  
          } 
          Fi.se[g] <- Fi.se.1^2 + Fi.se.2   
        } else {
          ## CASE 3 : There is no observation with same cell median, nor censored obs.

          min.dff <- min(abs((theta.g - theta.med[g])))
          time.dff <- abs(theta.g - theta.med[g])

          # Index of the nearlist theta within gth group
          index.2 <- which(time.dff == min.dff)[1] 
          
          # Index of the group median
          index.3 <- which(theta.g == theta.med[g])[1] 
          Fi.se.2 <- (F.i[index.3]-F.i[index.2])^2 / 2  
          Fi.se[g] <- Fi.se.1^2 + Fi.se.2   
        }  
      }
    }   

    ## Update Fi.se... Do not consider cells that doesn't have median
    Fi.se <- na.omit(Fi.se)
    weight <- 1 / Fi.se
    hi <- weight / sum(weight)

    ## Calculate Statistics
    C <- sum(weight * (eta.hat - sum(hi*eta.hat))^2 ) 
    return(C)

  } else {
    cat('Overall median survival time does not exist !!','\n')
    return('PASS')
  } 
}

###########################################################################
# For given snp.dat, Implement MDR approach and return Median.statistics. #
###########################################################################

Median.MDR <- function(snp.dat.part, time, death){
  death <- death
  time <- time
  snp.dat <- as.matrix(snp.dat.part)
  k <- ncol(snp.dat.part) # Number of SNP factors... likely to be 2 in this case.
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

  ## Find Median Survival time
  total.med <- SurvMed(time, death)[1] 
  
  ## If median exist, then implement MDR. If not, PASS
  if ( is.na(total.med) == FALSE ){
    # Number of observation without missing
    n <- length(fids)   

    ## split data into 9-cells
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

    ## A binary variable which has 2 if obs is in high risk, 1 otherwise. 
    G <- rep(1, n)
  
    ## Find survival times within each cells based on KM method. 
    Cell.membership <- rep(NA, n)
    for(cel in 1:cells.num){
      temp.ids <- cells[[cel]][,1] # cell ids in cel-th cell
      Cell.membership[temp.ids] <- cel 
      if(length(temp.ids) != 0){
        # Fit KM survival model within cel-th cell
        cell.fit <- survfit(Surv(time[temp.ids], death[temp.ids]) ~ 1) 
        # The median survival time in cel-th cell
        cell.med <- unname(as.matrix(summary(cell.fit)$table)["median",1])
        # IF cell median is not NA  
        if ( is.na(cell.med) == FALSE ){ 
          if(total.med > cell.med){        
            # Observations in cel-th cell are listed as high.risk members!
            G[temp.ids] <- 2  
          }
        }
      }
    }

    ## Empty G should be deleted... because G contains IDs that are not complete
    G <- G[is.na(Cell.membership) == FALSE]
   
    if (length(levels(as.factor(G))) == 1){ 
      # If observations in training sets are all in 'High' or all in 'Low'
      median.stat <- NA
    } else {	
      # Calculate Median test  
      median.stat <- CZ.Median.test(time, death, G)
    }   
    result$stat <- median.stat
    result$high.id <- which(G == 2)
    result$low.id <- which(G == 1)
    result$Cell.id <- Cell.membership
  } else {
    result$stat <- 'PASS'
  }
  return(result)
}

#######################################
# Evaluate k=2 model using 10-fold CV #
#######################################

CV.Median.MDR <- function(snp.dat.part, time, death, nfolds){
  snp.dat.part <- snp.dat.part
  time <- time
  death <- death
  K <- ncol(snp.dat.part)
  nfolds <- nfolds
  N <- length(time)
  
  ## split the whole data into folds
  cv.len <- floor(N/nfolds) 
  obs.id <- 1:N
  CV.result <- numeric(nfolds)

  ## Shuffle the data before assign the observations... DO NOT Replace the sample 
  shuffle.id <- sample(1:N, N, replace=FALSE)   
  time <- time[shuffle.id] 
  death <- death[shuffle.id] 
  snp.dat.part<- snp.dat.part[shuffle.id,]
  
  ## start CV procedure
  CV.result <- numeric(nfolds)
  for (cv in 1:nfolds){
    test.id <- ((cv - 1) * cv.len + 1) : (cv * cv.len) 
    train.time <- time[-test.id]
    train.death <- death[-test.id]
    train.snp.dat <- snp.dat.part[-test.id,]
    tmp <- Median.MDR(train.snp.dat, train.time, train.death)
    CV.result[cv] <- tmp$stat
  }
  index <- which(CV.result == NA)
  index <- c(index, which(CV.result == 'PASS'))
  
  if ( length(index) == nfolds ){
    cat('All results are not available.. this combination is useless','\n')
    result <- NA
  } else if ( length(index) != 0 ){
    # Delete Results that are not available
    CV.result <- CV.result[-index]
    result <- mean(CV.result)
  } else {
    result <- mean(CV.result)     
  }
  return(result)
}

#########################################################
# Find best combination of k=2 model among all SNP data #
#########################################################

Median.MDR.Model <- function(snp.dat, time, death, nfolds, k=2){
  snp.dat <- snp.dat # Enter total matrix
  time <- time
  death <- death
  k <- k
  numK <- ncol(snp.dat)
  nfolds <- nfolds
  comb.mat <- combn(numK,2)
  ncombi <- choose(numK,2)
  combi.vec <- numeric(ncombi)
  result <- list() 
  
  ## Save CV result from all combination
  for (combi in 1:ncombi){
    tmp.col <- comb.mat[ ,combi]
    tmp.mat <- snp.dat[ ,tmp.col]
    tmp <- CV.Median.MDR(tmp.mat, time, death, nfolds)
    combi.vec[combi] <- tmp
  }
  index <- is.na(combi.vec) == TRUE
  if (sum(index) != 0 ){
    if ( sum(index) == ncombi ){
      cat('All results are not available.. this combination is useless','\n')
      combi.vec <- NA
    } else {
      # Delete Results that are not available
      combi.vec <- combi.vec[-which(index)]
    }
  }

  ## Find best model
  if ( length(combi.vec) == 0 ){
    best.snp <- NA
    best.result <- NA
  } else {
    best.index <- which(combi.vec == max(combi.vec))
    best.snp <- snp.dat[,comb.mat[,best.index]]
    best.result <- Median.MDR(best.snp, time, death)
  }
  result$best.snp <- best.snp
  result$best.med.mdr <- best.result
  return(result)
}  