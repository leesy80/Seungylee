##################################
# Numorator calcuator for C.Stat #
##################################

C.stat <- function(vec){
  vec <- vec
  tmp <- vec[5] - vec[2]/vec[1]*vec[4]
  return(unname(tmp))
}

####################################
# Denomorator calcuator for C.Stat #
####################################

V.stat <- function(vec){
  vec <- vec
  tmp <- 0 
  if ( vec[1] != 0 & vec[1] != 1 ){
    tmp <- vec[2]*vec[3]*vec[4]*(vec[1]-vec[4]) / ( (vec[1]-1)*vec[1]^2 )
  }
  return(unname(tmp))
}

######################################################
# LogRank Statistic Calculator.. Returns Chisq value #
######################################################

SurvDif.LogRank <- function(time, death, G){
  time <- time
  death <- death  
  G <- G 

  ## Only uncensored times will be considered as event time
  complete.time <- time[which(death==1)]
  complete.time.1 <- time[intersect(which(death==1), which(G==1))] 

  ## Create time coordinates
  time.ref <- unique(sort(complete.time))
  time.1 <- na.omit( time[which(G == 1)] )
  time.2 <- na.omit( time[which(G == 2)] )
  J <- length(time.ref)

  ## Create TimeTable
  datamat <- matrix(0,J,6)
  for (j in 1:J){
    datamat[j,1] <- sum( time >= time.ref[j] )
    datamat[j,2] <- sum( time.1 >= time.ref[j] )

    # Only uncensored events should be counted!
    datamat[j,4] <- sum( complete.time == time.ref[j] )
    datamat[j,5] <- sum( complete.time.1 == time.ref[j] )
  } 
  datamat[,3] <- datamat[,1] - datamat[,2]
  datamat[,6] <- datamat[,4] - datamat[,5]

  ## Calculate LogRank Statistics
  num <- apply(datamat, 1, C.stat)
  denom <- apply(datamat, 1, V.stat)
  result <- (sum(num))^2 / sum(denom)
  return(result)
}

############################################
# C-Statistic Calculator.. Returns C-value #
############################################

SurvDif.C <- function(time, death, G){
  time <- time 
  death <- death 
  G <- G

  # Only uncensored times will be considered as event time #
  complete.time <- time[which(death==1)]
  complete.time.1 <- time[intersect(which(death==1), which(G==1))] 

  # Create time coordinates #
  time.ref <- unique(sort(complete.time))
  time.1 <- time[which(G == 1)]
  time.2 <- time[which(G == 2)]
  J <- length(time.ref)

  # Create TimeTable #
  datamat <- matrix(0,J,6)
  for (j in 1:J){
    datamat[j,1] <- sum( time >= time.ref[j] )
    datamat[j,2] <- sum( time.1 >= time.ref[j] )

    # Only uncensored events should be counted! #
    datamat[j,4] <- sum( complete.time == time.ref[j] )
    datamat[j,5] <- sum( complete.time.1 == time.ref[j] )
  } 
  datamat[,3] <- datamat[,1] - datamat[,2]
  datamat[,6] <- datamat[,4] - datamat[,5]

  # Calculate C-Statistics #
  result <- apply(datamat, 1, C.stat)
  return(sum(result))
}

Renyi.Stat <- function(time, death, G){
  time <- time
  death <- death
  G <- G 

  # Only uncensored times will be considered as event time #
  complete.time <- time[which(death==1)]
  complete.time.1 <- time[intersect(which(death==1), which(G==1))] 

  # Create time coordinates #
  time.ref <- unique(sort(complete.time))
  time.1 <- na.omit( time[which(G == 1)] )
  time.2 <- na.omit( time[which(G == 2)] )
  J <- length(time.ref)

  # Create TimeTable #
  datamat <- matrix(0,J,6)
  colnames(datamat) <- c('Y','Y1','Y2','D','D1','D2')
  for (j in 1:J){
    datamat[j,1] <- sum( time >= time.ref[j] )
    datamat[j,2] <- sum( time.1 >= time.ref[j] )

    # Only uncensored events should be counted #
    datamat[j,4] <- sum( complete.time == time.ref[j] )
    datamat[j,5] <- sum( complete.time.1 == time.ref[j] )
  } 
  datamat[,3] <- datamat[,1] - datamat[,2]
  datamat[,6] <- datamat[,4] - datamat[,5]
  Zt <- abs(cumsum( datamat[,5] - datamat[,2]*(datamat[,4]/datamat[,1])))
  sigsq <- sum( datamat[,2]/datamat[,1] * (datamat[,3]/datamat[,1]) * (datamat[,1] - datamat[,4])/(datamat[,1] - 1)*datamat[,4] )
  Q <- max(abs(Zt)) / sqrt(sigsq)
  return(Q)
}


most.freq <- function(x, n, freq){
  x <- x
  tmp <- sort(table(x),decreasing=TRUE)[1:n] 
  output <- names(tmp)
  if (freq == TRUE){
    output <- list()
    output$name <- names(tmp)
    output$freq <- unname(tmp)
  }
  return(output)
}