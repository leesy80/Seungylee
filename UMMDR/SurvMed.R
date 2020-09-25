##########################################
# SurvMed : Returns median time and S(t) #
##########################################
 
SurvMed <- function(time, death){
  time <- time 
  death <- death
  ID <- 1:length(death)
  St <- 1 ; d <- 0 ; c <- 0
  min.time <- NA

  while (St >= 0.5){
    if (length(time) == 0 ){
      break
    }
    min.time <- min(time)

    # Observe the smallest time
    del.index <- which(time == min.time)
    if (length(del.index) == 0 ){
      min.time <- St <- NA
      break
    }

    # Calculate Y : Number of obs who survive longer than min.time
    Y <- sum(time >= min.time)

    # Calculate d and c
    d.mat <- death[del.index]
    if (length(del.index) == 1){
      d <- sum((d.mat == 1)) # Number of events
      c <- sum((d.mat == 0)) # Number of censoring
    } else {
      d <- length(which(d.mat == 1)) # Number of events
      c <- length(which(d.mat == 0)) # Number of censoring
    }

    # Update datamat 
    St <- St*(1-d/Y)
    time <- time[-del.index]
    death <- death[-del.index]
  }
  return(c(min.time, St))
}