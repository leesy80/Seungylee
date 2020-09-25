CalType1 <- function(N, nsnp=2, trait.type = 'time', time.type='PH', cls.method, model,  
            adj.main, cova, maf, niter, alpha = 0.05, cp, gamma, beta = 0.0, delta, sigma = 1){
  ############################################
  ## calculate type I error -- no causal SNPs
  ## input: --  N : sample size
  ## ---------- nsnp : number of snps
  ## ---------- trait.type : type of trait, either "quantitative" or "binary"
  ## ---------- cls.method : classification of H/L, either "mean" or "obs-ratio"
  ## ---------- cova : covariates
  ## ---------- adj.main : adjust main effect of not
  ## ---------- maf : 0.2 or 0.4
  ## ---------- niter : total number of run for each penatrace model
  ## ---------- alpha : significant level
  ## ---------- nperm : small number of permutation times to estimate non-#central parameter
  ## output : the type I error and the raw type I error (without non-central correction)
  ############################################

  K <- 2 # 2 way interaction, first 2 are causal
  output <- test <- test2 <- cf <- NULL
  if ( maf %in% c(0.2, 0.4) ){
    nt <- sum(MAF == maf)
    ntt <- sample(1:nt, 1, replace=T)
    model.no <- which(MAF == maf)[ntt] # Randomly choose one among MAFs
  } else {
    model.no <- sample(1:dim(pen)[2],1,replace=T) # Randomly choose one
  }
  error.raw <- error <- stats <- rep(0, niter)
  for(run in 1:niter){
    ## Generate data under the null
    snp.mats <- dat[,3:4] # Only 2 Causal SNPs were generated but not correlated with time 
    if (gamma == 0){
      cova <- NULL
    } else {
      cova <- dat[,5] # fifth column of dat is Z
    }
    p <- ncol(snp.mats)
    snp.combs <- combn(p, K)   ## all SNP pairs
    ns <- ncol(snp.combs)
    phe <- dat[ ,2:1] ## dat[,1] : censoring indicator, dat[,2] : survival time 
    cf[run] <- (N-sum(phe[,2]))/N
    if ( model %in% c('umcox','umaft') ){
      res <- UM.MDR(snp.mats, phe, cova, cls.method, model, adj.main, nperm)
    } 
    test[run] <- res$p.values[1]
    test2[run] <- res$rpv[1]
    error[run] <- ifelse(res$p.values[1] < alpha/ns, 1, 0) # error under bonferoni correction
    ## raw chi-square statistics ##
    stats[run] <- res$stat[1]
    error.raw[run] <- ifelse(res$rpv[1] < alpha/ns, 1, 0)    ## error (without non-central correction)
  }
  merror <- mean(error) ; merror.raw <- mean(error.raw) ; mcf <- mean(cf)  
  output <- rbind(output, c(maf, merror,  merror.raw, mcf))
  colnames(output) <- c("maf", "error", "error.raw","cf")
  out <- list()
  out$output <- output
  out$test <- test
  out$test2 <- test2
  return(out)
}
CalType1 <- cmpfun(CalType1)

#------------------------------------------

CalPower2 <- function(N, nsnp, nway, trait.type = 'time', time.type, cls.method, model, adj.main, niter, alpha, 
                        cand.pmodels = i:i, cp, mu = 1, gamma, beta, delta, sigma, Balanced){
  set.seed(100)
  K <- nway # K-way interaction, first K are causal  
  output <- NULL
  cf <- NULL
  ##########################################################
  #           # Available options for 'cls.method' (Step 1):
  #
  # 1. cls.method == 'aft' : Use Aft residuals from survreg WITHOUT X to define H/L
  # 2. cls.method == 'aft2' : Use Aft residuals from survreg WITH X to define H/L
  # 3. cls.method == 'mg' : Use martingale residuals from coxph WITHOUT X to define H/L
  # 4. cls.method == 'mg2' : Use martingale residuals from coxph WITH X to define H/L
  # 5. cls.method == 'km' : Use Median Survival time (1 vs 8) to define H/L
  #
  #           # Available options for 'method' (Step 2):  
  # 
  # 1. model == 'umaft' : UM.MDR with Aft standardized residuals
  # 2. model == 'umcox' : UM.MDR with Martingale residuals
  # 3. model == 'aftmdr' : MDR with Aft standardized residuals
  # 4. model == 'coxmdr' : MDR with Martingale residuals
  # 5. model == 'kmmdr' : MDR with Cell Median
  #
  #           # Available options for 'adj.main' (Marginal Effect of SNP1, SNP2) :  
  # 
  # 1. adj.main == TRUE : SNP3 is included in data generating step
  # 2. adj.main == FALSE : SNP3 is not included in data generating step
  #
  #           # Available options for 'time.type' (Data Generating):  
  #
  # 1. time.type = 'PH' : Use coxph model to generate time
  # 2. time.type = 'LogN' : Use AFT model with error term following lognormal
  # 3. time.type == 'Weibull' : Use AFT model with error term following Weibull
  #
  #           # Available options for 'parameters' (Data Generating):  
  #
  # 1. gamma : Coefficient for covariate Z
  # 2. beta : Coefficient for H/L Indicator S
  # 3. delta : Coefficient for marginal SNP3
  # 4. sigma : Standard deviation of error term
  #
  ##########################################################

  ## Implement power Caculation for cand.pmodels = 70 models ##
  for(number in cand.pmodels){
    maf <- MAF[number]
    Pow.mdr <- Pow.rank <- Pow.Bonf <- pvalues <- smallest.combi <- rep(0, niter)
    cat('Starting',number,'th Penetraing function...','\n') 

    ## Generate data under the null hypothesis ##
    for (iter in 1:niter){
      ## Simulation Dataset is designed to have SNP1, SNP2 as a significant pairs ##  
      if (Balanced == TRUE){
        tmp <- simu_popu(N, M1 = K, M0 = nsnp-K, p.f1 = rep(maf, K), p.f0 = maf,
                         pen1 = pen[, number], mu1 = 1, trait.type, time.type='PH', cp, gamma, beta, delta, sigma = 1)
        prop <- sum(tmp$HL) / N
        if (prop == 0){
          prop <- 0.001
        }
        ## Start to generate the psuedo dataset ##
        NN <- floor(N/(2*prop) + 1)
        tmp <- simu_popu(NN, M1 = K, M0 = nsnp-K, p.f1 = rep(maf, K), p.f0 = maf,
                         pen1 = pen[, number], mu1 = 1, trait.type, time.type='PH', cp, gamma, beta, delta, sigma = 1)
        tmp.id <- which(tmp$HL == 1) ; tmp.L <- which(tmp$HL == 0)
        H.N <- sum(tmp$HL)
        if (H.N > N/2){
          tmp.id <- tmp.id[1:floor(N/2)]
        } else if (H.N < N/2) {
          d <- floor(N/2 - H.N + 1)
          tmp.id <- c(tmp.id, sample(tmp.id, d, replace=TRUE)) 
        }
        tmp.dat.H <- tmp$data[tmp.id,]
        L.N <- N - length(tmp.id)  
        if (length(tmp.L) >= L.N){    
          tmp.LL <- sample(tmp.L, L.N, replace = FALSE) 
        } else {
          tmp.LL <- sample(tmp.L, L.N, replace = TRUE)
        }
        tmp.dat.L <- tmp$data[tmp.LL,]
        dat <- rbind(tmp.dat.H, tmp.dat.L)
      } else {
        dat <- simu_popu(N, M1 = K, M0 = nsnp-K, p.f1 = rep(maf, K), p.f0 = maf,
                         pen1 = pen[, number], mu1 = 1, trait.type, time.type='PH', cp, gamma, beta, delta, sigma = 1)$data
      }
      phe <- as.matrix(dat[, 2:1]) # phe[,1] : Time, phe[,2] : Censoring Indicator
      cf[iter] <- (N - sum(phe[,2])) / N # Empirical Censoring fraction 
      pp <- dim(dat)[2] # Number of columns of dat
      snp.mats <- dat[, 3:(pp-1)] # Extract SNP
      if ( gamma != 0 ){
        cova <- dat[,pp] # Covariate : Z 
      } else {  
        cova <- rep(1, N)
      }
      p <- dim(snp.mats)[2] # Number of SNPs to be considered
      snp.combs <- combn(p, K) # Specify the all SNP pairs
      ncombi <- ncol(snp.combs) # Number of candidate combinations
      if ( model %in% c('umcox', 'umaft') ){
        ## Implement corresponding UM.MDR ##
        i.result <- UM.MDR(snp.mats, phe, cova, cls.method, nway = K, model, adj.main, nperm = 20)
        pvalues[iter] <- i.result$p.values[1] # i.result$p.values[1] represents the raw p-value of S created by SNP1 and SNP2
                                              # This implies correct selection 
        smallest.combi[iter] <- which.min(i.result$p.values)
        cat('   Operating Simulation..',iter,'th iteration is done','\n')
      } else { 
        ## Implement MDR procedure ##
        if ( model == 'aftmdr' ){
          i.result <- CV.Aft.MDR(snp.mats, phe, X = cova, nway = K, method = 'BA', nfold = 10)
        }
        if ( model == 'survmdr' ){
          i.result <- CV.Surv.MDR(snp.mats, phe, nway = K, nfold = 10)
        } 
        if ( model == 'coxmdr' ){
          i.result <- CV.Cox.MDR(snp.mats, phe, X = cova, nway = K, method = 'BA', nfold = 10)
        } 
        if ( model == 'kmmdr' ){
          i.result <- CV.Median.MDR(snp.mats, phe, nway = K, nfold = 10)
        }
        best.pair <- as.numeric(i.result$best.combi.No)      
        Pow.mdr[iter] <- ifelse(best.pair == 1, 1, 0) # Correct if the best.combi.No is '1' (means SNP1 and SNP2) 
        cat('   Operating Simulation..',iter,'th iteration is done',"best.pair is", best.pair,'\n')
      }
    }

    ## Summarize Powers ##
    if ( model %in% c('umcox', 'umaft') ){
      Pow.Bonf[which(pvalues <= alpha/ncombi)] <- 1 # Bonferroni correction is required  
      Pow.rank[which(smallest.combi == 1)] <- 1 # If Smallest p.value happens at 1, then success
    }
    cat(number,' th Penetrating function is done','\n')
    mcf <- mean(cf) ; mPow.mdr <- mean(Pow.mdr) ; mPow.Bonf <- mean(Pow.Bonf) ; mPow.rank <- mean(Pow.rank)
    tables <- rbind(output, c(number, mPow.Bonf, mPow.rank, mPow.mdr, mcf))
  }
  colnames(tables) <- c("number", "mPow.Bonf", "mPow.rank", "mPow.mdr", "cf")
  output <- list()
  output$tables <- tables
  return(output)
}
CalPower2 <- cmpfun(CalPower2)

#---------------------------------------------------


CalPower3 <- function(N, nsnp, nway, trait.type = 'time', time.type, cls.method, model, adj.main, niter, alpha, 
                        cand.pmodels = i:i, cp, mu = 1, gamma, beta, delta, sigma, Balanced = TRUE){
  set.seed(1)
  K <- nway # K-way interaction, first K are causal  
  output <- NULL
  cf <- NULL
  ##########################################################
  #           # Available options for 'cls.method' (Step 1):
  #
  # 1. cls.method == 'aft' : Use Aft residuals from survreg WITHOUT X to define H/L
  # 2. cls.method == 'aft2' : Use Aft residuals from survreg WITH X to define H/L
  # 3. cls.method == 'mg' : Use martingale residuals from coxph WITHOUT X to define H/L
  # 4. cls.method == 'mg2' : Use martingale residuals from coxph WITH X to define H/L
  # 5. cls.method == 'km' : Use Median Survival time (1 vs 8) to define H/L
  #
  #           # Available options for 'method' (Step 2):  
  # 
  # 1. model == 'umaft' : UM.MDR with Aft standardized residuals
  # 2. model == 'umcox' : UM.MDR with Martingale residuals
  # 3. model == 'aftmdr' : MDR with Aft standardized residuals
  # 4. model == 'coxmdr' : MDR with Martingale residuals
  # 5. model == 'kmmdr' : MDR with Cell Median
  #
  #           # Available options for 'adj.main' (Marginal Effect of SNP1, SNP2) :  
  # 
  # 1. adj.main == TRUE : SNP3 is included in data generating step
  # 2. adj.main == FALSE : SNP3 is not included in data generating step
  #
  #           # Available options for 'time.type' (Data Generating):  
  #
  # 1. time.type = 'PH' : Use coxph model to generate time
  # 2. time.type = 'LogN' : Use AFT model with error term following lognormal
  # 3. time.type == 'Weibull' : Use AFT model with error term following Weibull
  #
  #           # Available options for 'parameters' (Data Generating):  
  #
  # 1. gamma : Coefficient for covariate Z
  # 2. beta : Coefficient for H/L Indicator S
  # 3. delta : Coefficient for marginal SNP3
  # 4. sigma : Standard deviation of error term
  #
  ##########################################################

  ## Implement power Caculation for cand.pmodels = 70 models ##
  for(number in cand.pmodels){
    maf <- MAF[number]
    Pow.10.mdr <- Pow.30.mdr <- Pow.mdr <- Pow.rank <- Pow.Bonf <- pvalues <- smallest.combi <- rep(0, niter)
    cat('Starting',number,'th Penetraing function...','\n') 

    ## Generate data under the null hypothesis ##
    for (iter in 1:niter){
      ## Simulation Dataset is designed to have SNP1, SNP2 as a significant pairs ##  
      if (Balanced == TRUE){
        tmp <- simu_popu(N, M1 = K, M0 = 0, p.f1 = rep(maf, K), p.f0 = maf,
                         pen1 = pen[, model.no], mu1 = 1, trait.type, time.type='PH', cp, gamma, beta, delta, sigma = 1)
        prop <- sum(tmp$HL) / N
        if (prop == 0){
          prop <- 0.001
        }
        ## Start to generate the psuedo dataset ##
        NN <- floor(N/(2*prop) + 1)
        tmp <- simu_popu(NN, M1 = K, M0 = 0, p.f1 = rep(maf, K), p.f0 = maf,
                         pen1 = pen[, model.no], mu1 = 1, trait.type, time.type='PH', cp, gamma, beta, delta, sigma = 1)
        tmp.id <- which(tmp$HL == 1) ; tmp.L <- which(tmp$HL == 0)
        H.N <- sum(tmp$HL)
        if (H.N > N/2){
          tmp.id <- tmp.id[1:floor(N/2)]
        } else if (H.N < N/2) {
          d <- floor(H.N - N/2 + 1)
          tmp.id <- c(tmp.id, sample(tmp.id, d, replace=TRUE)) 
        }
        tmp.dat.H <- tmp$data[tmp.id,]
        L.N <- N - length(tmp.id)
        if (length(tmp.L) >= L.N){    
          tmp.LL <- sample(tmp.L, L.N, replace = FALSE) 
        } else {
          tmp.LL <- sample(tmp.L, L.N, replace = TRUE)
        }
        tmp.dat.L <- tmp$data[tmp.LL,]
        dat <- rbind(tmp.dat.H, tmp.dat.L)
      } else {
        dat <- simu_popu(N, M1 = K, M0 = 0, p.f1 = rep(maf, K), p.f0 = maf,
                         pen1 = pen[, model.no], mu1 = 1, trait.type, time.type='PH', cp, gamma, beta, delta, sigma = 1)$data
      }
      phe <- as.matrix(dat[, 2:1]) # phe[,1] : Time, phe[,2] : Censoring Indicator
      cf[iter] <- (N - sum(phe[,2])) / N # Empirical Censoring fraction 
      pp <- dim(dat)[2] # Number of columns of dat
      snp.mats <- dat[, 3:(pp-1)] # Extract SNP
      if ( gamma != 0 ){
        cova <- dat[,pp] # Covariate : Z 
      } else {  
        cova <- NULL
      }
      p <- dim(snp.mats)[2] # Number of SNPs to be considered
      snp.combs <- combn(p, K) # Specify the all SNP pairs
      ncombi <- ncol(snp.combs) # Number of candidate combinations
      if ( model %in% c('umcox', 'umaft') ){
        ## Implement corresponding UM.MDR ##
        i.result <- UM.MDR(snp.mats, phe, cova, cls.method, model, adj.main, nperm = 20)
        pvalues[iter] <- i.result$p.values[1] # i.result$p.values[1] represents the raw p-value of S created by SNP1 and SNP2
                                              # This implies correct selection 
        smallest.combi[iter] <- which.min(i.result$p.values)
        cat('   Operating Simulation..',iter,'th iteration is done','\n')
      } else { 
        ## Implement MDR procedure ##
        if ( model == 'aftmdr' ){
          i.result <- CV.Aft.MDR(snp.mats, phe, X = cova, nway = K, method = 'BA', nfold = 10)
        }
        if ( model == 'survmdr' ){
          i.result <- CV.Surv.MDR(snp.mats, phe, nway = K, nfold = 10)
        } 
        if ( model == 'coxmdr' ){
          i.result <- CV.Cox.MDR(snp.mats, phe, X = cova, nway = K, method = 'BA', nfold = 10)
        } 
        if ( model == 'kmmdr' ){
          i.result <- CV.Median.MDR(snp.mats, phe, nway = K, nfold = 10)
        }
        ## Extract ncombi-numbers of train scores ##
        score.list <- i.result$train.score
        Pow.10.mdr[iter] <- ifelse(rank(-score.list)[1] <= 10, 1, 0)
        Pow.30.mdr[iter] <- ifelse(rank(-score.list)[1] <= 30, 1, 0)
        best.pair <- as.numeric(i.result$best.combi.No)      
        Pow.mdr[iter] <- ifelse(best.pair == 1, 1, 0) # Correct if the best.combi.No is '1' (means SNP1 and SNP2) 
        cat('   Operating Simulation..',iter,'th iteration is done',"best.pair is", best.pair,'\n')
      }
    }

    ## Summarize Powers ##
    if ( model %in% c('umcox', 'umaft') ){
      Pow.Bonf[which(pvalues <= alpha/ncombi)] <- 1 # Bonferroni correction is required  
      Pow.rank[which(smallest.combi == 1)] <- 1 # If Smallest p.value happens at 1, then success
      mcf <- mean(cf) ; mPow.mdr <- mean(Pow.mdr) ; mPow.Bonf <- mean(Pow.Bonf) ; mPow.rank <- mean(Pow.rank)
      tables <- rbind(output, c(number, mPow.Bonf, mPow.rank, mPow.mdr, mcf))
      colnames(tables) <- c("number", "mPow.Bonf", "mPow.rank", "mPow.mdr", "cf")
    } else {
      mcf <- mean(cf) ; mPow.mdr <- mean(Pow.mdr) ; mPow.10.mdr <- mean(Pow.10.mdr) ; mPow.30.mdr <- mean(Pow.10.mdr)
      tables <- rbind(output, c(number, mPow.mdr, mPow.10.mdr, mPow.30.mdr, mcf))
      colnames(tables) <- c("number", "mPow.mdr", "mPow.10.mdr", "mPow.30.mdr", "cf")
    } 
    cat(number,' th Penetrating function is done','\n')
  }
  output <- list()
  output$tables <- tables
  return(output)
}
CalPower3 <- cmpfun(CalPower3)

#---------------------------------------------------