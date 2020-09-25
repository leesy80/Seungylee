## p = Pr(T > C) : Censoring fraction
lam <- 2 ; alpha <- 5 
N <- 10000
time <- rweibull(N, shape = alpha, scale = lam)
# theta <- 17.9 # 0.1
# theta <- 9.0 # 0.2 
theta <- 5.9 # 0.3
C <- runif(N,0,theta)
sum(time > C) / N ; cat('theta =',theta,'\n')

## Log - Normal Distribution ##
time <- rlnorm(N, meanlog = 0, sdlog = 1)
# theta <- 15.7 # 0.1
# theta <- 4.7 # 0.3 
theta <- 2.35 # 0.5

C <- runif(N,0,theta)
sum(time > C) / N ; cat('theta =',theta,'\n')

## Cox with Weibull(5,2) Distribution ##
u <- runif(N, min=0, max=1)
dim(u) <- c(N,1)
time <- ( -log(1-u) / (scale * exp(dTrt[,1]*beta + z*gamma + delta*nG[,1])) )^(1/shape) # nG[,1] : first non-causal SNP (SNP3)
theta <- 2.35 # 0.5
C <- runif(N,0,theta)
sum(time > C) / N ; cat('theta =',theta,'\n')
