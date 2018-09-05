
## classification step, assign high/low for each cell and test 
class_HL <- function(snp.dat, SS, classm = 'mg',cova=NULL){
  snp.dat = as.matrix(snp.dat)
  
  # remove missing values
  fids = which(complete.cases(snp.dat))
  snp.dat = as.matrix(snp.dat[fids, ])
  
  k = ncol(snp.dat) 
  n = nrow(snp.dat)
  
  # split completed data into cells
  tlist = vector('list', k)
  for(i in 1:k) tlist[[i]] = snp.dat[, i]
  cells = split(data.frame(cbind(fids, snp.dat)), tlist)
  
  # delete NULL cells
  obs.cell = sapply(cells, function(x) nrow(x))
  cell.null = which(obs.cell == 0)
  if(length(cell.null) > 0 ) cells = cells[-cell.null]
  
  ncells = length(cells)
  high.ids  = NULL
  weights = rep(1, n)  
  
  # classify each cell into H/L by various methods
  if(classm == 'mg'){
    # use mg to define H/L, corresponding to time trait
    library(splines)
    library(survival)

     Ycen <- Surv(SS[,2], SS[,1]==1)
     
     cox <- coxph(Ycen ~ 1)
     mtresid <- residuals(cox, type="martingale")

     for(i in 1:ncells){
      temp.ids = cells[[i]][, 1]
      if(length(temp.ids) == 0) next
      mSS = sum(mtresid[temp.ids])
      weights[temp.ids] = (length(temp.ids))
      if (mSS >= 0){
        high.ids = c(high.ids, temp.ids)
        }
     }
  }

  low.ids = setdiff(fids, high.ids)
  
  # labels indicate high/low
  labels = rep(NA, n)
  labels[high.ids] = 1
  labels[low.ids] = 0

  # note that labels are NA, corresponding the obs that having NA for SNP
  return(list('labels' = labels, 'weights' = weights)) 
}
class_HL <- cmpfun(class_HL)
#--------------------------------------------


## estimat non-center parameter of the null distribution due to classification
est_nullcenter <- function(phe, cova=NULL, snp.pair, classm = 'mg', lbd = 0, nperm = 10, adj.main ='FALSE'){
  library(splines)
  library(survival)
  n = nrow(phe)
  center = rep(0, nperm)
  
  for(i in 1:nperm){
    perm=sample(1:n, n)
    perm.phe = phe[perm,]
    if(length(cova)==n){perm.cova = cova[perm]}
    else{perm.cova = cova[perm,]
         k=ncol(cova)}
    res = class_HL(snp.pair, perm.phe, classm='mg')
    labels = res$labels    
    Ycen <- Surv(perm.phe[,2], perm.phe[,1]==1)

    if(adj.main == 'FALSE'){
     xx = cbind(labels,perm.cova)
     sumstat = summary(coxph(Ycen ~ xx))$coef
     center[i] = (sumstat[1, 1]/sumstat[1, 3])^2
    } 
 
    else{
     xx = cbind(labels,perm.cova[,1:k-2],cova[,k-1:k])
     sumstat = summary(coxph(Ycen ~ xx))$coef
     center[i] = (sumstat[1, 1]/sumstat[1, 3])^2
    }
  }
  return(mean(center))
}
est_nullcenter = cmpfun(est_nullcenter)

#----------------------------------


## the non-zero center parameter is estimated by a few permutation 
## snp.all ----snp matrix, n by p
## snp.combs ---all snp pairs 2 by (p choose 2)
## phe --- phenotype
## cova -- covariate
## classm -- classification rule of H/L
## adj.main -- adjust marginal effect or not
## nperm ----- small number of permutation time to estimate non-central 
#parameters
## Output: --- a list including the coefficient (beta), pvalue (pv),
## ----------- raw pvalue (rpv) and statistics(stat) for all snp pairs
coxumMDR <- function(snp.all, snp.combs, phe, cova = NULL, classm = 'mg', adj.main = 'TRUE', nperm = 10 ){
  set.seed(1234)
  ns = ncol(snp.combs)
  n = nrow(phe)
  d = ncol(phe)
  pv<-NULL
  rpvs = pvs = coefs =stat= rep(0, ns)
  cova<-as.matrix(cova)
  # save all coefs and pvs for S
  
  SS = phe
 
  
  # select best model(i.e. snp combination)
  # use glm
  
   for(j in 1:ns){
    res = class_HL(snp.all[, snp.combs[, j]], SS, classm)
    labels = res$labels
   
    if(adj.main == 'TRUE'){
      snp1 = snp.all[, snp.combs[1, j]]
      snp2 = snp.all[, snp.combs[2, j]]
      xx = cbind(cova,snp1,snp2)
      
      # use coxph
      
      library(splines)
      library(survival)
      Ycen <- Surv(SS[,2], SS[,1]==1)
      sumstat= summary(coxph(Ycen ~ labels+xx))$coef
      print(sumstat)
      coefs[j]=sumstat[1,1]
      stat[j] = ((sumstat[1,1])^2)/((sumstat[1,3])^2)
     
      rpvs[j] =pchisq(stat[j], df=1,ncp=0, lower.tail = FALSE )
      
      null.center = est_nullcenter(SS, xx, snp.all[, snp.combs[, j]], classm,lbd=0, nperm, adj.main = 'TRUE')
      
      #lambda = null.center
      lambda = max(null.center-1,0)

      pvs[j] =pchisq(stat, df=1,ncp= lambda, lower.tail = FALSE )
    
    }else{
       
      xx  = cbind(labels,cova)
      library(splines)
      library(survival)
      Ycen <- Surv(SS[,2], SS[,1]==1)
      sumstat= summary(coxph(Ycen ~ xx))$coef
      coefs[j]=sumstat[1,1]
      stat[j] = ((sumstat[1,1])^2)/((sumstat[1,3])^2)
   
      rpvs[j] = pchisq(stat[j], df=1,ncp=0, lower.tail = FALSE )
      
      null.center = est_nullcenter(SS, cova, snp.all[, snp.combs[, j]],classm, lbd = 0, nperm,adj.main = 'FALSE')
      
      lambda = max(null.center-1,0)
    
      pvs[j] = pchisq(stat, df=1,ncp= lambda, lower.tail = FALSE )
      
    }
  }  
  return(list( 'coef'=coefs,'stat'=stat,'rpv' = rpvs,'pv' = pvs))

}
coxumMDR <- cmpfun(coxumMDR)




