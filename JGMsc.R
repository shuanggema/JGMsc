########################################################################################################################
# JGMsc: Joint estimation for precision matrices with sign consistency
# Manuscript
#  Promote Sign Consistency in the Joint Estimation of Precision Matrices
#-----------------------------------------------------------------------------------------------------------------------
# Inputs
#   covarianceList: a list of sample covariance matrices 
#   nM: a vector whose elements are sample size for each datasets
#   p: dimension
#   lambda1_vector: a vector of lambda 1 values
#   lambda2: a single elements for lambda 2
#   a: value to be specified for the composite MCP
#   tol and maxIter: specifications on stoppoing criteria
#   weightsType: "weighted" for weighted by sample sizes and "equal" otherwise
#-----------------------------------------------------------------------------------------------------------------------
# Outputs
#   est:   a list of estimated precision matrices at the end of iterations
#          (length=number of lambda 1 tunings) 
#          each element of the list is another list whose elements are estimated matrices for indiviudal datasets
#   iters: a vector of number of iterations
#          (length=number of lambda 1 tunings)
#   time:  a vector of computation time 
#          (length=number of lambda 1 tunings)
#   NP:    a matrix of number of selected edges at the end of iterations
#          (row=number of lambda 1 tunings, col = number of datasets)
#   Lik:   a matrix of log-likelihood functions at the end of iterations
#          (row=number of lambda 1 tunings, col = number of datasets)
#   BIC:   a matrix of BIC at the end of iterations
#          (row=number of lambda 1 tunings, col = number of datasets)
########################################################################################################################


library(MASS)
library(Matrix)


JGMsc <- function(covarianceList,
                        nM,p,lambda1_vector, lambda2,
                        a=3,tau=1,u=1, weightsType="weighted",
                        tol=1E-7, maxIter =1000){
  
  #----------- * predefined functions * -------------------------#
  
  penSpMcp <- function(x,lambda,gamma){
    val <- (lambda*abs(x) -x^2/(2*gamma))*(abs(x)<(lambda*gamma)) +
      (1-(abs(x)<(lambda*gamma)))*(lambda*gamma)^2 /(2*gamma)
    return(val)
  }
  
  penSpMcp.Deriv <- function(x,lambda,gamma){
    val   <- (lambda -abs(x)/gamma)*(abs(x)<(lambda*gamma))
    return(val)
  }
  
  NormDiff <- function(xnew,xold){
    val <- sum((xnew-xold)^2)/max(c(1, sum(xnew^2), sum(xold^2)))
    return(val)
  }
  
  #----------- * import * -------------------------#
  
  M  = length(nM)
  if (weightsType=="equal")    weightsVal = rep(1,M)
  if (weightsType=="weighted") weightsVal= nM/sum(nM)
  
  n.tuning_1 = length(lambda1_vector)
  
  #----------- * initialize * -------------------------#
  
  time_cost=iters = rep(NA,times=n.tuning_1)
  NP = BIC = Lik= matrix(NA, nrow=maxIter,ncol=M)
  est=vector(mode = "list", length = n.tuning_1)
  
  precisionList  = zList = lapply(1:M,function(m) matrix(0,p,p))
  i_lam=n.tuning_1
  
  #----------- * iterations * -------------------------#
  while(i_lam>0){
    
    lambda1 = lambda1_vector[i_lam] 
    b = 0.5*M*lambda1^2*a
    
    iter       = 0
    diffVal    = 1000  
    LambdaList = lapply(1:M,function(m) matrix(0,p,p))
    precisionList  = zList = lapply(1:M,function(m) matrix(0,p,p))
    
    pt = proc.time()
    
    #----------- * update   * -------------------------#
    
    while((iter<10)||(iter <= (maxIter-1)) & (diffVal > tol)){
      
      iter             = iter + 1
      precisionList.pre = precisionList
      zList.pre         = zList
      
      # update precisionList
      
      for(m in 1:M){
        eigenComp          = eigen(covarianceList[[m]] + LambdaList[[m]]- u*zList[[m]])
        D                  = eigenComp$values
        V                  = eigenComp$vectors
        D2                 = (-D + sqrt(D^2 + 4*u))/(2*u)
        precisionList[[m]] = round((V %*% diag(D2) %*% t(V)),digits=10)
      }
      
      diffVal                      = max(sapply(1:M,function(m)  NormDiff(precisionList.pre[[m]], precisionList[[m]])))
      
      # update zList and  LambdaList
      ##### use (1-diag(1,p)) such that lambda 2 penalty has no effect on the diagonal.
      
      wList       = lapply(1:M, function(m) 1/sqrt(tau^2+precisionList[[m]]^2)*(1-diag(1,p)))
      wzList      = lapply(1:M, function(m) wList[[m]]*zList[[m]])
      penaltyList = lapply(1:M, function(m) penSpMcp(zList[[m]],lambda=lambda1,gamma = a))
      penalty_sum = Reduce('+',penaltyList)
      
      for (m in 1:M){
        a.m             = u * weightsVal[m] + 2* lambda2* (M-1) *(wList[[m]]^2)
        b.m             = u * weightsVal[m] * precisionList[[m]] + weightsVal[m] * LambdaList[[m]] + 2 * lambda2 * Reduce('+',wzList[-m]) * wList[[m]]
        ###### use (1-diag(1,p)) such that lambda 1 penalty has no effect on the diagonal.
        t.m             = penSpMcp.Deriv(x=penalty_sum,lambda=1,gamma=b)*
          penSpMcp.Deriv(x=zList[[m]],lambda=lambda1,gamma=a)*(1-diag(1,p))  
        zList[[m]]      = round((sign(b.m)*(abs(b.m)-t.m)*(abs(b.m)>t.m)/a.m),digits=10)
        LambdaList[[m]] = round((LambdaList[[m]] + u * (precisionList[[m]] - zList[[m]])),digits=10)
      } 
      
    }
    
    # ---------------------- record and evaluation ---------------------  
    
    time_cost[i_lam] = (proc.time() - pt)[3]
    iters[i_lam]     = iter
    est[[i_lam]]     = zList
    
    NP_tmp =BIC_tmp = Lik_tmp = rep(0,M)
    for (mm in 1:M){
      NP_tmp[mm] = sum(abs(zList[[mm]])>0) - p
      Lik_tmp[mm] =   nM[mm]*(-log(det(precisionList[[mm]])) + sum(covarianceList[[mm]]*precisionList[[mm]]))
      BIC_tmp[mm] = Lik_tmp[mm] +log(nM[mm])*(NP_tmp[mm]/2)
    }
    NP[i_lam,]  = NP_tmp
    BIC[i_lam,] = BIC_tmp
    Lik[i_lam,] = Lik_tmp
    
    i_lam = i_lam -1
  }
  
  #----------- * end: iterations * -------------------------#
  
  return(list(
    est= est,
    iters = iters,
    time=time_cost,
    NP = NP, 
    Lik = Lik,
    BIC = BIC))
}
