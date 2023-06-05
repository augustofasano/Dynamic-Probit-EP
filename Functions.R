# Functions
#-----------------------------------------------------------------------------#
# Get parameters for multivariate smoothing
getSmooothParams = function(y, mu0 = 0, Sigma0 = as.matrix(5), SigmaEps = as.matrix(0.1), matFF, GG = as.matrix(1)) {
  p = length(mu0)
  n = length(y)
  
  # compute the parameters
  Om = matrix(nrow = p*n, ncol = p*n)
  xi = matrix(nrow = p*n, ncol = 1)
  
  # create the variances
  varOld = Sigma0
  muOld = mu0
  for(tt in 1:n) {
    varNew = GG%*%varOld%*%t(GG)+SigmaEps
    Om[(p*(tt-1)+1):(p*tt),(p*(tt-1)+1):(p*tt)] = varNew
    varOld = varNew
    
    muNew = GG%*%muOld
    xi[(p*(tt-1)+1):(p*tt),1] = GG%*%muNew
    muOld = muNew
  }
  
  if(n > 1) {
    # create the covariances
    for(tt in 1:(n-1)) {
      covOld = Om[(p*(tt-1)+1):(p*tt),(p*(tt-1)+1):(p*tt)]
      for(ss in (tt+1):n) {
        covNew = covOld%*%t(GG)
        Om[(p*(tt-1)+1):(p*tt),(p*(ss-1)+1):(p*ss)] = covNew
        Om[(p*(ss-1)+1):(p*ss),(p*(tt-1)+1):(p*tt)] = t(covNew)
        covOld = covNew
      }
    }
  }
  
  om = sqrt(diag(Om[cbind(1:(p*n),1:(p*n))],p*n,p*n))
  bar_Om <- solve(om)%*%Om%*%solve(om)
  
  listFF = list()
  for(i in 1:n) {listFF[[i]] = matrix(matFF[i,],nrow = 1)}
  bigFF = bdiag(listFF)
  
  D = Diagonal(n=n, x=2*y-1)%*%bigFF
  s = Diagonal(n=n,x=sqrt((D%*%Om%*%t(D)+diag(1,n,n))[cbind(1:n,1:n)]))
  sSolve = Diagonal(n=n,x=1/diag(s))
  
  Delta = bar_Om%*%om%*%t(D)%*%sSolve
  gamma = sSolve%*%D%*%xi
  Gamma = sSolve%*%(D%*%Om%*%t(D)+diag(1,n,n))%*%sSolve
  
  Output <- list(xi=xi,omega=om,Omega=Om,Delta=as.matrix(Delta),gamma=as.matrix(gamma),Gamma=as.matrix(Gamma),GammaUnnormalized=as.matrix(D%*%Om%*%t(D)+diag(1,n,n)))
  return(Output)
}

# Get parameters for PFM approximation
getParamsPFM = function(X,y,Omega,moments = TRUE,tolerance = 1e-2, maxIter = 1e4) {
  ######################################################
  # PRECOMPUTATION
  ######################################################
  # define model dimensions
  n = dim(X)[1]
  p = dim(X)[2]
  
  # modify this part with Woodbury
  invOmega = solve(Omega)
  V = solve(t(X)%*%X+invOmega)
  H = X%*%V%*%t(X)
  invOmZ = solve(diag(1,nrow=n,ncol=n)+X%*%Omega%*%t(X)) # needed for ELBO
  
  # compute optimal sigma2
  h = diag(diag(H))
  sigma2 = matrix(1/(1-diag(H)), ncol = 1)
  sigma = sqrt(sigma2)
  
  # compute matrix to write the CAVI update in a vectorized form
  A = diag(as.double(sigma2), nrow = n, ncol = n)%*%(H - h)
  
  # alternative (better) writing
  # A = diag(as.double(sigma2), nrow = n, ncol = n)%*%H
  # A[cbind(1:n,1:n)] = 0
  
  # other useful quantities needed for ELBO
  diagInvOmZ = diag(invOmZ)
  coeffMean_Z2 = diagInvOmZ-1/sigma2
  
  # initialization of variables
  mean_Z2 = matrix(0,n,1)
  mu = matrix(0,n,1)
  
  musiRatio = as.double(mu/sigma)
  phiPhiRatio = exp(dnorm(musiRatio,log = T)-pnorm((2*y-1)*musiRatio,log.p = T))
  meanZ = mu + (2*y-1)*sigma*phiPhiRatio
  
  elbo = -Inf
  diff = 1
  nIter=0
  
  ######################################################
  # CAVI ALGORITHM
  ######################################################
  
  while(diff > tolerance & nIter < maxIter) {
    elboOld = elbo
    sumLogPhi = 0
    
    for(i in 1:n) {
      mu[i] = A[i,]%*%meanZ
      
      # compute first (needed for algorithm) and second (needed for ELBO) moments
      musiRatio = mu[i]/sigma[i]
      phiPhiRatio = exp(dnorm(musiRatio, log = T) - pnorm((2*y[i]-1)*musiRatio, log.p = T))
      meanZ[i] = mu[i] + (2*y[i]-1)*sigma[i]*phiPhiRatio
      mean_Z2[i] = mu[i]^2+sigma2[i]+(2*y[i]-1)*mu[i]*sigma[i]*phiPhiRatio # needed for ELBO
      sumLogPhi = sumLogPhi + pnorm((2*y[i]-1)*musiRatio, log.p = T)
    }
    
    # computation of ELBO (up to an additive constant not depending on mu)
    elbo = -(t(meanZ)%*%invOmZ%*%meanZ -
               sum((meanZ^2)*diagInvOmZ) +
               sum(mean_Z2*coeffMean_Z2))/2 -
      sum(meanZ*mu/sigma2) + sum((mu^2)/sigma2)/2 + sumLogPhi
    
    diff = abs(elbo-elboOld)
    nIter = nIter+1
    
    if(nIter%%100==0) {print(paste0("iter: ", nIter, ", ELBO: ", elbo))}
  }
  
  # get the optimal parameters of the normals before truncation, now that convergence has been reached
  mu = A%*%meanZ
  
  results = list(mu = mu, sigma2 = sigma2, nIter = nIter)
  
  ######################################################
  # (OPTIONAL) CLOSED-FORM MOMENTS' COMPUTATION
  ######################################################
  
  if(moments == TRUE) {
    # TO BE MODIFIED: compute V and V%*%t(X), directly or with Woodbury
    diagV = diag(V) # V already computed
    VXt = V%*%t(X)
    
    musiRatio = mu/sigma
    phiPhiRatio = exp(dnorm(musiRatio, log = T) - pnorm((2*y-1)*musiRatio, log.p = T))
    
    meanZ = mu + (2*y-1)*sigma*phiPhiRatio
    postVarZ = as.double(sigma2*(1-(2*y-1)*musiRatio*phiPhiRatio - phiPhiRatio^2))
    
    W = apply(VXt,1,function(x) sum(x*x*postVarZ))
    
    meanBeta = VXt%*%meanZ
    varBeta = diagV + W
    
    moments_PFM = list(meanBeta=meanBeta,varBeta=matrix(varBeta,ncol = 1))
    
    results = c(results,postMoments=moments_PFM)
  }
  
  return(results)
}

# Get parameters for MF approximation
getParamsMF = function(X,y,Omega, tolerance = 1e-2, maxIter = 1e4){
  ######################################################
  # PRECOMPUTATION
  ######################################################
  # define model dimensions
  n = dim(X)[1]
  p = dim(X)[2]
  
  # compute V and other useful quantities
  invOmega = solve(Omega)
  V = solve(t(X)%*%X+invOmega)
  H = X%*%V%*%t(X)
  diagV = diag(V)
  VXt = V%*%t(X)
  
  # other useful quantites
  # XVVXt = t(VXt)%*%VXt # needed for ELBO
  signH = H
  signH[y==0,] = -H[y==0,]
  
  # initialization of variables
  meanZ = matrix(0,n,1)
  diff = 1
  elbo = -Inf
  nIter=0
  
  ######################################################
  # CAVI ALGORITHM
  ######################################################
  
  while(diff > tolerance & nIter < maxIter) {
    elboOld = elbo
    
    # update parameters
    mu = H%*%meanZ
    meanZ = mu +(2*y-1)*exp(dnorm(mu, log = T) - pnorm((2*y-1)*mu, log.p = T))
    
    # compute ELBO
    # elbo = (t(meanZ)%*%XVVXt%*%meanZ)/nu2 + sum(pnorm(signH%*%meanZ, log.p = T))
    VXTt_meanZ = VXt%*%meanZ
    elbo = -0.5*t(VXTt_meanZ)%*%invOmega%*%VXTt_meanZ + sum(pnorm(signH%*%meanZ, log.p = T))
    
    # compute change in ELBO
    diff = abs(elbo-elboOld)
    
    nIter = nIter+1
    if(nIter %% 100 == 0) {print(paste0("iter: ", nIter, ", ELBO: ", elbo))}
  }
  meanBeta = VXt%*%meanZ
  
  return(list(meanBeta = meanBeta, diagV = diagV, nIter = nIter))
}

# ------------------------------- EFFICIENT EP ---------------------------------
zeta1 = function(x){exp(dnorm(x, log = T) - pnorm(x,log.p = T))}
zeta2 = function(x,z1){-z1^2-x*z1}

getParamsEP = function(X,y,Omega,tolerance=1e-3,maxIter=1e3,fullVar=FALSE,nPrint=100){
  
  n = dim(X)[1]
  pn = dim(X)[2]
  
  threshold = 50 # threshold to decide whether to use algo 1 or 2 in the paper
  algo=1
  if(pn>threshold){
    algo=2
  }
  
  ### Initialization
  if(algo==1){
    invQ = Omega
  }else{
    V = tcrossprod(Omega,X) # this could be made more efficient as X is sparse
  }
  
  r = rep(0,pn)
  
  k = double(length = n)
  m = double(length = n)
  
  diff = 1
  nIter = 0
  
  ### Iterations
  
  if(algo==1){
    
    while(diff > tolerance && nIter < maxIter){
      
      diff = 0.
      count = 0
      
      for(i in 1:n){
        
        xi = X[i,]
        
        r_i = r - m[i]*xi
        
        Oxi = invQ%*%xi
        
        Oi = invQ + tcrossprod(Oxi)*k[i]/as.double(1.-k[i]*crossprod(xi,Oxi))
        
        Oixi = Oi%*%xi
        xiOixi = as.double(crossprod(xi,Oixi))
        
        if(xiOixi>0){
          
          r_iOixi = as.double(crossprod(r_i,Oixi))
          
          s = (2.*y[i]-1.)/sqrt(1.+xiOixi)
          tau = s*r_iOixi
          
          z1 = zeta1(tau)
          z2 = zeta2(tau,z1)
          
          kNew = - z2/(1.+xiOixi+z2*xiOixi)
          mNew = s*z1 + kNew*r_iOixi + kNew*s*z1*xiOixi
          
          maxDiff = max(abs(c(kNew - k[i], mNew - m[i])))
          if(maxDiff>diff){diff = maxDiff}
          
          k[i] = kNew
          m[i] = mNew
          
          r = r_i + m[i]*xi

          invQ = Oi - tcrossprod(Oixi)*k[i]/(1.+k[i]*xiOixi)
        }else{
          count = count+1
          print(paste0(count," units skipped"))
        }
        
      }
      
      nIter = nIter + 1
      if(nIter %% nPrint == 0) {print(paste0("iteration ",nIter))}
    }
    
  } else {
    
    while(diff > tolerance && nIter < maxIter){
      
      diff = 0.
      count = 0
      
      for(i in 1:n){
        
        v = V[,i]
        xi = X[i,]
        xTv = crossprod(xi,v)[1]
        
        d = 1-k[i]*xTv
        w = v/d
        xTw = xTv/d
        
        if(xTw>0){
          
          r_iTw = crossprod(r,w) - m[i]*xTw
          
          s = (2*y[i]-1)/sqrt(1+xTw)
          tau = s*r_iTw
          
          z1 = zeta1(tau)
          z2 = zeta2(tau,z1)
          
          kNew = as.double(-z2/(1 + xTw + z2*xTw))
          mNew = as.double(z1*s + kNew*r_iTw + kNew*z1*s*xTw)
          
          r = r + (mNew - m[i])*xi
          
          ratio = (k[i]-kNew)/(1.+(kNew-k[i])*xTv)
          V = V + (v*ratio)%*%crossprod(xi,V)
          
          maxDiff = max(abs(c(kNew - k[i], mNew - m[i])))
          if(maxDiff>diff){diff = maxDiff}
          
          k[i] = kNew
          m[i] = mNew
          
        }else{
          count = count+1
          print(paste0(count," units skipped"))
        }
      }
      
      nIter = nIter + 1
      if(nIter %% nPrint == 0) {print(paste0("iteration ",nIter))}
    }
    
  }
  
  ### Posterior Approximate Moments
  
  if(algo==1){
    
    meanBeta = invQ%*%r
    diagOmegaEP = diag(invQ)
    
  }else{
    XOm = X%*%Omega # this could be made more efficient as X is sparse
    diagOmegaEP = diag(Omega) - rowSums(V*t(k*XOm))
    Om_r = Omega%*%r
    meanBeta = Om_r - V%*%(k*X%*%Om_r) # this could be made more efficient as X is sparse
  }
  
  results = list(meanBeta = meanBeta, diagOmegaEP = diagOmegaEP, 
                 nIter = nIter, kEP = k, mEP = m)
  
  if(fullVar==TRUE){
    if(algo==2) {
      invQ = Omega - V%*%(k*XOm)
    }
    results = append(list(OmegaEP=invQ),results)
  }
  
  return(results)
}


# Multivariate smoothing for generic dynamic probit
smoothing = function(y, mu0 = 0, Sigma0 = as.matrix(5), SigmaEps = as.matrix(0.1), matFF, GG = as.matrix(1), seed = 123, nSim = 1e4) {
  params=getSmooothParams(y = y,
                          mu0 = mu0,
                          Sigma0 = Sigma0,
                          SigmaEps = SigmaEps,
                          matFF = matFF,
                          GG = GG)
  
  simS_prov = t(rSUN(xi = params$xi,
                     Omega = params$Omega,
                     Delta = params$Delta,
                     gamma = params$gamma,
                     Gamma = params$Gamma,
                     n = nSim,
                     seed = seed))
  
  p = length(mu0)
  n = length(y)
  simS = array(dim = c(nSim,n,p))
  
  for(i in 1:p){
    seqInd = seq.int(from = i, by = p, length.out = n)
    simS[,,i] = simS_prov[,seqInd]
  }
  
  Output <- list(smoothState=list(values=simS))
  return(Output)
}

# Generic function to sample from a SUN with given parameters (one sample per column in output, i.e. p x nSamples output)
rSUN = function(xi,Omega,Delta,gamma,Gamma,n=1, seed = 123) {
  set.seed(seed)
  d = length(xi)
  m = length(gamma)
  
  if(d==1) {
    omega=sqrt(Omega)
  } else {
    omega=diag(sqrt(diag(Omega)))
  }
  
  #OmegaBar = cov2cor(Omega)
  omDe = omega%*%Delta
  invGamma = solve(Gamma)
  
  PsiDel = Omega-omDe%*%invGamma%*%t(omDe)
  PsiDel = (PsiDel+t(PsiDel))/2
  
  sampleMultNorm=t(rmvnorm(n = n, mean = rep(0,d), sigma = PsiDel))
  
  if(m==1){ # when m == 1 rtmvnorm returns the samples as a row vector (i.e. one sample per column)
    sampleTruncNorm = rtmvnorm(n = n, mu = rep(0,m), sigma = Gamma, lb = -gamma, ub = rep(Inf,m))
  } else{ # when m > 1 rtmvnorm returns the samples in an n x m matrix (i.e. one sample per row) 
    sampleTruncNorm = t(rtmvnorm(n = n, mu = rep(0,m), sigma = Gamma, lb = -gamma, ub = rep(Inf,m)))
  }
  
  rep_xi = matrix(nrow = d, ncol = n)
  rep_xi[,1:n] = xi
  
  sampleSUN=rep_xi+sampleMultNorm+omDe%*%invGamma%*%sampleTruncNorm
}
