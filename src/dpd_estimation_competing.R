
library(optimx)
library(MASS)

survival.units <- function(N,n){
  
  # Compute survival units on each step of the experiment
  
  survival = c(N)
  for(step in 2:length(n)){survival[step] = survival[step-1]-n[step-1]}
  return(survival)
}

shifted.times <- function(tauvec, lambda){
  
  #compute vector of shited times for the CE model
  
  k = length(tauvec)

  #calculate h_i vector
  h = c(0, (-1+lambda[2]/lambda[1])*tauvec[1])
  return(h)
}

shifted.times.ast <- function(stressvec, tauvec, lambda){
  
  k = length(stressvec)
  
  #calculate h_i vector
  h = c(0, tauvec[1]*(lambda[2]/lambda[1])*(stressvec[2]-stressvec[1]))
  return(h)
  
}

relative.risks <- function(lambdas, C){
  
  K1 = 0
  K2 = 0
  for (c in c(1:C)){
    K1 = K1+ 1/lambdas[(2*c-1)]
    K2 = K2 + 1/lambdas[(2*c)]
  }
  
  rr = vector()
  for (c in c(1:C)){
    rr[(2*c-1)] = (1/lambdas[(2*c-1)])/K1
    rr[(2*c)]= (1/lambdas[(2*c)])/K2
  }
  
  return(rr)
}

model.distribution.exponential <-function(lambda, tauvec, stressvec, t){

  #compute lifetime for the step-stress model under exponential distributions

  k = length(tauvec) #number of stress level/number of time of stress change

  hi1 = shifted.times(tauvec, lambda)

  #create position vector of the times of stress change
  stress.position = (t >= tauvec[1])+1
  ttras = (t+hi1[stress.position])

  if (t <= tauvec[length(tauvec)]){
    DFsol = 1-exp(-ttras/lambda[stress.position])
  }else{
    stop("lifetime out of experimental limits")
  }

  return(DFsol)
}

model.distribution.competing.exponential <-function(lambdas, tauvec, stressvec, t, C=2){
  
  # compute lifetime for the step-stress model under exponential distributions
  
  # verify the number of competing risks
  if (length(lambdas)/length(stressvec) != C){
    stop("number of competing risk does not match with number of lambdas")
  }
  
  hs = list()
  for (j in c(1:C)){
    hs[[j]] = shifted.times(tauvec, lambdas[(2*j-1):(2*j)])
  }
  
  #create position vector of the times of stress change
  stress.position = (t >= tauvec[1])+1
  DFsol = 0
  if (t <= tauvec[length(tauvec)]){
    for (j in c(1:C)){
      ttras = (t+hs[[j]][stress.position])
      lambda = lambdas[(2*j-1):(2*j)][stress.position]
      DFsol = DFsol + ttras/lambda
    }
  }else{
    stop("lifetime out of experimental limits")
  }
  
  return(1-exp(-DFsol))
}

model.probability.joint.competing.exponential <-function(lambdas, tauvec, stressvec, t, j, C=2){
    
    rs = relative.risks(lambdas, C)
    survival = 1 - model.distribution.competing.exponential(lambdas, tauvec, stressvec, t, C)
    
    stress.position = (t > tauvec[1])+1
    r = rs[(2*j-1):(2*j)][stress.position]
    return(r*survival)
}

theoretical.probability.marginal <- function(lambdas, tauvec, stressvec, ITvec, C=2){
  
  th= c()
  th[1] = model.distribution.competing.exponential(lambdas, tauvec, stressvec, ITvec[1])

  
  #Probability of failure
  for(step in 2:length(ITvec)){
      th[step] = -model.distribution.competing.exponential(lambdas, tauvec, stressvec, ITvec[step-1])+model.distribution.competing.exponential(lambdas, tauvec, stressvec, ITvec[step])
  }
  
  #survival probability
  th[step+1] = 1 - model.distribution.competing.exponential(lambdas, tauvec, stressvec, ITvec[step])
  
  return(th)
}

theoretical.probability <- function(lambdas, tauvec, stressvec, ITvec, C=2){
  
  th= c()
  for(j in c(1:C)){
    th[j] = model.probability.joint.competing.exponential(lambdas, tauvec, stressvec, 0, j) - model.probability.joint.competing.exponential(lambdas, tauvec, stressvec, ITvec[1], j)
  }
  
  #Ahora las diferencias
  for(step in 2:length(ITvec)){
    for(j in c(1:C)){
      #s = stress.position = sum(ITvec[step] > tauvec)+1
      s= (step-1)*C+j
      th[s] = model.probability.joint.competing.exponential(lambdas, tauvec, stressvec, ITvec[step-1], j)-model.probability.joint.competing.exponential(lambdas, tauvec, stressvec, ITvec[step], j)
    }
  }
  
  th[s+1] = 1 - model.distribution.competing.exponential(lambdas, tauvec, stressvec, ITvec[step])
  
  #alliavte nans problem
  th[th<1e-8]=1e-5
  th[th>(1-1e-8)]= (1-1e-5)
  
  th = th/sum(th)
  
  return(th)
}

simulate.sample <- function(theta, cont, N, tauvec, stressvec, ITvec,  C=2, seed){
  

  lambdas = c()
  for (j in c(1:C)){
    lambdas[(2*j-1):(2*j)] = exp(theta[(2*j-1)]+ theta[2*j]*stressvec)
   }
  
  set.seed(seed)
  
  th = theoretical.probability(lambdas, tauvec, stressvec, ITvec)
  th.cont = rep(0, length(th)) #theoretical.probability(lambdas.cont, tauvec, stressvec, ITvec)
  th.cont[c(3,4,5,6)] = 1/4
  
  ncont = floor(N*cont) 
  th = th/sum(th)
  n1 = drop(rmultinom(1, (N-ncont), prob = th)) 
  n2 = drop(rmultinom(1, ncont, prob = th.cont)) 
  n= n1+n2
  
  return(n)
}

loglikelihood <- function(theta, N, tauvec, stressvec, ITvec, p, C=2){
  
  lambdas = c()
  for (j in c(1:C)){
    lambdas[(2*j-1):(2*j)] = exp(theta[(2*j-1)]+ theta[2*j]*stressvec)
  }
  
  th = theoretical.probability(lambdas, tauvec, stressvec, ITvec, C)
  
  like = 0
  for(step in 1:length(p)){ like = like + p[step]*log(th[step]) }
  
  return(-like)
}

dbeta <- function(hat.p, th.p, beta){
  return(th.p^(beta+1)-((beta+1)/beta)*(hat.p*(th.p^beta)))
}

DPDloss <- function(theta, N, tauvec, stressvec, ITvec, n, beta, C=2){
  
  lambdas = c()
  for (j in c(1:C)){
    lambdas[(2*j-1):(2*j)] = exp(theta[(2*j-1)]+ theta[2*j]*stressvec)
  }
  
  th = theoretical.probability(lambdas, tauvec, stressvec, ITvec, C)
  
  p = n/N
  p[p==0]<-0.01
  p[p==1]<-0.99
  
  if (beta == 0){
    DPD=loglikelihood(theta, N, tauvec, stressvec, ITvec, p, C)
  }else{
    DPD = 0
    for(step in 1:length(n)){ DPD = DPD + dbeta(p[step],th[step], beta)}
  }
  return(DPD)
}


estimate.function <-function(N, tauvec, stressvec, ITvec, n, C=2, initial = c(7, -0.01, 5, -0.02), beta.list = c(0.2,0.4,0.6,0.8,1)){
  
  estimates = list()
  
  #MLE
  MLE =  tryCatch(optimr(par= initial, DPDloss,  N = N, tauvec = tauvec, 
                         stressvec = stressvec, ITvec = ITvec, n = n, beta=0, C=C,
                         method = 'Nelder-Mead'),
                  error=function(sol){sol$code=3;return(NA)})
  if(!is.na( sum(MLE$par)) ){estimates[["MLE"]] = as.vector(MLE$par)}else{estimates[["MLE"]] = NA}
  MLE
  
  for(beta in beta.list){
    DPDE =  tryCatch(optimr(par= initial,DPDloss,  N = N, tauvec = tauvec, 
                            stressvec = stressvec, ITvec = ITvec, n = n, beta= beta, C=C,
                            method = "Nelder-Mead"), 
                     error=function(sol){sol$code=3;return(NA)})
    if(!is.na( sum(DPDE$par)) ){estimates[[paste0("DPD", beta)]] = as.vector(DPDE$par)}else{estimates[[paste0("DPD", beta)]] = NA}
  }
  
  
  return(estimates)
}

estimate.beta <-function(N, tauvec, stressvec, ITvec, n, C=2, initial = c(7, -0.01, 5, -0.02), beta=0){
  
    DPDE =  tryCatch(optimr(par= initial,DPDloss,  N = N, tauvec = tauvec, stressvec = stressvec, ITvec = ITvec, n = n, beta= beta, C=C), 
                     error=function(sol){sol$code=3;return(NA)})

  return(DPDE$par)
}



RMSE.function <- function(theta, theta.hat){return(mean(((theta-theta.hat)/theta)^2) )}

validility.check <- function(n, n1_){
  #only for two components
  valid1 = (sum(n[seq(1,2*n1_,2)])  != 0)
  valid2 = (sum(n[seq(2,2*n1_,2)]) != 0)
  
  valid3 = (sum(n[seq(2*n1_+1,(length(n)-1), by=2)])  != 0)
  valid4 = (sum(n[seq(2*n1_+2,(length(n)-1), by=2)]) != 0)
  
  return (valid1 & valid2 & valid3 & valid4)
}

simulate <-function(theta, cont.vec,  N, tauvec, stressvec, ITvec, C=2, B=50, beta.list = c(0.2,0.4,0.6,0.8,1)){

  nbetas = length(beta.list)+1
  res.mat = matrix(0, nrow =  length(cont.vec), ncol = nbetas)
  
  t1_ = sum(ITvec <= tauvec[1])
  
  for(cont in 1:length(cont.vec)){
    
    count = 0
    
    for(b in 1:B){
       
      valid = FALSE
      while (!valid){
        n = simulate.sample(theta, cont.vec[cont], N, tauvec, stressvec, ITvec,  seed = 10*b, C=C)
        if( validility.check(n, t1_)){valid = TRUE}else{b=7*b}
      }
      
      estimators = estimate.function(N, tauvec, stressvec, ITvec, n, initial = c(7, -0.01, 4, -0.01), C=C, beta.list = beta.list)
      
      if(sum(sapply(estimators, function(x) sum(is.na(x))))<1){ #none of the estimators is nan
        res.mat[cont,] = res.mat[cont,]+ unlist(lapply(estimators, RMSE.function, theta = theta))
        count = count+1
      }
      
    }
    res.mat[cont,] = res.mat[cont,]/count
    
  }
  return(res.mat)
}

simulate.estimate <-function(theta, cont.vec,  N, tauvec, stressvec, ITvec, C=2, B=50, beta.list = c(0.2,0.4,0.6,0.8,1)){
  
  nbetas = length(beta.list)+1
  #res.mat list with estimates (four parameters)
  res.mat = list(matrix(0, nrow =  length(cont.vec), ncol = nbetas),
                 matrix(0, nrow =  length(cont.vec), ncol = nbetas),
                 matrix(0, nrow =  length(cont.vec), ncol = nbetas),
                 matrix(0, nrow =  length(cont.vec), ncol = nbetas) )
  
  t1_ = sum(ITvec <= tauvec[1])
  
  for(cont in 1:length(cont.vec)){
    
    count = 0
    
    for(b in 1:B){
      
      valid = FALSE
      while (!valid){
        n = simulate.sample(theta, cont.vec[cont], N, tauvec, stressvec, ITvec,  seed = 10*b, C=C)
        if( validility.check(n, t1_)){valid = TRUE}else{b=7*b}
      }
      
     
      estimators = estimate.function(N, tauvec, stressvec, ITvec, n, initial = c(7, -0.01, 4, -0.01), C=C, beta.list = beta.list)
      
      if(sum(sapply(estimators, function(x) sum(is.na(x))))<1){ #none of the estimators is nan
        for (be in 1:nbetas){
          for (e in 1:4){res.mat[[e]][cont,be] = res.mat[[e]][cont,be]+ estimators[[be]][[e]]}
        }
        count = count+1
      }
    }
    for (e in 1:4){res.mat[[e]][cont,] = res.mat[[e]][cont,]/count}
  }
  
  return(res.mat)
}


######## ------------- MAIN MATRICES--------- #############
Wmatrix.entryjk.1 <- function(lambdas, tauvec, stressvec, C, ITvec, step, j, k){

  rs = relative.risks(lambdas, C)
  hs = list()
  for (c in c(1:C)){ hs[[c]] = shifted.times(tauvec, lambdas[(2*c-1):(2*c)]) }
    
  survival = 1 - model.distribution.competing.exponential(lambdas, tauvec, stressvec, ITvec[step], C)
  
  stress.position = (ITvec[step] >= tauvec[1])+1
  rj = rs[(2*j-1):(2*j)][stress.position] #pi_ij
  rk = rs[(2*k-1):(2*k)][stress.position] #pi_ik
  
  lambdak = lambdas[(2*k-1):(2*k)][stress.position]
  ttrask = (ITvec[step]+hs[[k]][stress.position])
  
  if (j==k){
    entry = rk*survival*(rk-1+ ttrask/lambdak)
  }else{
    entry = rj*survival*(rk + ttrask/lambdak)
  }
  return(entry)
}

Wmatrix.entryjk.2 <- function(lambdas, tauvec, stressvec, C, ITvec, step, j, k){
  
  rs = relative.risks(lambdas, C)
  hs = list()
  hs_ast = list()
  for (c in c(1:C)){ hs[[c]] = shifted.times(tauvec, lambdas[(2*c-1):(2*c)]) }
  for (c in c(1:C)){ hs_ast[[c]] = shifted.times.ast(tauvec, stressvec, lambdas[(2*c-1):(2*c)]) }
  
  survival = 1 - model.distribution.competing.exponential(lambdas, tauvec, stressvec, ITvec[step], C)
  
  stress.position = (ITvec[step] >= tauvec[1])+1
  rj = rs[(2*j-1):(2*j)][stress.position] #pi_ij
  rk = rs[(2*k-1):(2*k)][stress.position] #pi_ik
  
  lambdak = lambdas[(2*k-1):(2*k)][stress.position]
  ttrask = (ITvec[step]+hs[[k]][stress.position])
  
  if (j==k){
    entry = rj*survival*(rk*stressvec[stress.position]-stressvec[stress.position] + (-hs_ast[[k]][stress.position]+ttrask*stressvec[stress.position])/lambdak)
  }else{
    entry = rj*survival*(rk*stressvec[stress.position] + (-hs_ast[[k]][stress.position]+ttrask*stressvec[stress.position])/lambdak)
  }
  return(entry)
}

Wmatrix.entr01 <- function(lambdas, tauvec, stressvec, C, k){
  
  survival = 1- model.distribution.competing.exponential(lambdas, tauvec, stressvec, tauvec[-1], C)
  
  hs = list()
  for (c in c(1:C)){ hs[[c]] = shifted.times(tauvec, lambdas[(2*c-1):(2*c)]) }
  
  stress.position = length(stressvec)
  
  ttras = tauvec[-1]+hs[[k]][stress.position]
  lambda = lambdas[(2*k-1):(2*k)][stress.position]
  
  return(survival*ttras/lambda)
  
}

Wmatrix.entr02 <- function(lambdas, tauvec, stressvec, C, k){
  
  survival = 1- model.distribution.competing.exponential(lambdas, tauvec, stressvec, tauvec[-1], C)
  
  hs = list()
  hs_ast = list()
  for (c in c(1:C)){ hs[[c]] = shifted.times(tauvec, lambdas[(2*c-1):(2*c)]) }
  for (c in c(1:C)){ hs_ast[[c]] = shifted.times.ast(tauvec, stressvec, lambdas[(2*c-1):(2*c)]) }
  
  stress.position = length(stressvec)
  
  ttras = -hs_ast[[k]][stress.position] +(tauvec[-1]+hs[[k]][stress.position])*stressvec[stress.position]
  lambda = lambdas[(2*k-1):(2*k)][stress.position]
  
  return(survival*ttras/lambda)
  
}

Wmatrix <-function(lambdas, tauvec, stressvec, ITvec, C=2){
  
  # The matrix W has dimensions (RL+1)x2R and contains L row blocks of j x 2R
  # Each row block has the (ordered) derivatives of the probabilities plj
  # with respect to each component 

  W = matrix(NA, nrow=length(ITvec)*C+1, ncol = 2*C)
  
  
  for (step in 1:(length(ITvec))){
    for (j in 1:C){
      for (k in 1:C){
        rw = (step-1)*C+j
        W[rw, (2*k-1)] = Wmatrix.entryjk.1(lambdas, tauvec, stressvec, C, c(0,ITvec), step, j, k) - Wmatrix.entryjk.1(lambdas, tauvec, stressvec, C, ITvec, step, j, k)
        W[rw, (2*k)] = Wmatrix.entryjk.2(lambdas, tauvec, stressvec, C, c(0,ITvec), step, j, k) - Wmatrix.entryjk.2(lambdas, tauvec, stressvec, C, ITvec, step, j, k)
      }
    }
  }
  for (k in 1:C){
    W[rw+1, (2*k-1)] = Wmatrix.entr01(lambdas, tauvec, stressvec, C, k)
    W[rw+1, (2*k)] = Wmatrix.entr02(lambdas, tauvec, stressvec, C, k)
  }

  return(W)
  
}

Dmatrix  <-function(lambdas, tauvec, stressvec, ITvec, C, beta){
  th_beta = theoretical.probability(lambdas, tauvec, stressvec, ITvec, C)^beta
  return (diag(th_beta))
}

MDPDE_ec <- function(estimate, tauvec, stressvec, ITvec, C, beta, n){
  
  N = sum(n)
  lambdas = c()
  for (j in c(1:C)){lambdas[(2*j-1):(2*j)] = exp(estimate[(2*j-1)]+ estimate[2*j]*stressvec) }
  
  W = Wmatrix(lambdas, tauvec, stressvec, ITvec, C)
  D= Dmatrix(lambdas, tauvec, stressvec, ITvec, C, beta-1)
  p_est = n/N
  th = theoretical.probability(lambdas, tauvec, stressvec, ITvec, C)

  
  return(t(W)%*%D%*%(p_est-th))
  
}

Sigmamatrix<- function(lambdas, tauvec, stressvec, ITvec, C, beta){
  W = Wmatrix(lambdas, tauvec, stressvec, ITvec, C)
  D = Dmatrix(lambdas, tauvec, stressvec, ITvec, C, beta-1)
  J = t(W)%*%D%*%W
  th = theoretical.probability(lambdas, tauvec, stressvec, ITvec, C)^beta
  
  K=  t(W)%*%(Dmatrix(lambdas, tauvec, stressvec, ITvec, C, beta-2)-th%*%t(th))%*%W
  
  Sigma = solve(J)%*%K%*%solve(J)
  return(Sigma)
}


norm <-function(v){return(sqrt(sum(v^2)))}


######## ------------- MEAN LIFETIME ESTIMATION -------- #############

mean.lifetime <- function(theta, C, x0, t0=FALSE, alpha0=FALSE){
  
  K0 = 0
  for(j in 1:C){ K0 = K0 + exp(-theta[(2*j-1)]- theta[(2*j)]*x0) }
  return(1/K0)
  
}

mean.lifetime.grad <- function(theta, C, x0){
  
  K0 = 0
  for(j in 1:C){ K0 = K0 + exp(-theta[(2*j-1)]- theta[(2*j)]*x0) }
  
 # pi0 = c()
  grad = c()
  for(j in 1:C){
    pi0 = exp(-theta[(2*j-1)]- theta[(2*j)]*x0)/K0
    grad = c(grad, c(pi0/K0, pi0*x0/K0))
  }
  
  return(grad)
}

mean.lifetime.sd <- function(theta, tauvec, stressvec, ITvec, C, N, x0, beta){
 
  lambdas = c()
  for (c in c(1:C)){lambdas[(2*c-1):(2*c)] = exp(theta[(2*c-1)]+ theta[2*c]*stressvec) }
  
  grad= mean.lifetime.grad(theta, C,x0)
  
  Sigma = Sigmamatrix(lambdas, tauvec, stressvec, ITvec, C, beta)/N
  sd = drop(sqrt(grad%*%Sigma%*%grad))
  
  return(sd)
}

mean.lifetime.ci <- function(theta, tauvec, stressvec, ITvec, C, N, x0, t0, alpha0, beta, alpha=0.05){
  
  m = mean.lifetime(theta, C,x0)
  sd = mean.lifetime.sd(theta = -c(log(0.004),0.02, log(0.003),0.04), tauvec, stressvec, ITvec, C, N, x0, beta)
  z= qnorm(1-alpha/2, 0,1)
  
  return(c(m-z*sd, m+z*sd))
}

mean.lifetime.transformed.ci <- function(theta, tauvec, stressvec, ITvec, C, N, x0, t0, alpha0, beta, alpha=0.05){
  
  m = mean.lifetime(theta, C,x0)
  sd = mean.lifetime.sd(theta, tauvec, stressvec, ITvec, C, N, x0, beta)
  z= qnorm(1-alpha/2, 0,1)
  
  i = m*exp(-z*(sd/m))
  s = m*exp(z*(sd/m))
  
  return(c(i, s))
}


reliability <- function(theta, C,x0, t0, alpha0 = FALSE){
  
  K0 = 0
  for(j in 1:C){ K0 = K0 + exp(-theta[(2*j-1)]- theta[(2*j)]*x0) }
  return(exp(-t0*K0))
  
}

reliability.grad <- function(theta, C,x0, t0){
  
  rel= reliability(theta, C,x0, t0)
  
  grad = c()
  for(j in 1:C){
    pi0 = rel*t0*exp(-theta[(2*j-1)]- theta[(2*j)]*x0)
    grad = c(grad, c(pi0, pi0*x0))
  }
  
  return(grad)
}

reliability.sd <- function(theta, tauvec, stressvec, ITvec, C, N, x0, t0, beta){
  
  lambdas = c()
  for (c in c(1:C)){lambdas[(2*c-1):(2*c)] = exp(theta[(2*c-1)]+ theta[2*c]*stressvec) }
  
  grad= reliability.grad(theta, C, x0, t0)
  
  Sigma = Sigmamatrix(lambdas, tauvec, stressvec, ITvec, C, beta)/N
  sd = drop(sqrt(grad%*%Sigma%*%grad))
  
  return(sd)
}

reliability.ci <- function(theta, tauvec, stressvec, ITvec, C, N, x0, t0, alpha0, beta, alpha=0.05){
  
  r = reliability(theta, C,x0, t0)
  sd = reliability.sd(theta, tauvec, stressvec, ITvec, C, N, x0, t0,beta)
  z= qnorm(1-alpha/2, 0,1)
  
  return(c(r-z*sd, r+z*sd))
}

reliability.transformed.ci <- function(theta, tauvec, stressvec, ITvec, C, N, x0, t0, alpha0, beta, alpha=0.05){
  
  r = reliability(theta, C,x0, t0)
  sd = reliability.sd(theta, tauvec, stressvec, ITvec, C, N, x0, t0,beta)
  z= qnorm(1-alpha/2, 0,1)
  
  S = exp(z*sd/(r*(1-r)) )  
  
  i = r/(r+(1-r)*S)
  s = r/(r+(1-r)/S)
  return(c(i,s))
}


quantile <- function(theta, C,x0, t0=FALSE, alpha0){
  
  K0 = 0
  for(j in 1:C){ K0 = K0 + exp(-theta[(2*j-1)]- theta[(2*j)]*x0) }
  return(-log(1-alpha0)/K0)
  
}

quantile.grad <- function(theta, C,x0, alpha0){
  
  q = quantile(theta, C,x0, alpha0=alpha0)
  
  K0 = 0
  for(j in 1:C){ K0 = K0 + exp(-theta[(2*j-1)]- theta[(2*j)]*x0) }
  
  grad = c()
  for(j in 1:C){
    pi0 = -log(1-alpha0)*exp(-theta[(2*j-1)]- theta[(2*j)]*x0)/K0
    grad = c(grad, c(pi0/K0, pi0*x0/K0))
  }
  
  return(grad)
}

quantile.sd <- function(theta, tauvec, stressvec, ITvec, C, N, x0,  alpha0, beta){
  
  lambdas = c()
  for (c in c(1:C)){lambdas[(2*c-1):(2*c)] = exp(theta[(2*c-1)]+ theta[2*c]*stressvec) }
  
  grad= quantile.grad(theta, C, x0, alpha0)
  
  Sigma = Sigmamatrix(lambdas, tauvec, stressvec, ITvec, C, beta)/N
  sd = drop(sqrt(grad%*%Sigma%*%grad))
  
  return(sd)
}


quantile.ci <- function(theta, tauvec, stressvec, ITvec, C, N, x0, t0, alpha0, beta,  alpha=0.05){
  
  r = quantile(theta, C,x0, alpha0=alpha0)
  sd = quantile.sd(theta, tauvec, stressvec, ITvec, C, N, x0, alpha0, beta)
  z= qnorm(1-alpha/2, 0,1)
  
  return(c(r-z*sd, r+z*sd))
}

quantile.transformed.ci <- function(theta, tauvec, stressvec, ITvec, C, N, x0, t0, alpha0, beta,  alpha=0.05){
  
  r = quantile(theta, C,x0, alpha0=alpha0)
  sd = quantile.sd(theta, tauvec, stressvec, ITvec, C, N, x0, alpha0, beta)
  z= qnorm(1-alpha/2, 0,1)
  
  i = r*exp(-z*(sd/r))
  s = r*exp(z*(sd/r))
  return(c(i,s))
}


simulate.charasteristic <-function(theta,  cont.vec,  N, tauvec, stressvec, ITvec, x0=30, t0=50, alpha = 0.05,
                                   C=2, B=50, beta.list = c(0.2,0.4,0.6,0.8,1),
                                    charasteristic="mean.lifetime", alpha0=0.05){
    
    char.function = match.fun(charasteristic)
    
    f <- paste0(charasteristic, ".ci")
    ft <- paste0(charasteristic, ".transformed.ci")
    ci.function = match.fun(f)
    ci.function.transformed = match.fun(ft)
  
    nbetas = length(beta.list)+1
    res = list()
    rest = list()
    for(b in 1:nbetas){
      res[[b]] = matrix(0, nrow =  length(cont.vec), ncol = 6)
      rest[[b]] = matrix(0, nrow =  length(cont.vec), ncol = 6)
    }
    
    n1_ = sum(ITvec <= tauvec[1])
    m = char.function(theta=theta, C=C, x0=x0, t0=t0, alpha0=alpha0)
    
    for(cont in 1:length(cont.vec)){
      
      count = 0
      
      for(b in 1:B){
        
        valid = FALSE
        while (!valid){
          n = simulate.sample(theta, cont.vec[cont], N, tauvec, stressvec, ITvec,  seed = 10*b, C=C)
          if( validility.check(n, n1_)){valid = TRUE}else{b=7*b}
        }
        
        estimators = estimate.function(N, tauvec, stressvec, ITvec, n, initial = c(7, -0.01, 5, -0.02), C=C, beta.list = beta.list)
        
        if(sum(sapply(estimators, function(x) sum(is.na(x))))<1){ #none of the estimators is nan
          
          for (e in 1:nbetas){
            me = char.function(theta = estimators[[e]], C=C, x0=x0, t0=t0, alpha0=alpha0)
            ci = ci.function(estimators[[e]], tauvec, stressvec, ITvec, C, N, x0=x0, t0=t0, alpha0=alpha0, beta= c(0,beta.list)[e], alpha=alpha)
            ci.transformed =  ci.function.transformed(estimators[[e]], tauvec, stressvec, ITvec, C, N, x0=x0, t0=t0, alpha0=alpha0, beta= c(0,beta.list)[e], alpha=alpha)
            
            res[[e]][cont,1] = res[[e]][cont,1]+ (me - m)^2
            res[[e]][cont,2] = res[[e]][cont,2]+ me
            res[[e]][cont,3:4] = res[[e]][cont,3:4] + ci
            res[[e]][cont,5] = res[[e]][cont,5]+ ((ci[1] < m) & (ci[2] > m))
            res[[e]][cont,6] = res[[e]][cont,6]+ (ci[2] -ci[1])
            
            rest[[e]][cont,1] = rest[[e]][cont,1]+ (me - m)^2
            rest[[e]][cont,2] = rest[[e]][cont,2]+ me
            rest[[e]][cont,3:4] = rest[[e]][cont,3:4] + ci.transformed
            rest[[e]][cont,5] = rest[[e]][cont,5]+ ((ci.transformed[1] < m) & (ci.transformed[2] > m))
            rest[[e]][cont,6] = rest[[e]][cont,6]+ (ci.transformed[2] -ci.transformed[1])
            }
         
          count = count+1
        }
        
      }
      for (b in 1:nbetas){
        res[[b]][cont,] = res[[b]][cont,]/count
        rest[[b]][cont,] = rest[[b]][cont,]/count
        }
    }
    return(list("direct" = res, "transformed" = rest))
  }



### BOOTSTRAP CONFIDENCE INTERVALS

simulate.times <- function(theta, tauvec, stressvec, ITvec, C, N){
  
  #C uniform samples
  U=matrix(NA, nrow=C, ncol=N)
  for (c in 1:C){U[c,] = runif(N,0,1)}
  
  #inverting from marginals
  Tm = matrix(NA, nrow=C, ncol=N)
  for (c in 1:C){
    lambda = exp(theta[(2*c-1)]+ theta[2*c]*stressvec)
    lim = model.distribution.exponential(lambda, tauvec, stressvec, tauvec[1])
    h = shifted.times(tauvec, lambda)
    for (k in 1:N){
      if (U[c,k] < lim ){ Tm[c,k]= -log(1-U[c,k])*lambda[1] }else{ Tm[c,k] = -log(1-U[c,k])*lambda[2] - h[2]}
    }
  }
  
  times = apply(Tm,2, min)
  crisks = apply(Tm,2, which.min)
  
  return(cbind(times,crisks))
}

simulate.times2 <- function(theta, tauvec, stressvec, ITvec, C, N){
  
  #C uniform samples
  U=matrix(NA, nrow=C, ncol=N)
  for (c in 1:C){U[c,] = runif(N,0,1)}
  
  #inverting from marginals
  # Tm = matrix(NA, nrow=C, ncol=N)
  
  #sample_size=N
  ni1=1 
  times = c()
  crisks = c()
  
  
  for (i in 1:2){
    Tmi=c()
    Cmi=c()
    # for every U_k 
    for (k in 1:ncol(U)){

      # apply stress
      Tm = matrix(NA, nrow=C, ncol=ncol(U))
      
      for (c in 1:C){
        # for each competing risk
        lambda = exp(theta[(2*c-1)]+ theta[2*c]*stressvec)
        h = shifted.times(tauvec, lambda)
        Tm[c,k] = -log(1-U[c,k])*lambda[i]-h[i]
      }
    
    Tmi[k] = min(Tm[,k])
    Cmi[k] = which(min(Tm[,k])==Tm[,k])
    }
    
    #sort elements of Tmi
    s = sort(Tmi, index.return=TRUE)
    sTmi = s$x
    sCmi=s$ix
    #find ni
    ni = sum(sTmi<tauvec[i])
    delcol = which(Tmi<tauvec[i])
    
    times[ni1:(ni+ni1-1)] =  sTmi[1:ni]
    crisks[ni1:(ni+ni1-1)] =  sCmi[1:ni]
    
    ni1=ni+ni1
    
    U = U[,-delcol] 
}
  
  ni = N
  times[ni1:length(times)] = times[length(times)]+10
  crisks[ni1:length(times)] = 3
  
  return(cbind(times,crisks))
}

build.censored.sample <- function(simulated.times, ITvec, C){
  
  n = c()
  for(c in 1:C){n[c]=sum((simulated.times[,1] <= ITvec[1])&(simulated.times[,2] == c))}
  
  failures = n
  
  for(l in 2:length(ITvec)){
   
    failures[2*l-1] = sum((simulated.times[,1] <= ITvec[l])&(simulated.times[,2] ==1)) 
    n[2*l-1] = failures[2*l-1] - failures[2*l-3]
    
    failures[2*l] = sum((simulated.times[,1] <= ITvec[l])&(simulated.times[,2] ==2)) 
    n[2*l] = failures[2*l] - failures[2*l-2]
    }
  n[2*l+1] = length(simulated.times[,1])-sum(n)
  
  return(n)
}

z.function <- function(bestimates, estimate){
  m = mean(bestimates <= estimate)
  if (m==0){m=0.0001}
  if (m==1){m=0.9999}
  return(qnorm(m,0,1))
}

jacknife.function <-function(n, beta, char.function,C,x0,t0,alpha0){
  
  NL = length(n)-1
  jacknifematrix = c() #matrix(NA, nrow=NL, ncol=1)
  
  for(nij in 1:NL){
    
    nn = n
    #delete an observed failure and add it to the survival failures
    nn[nij] = nn[nij]-1
    nn[NL+1] = nn[NL+1]+1
    
    nestimate = estimate.beta(N, tauvec, stressvec, ITvec, nn, initial = c(7, -0.01, 5, -0.02), C=C, beta=beta)
    nchar =  char.function(nestimate, C=C,x0=x0,t0=t0, alpha0=alpha0)
    
    jacknifematrix[nij] = nchar
  }
  
  #weighted mean
  jack = sum(jacknifematrix*n[1:NL])/sum(n[1:NL])
  return(jack)
}

gamma.function <- function( n, beta, char.function,C,x0,t0,alpha0){
  
  NL = length(n)-1
  gammamatrix = c()
  jack = jacknife.function(n, beta, char.function,C,x0,t0,alpha0)
  
  for(nij in 1:NL){
    
    nn = n
    #delete an observed failure and add it to the survival failures
    nn[nij] = nn[nij]-1
    nn[NL+1] = nn[NL+1]+1
    
    nestimate = estimate.beta(N, tauvec, stressvec, ITvec, nn, initial = c(7, -0.01, 5, -0.02), C=C, beta=beta)
    #nchars = lapply(nestimators, char.function, C=C,x0=x0,t0=t0, alpha0=alpha0)
    gammamatrix[nij] = char.function(nestimate, C=C,x0=x0,t0=t0, alpha0=alpha0)
  }
  
  
  num =  sum(n[1:NL]*(gammamatrix-jack)^3)
  denom = sum(n[1:NL]*(gammamatrix-jack)^2)
  
  return((1/6)*num*denom^(-3/2))
  
}

gamma1.function <-function(z0,gamma,alpha){
  qn = qnorm(1-alpha,0,1)
  q = z0+(z0-qn)/(1-gamma*(z0-qn))
  return(pnorm(q,0,1))
}

gamma2.function <-function(z0, gamma, alpha){
  qn = qnorm(1-alpha,0,1)
  q = z0+(z0+qn)/(1-gamma*(z0+qn))
  return(pnorm(q,0,1))
}

bootstrap.estimates <- function(theta.estimate, tauvec, stressvec, ITvec, C, N, char.function, x0, t0, alpha0, beta, B,  alpha=0.05){
  
  
  t1_ = sum(ITvec <= tauvec[1])
  
  #nbetas = length(beta.list)+1
  bestimates = c() #matrix(NA, nrow = B, ncol = nbetas)
  
  for(b in 1:B){
    
    valid = F
    #generate bootstrap sample
    while (!valid){
      simulated_sample =  simulate.times(theta.estimate, tauvec, stressvec, ITvec, C, N)
      n = build.censored.sample(simulated_sample, ITvec, C)
      if( validility.check(n, t1_)){valid = TRUE}
    }
    
    #fit MDPDE
    estimate = estimate.beta(N, tauvec, stressvec, ITvec, n, initial = theta.estimate, C=C, beta=beta)
    bestimates[b] = char.function(theta = estimate, C=C, x0=x0, t0=t0, alpha0=alpha0)
  }
    #return bootstrap sample
    return(bestimates)
}

bootstrap.beta.ci <- function(n, estimate, tauvec, stressvec, ITvec, C, N, char.function, x0, t0, alpha0, beta, B,  alpha=0.05){
  
  
  bestimates = bootstrap.estimates(estimate, tauvec, stressvec, ITvec, C, N, char.function, x0, t0, alpha0, beta, B,  alpha=0.05)
  
  #compute z
  sestimate =  char.function(estimate, C=C, x0=x0, t0=t0, alpha0 = alpha0)
  z0 = z.function(bestimates, sestimate)
   
  #compute gammaest
  gamma = gamma.function(n, beta,  char.function,C,x0,t0,alpha0)
  
  #compute gamma1
  gamma1s = gamma1.function(z0,gamma,alpha)
  
  #compute gamma2
  gamma2s = gamma2.function(z0,gamma,alpha)
  
  q1 = max(floor(gamma1s*B),1)
  q2 =  min(floor(gamma2s*B), B)
  
  bestimates = sort(bestimates)
  return( c(bestimates[q1], bestimates[q2]) )
  
  }


simulate.bootstrap.charasteristic <-function(theta,  cont.vec,  N, tauvec, stressvec, ITvec, x0=30, t0=50, alpha = 0.05,
                                   C=2, R=50, B=20, beta.list = c(0.2,0.4,0.6,0.8,1),
                                   charasteristic="mean.lifetime", alpha0=0.05){
  
  char.function = match.fun(charasteristic)
  
  f <- paste0(charasteristic, ".ci")
  ci.function = match.fun(f)
  
  beta.list = c(0,beta.list)
  nbetas = length(beta.list)
  res = list()
  for(b in 1:nbetas){res[[b]] = matrix(0, nrow =  length(cont.vec), ncol = 6)}
  
  
  n1_ = sum(ITvec <= tauvec[1])
  m = char.function(theta=theta, C=C, x0=x0, t0=t0, alpha0=alpha0)
  
  for(cont in 1:length(cont.vec)){
    
    count = 0
    
    for(r in 1:R){
      
      valid = FALSE
      while (!valid){
        n = simulate.sample(theta, cont.vec[cont], N, tauvec, stressvec, ITvec,  seed = 10*r, C=C)
        if( validility.check(n, n1_)){valid = TRUE}else{r=7*r}
      }
      
      estimators = estimate.function(N, tauvec, stressvec, ITvec, n, initial = c(7, -0.01, 5, -0.02), C=C, beta.list = beta.list[2:nbetas])
      
      if(sum(sapply(estimators, function(x) sum(is.na(x))))<1){ #none of the estimators is nan
        
        cis = sapply(c(1:nbetas), function (e) bootstrap.beta.ci(n, estimators[[e]], tauvec, stressvec, ITvec, C, N, char.function, x0, t0, alpha0, beta.list[e], B,  alpha))
        
        for (e in 1:nbetas){
          
          me = char.function(theta = estimators[[e]], C=C, x0=x0, t0=t0, alpha0=alpha0)
          ci = cis[,e]
          
          res[[e]][cont,1] = res[[e]][cont,1]+ (me - m)^2
          res[[e]][cont,2] = res[[e]][cont,2]+ me
          res[[e]][cont,3:4] = res[[e]][cont,3:4] + ci
          res[[e]][cont,5] = res[[e]][cont,5]+ ((ci[1] < m) & (ci[2] > m))
          res[[e]][cont,6] = res[[e]][cont,6]+ (ci[2] -ci[1])
        }
        
        count = count+1
      }
    
    }
    for (b in 1:nbetas){res[[b]][cont,] = res[[b]][cont,]/count}
  }
  return(res)
}

