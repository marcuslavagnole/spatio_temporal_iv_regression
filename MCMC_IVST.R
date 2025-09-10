#' Spatio-Temporal Instrumental Variables Regression with Missing Data
#'
#' This function refers to a Bayesian approach for a Spatio-Temporal 
#' Instrumental Variables Regression with Missing Data.
#'
#' @param y (nxt)-dimensional vector of responses.
#' @param x Matrix - (nxt)x(q+1) - of exogenous variables (include the intercept).
#' @param x_end (nxt)-dimensional vector of endogenous variables.
#' @param z Matrix - (nxt)x(l+1) - of instruments (include the intercept).
#' @param W Adjacency matrix.
#' @param n Number of units.
#' @param t Number of periods.
#' @param n_mcmc Number of iterations.
#' @param burnin_mcmc Number of initial iterations to be discarded.
#' @param thin_mcmc Thinning parameter.
#' 
#' @return A list with the chains of all parameters of interest.

## Packages
require(mvtnorm); require(MCMCpack)

## MCMC
spatio_temporal_IV <- function(y,x,x_end,z,W,n,t,n_mcmc,burnin_mcmc,thin_mcmc){
  y_ind     <- which(is.na(y))
  x_end_ind <- which(is.na(x_end))
  indices   <- c(y_ind,(n*t+x_end_ind))
  resultado <- list()
  # Create auxiliary objects 
  sigma <- array(dim=c(2,2,n_mcmc))
  beta  <- matrix(NA,n_mcmc,(dim(x)[2]+1))
  delta <- matrix(NA,n_mcmc,(dim(z)[2]))
  phi   <- matrix(NA,n_mcmc,1)
  rho   <- matrix(NA,n_mcmc,1)
  impute<- matrix(NA,n_mcmc,length(indices))
  # Set the initial values
  sigma[,,1]<- diag(1,2)
  reg1 = lm(y ~ -1 + x + x_end)
  beta[1,]  <- summary(reg1)$coefficients[,1]
  reg2 = lm(x_end ~ -1 + z)
  delta[1,] <- summary(reg2)$coefficients[,1]
  phi[1,1]  <- runif(1)
  rho[1,1]  <- runif(1)
  omega <- formataOMEGA(phi[1,1],1,t)
  tau   <- solve(diag(rowSums(W)) - rho[1,1]*W)
  y[y_ind] <- rnorm(length(y_ind),mean(y[-y_ind]),0.1)
  x_end[x_end_ind] <- rnorm(length(x_end_ind),mean(x_end[-x_end_ind]),0.1)
  #MCMC
  for (k in 2:n_mcmc){
    sigma[,,k] <- atualizarSIGMA(3,diag(1,2),x_end,z,delta[k-1,],y,x,beta[k-1,],omega,n,t,tau,rho[k-1,1],W)
    phi[k,1]   <- atualizarPHI(0,100,phi[k-1,1],x_end,z,delta[k-1,],y,x,beta[k-1,],sigma[,,k],tau,rho[k-1,1],W,1,n,t,2)
    omega      <- formataOMEGA(phi[k,1],1,t)
    rho[k,1]   <- atualizarRHO(0,100,rho[k-1,1],W,x_end,z,delta[k-1,],y,x,beta[k-1,],sigma[,,k],omega,t,2)
    tau        <- solve(diag(rowSums(W)) - rho[k,1]*W)
    delta[k,]  <- atualizarDELTA(rep(0,dim(z)[2]),diag(100,dim(z)[2]),sigma[,,k],x_end,z,y,x,beta[k-1,],omega,n,tau,rho[k,1],W)
    beta[k,]   <- atualizarBETA(rep(0,(dim(x)[2]+1)),diag(100,(dim(x)[2]+1)),sigma[,,k],x_end,z,delta[k,],y,x,omega,n,tau,rho[k,1],W)
    # Imputing value
    impute[k,]  <- imputeMISSING(indices,y,x,beta[k,],x_end,z,delta[k,],sigma[,,k],tau,omega,n,t)
    y[y_ind]         <- impute[k,1:length(y_ind)]
    x_end[x_end_ind] <- impute[k,(length(y_ind)+1):(length(x_end_ind)+length(y_ind))]
  }
  resultado[['beta_star']]<- beta[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),1:dim(x)[2]]
  resultado[['beta']]     <- beta[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),(dim(x)[2]+1)]
  resultado[['delta']]    <- delta[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[['rho']]      <- rho[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),1]
  resultado[['phi']]      <- phi[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),1]
  resultado[['sigma']]    <- sigma[,,seq(burnin_mcmc+1,n_mcmc,thin_mcmc)]
  
  return(resultado)
}
  

#######################################################################
## Auxiliary functions: Sampling from full conditional distributions ##
#######################################################################

### Full conditional distribution for sigma
atualizarSIGMA<-function(a0,A0,x_end,z,delta,y,x,beta,omega,n,t,tau,rho,W){
  v1 <- y - cbind(x,x_end)%*%beta
  v2 <- x_end - z%*%delta
  v  <- cbind(v1,v2)
  
  omega_inv <- chol2inv(chol(omega))
  tau_inv   <- (diag(rowSums(W)) - rho*W)
    
  a1 <- a0 + n*t
  A1 <- A0 + t(v)%*%kronecker(tau_inv,omega_inv)%*%v
  sigma <- riwish(a1,A1)
  return(sigma)
} 

### Full conditional distribution for beta
atualizarBETA<-function(b0,B0,sigma,x_end,z,delta,y,x,omega,n,tau,rho,W){
  sigma11  <- sigma[1,1] ; sigma12 <- sigma[1,2]  
  sigma21  <- sigma[2,1] ; sigma22 <- sigma[2,2]
  
  omega_inv     <- chol2inv(chol(omega))
  tau_inv       <- (diag(rowSums(W)) - rho*W)
  spat_temp     <- kronecker(tau,omega) 
  spat_temp_inv <- kronecker(tau_inv,omega_inv)
  
  sigma.yx    <- sigma11*spat_temp - (sigma12*spat_temp)%*%((1/sigma22)*spat_temp_inv)%*%(sigma21*spat_temp)
  sigma.yx_inv<- chol2inv(chol(sigma.yx))
  v           <- y - (sigma12*spat_temp)%*%((1/sigma22)*spat_temp_inv)%*%(x_end-z%*%delta)
  matrix.x    <- cbind(x,x_end)
  
  B1       <- solve( solve(B0) + t(matrix.x)%*%sigma.yx_inv%*%matrix.x )
  b1       <- B1%*%( solve(B0)%*%b0 + t(matrix.x)%*%sigma.yx_inv%*%v )
  beta     <- rmvnorm(1,b1,B1)
  return(beta)
}

### Full conditional distribution for delta
atualizarDELTA<-function(d0,D0,sigma,x_end,z,y,x,beta,omega,n,tau,rho,W){
  sigma11  <- sigma[1,1] ; sigma12 <- sigma[1,2]  
  sigma21  <- sigma[2,1] ; sigma22 <- sigma[2,2]
  
  omega_inv     <- chol2inv(chol(omega))
  tau_inv       <- (diag(rowSums(W)) - rho*W)
  spat_temp     <- kronecker(tau,omega) 
  spat_temp_inv <- kronecker(tau_inv,omega_inv)
  
  sigma.xy    <- sigma22*spat_temp - (sigma21*spat_temp)%*%((1/sigma11)*spat_temp_inv)%*%(sigma12*spat_temp)
  sigma.xy_inv<- chol2inv(chol(sigma.xy))
  matrix.x    <- cbind(x,x_end)
  v           <- x_end - (sigma21*spat_temp)%*%((1/sigma11)*spat_temp_inv)%*%(y-matrix.x%*%beta)
  
  D1      <- solve( solve(D0) + t(z)%*%sigma.xy_inv%*%z )
  d1      <- D1%*%( solve(D0)%*%d0 + t(z)%*%sigma.xy_inv%*%v )
  delta   <- rmvnorm(1,d1,D1)
  return(delta)
}

### Metropolis-Hastings for phi
atualizarPHI<-function(a,b,phi,x_end,z,delta,y,x,beta,sigma,tau,rho,W,sigma2,n,t,p){
  valoratual    <- log((1+phi)/(1-phi))
  valorproposto <- rnorm(1,valoratual,sqrt(0.05))
  candidato     <- exp(condicionalPSI(a,b,valorproposto,x_end,z,delta,y,x,beta,sigma,tau,rho,W,sigma2,n,t,p)
                       -condicionalPSI(a,b,valoratual,x_end,z,delta,y,x,beta,sigma,tau,rho,W,sigma2,n,t,p))
  
  chanceaceitar <- min(1,candidato)
  PHIfinal      <- ifelse(runif(1)<chanceaceitar,(exp(valorproposto)-1)/(1+exp(valorproposto)),phi)
  return(PHIfinal)
}

### Full conditional for psi
condicionalPSI<-function(a,b,psi,x_end,z,delta,y,x,beta,sigma,tau,rho,W,sigma2,n,t,p){
  v1 <- y - cbind(x,x_end)%*%beta
  v2 <- x_end - z%*%delta
  v  <- c(v1,v2)
  
  ar1.aux   <- NULL
  for(j in 1:(t-1)){
    ar1.aux[j] <- ((exp(psi)-1)/(1+exp(psi)))^j 
  }
  m.ar1     <- diag(1,t)
  for(j in 1:(t-1)){
    m.ar1[j,(j+1):t] <- ar1.aux[1:(t-j)] 
  }
  m.ar1[lower.tri(m.ar1)] <- t(m.ar1)[lower.tri(m.ar1)]
  omega <- (sigma2/(1-((exp(psi)-1)/(1+exp(psi)))^2))*m.ar1
  
  omega_inv     <- chol2inv(chol(omega))
  tau_inv       <- (diag(rowSums(W)) - rho*W)
  sigma_inv     <- chol2inv(chol(sigma))
  
  #priori  <- -0.5*(((exp(psi)-1)/(1+exp(psi)))-a)^2/b
  priori  <- -0.5*(psi-a)^2/b
  verossi <- -(0.5*n*p)*log(det(omega)) - 0.5*(v%*%kronecker(sigma_inv,kronecker(tau_inv,omega_inv))%*%v) + log((2*exp(psi))/(1+exp(psi))^2)
  funcao  <- priori + verossi
  return(funcao)
}

### Metropolis-Hastings for rho
atualizarRHO<-function(a,b,rho,W,x_end,z,delta,y,x,beta,sigma,omega,t,p){
  valoratual    <- log(rho/(1-rho))
  valorproposto <- rnorm(1,valoratual,sqrt(0.01))
  candidato     <- exp(condicionalGAMMA(a,b,valorproposto,W,x_end,z,delta,y,x,beta,sigma,omega,t,p)
                       -condicionalGAMMA(a,b,valoratual,W,x_end,z,delta,y,x,beta,sigma,omega,t,p))
  
  chanceaceitar <- min(1,candidato)
  RHOfinal      <- ifelse(runif(1)<chanceaceitar,exp(valorproposto)/(1+exp(valorproposto)),rho)
  return(RHOfinal)
}

### Full conditional for GAMMA
condicionalGAMMA<-function(a,b,gama,W,x_end,z,delta,y,x,beta,sigma,omega,t,p){
  v1 <- y - cbind(x,x_end)%*%beta
  v2 <- x_end - z%*%delta
  v  <- c(v1,v2)
  tau<- chol2inv(chol(diag(rowSums(W)) - (exp(gama)/(1+exp(gama)))*W))
  
  omega_inv     <- chol2inv(chol(omega))
  tau_inv       <- (diag(rowSums(W)) - (exp(gama)/(1+exp(gama)))*W)
  sigma_inv     <- chol2inv(chol(sigma))
  
  priori  <- -0.5*(gama-a)^2/b
  verossi <- -(0.5*t*p)*log(det(tau)) - 0.5*(v%*%kronecker(sigma_inv,kronecker(tau_inv,omega_inv))%*%v) + log(exp(gama)/(1+exp(gama))^2)
  funcao  <- priori + verossi
  return(funcao)
}

### Missing imputation
imputeMISSING<-function(entries,y,x,beta,x_end,z,delta,sigma,tau,omega,n,t){
  sigma_bar   <- kronecker(sigma,kronecker(tau,omega))
  sigma_bar11 <- sigma_bar[entries,entries]
  sigma_bar12 <- sigma_bar[entries,c(1:(2*n*t))[-entries]]
  sigma_bar21 <- t(sigma_bar12)
  sigma_bar22 <- sigma_bar[c(1:(2*n*t))[-entries],c(1:(2*n*t))[-entries]]
  
  v           <- c(y,x_end)
  matrix.x    <- cbind(x,x_end)
  mu          <- c(matrix.x%*%beta,z%*%delta)
  
  sigma_bar22_inv <- chol2inv(chol(sigma_bar22))
  media       <- mu[entries] + sigma_bar12%*%sigma_bar22_inv%*%(v[-entries]-mu[-entries])
  covariancia <- sigma_bar11 + sigma_bar12%*%sigma_bar22_inv%*%sigma_bar21
  imputation  <- rmvnorm(1,media,covariancia)
  return(imputation)
}


##################################################
## Auxiliary functions: Supplementary functions ##
##################################################

formataOMEGA<-function(phi,sigma2,t){
  ar1.aux   <- NULL
  for(j in 1:(t-1)){
    ar1.aux[j] <- phi^j 
  }
  m.ar1     <- diag(1,t)
  for(j in 1:(t-1)){
    m.ar1[j,(j+1):t] <- ar1.aux[1:(t-j)] 
  }
  m.ar1[lower.tri(m.ar1)] <- t(m.ar1)[lower.tri(m.ar1)]
  omega <- (sigma2/(1-phi^2))*m.ar1
  return(omega)

}
