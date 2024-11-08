# Full conditional distribution for sigma
atualizarSIGMA<-function(a0,A0,endogena,instrumentos,delta,y,covariaveis,beta,omega,n,t,tau,rho,W){
  v1 <- y - cbind(covariaveis,endogena)%*%beta
  #v1 <- y - covariaveis%*%beta
  v2 <- endogena - instrumentos%*%delta
  v  <- cbind(v1,v2)
  
  omega_inv <- chol2inv(chol(omega))
  tau_inv   <- (diag(rowSums(W)) - rho*W)
  #omega_inv <- diag(1,t)
  #tau_inv   <- diag(1,n)
    
  a1 <- a0 + n*t
  A1 <- A0 + t(v)%*%kronecker(tau_inv,omega_inv)%*%v
  sigma <- riwish(a1,A1)
  return(sigma)
} 

# Full conditional distribution for beta
atualizarBETA<-function(b0,B0,sigma,endogena,instrumentos,delta,y,covariaveis,omega,n,tau,rho,W){
  sigma11  <- sigma[1,1] ; sigma12 <- sigma[1,2]  
  sigma21  <- sigma[2,1] ; sigma22 <- sigma[2,2]
  
  omega_inv     <- chol2inv(chol(omega))
  tau_inv       <- (diag(rowSums(W)) - rho*W)
  #omega_inv <- diag(1,t)
  #tau_inv   <- diag(1,n)
  spat_temp     <- kronecker(tau,omega) 
  spat_temp_inv <- kronecker(tau_inv,omega_inv)
  
  sigma.yx    <- sigma11*spat_temp - (sigma12*spat_temp)%*%((1/sigma22)*spat_temp_inv)%*%(sigma21*spat_temp)
  sigma.yx_inv<- chol2inv(chol(sigma.yx))
  v           <- y - (sigma12*spat_temp)%*%((1/sigma22)*spat_temp_inv)%*%(endogena-instrumentos%*%delta)
  matrix.x    <- cbind(covariaveis,endogena)
  #matrix.x    <- covariaveis
  
  B1       <- solve( solve(B0) + t(matrix.x)%*%sigma.yx_inv%*%matrix.x )
  b1       <- B1%*%( solve(B0)%*%b0 + t(matrix.x)%*%sigma.yx_inv%*%v )
  beta     <- rmvnorm(1,b1,B1)
  return(beta)
}

# Full conditional distribution for delta
atualizarDELTA<-function(d0,D0,sigma,endogena,instrumentos,dados,covariaveis,beta,omega,n,tau,rho,W){
  sigma11  <- sigma[1,1] ; sigma12 <- sigma[1,2]  
  sigma21  <- sigma[2,1] ; sigma22 <- sigma[2,2]
  
  omega_inv     <- chol2inv(chol(omega))
  tau_inv       <- (diag(rowSums(W)) - rho*W)
  #omega_inv <- diag(1,t)
  #tau_inv   <- diag(1,n)
  spat_temp     <- kronecker(tau,omega) 
  spat_temp_inv <- kronecker(tau_inv,omega_inv)
  
  sigma.xy    <- sigma22*spat_temp - (sigma21*spat_temp)%*%((1/sigma11)*spat_temp_inv)%*%(sigma12*spat_temp)
  sigma.xy_inv<- chol2inv(chol(sigma.xy))
  matrix.x    <- cbind(covariaveis,endogena)
  #matrix.x    <- covariaveis
  v           <- endogena - (sigma21*spat_temp)%*%((1/sigma11)*spat_temp_inv)%*%(y-matrix.x%*%beta)
  
  D1      <- solve( solve(D0) + t(instrumentos)%*%sigma.xy_inv%*%instrumentos )
  d1      <- D1%*%( solve(D0)%*%d0 + t(instrumentos)%*%sigma.xy_inv%*%v )
  delta   <- rmvnorm(1,d1,D1)
  return(delta)
}

# Full conditional for psi
condicionalPSI<-function(a,b,psi,endogena,instrumentos,delta,y,covariaveis,beta,sigma,tau,rho,W,sigma2,n,t,p){
  v1 <- y - cbind(covariaveis,endogena)%*%beta
  #v1 <- y - covariaveis%*%beta
  v2 <- endogena - instrumentos%*%delta
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
# Metropolis-Hasting for phi
atualizarPHI<-function(a,b,phi,endogena,instrumentos,delta,y,covariaveis,beta,sigma,tau,rho,W,sigma2,n,t,p){
  valoratual    <- log((1+phi)/(1-phi))
  valorproposto <- rnorm(1,valoratual,sqrt(0.05))
  candidato     <- exp(condicionalPSI(a,b,valorproposto,endogena,instrumentos,delta,y,covariaveis,beta,sigma,tau,rho,W,sigma2,n,t,p)
                       -condicionalPSI(a,b,valoratual,endogena,instrumentos,delta,y,covariaveis,beta,sigma,tau,rho,W,sigma2,n,t,p))
  
  chanceaceitar <- min(1,candidato)
  PHIfinal      <- ifelse(runif(1)<chanceaceitar,(exp(valorproposto)-1)/(1+exp(valorproposto)),phi)
  return(PHIfinal)
}

# Full conditional for GAMMA
condicionalGAMMA<-function(a,b,gama,W,endogena,instrumentos,delta,y,covariaveis,beta,sigma,t,p){
  v1 <- y - cbind(covariaveis,endogena)%*%beta
  #v1 <- y - covariaveis%*%beta
  v2 <- endogena - instrumentos%*%delta
  v  <- c(v1,v2)
  tau<- chol2inv(chol(diag(rowSums(W)) - (exp(gama)/(1+exp(gama)))*W))
  
  omega_inv     <- chol2inv(chol(omega))
  #omega_inv <- diag(1,t)
  tau_inv       <- (diag(rowSums(W)) - (exp(gama)/(1+exp(gama)))*W)
  sigma_inv     <- chol2inv(chol(sigma))
  
  #priori  <- -0.5*((exp(gama)/(1+exp(gama)))-a)^2/b
  priori  <- -0.5*(gama-a)^2/b
  verossi <- -(0.5*t*p)*log(det(tau)) - 0.5*(v%*%kronecker(sigma_inv,kronecker(tau_inv,omega_inv))%*%v) + log(exp(gama)/(1+exp(gama))^2)
  funcao  <- priori + verossi
  return(funcao)
}
# Metropolis-Hasting for rho
atualizarRHO<-function(a,b,rho,W,endogena,instrumentos,delta,y,covariaveis,beta,sigma,t,p){
  valoratual    <- log(rho/(1-rho))
  valorproposto <- rnorm(1,valoratual,sqrt(0.01))
  candidato     <- exp(condicionalGAMMA(a,b,valorproposto,W,endogena,instrumentos,delta,y,covariaveis,beta,sigma,t,p)
                       -condicionalGAMMA(a,b,valoratual,W,endogena,instrumentos,delta,y,covariaveis,beta,sigma,t,p))
  
  chanceaceitar <- min(1,candidato)
  RHOfinal      <- ifelse(runif(1)<chanceaceitar,exp(valorproposto)/(1+exp(valorproposto)),rho)
  return(RHOfinal)
}
