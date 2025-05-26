# Missing imputation
imputMISSING<-function(entries,y,covariaveis,beta,endogena,instrumentos,delta,sigma,tau,omega,n,t){
  sigma_bar   <- kronecker(sigma,kronecker(tau,omega))
  sigma_bar11 <- sigma_bar[entries,entries]
  sigma_bar12 <- sigma_bar[entries,c(1:(2*n*t))[-entries]]
  sigma_bar21 <- t(sigma_bar12)
  sigma_bar22 <- sigma_bar[c(1:(2*n*t))[-entries],c(1:(2*n*t))[-entries]]
  
  v           <- c(y,endogena)
  matrix.x    <- cbind(covariaveis,endogena)
  mu          <- c(matrix.x%*%beta,instrumentos%*%delta)
  
  sigma_bar22_inv <- chol2inv(chol(sigma_bar22))
  media       <- mu[entries] + sigma_bar12%*%sigma_bar22_inv%*%(v[-entries]-mu[-entries])
  covariancia <- sigma_bar11 + sigma_bar12%*%sigma_bar22_inv%*%sigma_bar21
  imputation  <- rmvnorm(1,media,covariancia)
  return(imputation)
}
