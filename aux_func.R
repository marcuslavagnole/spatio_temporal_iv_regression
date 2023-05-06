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