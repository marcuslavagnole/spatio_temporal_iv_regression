# Load libraries
library(mvtnorm)
library(MCMCpack)
library(Matrix)
library(msm)
library(abind)
library(dplyr)
#memory.limit(size=34440)

# Set working directory to file location and load utils
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Full_Conditionals_IVST.R")
source("Missing_Imputation.R")
source("aux_func.R")

# Load data ans set its configuration
load("data.RData")
load("W_MS.RData")

dados = dados[which(dados$uf=="50"),]
aux1  = dados[which(dados$nomemun=="SELVÍRIA"),]
aux2  = dados[which(dados$nomemun=="SETE QUEDAS"),]
dados[c(which(dados$nomemun=="SETE QUEDAS"),which(dados$nomemun=="SELVÍRIA")),] = rbind(aux1,aux2)
colSums(is.na(dados))

# Define panel dimensions
n     <- length(unique(dados$codmun))
t     <- length(unique(dados$ano))
p     <- 2

# Define dependent, indepent, and instrumental variables
attach(dados)
y     <- lpib_real_1000
x     <- cbind(rep(1,dim(dados)[1]),lvinculos31dez,lvalor_PBF,lreceitas_mun,lValor_BPC,lPrevidencia)
x_end <- lDENS_BANDA_LARGA_SCM
I     <- cbind(rep(1,dim(dados)[1]),l512kbps)

y_ind     <- which(is.na(y))
x_end_ind <- which(is.na(x_end))
indices   <- c(y_ind,(n*t+x_end_ind))

# Define number of iterations to run
NN <- 21000

# Create auxiliary objects 
sigma <- array(dim=c(p,p,NN))
beta  <- matrix(NA,NN,(dim(x)[2]+1))
delta <- matrix(NA,NN,(dim(I)[2]))
sigma2<- matrix(NA,NN,1)
phi   <- matrix(NA,NN,1)
rho   <- matrix(NA,NN,1)
imput <- matrix(NA,NN,length(indices))
  
# Set the initial values
sigma[,,1] <- diag(1,p)
reg1 = lm(y ~ lvinculos31dez+lvalor_PBF+lreceitas_mun+lValor_BPC+lPrevidencia+lDENS_BANDA_LARGA_SCM)
beta[1,]   <- summary(reg1)$coefficients[,1]
reg2 = lm(x_end ~ l512kbps)
delta[1,]  <- summary(reg2)$coefficients[,1]
sigma2[1,1]<- 1
phi[1,1]   <- 0.8
rho[1,1]   <- 0.9
omega <- formataOMEGA(phi[1,1],sigma2[1,1],t)
tau   <- solve(diag(rowSums(W)) - rho[1,1]*W)
x_end[x_end_ind] <- rnorm(length(x_end_ind),mean(x_end[-x_end_ind]),0.1)

#MCMC
for (k in 2:NN){
  sigma[,,k] <- atualizarSIGMA(3,diag(1,2),x_end,I,delta[k-1,],y,x,beta[k-1,],omega,n,t,tau,rho[k-1,1],W)
  sigma2[k,1]<- 1
  phi[k,1]   <- atualizarPHI(0,100,phi[k-1,1],x_end,I,delta[k-1,],y,x,beta[k-1,],sigma[,,k],tau,rho[k-1,1],W,sigma2[k,1],n,t,p)
  omega      <- formataOMEGA(phi[k,1],sigma2[k,1],t)
  rho[k,1]   <- atualizarRHO(0,100,rho[k-1,1],W,x_end,I,delta[k-1,],y,x,beta[k-1,],sigma[,,k],t,p)
  tau        <- solve(diag(rowSums(W)) - rho[k,1]*W)
  delta[k,]  <- atualizarDELTA(rep(0,2),diag(100,2),sigma[,,k],x_end,I,y,x,beta[k-1,],omega,n,tau,rho[k,1],W)
  beta[k,]   <- atualizarBETA(rep(0,7),diag(100,7),sigma[,,k],x_end,I,delta[k,],y,x,omega,n,tau,rho[k,1],W)
  
  # Imputing value
  imput[k,]  <- imputMISSING(indices,y,x,beta[k,],x_end,I,delta[k,],sigma[,,k],tau,omega,n,t)
  #y[y_ind]         <- imput[k,1:length(y_ind)]
  x_end[x_end_ind] <- imput[k,(length(y_ind)+1):(length(x_end_ind)+length(y_ind))]
}
