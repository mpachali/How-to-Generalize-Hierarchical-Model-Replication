##################################################################################################
#########Constrained hierarchical prior using the log-normal distribution#########################
#######################Marginal-conditional decomposition#########################################
###############Simulation study (Section 5) in Pachali, Kurz and Otter, QME 2020##################
##################################################################################################

rm(list=ls())

library(Rcpp)
library(devtools)
library(RcppArmadillo)
library(MASS)
library(lattice)
library(Matrix)
library(bayesm)
library(tikzDevice)
library(plyr)
library(numDeriv)

set.seed(77) 

###
### FUNCTION GENERATING CHOICE DATA FROM LOGIT MODEL
###
simmnlv2 = function(n,beta){
  
  #
  # p. rossi 2004
  # Modified by Max Pachali & Thomas Otter 2014/2015/2016
  #
  # Purpose: simulate from MNL (including X values)
  #
  # Arguments:
  # n is number of choice tasks
  # beta is true parm value
  #
  # Output:
  # list of X 
  # y (indicator of choice --> 1, ...,p)
  # prob is a n x p matrix of choice probs
  #
  
  k = length(beta)
  # MATRIX OF ALL POSSIBLE ATTRIBUTE LEVEL COMBINATIONS
  X_full_original = rbind(c(1,0,1,0),c(0,1,1,0),c(1,0,0,1),c(0,1,0,1),c(0,0,0,0))
  # DEDUCE p
  p = dim(X_full_original)[1]
  # REPLICATE ORIGINAL MATRIX n TIMES
  for(i in 1:n){
    if(i==1){
      X_full = X_full_original
    }else{
      X_full = rbind(X_full,X_full_original)
    }
  }
  # ADD PRICE
  X_full = cbind(X_full[,1:2],runif(n*p,min=0.5,max=3),X_full[,3:4])
  # DELETE PRICES FOR OUTSIDE OPTION
  for(i in 1:n){
    X_full[((i-1)*p+p),3] = 0
  }
  # CONSTRUCT PROBABILITIES
  Xbeta=X_full%*%beta
  p=nrow(Xbeta)/n
  Xbeta=matrix(Xbeta,byrow=TRUE,ncol=p) #reshape Xbeta matrix to be of (nxp)-dimension
  Prob=exp(Xbeta)
  iota=rep(1,p)
  denom=Prob%*%iota
  Prob=Prob/as.vector(denom)
  # DRAW CHOICES Y
  y=vector("double",n)
  ind=1:p
  for (i in 1:n){
    yvec=rmultinom(1,1,Prob[i,])
    y[i]=ind%*%yvec
  }
  
  return(list(y=y,X=X_full,beta=beta,prob=Prob))
}

###
### DEFINE DIMENSIONS OF ARTIFICIAL DATA
###
nunits = 1000 # number of units
cmax = 4 # maximum number of tasks per unit
# CREATE Vbetastar FOR MULTIVARIATE NORMAL DISTRIBUTION OF HETEROGENEITY (BETASTARS)
trueVbetastar = diag(c(0.4,0.2,0.4,2,4),nrow=5,ncol=5)
#correlation between first level & second level of quality attribute
trueVbetastar[1,2] = 0.1; trueVbetastar[2,1] = 0.1
#consumers prefer the high quality attribute are less price sensitive
trueVbetastar[2,3] = -0.15; trueVbetastar[3,2] = -0.15
#consumers prefer the "known" brand are less price sensitive
trueVbetastar[4,3] = -0.05; trueVbetastar[3,4] = -0.05
#consumers prefer the "unkown" brand are more price sensitive
trueVbetastar[5,3] = 0.05; trueVbetastar[3,5] = 0.05
#Check whether Cholesky-decomposition works for the given data generating process...
chol(trueVbetastar)
# DEFINE MEANS FOR MULTIVARIATE NORMAL DISTRIBUTION OF HETEROGENEITY (BETASTARS)
truebetabarstar = c(0.5,-0.5,0.8,2.5,2.5)
# DRAW BETASTARS FOR EACH UNIT
betastar = mvrnorm(n=nunits, truebetabarstar, trueVbetastar) #draw true betas from assumed distribution of heterogeneity
# TRANSFORM STARS TO BETA
beta = exp(as.matrix(betastar[,1]))
beta = cbind(beta,beta[,1]+exp(betastar[,2]))
beta = cbind(beta, -exp(betastar[,3]))
beta = cbind(beta,betastar[,4:5])

# CREATE MULTINOMIAL LOGIT y AND X FOR EACH UNIT USING BETA
lgtdata_sim=NULL
hdata=NULL
for (i in 1:nunits) {
  hdata[[i]] = simmnlv2(cmax,beta[i,])
}
# IDENTIFICATION 
for(i in 1:nunits){
  lgtdata_sim[[i]]=list(X=hdata[[i]]$X[,-1],y=hdata[[i]]$y) 
}

#View(lgtdata_sim[[1000]]$X)
#View(lgtdata_sim[[1000]]$y)

N = length(lgtdata_sim)
t = cmax
p = dim(lgtdata_sim[[1]]$X)[1]/cmax

# ESTIMATION DATA
E_Data=list(p=p,lgtdata=lgtdata_sim) #construct a list with all relevant data for estimation sample

# IDENTIFIED BETA AND BETASTAR (FOR COMPARISON LATER)
beta_id_DGP = beta[,-1]
beta_id_DGP[,1] = exp(betastar[,2])
beta_id_DGP[,2] = -exp(betastar[,3])
beta_id_DGP[,3] = betastar[,4] + beta[,1]
beta_id_DGP[,4] = betastar[,5] + beta[,1]

betastar_id_DGP = betastar[,-1]
betastar_id_DGP[,1] = betastar[,2]
betastar_id_DGP[,2] = betastar[,3]
betastar_id_DGP[,3] = betastar[,4] + beta[,1]
betastar_id_DGP[,4] = betastar[,5] + beta[,1]

###  POPULATION DISTRIBUTION --> See Table 2 in the paper
nunits = 100000
betastar_pop = mvrnorm(n=nunits, truebetabarstar, trueVbetastar) #draw true betas from assumed distribution of heterogeneity
beta_pop = exp(as.matrix(betastar_pop[,1]))
beta_pop = cbind(beta_pop,beta_pop[,1]+exp(betastar_pop[,2]))
beta_pop = cbind(beta_pop, -exp(betastar_pop[,3]))
beta_pop = cbind(beta_pop,betastar_pop[,4:5])
# Infer betastar_id_pop & beta_id_pop
#betastar
betastar_id_pop = betastar_pop[,-1]
betastar_id_pop[,1] = betastar_pop[,2]
betastar_id_pop[,2] = betastar_pop[,3]
betastar_id_pop[,3] = betastar_pop[,4] + beta_pop[,1]
betastar_id_pop[,4] = betastar_pop[,5] + beta_pop[,1]
#beta
beta_id_pop = beta_pop[,-1]
beta_id_pop[,1] = exp(betastar_pop[,2])
beta_id_pop[,2] = -exp(betastar_pop[,3])
beta_id_pop[,3] = betastar_pop[,4] + beta_pop[,1]
beta_id_pop[,4] = betastar_pop[,5] + beta_pop[,1]

# UPPER LEVEL DGP
truebetabarstar_id = apply(betastar_id_pop,2,mean) 
trueVbetastar_id = cov(betastar_id_pop)

trueGamma_id = chol2inv(chol(trueVbetastar_id[1:2,1:2]))%*%trueVbetastar_id[1:2,3:4]
trueSigma_id = trueVbetastar_id[3:4,3:4] - t(trueGamma_id) %*% trueVbetastar_id[1:2,1:2] %*% trueGamma_id
trueVstar_id = trueVbetastar_id[1:2,1:2]

###################################################
### ESTIMATION ANALYSIS BASED ON SIMULATED DATA ###
###################################################
N = dim(beta_id_DGP)[1]

### SOURCE
source('Functions/rhierMnlRwMixture_main_sim_correctHess.R')
Rcpp::sourceCpp('Functions/rhierMnlRwMixture_rcpp_loop_Sim_modHP.cpp',showOutput = FALSE)
Rcpp::sourceCpp('Functions/Speed++_Sim.cpp',showOutput = FALSE)

### PRIOR & MCMC 
nvar_c = 2
nu = 15 + nvar_c
Amu = diag(0.1, nrow = nvar_c, ncol = nvar_c)
Prior = list(ncomp=1,Amu=Amu,nu=nu,V=nu*diag(nvar_c)*0.5)
# MARGINAL-CONDITIONAL SAMPLER   
Mcmc = list(R=120000,keep=60,w=0.1)
out_marginal_cond = rhierMnlRwMixture_SR(Data=E_Data, Prior=Prior, Mcmc=Mcmc, nvar_c=nvar_c, flag="approx")

betastar_id_mc = out_marginal_cond$betadraw
compdraw_mc = out_marginal_cond$nmix$compdraw
probdraw_mc = out_marginal_cond$nmix$probdraw
rejection_mc = out_marginal_cond$rejection
  
#########################################################
### ILLUSTRATE THE NUMERICAL PROPERTIES OF POSTERIORS ###
#########################################################

### INDIVIDUAL LEVEL DRAWS
x = 1:1000
index_resp = sample(x, 12, replace = FALSE, prob = NULL)
index_resp

# CONVERGENCE FOR SELECTED INDIVIDUALS (Att 1)
# quartz() # for mac
windows() # for windows
par(mfrow=c(3,4))  # multiple plots are filled by rows!!!
plot(betastar_id_mc[index_resp[1],1,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[1],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[1],1], lwd=5, col = "red")
plot(betastar_id_mc[index_resp[2],1,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[2],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[2],1], lwd=5, col = "red")
plot(betastar_id_mc[index_resp[3],1,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[3],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[3],1], lwd=5, col = "red")
plot(betastar_id_mc[index_resp[4],1,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[4],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[4],1], lwd=5, col = "red")
plot(betastar_id_mc[index_resp[5],1,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[5],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[5],1], lwd=5, col = "red")
plot(betastar_id_mc[index_resp[6],1,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[6],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[6],1], lwd=5, col = "red")
plot(betastar_id_mc[index_resp[7],1,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[7],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[7],1], lwd=5, col = "red")
plot(betastar_id_mc[index_resp[8],1,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[8],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[8],1], lwd=5, col = "red")
plot(betastar_id_mc[index_resp[9],1,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[9],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[9],1], lwd=5, col = "red")
plot(betastar_id_mc[index_resp[10],1,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[10],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[10],1], lwd=5, col = "red")
plot(betastar_id_mc[index_resp[11],1,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[11],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[11],1], lwd=5, col = "red")
plot(betastar_id_mc[index_resp[12],1,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[12],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[12],1], lwd=5, col = "red")


# CONVERGENCE FOR SELECTED INDIVIDUALS (Att 3)
# quartz() # for mac
windows() # for windows
par(mfrow=c(3,4))  # multiple plots are filled by rows!!!
plot(betastar_id_mc[index_resp[1],3,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[1],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[1],3], lwd=5, col = "red")
plot(betastar_id_mc[index_resp[2],3,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[2],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[2],3], lwd=5, col = "red")
plot(betastar_id_mc[index_resp[3],3,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[3],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[3],3], lwd=5, col = "red")
plot(betastar_id_mc[index_resp[4],3,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[4],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[4],3], lwd=5, col = "red")
plot(betastar_id_mc[index_resp[5],3,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[5],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[5],3], lwd=5, col = "red")
plot(betastar_id_mc[index_resp[6],3,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[6],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[6],3], lwd=5, col = "red")
plot(betastar_id_mc[index_resp[7],3,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[7],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[7],3], lwd=5, col = "red")
plot(betastar_id_mc[index_resp[8],3,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[8],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[8],3], lwd=5, col = "red")
plot(betastar_id_mc[index_resp[9],3,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[9],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[9],3], lwd=5, col = "red")
plot(betastar_id_mc[index_resp[10],3,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[10],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[10],3], lwd=5, col = "red")
plot(betastar_id_mc[index_resp[11],3,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[11],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[11],3], lwd=5, col = "red")
plot(betastar_id_mc[index_resp[12],3,],xlab="R",ylab="Betastar",lwd = 1 ,main=index_resp[12],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = betastar_id_DGP[index_resp[12],3], lwd=5, col = "red")


### RECOVER UPPER LEVEL MODEL COEFFICIENTS
Gamma_id_draw = NULL
Sigma_id_draw = NULL
Vbetastar_id_draw = NULL
Vstar_id_draw = NULL
betabarstar_id_draw = NULL

nvar = 4
nvar_uc = nvar - nvar_c
R = dim(betastar_id_mc)[3]

for(i in 1:R){
  Vbetastar_id_draw[[i]] = chol2inv(chol((compdraw_mc[[i]][[1]][[2]] %*% t(compdraw_mc[[i]][[1]][[2]]))))
  Vstar_id_draw[[i]] = Vbetastar_id_draw[[i]][1:nvar_c,1:nvar_c] 
  Gamma_id_draw[[i]] = chol2inv(chol(Vbetastar_id_draw[[i]][1:nvar_c,1:nvar_c]))%*%Vbetastar_id_draw[[i]][1:nvar_c,(nvar_c+1):nvar]
  Sigma_id_draw[[i]] = Vbetastar_id_draw[[i]][(nvar_c+1):nvar,(nvar_c+1):nvar] - t(Gamma_id_draw[[i]]) %*% Vbetastar_id_draw[[i]][1:nvar_c,1:nvar_c] %*% Gamma_id_draw[[i]]
  betabarstar_id_draw[[i]] = compdraw_mc[[i]][[1]][[1]]
}

Gamma_id_draw_cube <- array(0,dim=c(nvar_c,nvar_uc,R))
Sigma_id_cube <- array(0,dim=c(nvar_uc,nvar_uc,R))
betabarstar_id_draw_mat <- array(0,dim=c(nvar,R))
Vstar_id_draw_cube <- array(0,dim=c(nvar_c,nvar_c,R))
Vbetastar_id_draw_cube <- array(0,dim=c(nvar,nvar,R))

for(i in 1:R){
  Gamma_id_draw_cube[,,i] = Gamma_id_draw[[i]]
  Sigma_id_cube[,,i] = Sigma_id_draw[[i]]
  betabarstar_id_draw_mat[,i] = betabarstar_id_draw[[i]]
  Vstar_id_draw_cube[,,i] = Vstar_id_draw[[i]]
  Vbetastar_id_draw_cube[,,i] = Vbetastar_id_draw[[i]]
}

### betabarstar_id 
# quartz() # for mac
windows() # for windows
par(mfrow=c(2,2))  # multiple plots are filled by rows!!!
plot(betabarstar_id_draw_mat[1,],xlab="R",ylab="betabarstar_1",type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h=truebetabarstar_id[1],col="red",lwd=5)
plot(betabarstar_id_draw_mat[2,],xlab="R",ylab="betabarstar_2",type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h=truebetabarstar_id[2],col="red",lwd=5)
plot(betabarstar_id_draw_mat[3,],xlab="R",ylab="betabarstar_3",type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h=truebetabarstar_id[3],col="red",lwd=5)
plot(betabarstar_id_draw_mat[4,],xlab="R",ylab="betabarstar_4",type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h=truebetabarstar_id[4],col="red",lwd=5)

###Mains Vbetastar_id 
# quartz() # for mac
windows() # for windows
par(mfrow=c(2,2))  # multiple plots are filled by rows!!!
plot(Vbetastar_id_draw_cube[1,1,],xlab="R",ylab="Vbetastar_11",type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h=trueVbetastar_id[1,1],col="red",lwd=5)
plot(Vbetastar_id_draw_cube[2,2,],xlab="R",ylab="Vbetastar_22",type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h=trueVbetastar_id[2,2],col="red",lwd=5)
plot(Vbetastar_id_draw_cube[3,3,],xlab="R",ylab="Vbetastar_33",type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h=trueVbetastar_id[3,3],col="red",lwd=5)
plot(Vbetastar_id_draw_cube[3,3,],xlab="R",ylab="Vbetastar_44",type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h=trueVbetastar_id[4,4],col="red",lwd=5)


### Check Gamma-draws
# quartz() # for mac
windows() # for windows
par(mfrow=c(2,2))  # multiple plots are filled by rows!!!
plot(Gamma_id_draw_cube[1,1,],xlab="R",ylab="Gamma_11",type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
plot(Gamma_id_draw_cube[2,1,],xlab="R",ylab="Gamma_21",type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
plot(Gamma_id_draw_cube[1,2,],xlab="R",ylab="Gamma_12",type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
plot(Gamma_id_draw_cube[2,2,],xlab="R",ylab="Gamma_22",type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()

