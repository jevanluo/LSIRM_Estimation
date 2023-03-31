rm(list=ls())  
setwd("Psych paper")

library(difR)

data(verbal)
IRdata<-verbal[,-c(25,26)]

##################JAGS
library(rjags)

set.seed(1)

N<-nrow(IRdata)
I<- ncol(IRdata)

X1 <- IRdata

mymodel <- "modelfile.bug"

myinits <- function(){list(beta = rnorm(I, 0, 1),
                           theta = rnorm(N, 0, 1), 
                           z = matrix(cbind(rnorm(N, 0, 1), rnorm(N, 0, 1)), ncol=2),#revised the format, but the simulate result is the same
                           w = matrix(cbind(rnorm(I, 0, 1), rnorm(I, 0, 1)), ncol=2),
                           sigma1 = 1/3, 
                           gamma = 3) 
}

myparams= c("beta", "theta", "sigma","gamma", "z", "w", "log.ll")

data <- data.frame(r = X1)
mydata <- list(r=data, I=I, N=N)

st0 <- Sys.time()

par_jagssamples<-jags.parallel(data=mydata, inits=myinits, parameters=myparams, model.file=mymodel, n.chains=2, 
     n.iter=90000,n.burnin=40000, n.thin=50, n.cluster=2)

st1 <- Sys.time()
st <- st1-st0
##################STAN
library(rstan)

mydata<-list(N=nrow(IRdata),
             I=ncol(IRdata ),
             Y=IRdata,
             mu=c(0,0),
             Sigma=diag(2),
             tau_beta=2,
             mu_gamma=0.5,
             tau_gamma=1,
             a_sigma=1,
             b_sigma=1)

myinits <- list(chain1=list(beta = rnorm(mydata$I, 0, 1),
                            theta = rnorm(mydata$N, 0, 1),
                            z = matrix(cbind(rnorm(mydata$N, 0, 1), rnorm(mydata$N, 0, 1)), ncol=2), #revised the format, but the simulate result is the same
                            w = matrix(cbind(rnorm(mydata$I, 0, 1), rnorm(mydata$I, 0, 1)), ncol=2),
                            sigma2 = 3,
                            gamma = 3),
                chain2=list(beta = rnorm(mydata$I, 0, 1),
                            theta = rnorm(mydata$N, 0, 1),
                            z = matrix(cbind(rnorm(mydata$N, 0, 1), rnorm(mydata$N, 0, 1)), ncol=2), #revised the format, but the simulate result is the same
                            w = matrix(cbind(rnorm(mydata$I, 0, 1), rnorm(mydata$I, 0, 1)), ncol=2),
                            sigma2 = 3,
                            gamma = 3))

myparams= c("beta", "theta", "sigma2","gamma", "z", "w", "log_ll")

st0 <- Sys.time()
par_stansamples <- stan(
              file = "modelfile.stan", 
              data = mydata,       # named list of data
              init =myinits, #list that contains initial values
              pars = myparams, # named list of parameters of interest
              chains = 2,             # number of Markov chains
              warmup = 10000,          #burn-in period per chain 
              iter = 15000,           # total number of iterations per chain
              cores = 2,              # number of cores 
              control = list(max_treedepth = 10), #For algorithm NUTS
              thin=5
              )
st1 <- Sys.time()
st <- st1-st0

##################NIMBLE 
library(parallel)
this_cluster <- makeCluster(2)

set.seed(1)
st0 <- Sys.time()
#create a function that includes all the modeling steps and run that function in parallel
run_MCMC_allcode <- function(seed, data) {
  library(nimble)
  
  mycode <- nimbleCode({
    for (i in 1:I) {
      for (j in 1:N) {
        r[j,i] ~ dbern(prob=p[j,i])
        
        logit(p[j,i]) <- theta[j] + beta[i] - gamma * sqrt((z[j,1] - w[i,1])^2 + (z[j,2] - w[i,2])^2)
        
      }
      beta[i] ~ dnorm(0, sd = tau_beta)
      w[i,] ~ dmnorm(pos_mu[1:2], pos_covmat[1:2,1:2])
    }
    for (j in 1:N) {
      theta[j] ~ dnorm(0, var = sigma) 
      z[j,] ~ dmnorm(pos_mu[1:2], pos_covmat[1:2,1:2])
      
    }
    
    
    log(gamma) ~ dlnorm(mu_gamma, sd = tau_gamma) 
    sigma ~ dinvgamma(shape=a_sigma,scale = b_sigma)
    
    for (i in 1:I) {
      for (j in 1:N) {
        log.ll[j,i] <- r[j,i] * log(p[j,i])  + (1-r[j,i])*log(1-p[j,i])
      }
    }
  })
  
  N<-nrow(data)
  I<- ncol(data)
  
  myconst <- list(I = ncol(data), 
                  N = nrow(data),
                  pos_mu = rep(0,2),
                  pos_covmat = diag(2),
                  tau_beta = 2,
                  mu_gamma = 0.5,
                  tau_gamma = 1,
                  a_sigma = 1,
                  b_sigma = 1)
  mydata <- list(r = data)
  
  set.seed(1)
  
  myinits <- list(beta    = rnorm(myconst$I, 0, 1),
                  theta   = rnorm(myconst$N, 0, 1), 
                  z    = matrix(cbind(rnorm(myconst$N, 0, 1), rnorm(myconst$N, 0, 1)), ncol=2), 
                  w    = matrix(cbind(rnorm(myconst$I, 0, 1), rnorm(myconst$I, 0, 1)), ncol=2), 
                  sigma = 3, 
                  gamma = 3,
                  log_gamma = log(3)) 
  
  myparams = c("beta", "theta", "sigma","gamma", "z", "w","log.ll")
  
  mymodel <- nimbleModel(code = mycode,
                         data = mydata,
                         constants = myconst,
                         inits = myinits)
  
  cmymodel <- compileNimble(mymodel)
  
  confMCMC <- configureMCMC(mymodel, monitors = myparams)
  
  myMCMC <- buildMCMC(cmymodel)
  cmyMCMC <- compileNimble(myMCMC, project = mymodel)
  
  results <- runMCMC(cmyMCMC, niter = 90000, nburnin = 40000, thin = 50, 
                     nchains = 1, 
                     samplesAsCodaMCMC = TRUE, 
                     setSeed = seed)
  
  return(results)
}

#execute the desired code using parLapply, with each process running in parallel.
par_nimblesamples <- parLapply(cl = this_cluster, X = 1:2, 
                          fun = run_MCMC_allcode, 
                          data = IRdata)

st1 <- Sys.time()
st <- st1-st0
#close the cluster when you're done with it.
stopCluster(this_cluster)

#We ran four independent MCMCs

#you can see each chain here
par_nimblesamples[[1]]
par_nimblesamples[[2]]