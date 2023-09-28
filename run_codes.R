##################################################
## Project: Estimating LSIRM with JAGS, Stan and NIMBLE
## Author: Jinwen Luo, Ludovica De Carolis, Biao Zeng, Minjeong Jeon
##################################################
rm(list=ls())  
setwd("Psych_paper")

##################################################
## Section 0: Load Data Example
##################################################
library(difR)
library(prolsirm) 
source("vislsirm.R")
library(coda)
library(ggplot2)


data(verbal)
IRdata<-verbal[,-c(25,26)]
N<-nrow(IRdata)
I<- ncol(IRdata)
r <- IRdata

##################################################
## Section 1A: Using JAGS to Estimate LSIRM
## Please stick to only one software at a time,
## the same procedure will be shown with Stan and Nimble as well
##################################################
##################JAGS
#load the required package
library(rjags)

#model definition in BUGS 
mymodel <- "modelfile.bug"

#data specification
mydata <- list(N=N,
               r=r, 
               I=I)

#initial values setting
myinits <- function(){list(beta = rnorm(I, 0, 1),
                           theta = rnorm(N, 0, 1), 
                           z = matrix(cbind(rnorm(N, 0, 1), rnorm(N, 0, 1)), ncol=2),
                           w = matrix(cbind(rnorm(I, 0, 1), rnorm(I, 0, 1)), ncol=2),
                           sigma1 = 1/3, 
                           gamma = 3)}

#parameters to save
myparams= c("beta", "theta", "sigma2","gamma", "z", "w", "log_ll")

#seed for reproducibility
set.seed(1)

#compilation
compiled_LSIRM <- jags.model(file = mymodel, 
                             data = mydata, 
                             inits = myinits, 
                             n.chains = 2)
#sampling
fitted_LSIRM <- coda.samples(model = compiled_LSIRM, 
                             inits=myinits, 
                             variable.names = myparams,
                             n.iter =90000, 
                             n.burnin=40000, 
                             thin =50, 
                             n.chains = 2)

#wrap+parallel
fitted_LSIRM_p <-jags(data=mydata, 
                    inits=myinits, 
                    parameters=myparams, 
                    model.file=mymodel, 
                    n.chains=2, 
                    n.iter=90000,
                    n.burnin=40000, 
                    n.thin=50, 
                    parallel=T, 
                    n.cores=2)

##################################################
## Section 1B: Using Stan to Estimate LSIRM
##################################################
##################STAN
#load the required package
library(rstan)

#model specification in Stan putting into an R object as a text 
mymodel <- "
data{
    int <lower=0> N; //respondents
    int <lower=0> I; //items
    int <lower=0,upper=1> r[N,I]; 
    
    real <lower=0> tau_beta; 
    real mu_gamma;
    real <lower=0> tau_gamma; 
    real a_sigma;
    real b_sigma;
    
    vector[2] mu;           
    cov_matrix[2] Sigma;  
    }
    
    parameters {
    vector[N] theta; 
    vector[I] beta; 
    
    real gamma;
    
    real <lower=0> sigma2; 
    
    matrix [N, 2] z;
    matrix [I, 2] w;
    }
    
    
    model{
    theta ~ normal(0, sqrt(sigma2)); 
    beta ~ normal(0, tau_beta);
    gamma ~ lognormal(mu_gamma,tau_gamma); 
    sigma2 ~ inv_gamma(a_sigma, b_sigma);   
    
    
    for (j in 1: N){
    z[j] ~ multi_normal(mu,Sigma); 
    }
    for (i in 1: I){
    w[i] ~ multi_normal(mu,Sigma);
    }
    
    for (j in 1: N){
    for (i in 1: I){
    r[j, i]~ bernoulli_logit(theta[j]+beta[i]-gamma*sqrt((z[j,1] - w[i,1])^2 + (z[j,2] - w[i,2])^2));
    }}
    }
    
    generated quantities {
    vector[I] log_ll[N];
    for (j in 1: N){
    for (i in 1: I){
    log_ll[j, i] = bernoulli_logit_lpmf(r[j, i]|theta[j]+beta[i]-gamma*sqrt((z[j,1] - w[i,1])^2 + (z[j,2] - w[i,2])^2));
    }}
}
"
#data specification
mydata<-list(N=N,
             I=I,
             r=r,
             mu=c(0,0),
             Sigma=diag(2),
             tau_beta=2,
             mu_gamma=0.5,
             tau_gamma=1,
             a_sigma=1,
             b_sigma=1)

#initial values setting
myinits <- function(){list(beta = rnorm(I, 0, 1),
                           theta = rnorm(N, 0, 1), 
                           z = matrix(cbind(rnorm(N, 0, 1), rnorm(N, 0, 1)), ncol=2),
                           w = matrix(cbind(rnorm(I, 0, 1), rnorm(I, 0, 1)), ncol=2),
                           sigma2 = 3,
                           gamma = 3)}

#parameters to save
myparams= c("beta", "theta", "sigma2","gamma", "z", "w", "log_ll")

#seed for reproducibility
set.seed(1)

#compilation
compiled_LSIRM <- stan_model(model_code = mymodel) 

#sampling
fitted_LSIRM <-sampling(compiled_LSIRM, 
                       data = mydata,
                       pars = myparams,       # named list of parameters of interest
                       chains = 2,            # number of Markov chains
                       warmup = 10000,        #burn-in period per chain 
                       iter = 15000,
                       init = myinits,         # total number of iterations per chain
                       cores = 2,              # number of cores (using 2 just for the vignette)
                       thin=5,
                       save_warmup=F)
#FOR VISUALIZATION:
fitted_LSIRM <-As.mcmc.list(fitted_LSIRM)
  
#wrap+parallel
fitted_LSIRM_p <- stan(file = "modelfile.stan", #file with the model specification in Stan 
                      data = mydata,      #named list of data
                      init =myinits,      #list that contains initial values
                      pars = myparams,    #named list of parameters of interest
                      chains = 2,        #number of Markov chains
                      warmup = 10000,    #burn-in period per chain 
                      iter = 15000,       #total number of iterations per chain
                      thin=5,            #thinning interval
                      cores = 2,         #number of cores 
                      seed = 1,          #seed for reproducibility
                      control = list(max_treedepth = 10)) #parameter for algorithm NUTS
                      
##################################################
## Section 1C: Using NIMBLE to Estimate LSIRM
##################################################
##################NIMBLE 
#load the required package
library(nimble)

#model definition in NIMBLE
mymodel <- nimbleCode({
  for (i in 1:I) {
    for (j in 1:N) {
      r[j,i] ~ dbern(prob=p[j,i])
      
      logit(p[j,i]) <- theta[j] + beta[i] - gamma * sqrt((z[j,1] - w[i,1])^2 + (z[j,2] - w[i,2])^2)
      
    }
    beta[i] ~ dnorm(0, sd = tau_beta)
    w[i,] ~ dmnorm(pos_mu[1:2], pos_covmat[1:2,1:2])
  }
  for (j in 1:N) {
    theta[j] ~ dnorm(0, var = sigma2) 
    z[j,] ~ dmnorm(pos_mu[1:2], pos_covmat[1:2,1:2])
    
  }
  
  log(gamma) ~ dlnorm(mu_gamma, sd = tau_gamma) 
  sigma2 ~ dinvgamma(shape=a_sigma,scale = b_sigma)
  
  for (i in 1:I) {
    for (j in 1:N) {
      log_ll[j,i] <- r[j,i] * log(p[j,i])  + (1-r[j,i])*log(1-p[j,i])
    }
  }
})

#data specification
mydata <- list(r=r)

#constants specification
myconst <- list(N=N,
                I=I,
                pos_mu = rep(0,2),
                pos_covmat = diag(2),
                tau_beta = 2,
                mu_gamma = 0.5,
                tau_gamma = 1,
                a_sigma = 1,
                b_sigma = 1)

#initial values setting
myinits <- list(beta = rnorm(I, 0, 1),
               theta = rnorm(N, 0, 1), 
               z = matrix(cbind(rnorm(N, 0, 1), rnorm(N, 0, 1)), ncol=2),
               w = matrix(cbind(rnorm(I, 0, 1), rnorm(I, 0, 1)), ncol=2),
               sigma2 = 3,
               gamma = 3,
               log_gamma = log(3))

#seed for reproducibility
myparams = c("beta", "theta", "sigma2","gamma", "z", "w","log_ll")

#seed for reproducibility
set.seed(1)

#compilation
model <- nimbleModel(code = mymodel,
                     data = mydata,
                     constants = myconst,
                     inits = myinits)

cmymodel <- compileNimble(model)
confMCMC <- configureMCMC(model, 
                          monitors = myparams)

myMCMC <- buildMCMC(confMCMC)
cmyMCMC <- compileNimble(myMCMC, project = model)

#sampling
fitted_LSIRM <- runMCMC(cmyMCMC, 
                       niter = 90000, 
                       nburnin = 40000, 
                       thin = 50, 
                       nchains = 2, 
                       setSeed = 1,
                       samplesAsCodaMCMC = T)

#wrapper
fitted_LSIRM_w <- nimbleMCMC(code = mymodel,
                          data = mydata,
                          constants = myconst, 
                          inits = myinits,
                          monitors = myparams,
                          niter = 90000,
                          nburnin = 40000,
                          thin = 50, 
                          nchains = 2, 
                          progressBar = TRUE,
                          samplesAsCodaMCMC = T)

##################################################
## Section 1D: Paralleling NIMBLE to Estimate LSIRM
##################################################
#parallel
library(difR)

data(verbal)
IRdata<-verbal[,-c(25,26)]

library(parallel)
this_cluster <- makeCluster(2)

#create a function that includes all the modeling steps and run that function in parallel
run_MCMC_allcode <- function(seed, data) {
  library(nimble)
  
  mymodel <- nimbleCode({
    for (i in 1:I) {
      for (j in 1:N) {
        r[j,i] ~ dbern(prob=p[j,i])
        
        logit(p[j,i]) <- theta[j] + beta[i] - gamma * sqrt((z[j,1] - w[i,1])^2 + (z[j,2] - w[i,2])^2)
        
      }
      beta[i] ~ dnorm(0, sd = tau_beta)
      w[i,] ~ dmnorm(pos_mu[1:2], pos_covmat[1:2,1:2])
    }
    for (j in 1:N) {
      theta[j] ~ dnorm(0, var = sigma2) 
      z[j,] ~ dmnorm(pos_mu[1:2], pos_covmat[1:2,1:2])
      
    }
    
    
    log(gamma) ~ dlnorm(mu_gamma, sd = tau_gamma) 
    sigma2 ~ dinvgamma(shape=a_sigma,scale = b_sigma)
    
    for (i in 1:I) {
      for (j in 1:N) {
        log_ll[j,i] <- r[j,i] * log(p[j,i])  + (1-r[j,i])*log(1-p[j,i])
      }
    }
  })
  
  N <-nrow(data)
  I <- ncol(data)
  r <- data
  
  myconst <- list(I = I, 
                  N = N,
                  pos_mu = rep(0,2),
                  pos_covmat = diag(2),
                  tau_beta = 2,
                  mu_gamma = 0.5,
                  tau_gamma = 1,
                  a_sigma = 1,
                  b_sigma = 1)
  
  mydata <- list(r = r)
  
  myinits <- list(beta    = rnorm(I, 0, 1),
                  theta   = rnorm(N, 0, 1), 
                  z    = matrix(cbind(rnorm(N, 0, 1), rnorm(N, 0, 1)), ncol=2), 
                  w    = matrix(cbind(rnorm(I, 0, 1), rnorm(I, 0, 1)), ncol=2), 
                  sigma2 = 3, 
                  gamma = 3,
                  log_gamma = log(3)) 
  
  myparams = c("beta", "theta", "sigma2","gamma", "z", "w","log_ll")
  
  set.seed(1)
  
  model <- nimbleModel(code = mymodel,
                       data = mydata,
                       constants = myconst,
                       inits = myinits)
  
  cmymodel <- compileNimble(model)
  confMCMC <- configureMCMC(model, 
                            monitors = myparams)
  
  myMCMC <- buildMCMC(confMCMC)
  cmyMCMC <- compileNimble(myMCMC, project = model)
  
  results <- runMCMC(cmyMCMC, niter = 90000, nburnin = 40000, thin = 50, 
                     samplesAsCodaMCMC = TRUE,
                     nchains = 1, 
                     setSeed = 1)
  
  return(results)
}

#execute the desired code using parLapply, with each process running in parallel.
fitted_LSIRM_p <- parLapply(cl = this_cluster, X = 1:2, 
                          fun = run_MCMC_allcode, 
                          data = IRdata)

#close the cluster when you're done with it.
stopCluster(this_cluster)


##################################################
## Section 2: Diagnostics
## Here we assume the fitted_LSIRM object has been obtained from 
## one of the above methods.
##################################################
############DIAGNOSTICS
library(ggmcmc)

postsample<-ggs(fitted_LSIRM)
ggs_density(postsample, family="gamma", greek = T)+theme(axis.ticks.y =element_blank(),
                                                  axis.text.y =element_blank(),
                                                  axis.title.y=element_blank(),
                                                  axis.title.x = element_blank(),
                                                  axis.text=element_text(size=20,angle = 90,face="bold"))

ggs_traceplot(postsample, family="gamma", original_thin = F, original_burnin =F, greek = T)+theme(axis.ticks.y =element_blank(),
                                                                                                  axis.text.y =element_blank(),
                                                                                                  axis.title.y=element_blank(),
                                                                                                  axis.title.x = element_blank(),
                                                                                                  axis.text=element_text(size=20,angle = 90,face="bold"))                                                                                           

ggs_autocorrelation(postsample, family="gamma", greek = T)+theme(axis.ticks.y =element_blank(),
                                                         axis.text.y =element_blank(),
                                                         axis.title.y=element_blank(),
                                                         axis.title.x = element_blank(),
                                                         axis.text=element_text(size=20,angle = 90,face="bold"))

##################################################
## Section 3: Visualizing Latent Space
## Here we assume the fitted_LSIRM object has been obtained from 
## one of the above methods.
##################################################
library(prolsirm) 
source("vislsirm.R")

fitted_samples <- as.matrix(fitted_LSIRM[[1]]) ### !!!! add in revision
llCols <- grep("log.ll", colnames(fitted_samples))
zCols <- grep("z", colnames(fitted_samples))
wCols <- grep("w", colnames(fitted_samples))
log_ll.samples <- fitted_samples[,llCols]
ll.samples <- fitted_samples[,c(llCols)]
zz.samples <- fitted_samples[,c(zCols)]
ww.samples <- fitted_samples[,c(wCols)]

M.keep = 1000
zz <- vector("list", M.keep)
for (i in 1:M.keep) {
        ### check the ordering of the parameters
        ### stan uses a different ordering so change the following code if stan is used
        
        ### for jags and nimble
        zz[[i]] <- cbind(zz.samples[i,1:N], zz.samples[i,(N+1):(2*N)])
        ### for stan
        #zz[[i]] <- cbind(zz0[i,2*(1:nz-1)+1], zz0[i,2*(1:nz)])
}

ww <- vector("list", M.keep)
for (i in 1:M.keep) {
        ### check the ordering of the parameters
        ### stan uses a different ordering so change the following code if stan is used
        
        ### for jags and nimble
        ww[[i]] <- cbind(zz.samples[i,1:I], zz.samples[i,(I+1):(2*I)])
        ### for stan
        #ww[[i]] <- cbind(ww0[i,2*(1:nw-1)+1], ww0[i,2*(1:nw)])
}

### use procrustes from prolsirm package to do the matching
matched <- procrustes(z=zz,w=ww,ref=1)
zz.matched <- matched$zz
ww.matched <- matched$ww

ItemLabel <- c ( rep ( c ( "Curse" ,"Scold" ,"Shout" ) ,8) )
plot_latent ( zz.matched , ww.matched , ItemGroup = ItemLabel ) 

