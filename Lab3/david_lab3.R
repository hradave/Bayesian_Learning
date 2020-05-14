#################### QUESTION 1 ####################

RNGversion('3.5.1')

library(mvtnorm)

#### a
data = read.table('rainfall.dat', col.names = 'rainfall')

#### (i)

n = length(data$rainfall)
avg = mean(data$rainfall)

#prior setup
mu_0 = 30 #mu prior from prior knowledge, or looking at the data if there is no knowledge
tau2_0 = 100 #std(mu_0) high, because we have no prior knowledge
sigma2_0 = var(data$rainfall) # should use var(data$rainfall for best guess)
nu_0 = 4 #df, small, because we have no prior knowledge

# scaled-inverse chi-square simulator (from lab1 and from NormalMixtureGibbs.R)
rScaledInvChi2 <- function(n, df, scale){
  return((df*scale)/rchisq(n,df=df))
}


# Gibbs sampling
steps = 1000
gibbsDraws <- matrix(0,steps+1,2) #1st column: mu, 2nd column: sigma2
gibbsDraws[1,1] = mu_0
gibbsDraws[1,2] = sigma2_0

for (i in 1:steps) {
  #simulate mu from conditional posterior (normal)
  w = (n/gibbsDraws[i,2])/(n/gibbsDraws[i,2] + 1/tau2_0)
  mu_n = w * avg + (1-w) * mu_0
  tau2_n = 1 / (n/gibbsDraws[i,2] + 1/tau2_0)
  
  #update mu
  gibbsDraws[i+1,1] = rnorm(1, mu_n, sqrt(tau2_n))
  
  #simulate sigma2 from conditional posterior (inv-chisq)
  nu_n = nu_0 + n
  scale = (nu_0 * sigma2_0 + sum((data$rainfall - gibbsDraws[i+1,1])^2)) / (n + nu_0)
  
  #update sigma2
  gibbsDraws[i+1,2] = rScaledInvChi2(1, nu_n, scale)
}

#### (ii)

burn_in_mu = 10
burn_in_sigma2 = 10
hist(gibbsDraws[burn_in_mu:length(gibbsDraws[,1]),1], breaks = 100)
hist(gibbsDraws[burn_in_sigma2:length(gibbsDraws[,2]),2], breaks = 100)

#Plot Markov chains and cumulative mean plots
cum_mean = numeric(steps + 1)
cum_var = numeric(steps + 1)
for (i in 1:(steps+1)) {
  cum_mean[i] = mean(gibbsDraws[1:i,1])
  cum_var[i] = mean(gibbsDraws[1:i,2])
}
plot(1:(steps+1), gibbsDraws[,1], type = 'l')
lines(1:(steps+1), cum_mean, col = 'blue', lwd = 5)
plot(1:(steps+1), gibbsDraws[,2], type = 'l')
lines(1:(steps+1), cum_var, col = 'blue', lwd = 5)

# result: mean of all posterior draws
mu_gibbs = mean(gibbsDraws[burn_in_mu:length(gibbsDraws[,1]),1])
sigma2_gibbs = mean(gibbsDraws[burn_in_sigma2:length(gibbsDraws[,2]),2])

#### b

#mixture model (normals)
# Code from NormalMixtureGibbs.R (by Mattias Villani)

# Estimating a simple mixture of normals
# Author: Mattias Villani, IDA, Linkoping University. http://mattiasvillani.com

##########    BEGIN USER INPUT #################
# Data options
x <- as.matrix(data$rainfall)

# Model options
nComp <- 2    # Number of mixture components

# Prior options
alpha <- 10*rep(1,nComp) # Dirichlet(alpha)
muPrior <- rep(30,nComp) # Prior mean of mu
tau2Prior <- rep(100,nComp) # Prior std of mu: high, because we have no prior knowledge
sigma2_0 <- rep(var(x),nComp) # s20 (best guess of sigma2)
nu0 <- rep(4,nComp) # degrees of freedom for prior on sigma2: small, because we have no prior knowledge

# MCMC options
nIter <- 100 # Number of Gibbs sampling draws

# Plotting options
plotFit <- TRUE
lineColors <- c("blue", "green", "magenta", 'yellow')
sleepTime <- 0.1 # Adding sleep time between iterations for plotting
################   END USER INPUT ###############

###### Defining a function that simulates from the 
rScaledInvChi2 <- function(n, df, scale){
  return((df*scale)/rchisq(n,df=df))
}

####### Defining a function that simulates from a Dirichlet distribution
rDirichlet <- function(param){
  nCat <- length(param)
  piDraws <- matrix(NA,nCat,1)
  for (j in 1:nCat){
    piDraws[j] <- rgamma(1,param[j],1)
  }
  piDraws = piDraws/sum(piDraws) # Diving every column of piDraws by the sum of the elements in that column.
  return(piDraws)
}

# Simple function that converts between two different representations of the mixture allocation
S2alloc <- function(S){
  n <- dim(S)[1]
  alloc <- rep(0,n)
  for (i in 1:n){
    alloc[i] <- which(S[i,] == 1)
  }
  return(alloc)
}
# Seed
set.seed(1234567890)

# Initial value for the MCMC
nObs <- length(x)
S <- t(rmultinom(nObs, size = 1 , prob = rep(1/nComp,nComp))) # nObs-by-nComp matrix with component allocations.
mu <- quantile(x, probs = seq(0,1,length = nComp))
sigma2 <- rep(var(x),nComp)
probObsInComp <- rep(NA, nComp)

# Setting up the plot
xGrid <- seq(min(x)-1*apply(x,2,sd),max(x)+1*apply(x,2,sd),length = 100)
xGridMin <- min(xGrid)
xGridMax <- max(xGrid)
mixDensMean <- rep(0,length(xGrid))
effIterCount <- 0
ylim <- c(0,2*max(hist(x)$density))


for (k in 1:nIter){
  message(paste('Iteration number:',k))
  alloc <- S2alloc(S) # Just a function that converts between different representations of the group allocations
  nAlloc <- colSums(S)
  print(nAlloc)
  # Update components probabilities
  pi <- rDirichlet(alpha + nAlloc)
  
  # Update mu's
  for (j in 1:nComp){
    precPrior <- 1/tau2Prior[j]
    precData <- nAlloc[j]/sigma2[j]
    precPost <- precPrior + precData
    wPrior <- precPrior/precPost
    muPost <- wPrior*muPrior + (1-wPrior)*mean(x[alloc == j])
    tau2Post <- 1/precPost
    mu[j] <- rnorm(1, mean = muPost, sd = sqrt(tau2Post))
  }
  
  # Update sigma2's
  for (j in 1:nComp){
    sigma2[j] <- rScaledInvChi2(1, df = nu0[j] + nAlloc[j], scale = (nu0[j]*sigma2_0[j] + sum((x[alloc == j] - mu[j])^2))/(nu0[j] + nAlloc[j]))
  }
  
  # Update allocation
  for (i in 1:nObs){
    for (j in 1:nComp){
      probObsInComp[j] <- pi[j]*dnorm(x[i], mean = mu[j], sd = sqrt(sigma2[j]))
    }
    S[i,] <- t(rmultinom(1, size = 1 , prob = probObsInComp/sum(probObsInComp)))
  }
  
  # Printing the fitted density against data histogram
  if (plotFit && (k%%1 ==0)){
    effIterCount <- effIterCount + 1
    hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = paste("Iteration number",k), ylim = ylim)
    mixDens <- rep(0,length(xGrid))
    components <- c()
    for (j in 1:nComp){
      compDens <- dnorm(xGrid,mu[j],sd = sqrt(sigma2[j]))
      mixDens <- mixDens + pi[j]*compDens
      lines(xGrid, compDens, type = "l", lwd = 2, col = lineColors[j])
      components[j] <- paste("Component ",j)
    }
    mixDensMean <- ((effIterCount-1)*mixDensMean + mixDens)/effIterCount
    
    lines(xGrid, mixDens, type = "l", lty = 2, lwd = 3, col = 'red')
    legend("topright", box.lty = 1, legend = c("Data histogram",components, 'Mixture'), 
           col = c("black",lineColors[1:nComp], 'red'), lwd = 2)
    Sys.sleep(sleepTime)
  }
  
}

hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = "Final fitted density")
lines(xGrid, mixDensMean, type = "l", lwd = 2, lty = 4, col = "red")
lines(xGrid, dnorm(xGrid, mean = mean(x), sd = apply(x,2,sd)), type = "l", lwd = 2, col = "blue")
legend("topright", box.lty = 1, legend = c("Data histogram","Mixture density","Normal density"), col=c("black","red","blue"), lwd = 2)

#########################    End of code from Mattias Villani    ##############################################




#### c

#histogram of data
hist(data$rainfall, freq = F, breaks = 30)

#normal density from a)
norm_gibbs = rnorm(1000, mean = mu_gibbs, sd = sqrt(sigma2_gibbs))
lines(density(norm_gibbs), col = 'blue', lwd = 2)

#mixture of normals from b)
lines(xGrid, mixDensMean, type = "l", lwd = 2, col = "red")




#################### QUESTION 2 ####################

data2 = read.table('eBayNumberOfBidderData.dat', header = T)


#a
model = glm(nBids ~ . - Const, data = data2, family = 'poisson')
model$coefficients

#b
# The following code was written using Mattias Villani's implementation as a template:
# https://github.com/mattiasvillani/BayesLearnCourse/raw/master/Code/MainOptimizeSpam.zip


y = as.vector(data2[,1])
x = as.matrix(data2[,2:10])
nfeatures = dim(x)[2]

# prior
mu = as.vector(rep(0,nfeatures))
sigma = 100*solve(t(x)%*%x)

LogPrior = function(theta, mu, sigma){
  return(dmvnorm(theta, as.matrix(mu), sigma, log = T))
}

# log likelihood
LogLikelihood = function(theta, y, x){
  sum(y * x%*%theta - exp(x%*%theta) - log(factorial(y))) #log-likelihood of poisson regression
}

# log posterior
LogPosterior = function(theta, y, x, mu, sigma){
  LogPrior(theta, mu, sigma) + LogLikelihood(theta, y, x) 
  # we can sum them instead of multiplying, because of the log
}

# initialize Beta vector randomly
# Seed
set.seed(1234567890)
Beta_init = as.vector(rnorm(dim(x)[2]))

# optimizing the log posterior by changing the Betas (maximize)
res = optim(Beta_init, LogPosterior, gr = NULL, y, x, mu, sigma, 
            method="BFGS", control=list(fnscale=-1), hessian=T)

Beta_hat = res$par # posterior mode
Hessian = res$hessian
post_sigma = solve(-Hessian) # posterior cov matrix
stdev = sqrt(diag(post_sigma))

### end of code from Mattias Villani's template


#c

RWMSampler = function(logPostFunc, c, iter, initial, proposal_sigma, ...){
  X = matrix(0,nrow = iter+1, ncol = length(initial))
  X[1,] = initial
  
  #acceptance probability stats
  alphas = numeric(iter)
  
  for (t in 1:iter) {
    print(t)
    Y = as.numeric(rmvnorm(1, mean = X[t,], sigma = c * proposal_sigma)) # Proposal needs to be a numeric vector
    U = runif(1, min = 0, max = 1) # Generate uniform
    
    alpha = min(1,
                (exp(logPostFunc(Y, ...) - logPostFunc(X[t,], ...))))
    
    if (U < alpha) {
      X[t+1,] = Y
    } else {
      X[t+1,] = X[t,]
    }
    alphas[t] = alpha # save acceptance probability
    t = t+1
  }
  return(list('mcmc' = X, 'alphas' = alphas))
}


# Seed
set.seed(1234567890)
nIter = 100
mc = RWMSampler(LogPosterior, c = 0.65, iter = nIter, initial = Beta_init, proposal_sigma = post_sigma, y, x, mu, sigma)

for (p in 1:9) {
  plot(mc$mcmc[,p], type = 'l', main = p)
}

burn_in_theta = 1000

for (p in 1:9) {
  hist(mc$mcmc[burn_in_theta:nIter,p], main = p, breaks = 30)
}

mean(mc$alphas) #0.27

Beta_mcmc = mc$mcmc[burn_in_theta:nIter,]


#d

Const = 1
PowerSeller = 1
VerifyID = 1
Sealed = 1
MinBlem = 0
MajBlem = 0
LargNeg = 0
LogBook = 1
MinBidShare = 0.5
x_pred = c(Const, PowerSeller, VerifyID, Sealed, MinBlem, MajBlem, LargNeg, LogBook, MinBidShare)

predictive_dist = numeric(dim(Beta_mcmc)[1])
# Seed
set.seed(1234567890)
for (i in 1:dim(Beta_mcmc)[1]) {
  predictive_dist[i] = rpois(1, exp(t(as.matrix(x_pred)) %*% as.matrix(Beta_mcmc[i,])))
}
hist(predictive_dist, freq = T)
sum(predictive_dist==0)/length(predictive_dist) #0.3548495

#normalized predictive distribution
values = unique(predictive_dist)
values = cbind(values,numeric(length(values)))
for (v in values[,1]) {
  count = sum(predictive_dist==v)
  values[values[,1]==v] = c(v, count)
}

values = data.frame(values)
colnames(values) = c('value', 'count')
values$ratio = values$count/sum(values$count)
values = values[order(values$value),]

barplot(values$ratio)
