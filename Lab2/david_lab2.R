# R version
RNGversion('3.5.1')


#libraries
library("mvtnorm")

# 2. Posterior approximation for classification with logistic regression

# 2.a

# The following code was written using Mattias Villani's implementation as a template:
# https://github.com/mattiasvillani/BayesLearnCourse/raw/master/Code/MainOptimizeSpam.zip

# read data
women = read.table('WomenWork.dat', header = T)

y = as.vector(women[,1])
x = as.matrix(women[,2:9])
tau = 10
nfeatures = dim(x)[2]


# prior
mu = as.vector(rep(0,nfeatures))
sigma = diag(tau^2, nfeatures, nfeatures)

LogPrior = function(Beta, mu, sigma){
  return(dmvnorm(Beta, as.matrix(mu), sigma, log = T))
}

# log likelihood
LogLikelihood = function(Beta, y, x){
  sum(y * x%*%Beta - log(1 + exp(x%*%Beta)))
}

# log posterior
LogPosterior = function(Beta, y, x, mu, sigma){
  LogPrior(Beta, mu, sigma) + LogLikelihood(Beta, y, x) 
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

# approximate 95% credible interval for NSmallChild
lower = Beta_hat[7] - 1.96 * stdev[7] #-2.121445
upper = Beta_hat[7] + 1.96 * stdev[7] #-0.5968567

# ANSWER:
# Yes, NSmallChild is an important determinant of whether a woman is working or not,
# because it has the highest absolute valued posterior coefficient (Beta_hat) of all features.
# Its coefficient is -1.36, which means that the more small children a woman has, the less likely
# it is that she's working, because the sign is negative.

# check results
glmModel <- glm(Work~0 + ., data = women, family = binomial)
glmModel$coefficients #roughly the same as Beta_hat

# 2.b

# The mode of a normal distribution is equal to the mean.

log_reg = function(Beta, x){
  return(exp(x%*%Beta) / (1 + exp(x%*%Beta)))
}


Constant = 1
HusbandInc = 10
EducYears = 8
ExpYears = 10
ExpYears2 = (ExpYears/10)^2
Age = 40
NSmallChild = 1
NBigChild = 1

x_pred = c(Constant, HusbandInc, EducYears, ExpYears, ExpYears2, Age, NSmallChild, NBigChild)

# simulate posterior draws from Beta
n = 10000
post_sim = rmvnorm(n, mean = Beta_hat, sigma = post_sigma) # should i add population variance to post_sigma?

# calculate target y
pred = apply(post_sim, 1, log_reg, x_pred)
hist(pred, freq = F)
lines(density(pred))


# 2.c





