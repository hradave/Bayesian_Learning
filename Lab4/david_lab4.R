# R version
RNGversion('3.5.1')

# libraries
library(rstan)

# Seed
set.seed(1234567890)



################# QUESTION 1 ###################


#### a
mu = 10
sigma2  = 2
t = 200

sim_AR_1 = function(start, mu, sigma2, phi, t){
  x = numeric(t)
  x[1] = start
  for (i in 2:t) {
    eps = rnorm(1, mean = 0, sd = sqrt(sigma2))
    x[i] = mu + phi*(x[i-1] - mu) + eps
  }
  return(x)
}


# Seed
set.seed(1234567890)
x1 = sim_AR_1(mu, mu, sigma2, -1, t)
plot(x1, type = 'l')

set.seed(1234567890)
x2 = sim_AR_1(mu, mu, sigma2, 0, t)
plot(x2, type = 'l')

set.seed(1234567890)
x3 = sim_AR_1(mu, mu, sigma2, 1, t)
plot(x3, type = 'l')

#### b

# Seed
set.seed(1234567890)
x = sim_AR_1(mu, mu, sigma2, 0.3, t)
y = sim_AR_1(mu, mu, sigma2, 0.95, t)

plot(x, type = 'l')
plot(y, type = 'l')

StanModel = '
data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real mu;
  real<lower=-1, upper=1> phi; // this is not recommended, but runs better than with a restrictive prior
  real<lower=0> sigma;
}
model {
  mu ~ normal(10,100); //non-informative, because of large sd
  phi ~ normal(0,10); //weakly informative, because we know it has to be between -1 and 1
  sigma ~ scaled_inv_chi_square(1,1.4); //non-informative, because of small df
  for (n in 2:N)
    y[n] ~ normal(mu+phi*(y[n-1]-mu),sigma);
}'

#stan uses uniform prior by default anyway (uniform(-infinity, +infinity)), but we can use that for a non-informative prior 

N = length(y)
data_y = list(N=N, y=y)
burnin = 1000
niter = 2000
fit_y = stan(model_code = StanModel, data = data_y, warmup = burnin, iter = niter, chains = 4)


# Print the fitted model
print(fit_y)
# Extract posterior samples
postDraws_y <- extract(fit_y)
# Do traceplots of the first chain
par(mfrow = c(1,1))
plot(postDraws_y$mu[1:(niter-burnin)],type="l",ylab="mu",main="Traceplot")
# Do automatic traceplots of all chains
traceplot(fit_y)
# Bivariate posterior plots
pairs(fit_y)




data_x = list(N=N, y=x)
fit_x = stan(model_code = StanModel, data = data_x, warmup = burnin, iter = niter, chains = 4)

# Print the fitted model
print(fit_x)
# Extract posterior samples
postDraws_x <- extract(fit_x)
# Do traceplots of the first chain
par(mfrow = c(1,1))
plot(postDraws_x$mu[1:(niter-burnin)],type="l",ylab="mu",main="Traceplot")
# Do automatic traceplots of all chains
traceplot(fit_x)
# Bivariate posterior plots
pairs(fit_x)


############# ANSWER #############

#It works well for process x (with phi=0.3), but not so well for process y (with phi=0.95) using non-informative priors.

##################################

#i


########### y ##########

#posterior means:
#mu = 6.6
#phi = 0.97
#sigma = 1.29

#95% credible intervals: (or just use the one from the output)
#mu = [-56.12, 69.32]
c(6.6 - 1.96 * 32, 6.6 + 1.96 * 32)
#phi = [0.9308, 1.0092]
c(0.97 - 1.96 * 0.02, 0.97 + 1.96 * 0.02)
#sigma = [1.1724, 1.4076]
c(1.29 - 1.96 * 0.06, 1.29 + 1.96 * 0.06)

#effective posterior samples
#mu = 687
#phi = 83
#sigma = 384




########### x ##########

#posterior means:
#mu = 10.25
#phi = 0.39
#sigma = 1.49

#95% credible intervals:
#mu = [9.91, 10.60]
#phi = [0.26, 0.44]
#sigma = [1.35, 1.65]

#effective posterior samples
#mu = 3159
#phi = 3847
#sigma = 3365



#YES, WE CAN ESTIMATE THE TRUE VALUES FROM X, BUT NOT FROM Y

#ii

######### y #########

traceplot(fit_y) #bad mixing, bad convergence, because the chains do not oscillate too much, and they don't cover the same areas.
plot(postDraws_y$mu,type="l",ylab="mu",main="Traceplot")
plot(postDraws_y$sigma,type="l",ylab="mu",main="Traceplot")
hist(postDraws_y$mu, breaks = 30)
hist(postDraws_y$sigma, breaks = 30)

# run only if permuted = FALSE in extract()
plot(postDraws_y[,1,1], type = 'l')
plot(postDraws_y[,2,1], type = 'l')
plot(postDraws_y[,3,1], type = 'l')
plot(postDraws_y[,4,1], type = 'l')



######### x #########

traceplot(fit_x) #good mixing, good convergence, because the chains oscillate a lot, and they cover the same areas, so they arrived at the same conclusion.
plot(postDraws_x$mu,type="l",xlab="mu",main="Traceplot")
plot(postDraws_x$sigma,type="l",xlab="mu",main="Traceplot")
hist(postDraws_x$mu, breaks = 30)
hist(postDraws_x$sigma, breaks = 30)

# run only if permuted = FALSE in extract()
plot(postDraws_x[,1,1], type = 'l')
plot(postDraws_x[,2,1], type = 'l')
plot(postDraws_x[,3,1], type = 'l')
plot(postDraws_x[,4,1], type = 'l')



#### c
data_campy = read.table('campy.dat', header = T)
plot(x=c(1:140), y=data_campy$c, type = 'l')

StanModel_poisson = '
data {
  int<lower=0> N;
  int<lower=0> c[N];
}
parameters {
  real mu;
  real<lower=-1, upper=1> phi; // this is not recommended, but runs better than with a restrictive prior
  real<lower=0> sigma;
  real x[N];
}
model {
  mu ~ normal(0,100); // non-informative prior, we know nothing about mu
  phi ~ normal(0,10); //weakly informative, because we know it has to be between -1 and 1
  sigma ~ scaled_inv_chi_square(1,2); //non-informative, because of small df, we know nothing about sigma
  for (n in 2:N){
    x[n] ~ normal(mu+phi*(x[n-1]-mu),sigma);
    c[n] ~ poisson(exp(x[n]));
  }
}'

N = length(data_campy$c)
data = list(N=N, c=data_campy$c)
burnin = 1000
niter = 2000
fit_poisson = stan(model_code = StanModel_poisson, data = data, warmup = burnin, iter = niter, chains = 4)

# Print the fitted model
print(fit_poisson)
# Extract posterior samples
postDraws <- extract(fit_poisson)
# Do traceplots of the first chain
par(mfrow = c(1,1))
plot(postDraws$mu[1:(niter-burnin)],type="l",ylab="mu",main="Traceplot")
# Do automatic traceplots of all chains
traceplot(fit_poisson, pars=c('mu', 'phi', 'sigma'))
# Bivariate posterior plots
pairs(fit_poisson)



#plot
intensity_posterior = data.frame(exp(postDraws[["x"]]))
intensity_posterior_means = colMeans(intensity_posterior)
intensity_posterior_sd = apply(intensity_posterior, 2, sd)

#first row = lower bound of 95% interval
#second row = upper bound of 95% interval
intensity_95_intervals = rbind(intensity_posterior_means - 1.96*intensity_posterior_sd,
                               intensity_posterior_means + 1.96*intensity_posterior_sd)

plot(x=c(1:140), y=data_campy$c, type = 'l', lwd = 2, xlab = 'time', ylab = 'infections')
lines(x=c(1:140), y=intensity_posterior_means, col = 'red', lwd = 2)
lines(x=c(1:140), y=intensity_95_intervals[1,], col = 'green')
lines(x=c(1:140), y=intensity_95_intervals[2,], col = 'orange')
legend(x='topleft', legend=c('Data', 'Posterior mean', '95% lower CI', '95% upper CI'), 
       col = c('black', 'red', 'green', 'orange'), lwd = 2)

#### d

StanModel_poisson2 = '
data {
  int<lower=0> N;
  int<lower=0> c[N];
}
parameters {
  real mu;
  real<lower=-1, upper=1> phi; // this is not recommended, but runs better than with a restrictive prior
  real<lower=0> sigma;
  real x[N];
}
model {
  mu ~ normal(0,100); // non-informative prior, we know nothing about mu
  phi ~ normal(0,10); //weakly informative, because we know it has to be between -1 and 1
  sigma ~ scaled_inv_chi_square(100,1); // informative, because of big df, we know that the error is small
  for (n in 2:N){
    x[n] ~ normal(mu+phi*(x[n-1]-mu),sigma);
    c[n] ~ poisson(exp(x[n]));
  }
}'

N = length(data_campy$c)
data = list(N=N, c=data_campy$c)
burnin = 1000
niter = 2000
fit_poisson2 = stan(model_code = StanModel_poisson2, data = data, warmup = burnin, iter = niter, chains = 4)

# Print the fitted model
print(fit_poisson2)
# Extract posterior samples
postDraws2 <- extract(fit_poisson2)
# Do traceplots of the first chain
par(mfrow = c(1,1))
plot(postDraws2$mu[1:(niter-burnin)],type="l",ylab="mu",main="Traceplot")
# Do automatic traceplots of all chains
traceplot(fit_poisson2, pars=c('mu', 'phi', 'sigma'))
# Bivariate posterior plots
pairs(fit_poisson2)



#plot
intensity_posterior2 = data.frame(exp(postDraws2[["x"]]))
intensity_posterior_means2 = colMeans(intensity_posterior2)
intensity_posterior_sd2 = apply(intensity_posterior2, 2, sd)

#first row = lower bound of 95% interval
#second row = upper bound of 95% interval
intensity_95_intervals2 = rbind(intensity_posterior_means2 - 1.96*intensity_posterior_sd2,
                               intensity_posterior_means2 + 1.96*intensity_posterior_sd2)

plot(x=c(1:140), y=data_campy$c, type = 'l', lwd = 2, xlab = 'time', ylab = 'infections')
lines(x=c(1:140), y=intensity_posterior_means2, col = 'red', lwd = 2)
lines(x=c(1:140), y=intensity_95_intervals2[1,], col = 'green')
lines(x=c(1:140), y=intensity_95_intervals2[2,], col = 'orange')
legend(x='topleft', legend=c('Data', 'Posterior mean', '95% lower CI', '95% upper CI'), 
       col = c('black', 'red', 'green', 'orange'), lwd = 2)


















# stan test
library(rstan)
library(rstan)
y = c(4,5,6,4,0,2,5,3,8,6,10,8)
N = length(y)
StanModel = '
data {
  int<lower=0> N; // Number of observations
  int<lower=0> y[N]; // Number of flowers
}
parameters {
  real mu;
  real<lower=0> sigma2;
}
model {
  mu ~ normal(0,100); // Normal with mean 0, st.dev. 100
  sigma2 ~ scaled_inv_chi_square(1,2); // Scaled-inv-chi2 with nu 1, sigma 2
  for(i in 1:N)
    y[i] ~ normal(mu,sqrt(sigma2));
}'


data = list(N=N, y=y)
burnin = 1000
niter = 2000
fit = stan(model_code=StanModel,data=data,
           warmup=burnin,iter=niter,chains=4)
# Print the fitted model
print(fit)
# Extract posterior samples
postDraws <- extract(fit)
# Do traceplots of the first chain
par(mfrow = c(1,1))
plot(postDraws$mu[1:(niter-burnin)],type="l",ylab="mu",main="Traceplot")
# Do automatic traceplots of all chains
traceplot(fit)
# Bivariate posterior plots
pairs(fit)
