# R version
RNGversion('3.5.1')

# libraries


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

x2 = sim_AR_1(mu, mu, sigma2, 0, t)
plot(x2, type = 'l')

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
  real phi;
  real<lower=0> sigma;
}
model {
  mu ~ normal(10,10);
  // phi ~ uniform(-1,1); // Setting the limit of beta (or phi)
  // sigma ~ uniform(0,100);
  for (n in 2:N)
    y[n] ~ normal(mu+phi*(y[n-1]-mu),sigma);
}'

#stan uses uniform prior by default anyway (uniform(-infinity, +infinity)), but we can use that for a non-informative prior 

N = length(y)
data = list(N=N, y=y)
burnin = 1000
niter = 2000
fit = stan(model_code = StanModel, data = data, warmup = burnin, iter = niter, chains = 4)


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

#i

#posterior means:
#alpha = 0.43
#beta = 0.96
#sigma = 1.29

#95% credible intervals:
#alpha = [-0.06, 0.92]
c(0.43 - 1.96 * 0.25, 0.43 + 1.96 * 0.25)
#beta = [0.9208, 0.9992]
c(0.96 - 1.96 * 0.02, 0.96 + 1.96 * 0.02)
#sigma = [1.1528, 1.4370]
c(1.29 - 1.96 * 0.07, 1.29 + 1.96 * 0.075)

#effective posterior samples
#alpha = 1775
#beta = 1774
#sigma = 1965

#YES, WE CAN ESTIMATE THE TRUE VALUES

#ii
traceplot(fit) #good mixing, good convergence, because all 4 chains arrive to the same conclusion

plot(postDraws$alpha,type="l",ylab="mu",main="Traceplot")
plot(postDraws$sigma,type="l",ylab="mu",main="Traceplot")
hist(postDraws$alpha, breaks = 30)
hist(postDraws$sigma, breaks = 30)

# run only if permuted = FALSE in extract()
plot(postDraws[,1,1], type = 'l')
plot(postDraws[,2,1], type = 'l')
plot(postDraws[,3,1], type = 'l')
plot(postDraws[,4,1], type = 'l')




#### c
data = read.table('campy.dat', header = T)




#### d




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
