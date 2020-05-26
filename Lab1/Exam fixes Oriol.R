
library(geoR)
library(HDInterval)
library(MESS)
library(ggplot2)

# 2a

Y = c(44, 25, 45, 52, 30, 63, 19, 50, 34, 67)
logY = log(Y)
n = length(Y)
mu = 3.7
tau2 = sum((logY - mu)^2) / n

v <- n #degrees of freedom
s <- tau2 #scale

phi <- seq(0, 1.5, by=0.001) #random variable 

p_inv_chi <- ((s*v/2)^(v/2)*exp((-v*s)/(2*phi)))/(gamma(v/2)*phi^(1+v/2))

# Seed
set.seed(1234567890)

nDraws = 10000
X = rchisq(nDraws,n)
sigma2_sim = n*tau2/X
hist(sigma2_sim, freq = F, ylim=c(0,6.5), breaks=30, xlim=c(0,1.5), 
     main = "Simulated Inv-Chi-squared", xlab="sigma-squared", col="lightgrey")


#plot(p_inv_chi, from=0, to=1.5, na.rm = T), col="red", lwd=2)
lines(phi, p_inv_chi, col="red", lwd=2, type="l")
  


### 3a

# 3 Bayesian inference for the concentration parameter in the von Mises distribution.

y = c(40, 303, 326, 285, 296, 314, 20, 308, 299, 296)
y_rad = c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)
n = length(y_rad)
mu = 2.39

#a
k = seq(0.01,7,by = 0.01)

posterior_k = 1 / ((2*pi*besselI(k, nu = 0))^n) * 
  exp(k * (sum(cos(y_rad-mu))-1))

normalized_posterior <- posterior_k/(sum(posterior_k)*0.01)

plot(k,normalized_posterior, type = 'l', main = 'Posterior of k', xlab="K values", ylab="Posterior")

