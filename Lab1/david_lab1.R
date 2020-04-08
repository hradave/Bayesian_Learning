# libraries
library(geoR)
library(HDInterval)
library(MESS)

####################################### 1 #######################################

#a

n = 20
s = 5
f = n-s
alpha_0 = 2
beta_0 = 2

posterior_draw = rbeta(n = 10000, shape1 = alpha_0 + s, shape2 = beta_0 + f)
hist(posterior_draw, freq = F)
mean(posterior_draw)
sd(posterior_draw)

# [2:31 PM] Hariprasath Govindarajan
# Histogram is probably not appropriate because you want to check for different number of random draws.
# You could calculate the mean and SD for different number of draws and make plots of mean vs number of draws, SD vs number of draws to analyze their convergence.
# you cannot conclude about convergence from just 3 data points. 


#b
nDraws = 10000
posterior_draw = rbeta(n = nDraws, shape1 = alpha_0 + s, shape2 = beta_0 + f)
larger_than_03_sim = sum(posterior_draw > 0.3) / length(posterior_draw)
larger_than_03_true = 1 - pbeta(q = 0.3, shape1 = alpha_0 + s, shape2 = beta_0 + f)


#c
nDraws = 10000
log_odds = log(posterior_draw/(1-posterior_draw))
hist(log_odds, freq = F)
plot(density(log_odds))


####################################### 2#######################################

#a
Y = c(44, 25, 45, 52, 30, 63, 19, 50, 34, 67)
logY = log(Y)
plot(density(Y))
plot(density(logY))

n = length(Y)
mu = 3.7
tau2 = sum((logY - mu)^2) / n

nDraws = 10000

# simulate from posterior for sigma^2
X = rchisq(nDraws,n)
sigma2_sim = n*tau2/X
plot(density(sigma2_sim))


# theoretical posterior for sigma^2?

sigma_posterior_draw = rinvchisq(n = nDraws, df = n, scale = tau2) # SHOULD WE USE THE RINVCHISQ?

plot(density(sigma_posterior_draw))
hist(sigma_posterior_draw, freq = F)



#b

G = 2 * pnorm(sqrt(sigma2_sim / 2)) - 1
hist(G)
plot(density(G))

# Gini index mean is closer to 0 than 0 -> more equal income than not



#c

#90% equal tail
equal_tail = quantile(G, probs = c(0.05, 0.95)) #0.1606429, 0.3368407

#HPD
densx = density(G)$x
densy = density(G)$y
sorted = cbind(seq(1,length(densy)), sort(densy, decreasing = T))
plot(sorted, type = 'l')
total_area = auc(sorted[,1], sorted[,2])

area = 0
i = 1
while (area/total_area < 0.9) {
  area = auc(sorted[1:(i+1),1], sorted[1:(i+1),2])
  i = i+1
}

cutoff = sorted[(i-1),2] # 1.773939 (subtract 1, because of the last unwanted loop iteration)
# the closest G-s from both tails from densx are 0.14747841, 0.31634166
HPD = c(0.14747841,0.31634166)

# alternative using package HDInterval
hdi(G, credMass = 0.9) #0.1513655, 0.3413572
hdi = c(0.1513655, 0.3413572)

# plotting
plot(density(G), main = "Kernel density estimate of G")
abline(v = equal_tail, lwd = 1, col = "blue")
abline(v = HPD, lwd = 1, col = "red")
abline(v = hdi, lwd = 1, col = "green")
legend(x = "topright", legend=c("equal tail 90%","HPD analytical 90%", "HPD package 90%"),
       col = c("blue", "red", "green"), lwd = 1)

# normalize PDF by dividing it by the area under the curve (general comment)




####################################### 3 #######################################


#a
y = c(40, 303, 326, 285, 296, 314, 20, 308, 299, 296)
y_rad = c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)
mu = 2.39
k = seq(0.01,5,by = 0.01)

# posterior_k = numeric()
# counter = 1
# for (current_k in k) {
#   posterior_k[counter] = 1 / ((2*pi*besselI(current_k, nu = 0))^10) * 
#                 exp(current_k * (sum(cos(y_rad-mu))-1))
#   counter = counter + 1
#   
# }
# plot(density(posterior_k))

posterior_k = 1 / ((2*pi*besselI(k, nu = 0))^10) * 
  exp(k * (sum(cos(y_rad-mu))-1))

plot(k,posterior_k, type = 'l', main = 'Posterior of k')

#b
k[which.max(posterior_k)] #2.12





# # test
# k = 4
# ytest = seq(-pi,pi,by=0.1)
# test = exp(k * (cos(ytest-pi))) / (2*pi*besselI(k, nu = 0))
#   
# plot(ytest,test)
# 
# xnorm = seq(1,2,by=0.01)
# pn = 1/(0.0741^2 * sqrt(2*pi)) * exp(-1*(1/(2*0.0741^2))*(xnorm-1.512)^2)
# plot(xnorm, pn)
