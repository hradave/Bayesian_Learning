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







#b









#c











####################################### 3 #######################################






#a








#b













