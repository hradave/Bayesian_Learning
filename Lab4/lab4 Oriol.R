# Time series models in Stan

### (a)

mu <- 10
sigma2 <- 2
T <- 200
phi <- 0

simulation_AR <- function(mu, phi, sigma2, T){
  x_t <- numeric(T)
  x_t[1] <- mu
  for (i in 2:T) {
    x_t[i] <- mu + phi*(x_t[i-1] - mu)  + rnorm(1, 0, sqrt(sigma2)) 
  }
  return(x_t)
}

simulation <- simulation_AR(mu,phi,sigma2,T)
plot(1:T, simulation, type='l')


### (b)

# Synthetic data
phi_1 <- 0.3
x1_t <- simulation_AR(mu, phi_1, sigma2, T)

phi_2 <- 0.95
y1_t <- simulation_AR(mu, phi_2, sigma2, T)

