library("mvtnorm")

# Linear and polynomial regression

data <- read.table("TempLinkoping.txt", header = TRUE)

y <- data$temp
time <- data$time

n <- length(y)
X <- cbind(matrix(1, nrow = length(y)),time,time^2)

# Initial prior hyperparameters
mu_0 <- c(-10,100,-100)
omega_0 <- 0.1*diag(3)
nu_0 <- 4
sigma2_0 <- 1

N_draws=1000

# Conjugate prior

draw_sigma2 <- function(N_draws, nu, sigma2){
  sigma2 <- nu*sigma2/rchisq(Ndraws, nu)
}

draw_betas <- function(sigma2_draws, mu, omega){
  draw_beta <- matrix(0, nrow=length(sigma2_draws), ncol=3 )
  for (sigma in sigma2_draws) {
    i <- which(sigma2_draws==sigma)
    draw_beta[i,] <- rmvnorm(1, mu, sigma*omega)
  }
  return(draw_beta)
}

plot_regression_draws <- function(X, betas, time){
  for (i in 1:length(betas[,1])) {
    if (i==1) plot(time, X%*%betas[i,], type = 'l', ylim=c(-15, 25), main = 'Temperature over time', ylab='temperature')
    else points(time, X%*%betas[i,], type='l')
    points(time,X%*%colMeans(betas),type = 'l',lwd=2,col='yellow')
  }
}


# Tune hyperparameters

# Initial hyperparameters
beta_0 <- compute_betas(Ndraws, mu_0, omega_0, nu_0, sigma2_0)

# 1st iteration
mu_1 <- c(-5,90,-90)
beta_1 <- compute_betas(Ndraws, mu_1, omega_0, nu_0, sigma2_0)

# 2nd iteration
mu_2 <- c(0,100,-100)
beta_2 <- compute_betas(Ndraws, mu_2, omega_0, nu_0, sigma2_0)

# 3rd iteration
mu_3 <- c(0,85,-85)
beta_3 <- compute_betas(Ndraws, mu_3, omega_0, nu_0, sigma2_0)

# 4th iteration
mu_4 <- c(-5, 85, -85)
beta_4 <- compute_betas(Ndraws, mu_4, omega_0, nu_2, sigma2_0)

plot(time, X %*% beta_0, type = 'l', lwd=2, col='red', ylim=c(-10, 25),
     ylab='Temperature', xlab='Time')
points(time, X %*% beta_1, lwd=2, col='blue', type='l')
points(time, X %*% beta_2, lwd=2, col='pink', type='l')
points(time, X %*% beta_3, lwd=2, col='green', type='l')
points(time, X %*% beta_4, lwd=2, col='yellow', type='l')


# We like blue

beta_hat <- solve(t(X)%*%X)%*%t(X)%*%y
mu_n <- solve(t(X)%*%X + omega_0) %*% (t(X)%*%X%*%beta_hat + omega_0%*%mu_0)
omega_n <- t(X) %*% X + omega_0
nu_n <- nu_0 + n # what is n?????????????????????
nu_sigma2_n <- nu_0*sigma2_0+(t(y)%*%y + t(mu_0)%*%omega_0%*%mu_0 - t(mu_n)%*%omega_n%*%mu_n)
sigma2_n <- nu_sigma2_n[1] / nu_n

draw_sigma2 <- nu_n*sigma2_n/rchisq(Ndraws, nu_n)
draw_beta <- matrix(0, nrow=length(draw_sigma2), ncol=3 )
for (sigma in draw_sigma2) {
  i <- which(draw_sigma2==sigma)
  draw_beta[i,] <- rmvnorm(1, mu_n, sigma*omega_n)
}






#### X = rchisq(nDraws,n)
#### sigma2_sim = n*tau2/X
