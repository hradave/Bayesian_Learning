# Linear and polynomial regression

data <- read.table("TempLinkoping.txt", header = TRUE)
View(data)

y <- data$temp
time <- data$time

n <- length(y)
X <- cbind(matrix(1, nrow = length(y)),time,time^2)
View(X)

# Initial prior hyperparameters
mu_0 <- c(-10,100,-100)
omega_0 <- 0.1*diag(3)
nu_0 <- 4
sigma2_0 <- 1

# Conjugate prior
draw_sigma2 <- nu_0*sigma2_0/rchisq(10, nu_0)
draw_beta <- matrix(0, nrow=length(draw_sigma2), ncol=3 )
for (sigma in draw_sigma2) {
  i <- which(draw_sigma2==sigma)
  draw_beta[i,] <- rmvnorm(1, mu_0, sigma*omega_0)
}


## ???? How many?? ###

#### X = rchisq(nDraws,n)
#### sigma2_sim = n*tau2/X
