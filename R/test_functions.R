rm(list=ls())
# Creating the dataset
n_train <- n_test <- 500
sd_val <- 0.1
n_gp <- 1
x_train <- matrix(seq(-pi,pi,length.out = n_train))
y <- sin(x_train) + rnorm(n = n_train,sd = sd_val)
y_max <- max(y)
y_min <- min(y)
y <- (y - min(y))/ (max(y) - min(y)) - 0.5
x_test <- matrix(seq(-pi,pi, length.out = n_test))

# Setting other values as default
phi <- rep(1,n_gp)
nu <- 16*n_gp
tau <- 100*(y_max-y_min)^2
n_mcmc <- 1000
a_tau <- 0.001
d_tau <- 0.001


# Sum of GP's
gpsum_obj <- gp_sum(x_train = x_train,
       y = y,
       x_test = x_test,phi = phi,nu = nu,n_gp = n_gp,n_mcmc = n_mcmc,n_burn = 100)

# plotting the data
pi_coverage(y = y,
            y_hat_post = gpsum_obj$y_hat,
            sd_post = gpsum_obj$tau^(-1/2),
            prob = 0.5,n_mcmc_replications = 100)


# Generating the posterior matrix
posterior_y <- matrix(0, nrow = n_mcmc,ncol = length(y))

# Generating the quantile posterior for y
for(i in 1:n_mcmc){
  for(j in 1:n_train){
    posterior_y[i,j] <- rnorm(n = 1,mean = gpsum_obj$y_hat[i,j],sd = gpsum_obj$tau[i]^(-1/2))
    # posterior_y[i,j] <- rnorm(n = 1,mean = gpsum_obj$y_hat[i,j],sd = tau^(-1/2))

  }
}

# Compare the PI intervals calculated from the mean values on

plot(x_train,y,pch = 20 )
lines(x_train,colMeans(gpsum_obj$y_hat),col = "blue")
lines(x_train,apply(posterior_y,2,function(x){quantile(x,probs = 0.25)}),lty = "dashed", col = "blue")
lines(x_train,apply(posterior_y,2,function(x){quantile(x,probs = 0.75)}),lty = "dashed", col = "blue")

ci_low <- apply(posterior_y,2,function(x){quantile(x,probs = 0.25)})
ci_up <- apply(posterior_y,2,function(x){quantile(x,probs = 0.75)})

mean((ci_low <= y) & (ci_up>=y))


for(i in 1:n_gp){
  lines(x_train,gpsum_obj$gps[[1000]][i,], col = i)
}

