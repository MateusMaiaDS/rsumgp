# Creating the omega matrix
omega <- function(x,nu,phi){

  # Returning the omega matrix
  omega_matrix <- exp(-as.matrix(stats::dist(x))^(2)/(2*phi^2) )/nu

  return(omega_matrix)
}


# GP sample
gp_sample <- function(x_train, resid, x_test,
                      phi, tau, nu){

  # Getting the main components
  omega_matrix <- omega(x = x_train,nu = nu,phi = phi)
  K_y <- omega_matrix + diag(1/tau, nrow = nrow(x_train))

  # Calculating the mean

  sample_mean <- crossprod(omega_matrix,solve(K_y,resid))
  sample_var <- omega_matrix  - crossprod(omega_matrix,solve(K_y,omega_matrix))

  # Sampling from a normal
  sampled_residual <- mvtnorm::rmvnorm(n = 1,mean = sample_mean,sigma = sample_var)

  return(sampled_residual)

}

# Sample tau
sample_tau <- function(y,predict_matrix,a_tau, d_tau){

  stats::rgamma(n = 1,shape = 0.5*length(y) + a_tau, rate = crossprod((y-colSums(predict_matrix))))

}


pi_coverage <- function(y, y_hat_post, sd_post,only_post = FALSE, prob = 0.5,n_mcmc_replications = 1000){

  # Getting the number of posterior samples and columns, respect.
  np <- nrow(y_hat_post)
  nobs <- ncol(y_hat_post)

  full_post_draw <- list()

  # Setting the progress bar
  progress_bar <- utils::txtProgressBar(
    min = 1, max = n_mcmc_replications,
    style = 3, width = 50 )

  # Only post matrix
  if(only_post){
    post_draw <- y_hat_post
  } else {
    for(i in 1:n_mcmc_replications){
      utils::setTxtProgressBar(progress_bar, i)

      full_post_draw[[i]] <-(y_hat_post + replicate(sd_post,n = nobs)*matrix(stats::rnorm(n = np*nobs),
                                                                             nrow = np))
    }
  }

  if(!only_post){
    post_draw<- do.call(rbind,full_post_draw)
  }

  # CI boundaries
  low_ci <- apply(post_draw,2,function(x){stats::quantile(x,probs = prob/2)})
  up_ci <- apply(post_draw,2,function(x){stats::quantile(x,probs = 1-prob/2)})

  pi_cov <- sum((y<=up_ci) & (y>=low_ci))/length(y)

  return(pi_cov)
}
