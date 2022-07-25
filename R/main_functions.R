# Creating the main function
gp_sum <- function(x_train,
                   y,
                   x_test,
                   phi,
                   nu,
                   n_gp,
                   n_mcmc,
                   n_burn,
                   a_tau = 0.001,
                   d_tau = 0.001){

    # Creating the auxiliar functions
    current_partial_residuals <- numeric(length(y))

    # Getting tree predictions
    gp_predictions <- matrix(0, nrow = n_gp,ncol = length(y))
    tau_post <- numeric(n_mcmc)
    y_hat_post <- matrix(0,nrow = n_mcmc,ncol = length(y))
    gp_pred_list = list()

    # Iterating all observations
    for(i in 1:n_mcmc){

      for(t in 1:n_gp){
        current_partial_residuals <- y - colSums(gp_predictions[-t,,drop = FALSE])

        gp_predictions[t,] <- gp_sample(x_train = x_train,
                  resid = current_partial_residuals,nu = nu,
                  x_test = x_test,phi = phi[t],tau = tau)

        # Plotting and observing the current GP prediction
        # points(x_train,gp_predictions[t,])
      }

      # plot(x_train,y,pch=20)
      # for(i in 1:n_gp){
      #   points(x_train,gp_predictions[i,],col = i)
      # }

      tau <- tau_post[i] <- sample_tau(y = y,predict_matrix = gp_predictions,
                                       a_tau = a_tau,d_tau = d_tau)

      y_hat_post[i,] <- colSums(gp_predictions)


      gp_pred_list[[i]] <- gp_predictions
    }

    return(list(y_hat = y_hat_post,
                tau = tau_post,
                gps = gp_pred_list))

}
