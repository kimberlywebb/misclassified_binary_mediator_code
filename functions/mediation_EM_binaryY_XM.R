mediation_EM_binaryY_XM <- function(Mstar, # Observed mediator vector
                                    outcome, # Outcome vector
                                    # Predictor matrices
                                    x_matrix, z_matrix, c_matrix,
                                    # Start values for parameters
                                    beta_start, gamma_start,
                                    theta_start, sigma_start,
                                    # EM settings
                                    tolerance = 1e-7, max_em_iterations = 1500,
                                    em_method = "squarem"){

  n_cat = 2 # Number of categories in mediator
  sample_size = length(Mstar) # Sample size
  
  # Create design matrices
  X = matrix(c(rep(1, sample_size), c(x_matrix)),
             byrow = FALSE, nrow = sample_size)
  Z = matrix(c(rep(1, sample_size), c(z_matrix)),
             byrow = FALSE, nrow = sample_size)
  
  # Create a matrix of observed mediator variables using dummy coding
  obs_M_reps = matrix(rep(Mstar, n_cat), nrow = sample_size, byrow = FALSE)
  category_matrix = matrix(rep(1:n_cat, each = sample_size), nrow = sample_size,
                           byrow = FALSE)
  obs_M_matrix = 1 * (obs_M_reps == category_matrix)
  
  # EM algorithm settings
  control_settings = list(convtype = "parameter", tol = tolerance,
                          stoptype = "maxiter", maxiter = max_em_iterations)
  
  # Run EM algorithm using turboEM
  results = turboEM::turboem(par = c(c(beta_start), c(gamma_start),
                                     c(theta_start)),
                             fixptfn = EM_function_binaryY_XM,
                             method = c(em_method),
                             obs_mediator = Mstar,
                             obs_outcome = outcome,
                             X = X, Z = Z, c_matrix = c_matrix,
                             sample_size = sample_size, n_cat = n_cat,
                             control.run = control_settings)
  
  # Recode observed mediator from 1/2, make 1/0
  Mstar01 = ifelse(Mstar == 1, 1, ifelse(Mstar == 2, 0, NA))
  
  # Do label switching correction within the EM algorithm simulation
  results_i_gamma <- matrix(turboEM::pars(results)[(ncol(X) + 2):(ncol(X) + 1 + (n_cat * ncol(Z)))],
                            ncol = n_cat, byrow = FALSE)
  results_i_pistar_v <- pistar_compute(results_i_gamma, Z, sample_size, n_cat)
  pistar_11 <- mean(results_i_pistar_v[1:sample_size, 1])
  pistar_22 <- mean(results_i_pistar_v[(sample_size + 1):(2*sample_size), 2])
  
  J <- pistar_11 + pistar_22 - 1
  J_flip <- (1 - pistar_22) + (1 - pistar_11) - 1
  estimates_i <- if ((J >= J_flip) |
                     (is.na(pistar_11) | is.na(pistar_22))) {
    # If turboem cannot estimate the parameters they will be NA.
    turboEM::pars(results)
  } else {
    gamma_index = (ncol(X) + 2):(ncol(X) + 1 + (n_cat * ncol(Z)))
    n_gamma_param = length(gamma_index) / n_cat
    gamma_flip_index = ncol(X) + 1 + c((n_gamma_param + 1):length(gamma_index), 1:n_gamma_param)
    c(-1*turboEM::pars(results)[1:(ncol(X) + 1)], turboEM::pars(results)[gamma_flip_index],
      turboEM::pars(results)[8] + turboEM::pars(results)[10],
      turboEM::pars(results)[9] + turboEM::pars(results)[12],
      -1 * turboEM::pars(results)[10],
      turboEM::pars(results)[11],
      -1 * turboEM::pars(results)[12])
  }
  
  # Set parameter names
  beta_param_names <- paste0(rep("beta", ncol(X)), 1:3)
  gamma_param_names <- paste0(rep("gamma", (n_cat * ncol(Z))),
                              rep(1:ncol(Z), n_cat),
                              rep(1:n_cat, each = ncol(Z)))
  theta_param_names <- c("theta_0", "theta_x", "theta_m", "theta_c",
                         "theta_xm")
  
  # Return data frame of estimates
  estimates <- data.frame(Parameter = c(beta_param_names,
                                        gamma_param_names,
                                        theta_param_names),
                          Estimates = c(estimates_i),
                          Convergence = c(rep(results$convergence,
                                              length(c(beta_param_names,
                                                       gamma_param_names,
                                                       theta_param_names)))))

  return(estimates)
}
