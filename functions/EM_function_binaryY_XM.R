EM_function_binaryY_XM <- function(param_current,
                                   obs_mediator, obs_outcome,
                                   X, Z, c_matrix,
                                   sample_size, n_cat){

  # Separate current parameters into beta, gamma, theta, sigma vectors
  beta_current = matrix(param_current[1:3], ncol = 1)
  gamma_current = matrix(c(param_current[4:7]),
                         ncol = n_cat, byrow = FALSE)
  theta_current = matrix(c(param_current[8:12]),
                         ncol = 1)

  # Create design matrix for true mediator model
  design_matrix = cbind(X, c_matrix)
  
  # Compute probability of each latent mediator value
  probabilities = pi_compute(beta_current, design_matrix, sample_size, n_cat)
  
  # Compute probability of observed mediator, given latent mediator
  conditional_probabilities = pistar_compute(gamma_current, Z, sample_size, n_cat)
  
  # Compute likelihood value of Y based on x, m, c, theta, and sigma
  p_yi_m0 = p_binary_yi_XM(theta_v = theta_current,
                           m = 0, x = X[,2], c = c_matrix,
                           outcome = obs_outcome,
                           sample_size, n_cat)
  p_yi_m1 = p_binary_yi_XM(theta_v = theta_current,
                           m = 1, x = X[,2], c = c_matrix,
                           outcome = obs_outcome,
                           sample_size, n_cat)

  # Create a matrix of observed mediator variables using dummy coding
  mstar_matrix = matrix(c(ifelse(obs_mediator == 1, 1, 0), 
                          ifelse(obs_mediator == 2, 1, 0)),
                        nrow = sample_size, byrow = FALSE)
  
  # Create a matrix of outcomes
  outcome_matrix = matrix(c(obs_outcome,
                            1 - obs_outcome),
                          nrow = sample_size, byrow = FALSE)
  
  # Compute E-Step weights
  weights = w_m_binaryY(mstar_matrix, outcome_matrix,
                        pistar_matrix = conditional_probabilities,
                        pi_matrix = probabilities,
                        p_yi_m0, p_yi_m1,
                        sample_size, n_cat)

  # Estimate gamma parameters using weighted logistic regression
  ## Weights from E-Step (split by value of latent mediator, m)
  ## Outcome is the observed mediator
  Mstar01 = mstar_matrix[,1]
  fit.gamma1 <- suppressWarnings( stats::glm(Mstar01 ~ . + 0, as.data.frame(Z),
                           weights = weights[,1],
                           family = "binomial"(link = "logit")) )
  gamma1_new <- unname(coefficients(fit.gamma1))

  fit.gamma2 <- suppressWarnings( stats::glm(Mstar01 ~ . + 0, as.data.frame(Z),
                           weights = weights[,2],
                           family = "binomial"(link = "logit")) )
  gamma2_new <- unname(coefficients(fit.gamma2))

  # Estimate beta parameters using logistic regression
  ## Outcome is the E-Step weight
  fit.beta <- suppressWarnings( stats::glm(weights[,1] ~ . + 0, as.data.frame(design_matrix),
                         family = stats::binomial()) )
  beta_new <- unname(coefficients(fit.beta))

  # Estimate theta parameters using a weighted logistic regression
  ## Duplicate the data, half has m = 0 and half has m = 1
  ## Weights from E-Step (split by value of latent mediator, m)
  ## Outcome is y
  x_vector = X[,2]
  data1 = data.frame(x = x_vector, c = c_matrix, m = 0,
                     w = weights[,2], y = obs_outcome)
  data2 = data.frame(x = x_vector, c = c_matrix, m = 1,
                     w = weights[,1], y = obs_outcome)
  doubled_data_theta = rbind(data1, data2)
  
  theta_update = glm(y ~ x + m + c + x*m, weights = w, data = doubled_data_theta,
                     family = "binomial"(link = "logit"))

  theta_new <- unname(coef(theta_update))
  
  # Save new parameters
  param_new = c(beta_new, gamma1_new, gamma2_new, theta_new)
  param_current = param_new
  return(param_new)

}

