EM_function_normalY_XM <- function(param_current,
                                   obs_mediator, obs_outcome,
                                   X, Z, c_matrix,
                                   sample_size, n_cat){

  # Separate current parameters into beta, gamma, theta, sigma vectors
  beta_current = matrix(param_current[1:3], ncol = 1)
  gamma_current = matrix(c(param_current[4:7]),
                         ncol = n_cat, byrow = FALSE)
  theta_current = matrix(c(param_current[8:12]),
                         ncol = 1)
  sigma_current = matrix(param_current[13])

  # Create design matrix for true mediator model
  design_matrix = cbind(X, c_matrix)
  
  # Compute probability of each latent mediator value
  probabilities = pi_compute(beta_current, design_matrix, sample_size, n_cat)
  
  # Compute probability of observed mediator, given latent mediator
  conditional_probabilities = pistar_compute(gamma_current, Z, sample_size, n_cat)
  
  # Compute likelihood value of Y based on x, m, c, theta, and sigma
  p_yi_m0 = p_normal_yi_XM(theta_v = theta_current,
                           sigma_v = sigma_current,
                           m = 0, x = X[,2], c = c_matrix,
                           outcome = obs_outcome,
                           sample_size, n_cat)
  p_yi_m1 = p_normal_yi_XM(theta_v = theta_current,
                           sigma_v = sigma_current,
                           m = 1, x = X[,2], c = c_matrix,
                           outcome = obs_outcome,
                           sample_size, n_cat)

  # Create a matrix of observed mediator variables using dummy coding
  mstar_matrix = matrix(c(ifelse(obs_mediator == 1, 1, 0), 
                          ifelse(obs_mediator == 2, 1, 0)),
                        nrow = sample_size, byrow = FALSE)
  
  # Compute E-Step weights
  weights = w_m_normalY(mstar_matrix,
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

  # Solve for theta parameters using a system of equations
  a_row1 <- c(1,
              mean(X[,2]),
              mean(c_matrix),
              mean(weights[,1]),
              mean(X[,2] * weights[,1]))
  
  a_row2 <- c(sum(X[,2]) / sum(X[,2]^2),
              1,
              sum(X[,2] * c_matrix) / sum(X[,2]^2),
              sum(X[,2] * weights[,1]) / sum(X[,2]^2),
              sum(X[,2] * X[,2] * weights[,1]) / sum(X[,2]^2))
  
  a_row3 <- c(sum(c_matrix) / sum(c_matrix^2),
              sum(c_matrix * X[,2]) / sum(c_matrix^2),
              1,
              sum(c_matrix * weights[,1]) / sum(c_matrix^2),
              sum(c_matrix * X[,2] * weights[,1]) / sum(c_matrix^2))
  
  a_row4 <- c(1,
              sum(X[,2] * weights[,1]) / sum(weights[,1]),
              sum(c_matrix * weights[,1]) / sum(weights[,1]),
              1,
              sum(X[,2] * weights[,1] * weights[,1]) / sum(weights[,1]))
  
  a_row5 <- c(sum(X[,2] * weights[,1]) / sum((X[,2]^2 * weights[,1])),
              sum(X[,2] * weights[,1] * X[,2]) / sum((X[,2]^2 * weights[,1])),
              sum(X[,2] * weights[,1] * c_matrix) / sum((X[,2]^2 * weights[,1])),
              sum(X[,2] * weights[,1]) / sum((X[,2]^2 * weights[,1])),
              1)
  
  A = matrix(c(a_row1, a_row2, a_row3, a_row4, a_row5), byrow = TRUE, nrow = 5)
  B = matrix(c(mean(obs_outcome),
               sum(X[,2] * obs_outcome) / sum(X[,2]^2),
               sum(c_matrix * obs_outcome) / sum(c_matrix^2),
               sum(weights[,1] * obs_outcome) / sum(weights[,1]),
               sum(X[,2] * weights[,1] * obs_outcome) / sum((X[,2]^2 * weights[,1]))),
             ncol = 1)
  
  theta_update <- solve(A, B)
  
  # Compute sigma estimate
  sigma_update <- (1 / sample_size) * sum(
    ((obs_outcome -
        theta_update[1] - theta_update[2] * X[,2]
      - theta_update[3] * c_matrix) ^ 2)
    - 2 * theta_update[4] * weights[,1] *
      (obs_outcome - theta_update[1] - theta_update[2] * X[,2]
                                            - theta_update[3] * c_matrix)
     + (theta_update[4] ^ 2 * weights[,1]) +
      (-2 * theta_update[5] * weights[,1] * X[,2] *
         (obs_outcome - theta_update[1] - theta_update[2] * X[,2]
          - theta_update[3] * c_matrix)) +
      (2 * theta_update[4] * weights[,1] * theta_update[5] * X[,2]) +
      (weights[,1] * theta_update[5]^2 * X[,2]^2)
  )
  
  # Reorder theta estimates
  theta_new <- theta_update[c(1,2,4,3,5)]
  
  sigma_new <- sigma_update
  
  # Save new parameters
  param_new = c(beta_new, gamma1_new, gamma2_new, theta_new, sigma_new)
  param_current = param_new
  return(param_new)

}




