# Function to generate outcomes Y from a Normal distribution based 
## on regression parameters theta and predictor matrices X, C, and M
generate_normalY <- function(theta, # Parameter values
                             true_M, # True mediator
                             x_vector, # Predictor of interest
                             c_vector, # Covariates
                             n, # Sample size
                             n_cat){

  m_indicator <- ifelse(true_M == 1, 1, 0) # M is coded 1/2, make 1/0
  
  # Combine all predictor variables into single matrix
  data_matrix <- matrix(c(rep(1, n), x_vector, m_indicator, c_vector),
                        ncol = 4, byrow = FALSE)
  
  # Generate mean and Normal errors
  additive_term <- data_matrix %*% theta
  additive_term_error <- rnorm(n) # Errors generated with SD, sigma = 1
  
  # Return value of Y
  return_normalY <- additive_term + additive_term_error

  return(return_normalY)
}
