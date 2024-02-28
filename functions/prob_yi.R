prob_yi <- function(theta, # Parameter vectors
                    true_M, x_vector, c_vector, # Predictors
                    n, n_cat){
  
  # compute response probability based on m, x, c, and theta 
  m_indicator <- ifelse(true_M == 1, 1, 0)
  data_matrix <- matrix(c(rep(1, n), x_vector, m_indicator, c_vector),
                        ncol = 4, byrow = FALSE)
  exp_result <- exp(data_matrix %*% theta)

  exp_denominator <- 1 + exp_result
  
  rho_matrix <- matrix(c(exp_result / exp_denominator, 1 / exp_denominator),
                       ncol = 2, byrow = FALSE)

  return(rho_matrix)
}
