p_normal_yi <- function(theta_v, sigma_v, # Parameter vectors
                        m, x, c, outcome, # Predictors and outcome
                        sample_size, n_cat){
  
  # Compute likelihood value of Y based on x, m, c, theta, and sigma
  model_y = theta_v[1] + theta_v[2] * x + theta_v[3] * m + theta_v[4] * c
  
  residual_term = outcome - model_y
  
  term1 = 1 / sqrt(2 * pi * c(sigma_v ^ 2))
  
  exp_term = exp(-1 * residual_term^2 * (1 / c(2 * sigma_v^2)))
  
  result = term1 * exp_term
  
  return(result)
}
