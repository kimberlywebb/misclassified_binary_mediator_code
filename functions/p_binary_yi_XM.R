p_binary_yi_XM <- function(theta_v, m, x, c, outcome,
                        sample_size, n_cat){
  
  # Compute likelihood value of Y based on x, m, c, and theta
  model_y = theta_v[1] + theta_v[2] * x + theta_v[3] * m + theta_v[4] * c +
    theta_v[5] * x * m
  exp_model_y = exp(model_y)
  p_result = exp_model_y / (1 + exp_model_y)
  p_matrix = matrix(c(p_result, 1 - p_result),
                     ncol = n_cat, nrow = sample_size,
                     byrow = FALSE)
  
  return(p_matrix)
}
