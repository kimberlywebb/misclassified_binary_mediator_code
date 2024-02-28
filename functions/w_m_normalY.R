w_m_normalY <- function(mstar_matrix, # Observed mediator matrix
                        pistar_matrix, # Probability of observed mediator given latent mediator
                        pi_matrix, # Mediator probabilities
                        p_yi_m0, p_yi_m1, # Likelihood value of Y for given M
                        sample_size, n_cat){
  
  y_m1_mstar1_frac1 = p_yi_m1 * mstar_matrix[,1] * pistar_matrix[1:sample_size, 1] *
    pi_matrix[,1]
  y_m2_mstar1_frac1 = p_yi_m0 * mstar_matrix[,1] * pistar_matrix[1:sample_size, 2] *
    pi_matrix[,2]
  
  y_m1_mstar2_frac1 = p_yi_m1 * mstar_matrix[,2] *
    pistar_matrix[(sample_size + 1):(n_cat * sample_size), 1] * pi_matrix[,1]
  y_m2_mstar2_frac1 = p_yi_m0 * mstar_matrix[,2] *
    pistar_matrix[(sample_size + 1):(n_cat * sample_size), 2] * pi_matrix[,2]
  
  mstar1_frac2 = (p_yi_m1 * pistar_matrix[1:sample_size, 1] * pi_matrix[,1]) +
    (p_yi_m0 * pistar_matrix[1:sample_size, 2] * pi_matrix[,2])
  
  mstar2_frac2 = (p_yi_m0 *
                    pistar_matrix[(sample_size + 1):(n_cat * sample_size), 2]
                  * pi_matrix[,2]) +
    (p_yi_m1 * pistar_matrix[(sample_size + 1):(n_cat * sample_size), 1] * pi_matrix[,1])
  
  m1_term = (y_m1_mstar1_frac1 / mstar1_frac2) + (y_m1_mstar2_frac1 / mstar2_frac2)
  m2_term = (y_m2_mstar1_frac1 / mstar1_frac2) + (y_m2_mstar2_frac1 / mstar2_frac2)
  
  weight_m = matrix(c(m1_term, m2_term), nrow = sample_size, byrow = FALSE)
  
  return(weight_m)
}
