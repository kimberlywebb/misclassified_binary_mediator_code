w_m_binaryY <- function(mstar_matrix, # Observed mediator matrix
                        outcome_matrix, # Outcome matrix
                        pistar_matrix, # Probability of observed mediator given latent mediator
                        pi_matrix, # Mediator probabilities
                        p_yi_m0, p_yi_m1, # Likelihood value of Y for given M
                        sample_size, n_cat){
  
  y1_m1_mstar1_frac1 = p_yi_m1[,1] * outcome_matrix[,1] * mstar_matrix[,1] *
    pistar_matrix[1:sample_size, 1] * pi_matrix[,1]
  y1_m2_mstar1_frac1 = p_yi_m0[,1] * outcome_matrix[,1] * mstar_matrix[,1] *
    pistar_matrix[1:sample_size, 2] * pi_matrix[,2]
  
  y1_m1_mstar2_frac1 = p_yi_m1[,1] * outcome_matrix[,1] * mstar_matrix[,2] *
    pistar_matrix[(sample_size + 1):(n_cat * sample_size), 1] * pi_matrix[,1]
  y1_m2_mstar2_frac1 = p_yi_m0[,1] * outcome_matrix[,1] * mstar_matrix[,2] *
    pistar_matrix[(sample_size + 1):(n_cat * sample_size), 2] * pi_matrix[,2]
  
  y2_m1_mstar1_frac1 = p_yi_m1[,2] * outcome_matrix[,2] * mstar_matrix[,1] *
    pistar_matrix[1:sample_size, 1] * pi_matrix[,1]
  y2_m2_mstar1_frac1 = p_yi_m0[,2] * outcome_matrix[,2] * mstar_matrix[,1] *
    pistar_matrix[1:sample_size, 2] * pi_matrix[,2]
  
  y2_m1_mstar2_frac1 = p_yi_m1[,2] * outcome_matrix[,2] * mstar_matrix[,2] *
    pistar_matrix[(sample_size + 1):(n_cat * sample_size), 1] * pi_matrix[,1]
  y2_m2_mstar2_frac1 = p_yi_m0[,2] * outcome_matrix[,2] * mstar_matrix[,2] *
    pistar_matrix[(sample_size + 1):(n_cat * sample_size), 2] * pi_matrix[,2]
  
  y1_mstar1_frac2 = (p_yi_m1[,1] * pistar_matrix[1:sample_size, 1] * pi_matrix[,1]) +
    (p_yi_m0[,1] * pistar_matrix[1:sample_size, 2] * pi_matrix[,2])
  y2_mstar1_frac2 = (p_yi_m1[,2] * pistar_matrix[1:sample_size, 1] * pi_matrix[,1]) +
    (p_yi_m0[,2] * pistar_matrix[1:sample_size, 2] * pi_matrix[,2])
  
  y1_mstar2_frac2 = (p_yi_m0[,1] *
                    pistar_matrix[(sample_size + 1):(n_cat * sample_size), 2]
                  * pi_matrix[,2]) +
    (p_yi_m1[,1] * pistar_matrix[(sample_size + 1):(n_cat * sample_size), 1] * pi_matrix[,1])
  y2_mstar2_frac2 = (p_yi_m0[,2] *
                       pistar_matrix[(sample_size + 1):(n_cat * sample_size), 2]
                     * pi_matrix[,2]) +
    (p_yi_m1[,2] * pistar_matrix[(sample_size + 1):(n_cat * sample_size), 1] * pi_matrix[,1])
  
  
  m1_mstar1_term = (y1_m1_mstar1_frac1 / y1_mstar1_frac2) + (y2_m1_mstar1_frac1 / y2_mstar1_frac2)
  m1_mstar2_term = (y1_m1_mstar2_frac1 / y1_mstar2_frac2) + (y2_m1_mstar2_frac1 / y2_mstar2_frac2)
  
  m2_mstar1_term = (y1_m2_mstar1_frac1 / y1_mstar1_frac2) + (y2_m2_mstar1_frac1 / y2_mstar1_frac2)
  m2_mstar2_term = (y1_m2_mstar2_frac1 / y1_mstar2_frac2) + (y2_m2_mstar2_frac1 / y2_mstar2_frac2)
  
  m1_term = m1_mstar1_term + m1_mstar2_term
  m2_term = m2_mstar1_term + m2_mstar2_term
  
  weight_m = matrix(c(m1_term, m2_term), nrow = sample_size, byrow = FALSE)
  
  return(weight_m)
}
