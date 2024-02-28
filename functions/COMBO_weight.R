COMBO_weight <- function(ystar_matrix, pistar_matrix, pi_matrix, sample_size, n_cat){

  pi_ij_1_allj_repped = do.call(rbind,
                                list(pi_matrix, pi_matrix))
  pistar_pi_1 = pistar_matrix * pi_ij_1_allj_repped
  suml_pistar_pi_1 = rowSums(pistar_pi_1)

  suml_pistar_pi_denominator_1 <- matrix(rep(suml_pistar_pi_1, n_cat),
                                         nrow = n_cat * sample_size,
                                         byrow = FALSE)
  obs_Y_matrix_repped_1 <- matrix(rep(c(ystar_matrix), each = n_cat),
                                  nrow = n_cat * sample_size, byrow = TRUE)
  weight_not_summed_1 <- obs_Y_matrix_repped_1 * (pistar_pi_1 / ifelse(suml_pistar_pi_denominator_1 == 0, .00000001, suml_pistar_pi_denominator_1))

  weight_1 <- matrix(NA, nrow = sample_size, ncol = n_cat)
  for(i in 1:sample_size){
    for(j in 1:n_cat){
      k_set = c(i, sample_size + i)
      sum_terms = weight_not_summed_1[c(k_set), j]
      weight_1[i, j] = sum(sum_terms)
    }
  }

  return(weight_1)
}
