
trait_distribution <- function(relative_abundance, trait_values) {
  # cwm: mean
  cwm_i <- sum(trait_values * relative_abundance)
  # cwv: variance
  cwv_i <- sum(relative_abundance * (trait_values - cwm_i)^2)
  # cws: skewness
  cws_i <- sum(relative_abundance * (trait_values - cwm_i)^3) / sqrt(cwv_i)^3
  # cwk: kurtosis
  cwk_i <- sum(relative_abundance * (trait_values - cwm_i)^4) / cwv_i^2 -3
  
  results <- c(cwm_i,cwv_i,cws_i,cwk_i)
  names(results) <- c('CWM','CWV','CWS','CWK')
  
  return(results)
}
