
# function to compute 4 trait moment indexes informative of distribution patterns
# abundance and trait are vectors of equal length and should have similar orders

trait_moment <- function(abundance, trait) {

  # CWM: community weighted mean using raw values
  cwm_i <- sum(trait * abundance) / sum(abundance)
  # CWV: community variance
  cwv_i <- sum(abundance * (trait - cwm_i)^2)  / sum(abundance)
  # CWS: community skewness
  cws_i <- sum(abundance * (trait - cwm_i)^3) / sum(abundance) / sqrt(cwv_i)^3
  # CWK: community kurtosis
  cwk_i <- sum(abundance * (trait - cwm_i)^4) / sum(abundance) / cwv_i^2 - 3
  
  results <- c(cwm_i,cwv_i,cws_i,cwk_i)
  names(results) <- c('CWM','CWV','CWS','CWK')
  
  return(results)
}

