
source('C:/Users/user/OneDrive/TESIS Y PUBLICACIONES/COTTBUS/montejo/functional_structure/scripts/import data.R')
source('C:/Users/user/OneDrive/TESIS Y PUBLICACIONES/COTTBUS/montejo/functional_structure/scripts/trait_distribution function.R')


v_traits <- colnames(traits)[3:19]
trait_dist <- list()

for (t in 1:length(v_traits)) {
  
  sites$CWM <- NA
  sites$CWV <- NA
  sites$CWS <- NA
  sites$CWK <- NA
  
  for (s in 1:nrow(sites)) {
  
  pl1 <- as.character(sites[s,'plot']) # plot
  ab1 <- abundances[pl1,] # abundances
  sp1 <- which(ab1>0) %>% names() # species names
  ab1 <- ab1[sp1] # remove missing species

  # filter trait database
  df1 <- traits %>% filter(forest==sites[s,'forest'] & species%in%sp1) %>%
    dplyr::select(species,v_traits[t]) %>% na.omit()
  colnames(df1) <- c("species","trait")
  
  # remove species without trait data
  ab1 <- ab1[unique(df1$species)]
  # recalculate abundances
  ab1 <- ab1/sum(ab1)
  
  # weight observations by missing data
  ab2 <- ab1 / table(df1$species)[names(ab1)]
  ab2 <- as.data.frame(ab2); colnames(ab2) <- c('species','abundance')
  
  # merge dataframes and run function
  df1 <- merge(df1, ab2, by='species')
  sites[s,c('CWM','CWV','CWS','CWK')] <- trait_distribution(relative_abundance=df1$abundance,
                                                            trait_values=df1$trait)
  
  }
  
  # save results
  trait_dist[[v_traits[t]]] <- sites
  # progress
  print(v_traits[t])
}
