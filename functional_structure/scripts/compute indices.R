
source('C:/Users/user/OneDrive/TESIS Y PUBLICACIONES/COTTBUS/montejo/functional_structure/scripts/import data.R')
source('C:/Users/user/OneDrive/TESIS Y PUBLICACIONES/COTTBUS/montejo/functional_structure/scripts/trait_moment function.R')

# scale traits
traits[,v_traits] <- apply(traits[,v_traits], 2, scale)

# I am not going to use speceis x site averages, instead I will use the raw data
# and weight NAs by abundances

for (t in 1:length(v_traits)) {
  
  sites$trait <- v_traits[t]
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
  sites[s,c('CWM','CWV','CWS','CWK')] <- trait_moment(abundance=df1$abundance,
                                                      trait=df1$trait,
                                                      presence=F)
  
  }
  
  # save results
  trait_dist[[v_traits[t]]] <- sites
  # progress
  print(v_traits[t])
  
}

# trait_dist_weighted <- trait_dist
# trait_dist_unweighted <- trait_dist

trait_dist_weighted <- do.call('rbind', trait_dist_weighted)
rownames(trait_dist_weighted) <- 1:nrow(trait_dist_weighted)
write.table(trait_dist_weighted, 'results/trait_dist_weighted.txt')

trait_dist_unweighted <- do.call('rbind', trait_dist_unweighted)
rownames(trait_dist_unweighted) <- 1:nrow(trait_dist_unweighted)
write.table(trait_dist_unweighted, 'results/trait_dist_unweighted.txt')


par(mfrow=c(1,1))
plot( trait_dist_weighted$CWM[trait_dist_weighted$trait=='SLA'],
      trait_dist_unweighted$CWM[trait_dist_unweighted$trait=='SLA'])
abline(0,1)





