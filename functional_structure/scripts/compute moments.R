
# library(ggpubr)

source('scripts/PCA_traits.R')
source('scripts/trait_moment function.R')

# I will use the raw data instead of species x site averages and weight NAs by abundances
# Kurtosis require at least four values (Hulme & Bernard-Verdier, 2018 JVS)

trait_dist_weighted <- list()

for (t in c('Leaf_Economics','Root_Economics','HyArq')) {
  
  sites_ab <- sites
  sites_ab$trait <- t

  for (s in 1:nrow(sites)) {
  
  pl1 <- as.character(sites_ab[s,'plot']) # plot
  ab1 <- abundances[pl1,]
  
  sp1 <- which(ab1>0) %>% names() # species names
  ab1 <- ab1[sp1] # abundances

  # filter trait database
  df1 <- pca1 %>% filter(forest==sites_ab[s,'forest'] & species%in%sp1) %>%
    dplyr::select(species,t) %>% na.omit()
  colnames(df1) <- c("species","trait")
  
  # remove species without trait data
  ab1 <- ab1[unique(df1$species)]

  # weight observations by missing data
  ab2 <- ab1 / table(df1$species)[names(ab1)]
  ab2 <- as.data.frame(ab2); colnames(ab2) <- c('species','abundance')
  ab2 <- merge(df1, ab2, by='species')

  # run function
  sites_ab[s,c('CWM','CWV','CWS','CWK')] <- trait_moment(abundance=ab2$abundance, trait=ab2$trait)

  }
  
  # save results
  trait_dist_weighted[[t]] <- sites_ab

  # progress
  print(t)
  
}

trait_dist_weighted <- do.call('rbind', trait_dist_weighted)
rownames(trait_dist_weighted) <- 1:nrow(trait_dist_weighted)

# long format for plotting
lg_moments_weighted <- trait_dist_weighted %>% pivot_longer(4:7, names_to='moment')

lg_moments_weighted <- lg_moments_weighted %>% subset(moment!='CWV') # remove CWV

lg_moments_weighted$trait <- as.factor(lg_moments_weighted$trait)
levels(lg_moments_weighted$trait)
lg_moments_weighted$trait <- factor(lg_moments_weighted$trait, levels=c("Leaf_Economics","Root_Economics","HyArq"))

lg_moments_weighted$moment <- as.factor(lg_moments_weighted$moment)
levels(lg_moments_weighted$moment)
lg_moments_weighted$moment <- factor(lg_moments_weighted$moment, levels=c("CWM","CWS","CWK"))

lg_moments_weighted$forest <- as.factor(lg_moments_weighted$forest)
levels(lg_moments_weighted$forest)
levels(lg_moments_weighted$forest)[levels(lg_moments_weighted$forest)=='Q. pyrenaica-F. sylvatica'] <- 'Mixed 1'
levels(lg_moments_weighted$forest)[levels(lg_moments_weighted$forest)=='Mixed forest'] <- 'Mixed 2'
lg_moments_weighted$forest <- factor(lg_moments_weighted$forest, levels=c("Shrubland","Quercus pyrenaica","Mixed 2","Mixed 1","Quercus petraea","Fagus sylvatica"))

# save
write.table(lg_moments_weighted, 'results/lg_moments_weighted.txt')

# plots
lg_moments_weighted <- lg_moments_weighted %>% subset(moment!='CWV')
ggplot(aes(x=forest, y=value, color=forest), data=lg_moments_weighted) +
  geom_boxplot() +
    theme_classic() + xlab('') + ylab('') +
  facet_grid(moment~trait, scales='free_y') +
  geom_hline(aes(yintercept=0), linetype='dashed', data=lg_moments_weighted %>% subset(moment!='CWM')) +
  theme(legend.position='top', legend.title=element_blank(), legend.text=element_text(face="italic"),
        axis.text.x=element_blank())
