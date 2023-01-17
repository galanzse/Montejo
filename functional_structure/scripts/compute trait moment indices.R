
library(factoextra)
library(vegan)
# library(Weighted.Desc.Stat)
source('C:/Users/user/OneDrive/TESIS Y PUBLICACIONES/COTTBUS/montejo/functional_structure/scripts/import data.R')
source('C:/Users/user/OneDrive/TESIS Y PUBLICACIONES/COTTBUS/montejo/functional_structure/scripts/trait_moment function.R')


# IMPUTATION ####
# imputation of species x trait combinations in sites with 0 obs

traits[traits$species=='a.his' & traits$forest=='Quercus pyrenaica',c('RDMC','SRL','TMDr')] <- traits[traits$species=='a.his' & traits$forest=='Shrubland',c('RDMC','SRL','TMDr')]
traits[traits$species=='c.sco' & traits$forest=='Q. pyrenaica-F. sylvatica',c('RDMC','SRL','TMDr')] <- traits[traits$species=='c.sco' & traits$forest=='Quercus pyrenaica',c('RDMC','SRL','TMDr')]
traits[traits$species=='q.pyr' & traits$forest=='Shrubland','LDMC'][1] <- traits %>% filter(species=='q.pyr') %>% dplyr::select(LDMC) %>% colMeans(na.rm=T)
# check
# aggregate(.~ species + forest, traits, mean)  %>% view()

root_traits$root_C[root_traits$species=='c.mon' & root_traits$forest=='Quercus pyrenaica'] <- root_traits %>% filter(species=='c.mon') %>% dplyr::select(root_C) %>% colMeans(na.rm=T)
root_traits[root_traits$species=='c.sco' & root_traits$forest=='Quercus pyrenaica',c('root_C','root_N')] <- root_traits %>% filter(species=='c.sco') %>% dplyr::select(root_C, root_N) %>% colMeans(na.rm=T)
root_traits[root_traits$species=='e.aus' & root_traits$forest=='Quercus pyrenaica',c('root_C','root_N')] <- root_traits %>% filter(species=='e.arb') %>% dplyr::select(root_C, root_N) %>% colMeans(na.rm=T) # e.aus needs to be imputed from a congeneric species
root_traits$root_C[root_traits$species=='g.flo' & root_traits$forest=='Quercus pyrenaica'] <- root_traits %>% filter(species=='g.flo') %>% dplyr::select(root_C) %>% colMeans(na.rm=T)
root_traits$root_C[root_traits$species=='i.aqu' & root_traits$forest=='Quercus pyrenaica'] <- root_traits %>% filter(species=='i.aqu') %>% dplyr::select(root_C) %>% colMeans(na.rm=T)
root_traits$root_C[root_traits$species=='l.sto' & root_traits$forest=='Shrubland'] <- root_traits %>% filter(species=='l.sto') %>% dplyr::select(root_C) %>% colMeans(na.rm=T)
root_traits$root_C[root_traits$species=='p.avi' & root_traits$forest=='Quercus pyrenaica'] <- root_traits %>% filter(species=='p.avi') %>% dplyr::select(root_C) %>% colMeans(na.rm=T)


# TRAITS ####
#  measures of kurtosis require at least four values (Hulme & Bernard-Verdier, 2018 JVS)

# scale traits
v_traits <- colnames(traits)[3:11]
s_traits <- traits
s_traits[,v_traits] <- apply(s_traits[,v_traits], 2, scale)
trait_dist_weighted <- list() # abundances
trait_dist_unweighted <- list() # presences

# I am not going to use species x site averages, instead I will use the raw data
# and weight NAs by abundances

for (t in 1:length(v_traits)) {
  
  sites$trait <- v_traits[t]
  sites$CWM <- NA
  sites$CWV <- NA
  sites$CWS <- NA
  sites$CWK <- NA
  
  sites_ab <- sites
  sites_pr <- sites_ab
  
  for (s in 1:nrow(sites)) {
  
  pl1 <- as.character(sites[s,'plot']) # plot
  ab1 <- abundances[pl1,]
  
  sp1 <- which(ab1>0) %>% names() # species names
  ab1 <- ab1[sp1] # abundaces
  pr1 <- ab1; pr1[] <- 1 # presences

  # filter trait database
  df1 <- s_traits %>% filter(forest==sites[s,'forest'] & species%in%sp1) %>%
    dplyr::select(species,v_traits[t]) %>% na.omit()
  colnames(df1) <- c("species","trait")
  
  # remove species without trait data
  ab1 <- ab1[unique(df1$species)]

  # weight observations by missing data
  ab2 <- ab1 / table(df1$species)[names(ab1)]
  ab2 <- as.data.frame(ab2); colnames(ab2) <- c('species','abundance')
  ab2 <- merge(df1, ab2, by='species')
  
  pr2 <- pr1 / table(df1$species)[names(pr1)]
  pr2 <- as.data.frame(pr2); colnames(pr2) <- c('species','presence')
  pr2 <- merge(df1, pr2, by='species')
  
  # run function
  sites_ab[s,c('CWM','CWV','CWS','CWK')] <- trait_moment(abundance=ab2$abundance, trait=ab2$trait)
  sites_pr[s,c('CWM','CWV','CWS','CWK')] <- trait_moment(abundance=pr2$presence, trait=pr2$trait)
  
  }
  
  # save results
  trait_dist_weighted[[v_traits[t]]] <- sites_ab
  trait_dist_unweighted[[v_traits[t]]] <- sites_pr
  
  # progress
  print(v_traits[t])
  
}

trait_dist_weighted <- do.call('rbind', trait_dist_weighted)
rownames(trait_dist_weighted) <- 1:nrow(trait_dist_weighted)

trait_dist_unweighted <- do.call('rbind', trait_dist_unweighted)
rownames(trait_dist_unweighted) <- 1:nrow(trait_dist_unweighted)


# HV_thick ####

# for HubVal and Thickness there is data

# scale traits
v_traits <- colnames(HV_thick)[3:4]
s_HV_thick <- HV_thick
s_HV_thick[,v_traits] <- apply(HV_thick[,v_traits], 2, scale)
HVT_dist_weighted <- list() # abundances
HVT_dist_unweighted <- list() # presences

for (t in 1:length(v_traits)) {
  
  sites$trait <- v_traits[t]
  sites$CWM <- NA
  sites$CWV <- NA
  sites$CWS <- NA
  sites$CWK <- NA
  
  sites_ab <- sites
  sites_pr <- sites_ab
  
  c(113,83,84,96,107,117,118)
  
  for (s in 1:nrow(sites)) {
    
    pl1 <- as.character(sites[s,'plot']) # plot
    ab1 <- abundances[pl1,]
    
    sp1 <- which(ab1>0) %>% names() # species names
    ab1 <- ab1[sp1] # abundances
    pr1 <- ab1; pr1[] <- 1 # presences
    
    # filter trait database
    df1 <- s_HV_thick %>% filter(forest==sites[s,'forest'] & species%in%sp1) %>% dplyr::select(species,v_traits[t]) %>% na.omit()
    df2 <- s_HV_thick %>% filter(species %in% setdiff(sp1, df1$species)) %>% dplyr::select(species,v_traits[t]) %>% na.omit()
    df1 <- rbind(df1, df2); colnames(df1) <- c("species","trait")
    
    # remove species without trait data
    ab1 <- ab1[unique(df1$species)]
    
    # weight observations by missing data
    ab2 <- ab1 / table(df1$species)[names(ab1)]
    ab2 <- as.data.frame(ab2); colnames(ab2) <- c('species','abundance')
    ab2 <- merge(df1, ab2, by='species')
    
    pr2 <- pr1 / table(df1$species)[names(pr1)]
    pr2 <- as.data.frame(pr2); colnames(pr2) <- c('species','presence')
    pr2 <- merge(df1, pr2, by='species')
    
    # run function
    sites_ab[s,c('CWM','CWV','CWS','CWK')] <- trait_moment(abundance=ab2$abundance, trait=ab2$trait)
    sites_pr[s,c('CWM','CWV','CWS','CWK')] <- trait_moment(abundance=pr2$presence, trait=pr2$trait)
    
  }
  
  # save results
  HVT_dist_weighted[[v_traits[t]]] <- sites_ab
  HVT_dist_unweighted[[v_traits[t]]] <- sites_pr
  
  # progress
  print(v_traits[t])
  
}

HVT_dist_weighted <- do.call('rbind', HVT_dist_weighted)
rownames(HVT_dist_weighted) <- 1:nrow(HVT_dist_weighted)

HVT_dist_unweighted <- do.call('rbind', HVT_dist_unweighted)
rownames(HVT_dist_unweighted) <- 1:nrow(HVT_dist_unweighted)


# ROOT TRAITS ####

# for roots there is only one mean value of species x site

# scale traits
v_traits <- colnames(root_traits)[3:4]
s_root_traits <- root_traits
s_root_traits[,v_traits] <- apply(s_root_traits[,v_traits], 2, scale)
root_dist_weighted <- list() # abundances
root_dist_unweighted <- list() # presences

for (t in 1:length(v_traits)) {
  
  sites$trait <- v_traits[t]
  sites$CWM <- NA
  sites$CWV <- NA
  sites$CWS <- NA
  sites$CWK <- NA
  
  sites_ab <- sites
  sites_pr <- sites_ab
  
  for (s in 1:nrow(sites)) {
    
    pl1 <- as.character(sites[s,'plot']) # plot
    ab1 <- abundances[pl1,]
    
    sp1 <- which(ab1>0) %>% names() # species names
    ab1 <- ab1[sp1] # abundaces
    pr1 <- ab1; pr1[] <- 1 # presences
    
    # filter trait database
    df1 <- s_root_traits %>% filter(forest==sites[s,'forest'] & species%in%sp1) %>%
      dplyr::select(species,v_traits[t]) %>% na.omit()
    colnames(df1) <- c("species","trait")
    
    # remove species without trait data
    ab1 <- ab1[unique(df1$species)]
    
    # weight observations by missing data
    ab1 <- as.data.frame(ab1); ab1$species <- rownames(ab1)
    colnames(ab1) <- c('abundance','species')
    ab1 <- merge(df1, ab1, by='species')
    
    pr1 <- as.data.frame(pr1); pr1$species <- rownames(pr1)
    colnames(pr1) <- c('presence','species')
    pr1 <- merge(df1, pr1, by='species')
    
    # run function
    sites_ab[s,c('CWM','CWV','CWS','CWK')] <- trait_moment(abundance=ab1$abundance, trait=ab1$trait)
    sites_pr[s,c('CWM','CWV','CWS','CWK')] <- trait_moment(abundance=pr1$presence, trait=pr1$trait)
    
  }
  
  # save results
  root_dist_weighted[[v_traits[t]]] <- sites_ab
  root_dist_unweighted[[v_traits[t]]] <- sites_pr
  
  # progress
  print(v_traits[t])
  
}

root_dist_weighted <- do.call('rbind', root_dist_weighted)
rownames(root_dist_weighted) <- 1:nrow(root_dist_weighted)

root_dist_unweighted <- do.call('rbind', root_dist_unweighted)
rownames(root_dist_unweighted) <- 1:nrow(root_dist_unweighted)


# EXPLORATORY ####

moments_weighted <- rbind(trait_dist_weighted, HVT_dist_weighted, root_dist_weighted)
write.table(moments_weighted, 'results/moments_weighted.txt')
moments_unweighted <- rbind(trait_dist_unweighted, HVT_dist_unweighted, root_dist_unweighted)
write.table(moments_unweighted, 'results/moments_unweighted.txt')

# correlation between weighted and unweighted
par(mfrow=c(2,2), mar=c(5,5,5,5))
plot(moments_weighted$CWM ~ moments_unweighted$CWM, ylab='Weighted', xlab='Unweighted', main='CWM'); abline(a=0, b=1, col='red', lwd=2)
plot(moments_weighted$CWV ~ moments_unweighted$CWV, ylab='Weighted', xlab='Unweighted', main='CWV'); abline(a=0, b=1, col='red', lwd=2)
plot(moments_weighted$CWS ~ moments_unweighted$CWS, ylab='Weighted', xlab='Unweighted', main='CWS'); abline(a=0, b=1, col='red', lwd=2)
plot(moments_weighted$CWK ~ moments_unweighted$CWK, ylab='Weighted', xlab='Unweighted', main='CWK'); abline(a=0, b=1, col='red', lwd=2)

# plots
moments_weighted$forest <- as.factor(moments_weighted$forest) # organize factor variables
moments_weighted$forest <- factor(moments_weighted$forest, levels=c("Fagus sylvatica","Quercus petraea","Q. pyrenaica-F. sylvatica","Mixed forest","Quercus pyrenaica","Shrubland"))
moments_weighted <- na.omit(moments_weighted)

ggplot(aes(y=CWM, x=forest), data=moments_weighted) +
  geom_boxplot() + facet_wrap(~trait, scales='free_y') +
  theme_classic() + xlab('') +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5))

ggplot(aes(y=CWV, x=forest), data=moments_weighted) +
  geom_boxplot() + facet_wrap(~trait, scales='free_y') +
  theme_classic() + xlab('') +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5))

ggplot(aes(y=CWS, x=forest), data=moments_weighted) +
  geom_boxplot() + facet_wrap(~trait, scales='free_y') +
  theme_classic() + xlab('') +
  geom_hline(yintercept=0) +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5))

ggplot(aes(y=CWK, x=forest), data=moments_weighted) +
  geom_boxplot() + facet_wrap(~trait, scales='free_y') +
  theme_classic() + xlab('') +
  geom_hline(yintercept=0) +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5))


moments_unweighted$forest <- as.factor(moments_unweighted$forest) # organize factor variables
moments_unweighted$forest <- factor(moments_unweighted$forest, levels=c("Fagus sylvatica","Quercus petraea","Q. pyrenaica-F. sylvatica","Mixed forest","Quercus pyrenaica","Shrubland"))
moments_unweighted <- na.omit(moments_unweighted)

ggplot(aes(y=CWM, x=forest), data=moments_unweighted) +
  geom_boxplot() + facet_wrap(~trait, scales='free_y') +
  theme_classic() + xlab('') +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5))

ggplot(aes(y=CWV, x=forest), data=moments_unweighted) +
  geom_boxplot() + facet_wrap(~trait, scales='free_y') +
  theme_classic() + xlab('') +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5))

ggplot(aes(y=CWS, x=forest), data=moments_unweighted) +
  geom_boxplot() + facet_wrap(~trait, scales='free_y') +
  theme_classic() + xlab('') +
  geom_hline(yintercept=0) +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5))

ggplot(aes(y=CWK, x=forest), data=moments_unweighted) +
  geom_boxplot() + facet_wrap(~trait, scales='free_y') +
  theme_classic() + xlab('') +
  geom_hline(yintercept=0) +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5))
