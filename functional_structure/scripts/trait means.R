
# generate a dataset of trait means per site
# for each dataset, I will impute species x trait combinations in sites with 0 obs, and average observations per species and site

source('C:/Users/user/OneDrive/TESIS Y PUBLICACIONES/COTTBUS/montejo/functional_structure/scripts/import data.R')


# TRAITS
traits[traits$species=='a.his' & traits$forest=='Quercus pyrenaica',c('RDMC','SRL','TMDr')] <- traits[traits$species=='a.his' & traits$forest=='Shrubland',c('RDMC','SRL','TMDr')]
traits[traits$species=='c.sco' & traits$forest=='Q. pyrenaica-F. sylvatica',c('RDMC','SRL','TMDr')] <- traits[traits$species=='c.sco' & traits$forest=='Quercus pyrenaica',c('RDMC','SRL','TMDr')]
traits[traits$species=='q.pyr' & traits$forest=='Shrubland','LDMC'][4] <- traits %>% filter(species=='q.pyr') %>% dplyr::select(LDMC) %>% colMeans(na.rm=T)

# averages traits x species using a loop because aggregate and summarise miss some combinations
v_traits <- colnames(traits)[3:11]
trait_means <- unique(traits[,c('species','forest')])
trait_means$LDMC <- NA
trait_means$SLA <- NA
trait_means$SDMC <- NA
trait_means$RDMC <- NA
trait_means$SRL <- NA
trait_means$TMDr <- NA
trait_means$leaf_d13C <- NA
trait_means$leaf_C <- NA
trait_means$leaf_N <- NA

for (i in 1:nrow(trait_means)) {
  trait_means[i,v_traits] <- traits %>% filter(species==trait_means$species[i] & forest==trait_means$forest[i]) %>% dplyr::select(all_of(v_traits)) %>% colMeans(na.rm=T)
}
trait_means <- na.omit(trait_means) # remove malus


# ROOT
root_traits$root_C[root_traits$species=='c.mon' & root_traits$forest=='Quercus pyrenaica'] <- root_traits %>% filter(species=='c.mon') %>% dplyr::select(root_C) %>% colMeans(na.rm=T)
root_traits[root_traits$species=='c.sco' & root_traits$forest=='Quercus pyrenaica',c('root_C','root_N')] <- root_traits %>% filter(species=='c.sco') %>% dplyr::select(root_C, root_N) %>% colMeans(na.rm=T)
root_traits[root_traits$species=='e.aus' & root_traits$forest=='Quercus pyrenaica',c('root_C','root_N')] <- root_traits %>% filter(species=='e.arb') %>% dplyr::select(root_C, root_N) %>% colMeans(na.rm=T) # e.aus needs to be imputed from a congeneric species
root_traits$root_C[root_traits$species=='g.flo' & root_traits$forest=='Quercus pyrenaica'] <- root_traits %>% filter(species=='g.flo') %>% dplyr::select(root_C) %>% colMeans(na.rm=T)
root_traits$root_C[root_traits$species=='i.aqu' & root_traits$forest=='Quercus pyrenaica'] <- root_traits %>% filter(species=='i.aqu') %>% dplyr::select(root_C) %>% colMeans(na.rm=T)
root_traits$root_C[root_traits$species=='l.sto' & root_traits$forest=='Shrubland'] <- root_traits %>% filter(species=='l.sto') %>% dplyr::select(root_C) %>% colMeans(na.rm=T)
root_traits$root_C[root_traits$species=='p.avi' & root_traits$forest=='Quercus pyrenaica'] <- root_traits %>% filter(species=='p.avi') %>% dplyr::select(root_C) %>% colMeans(na.rm=T)

# averages traits x species using a loop because aggregate and summarise miss some combinations
v_traits <- colnames(root_traits)[3:4]
root_means <- unique(root_traits[,c('species','forest')])
root_means$root_C <- NA
root_means$root_N <- NA

for (i in 1:nrow(root_means)) {
  root_means[i,v_traits] <- root_traits %>% filter(species==root_means$species[i] & forest==root_means$forest[i]) %>% dplyr::select(all_of(v_traits)) %>% colMeans(na.rm=T)
}
root_means <- na.omit(root_means) # remove malus


# MERGE AND ADD HUBVAL

trait_means
root_means

# setdiff(paste(trait_means$species, trait_means$forest), paste(root_means$species, root_means$forest))
trait_means <- merge(trait_means, root_means, by=c('species','forest'), all.x=T)

trait_means$Thickness <- NA
trait_means$HubVal <- NA
v_traits <- colnames(HV_thick)[3:4]

for (i in 1:nrow(trait_means)) {
  sp <- trait_means$species[i]; fo <- trait_means$forest[i]
  if (paste(sp,fo) %in% paste(HV_thick$species,HV_thick$forest)) {
      trait_means[i,v_traits] <- HV_thick %>% filter(species%in%sp & forest%in%fo) %>% dplyr::select(all_of(v_traits)) %>% colMeans(na.rm=T)
  } else {
      trait_means[i,v_traits] <- HV_thick %>% filter(species==trait_means$species[i]) %>% dplyr::select(all_of(v_traits)) %>% colMeans(na.rm=T)
  }
}


# IMPUTE root_C and root_N USING SPECIES MEANS
colSums(is.na(trait_means))

trait_means[trait_means$species=='a.his' & trait_means$forest=='Quercus pyrenaica',c('root_C','root_N')] <- trait_means %>% filter(species=='a.his') %>% dplyr::select(root_C, root_N) %>% colMeans(na.rm=T)

trait_means[trait_means$species=='c.sco' & trait_means$forest=='Q. pyrenaica-F. sylvatica',c('root_C','root_N')] <- trait_means %>% filter(species=='c.sco') %>% dplyr::select(root_C, root_N) %>% colMeans(na.rm=T)

trait_means[trait_means$species=='i.aqu' & trait_means$forest=='Quercus petraea',c('root_C','root_N')] <- trait_means %>% filter(species=='i.aqu') %>% dplyr::select(root_C, root_N) %>% colMeans(na.rm=T)

trait_means[trait_means$species=='s.auc' & trait_means$forest=='Quercus petraea',c('root_C','root_N')] <- trait_means %>% filter(species=='s.auc') %>% dplyr::select(root_C, root_N) %>% colMeans(na.rm=T)


# SAVE
write.table(trait_means, 'results/trait_means.txt')
