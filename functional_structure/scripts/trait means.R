
# generate a dataset of trait means per site
# for each dataset, I will impute species x trait combinations in sites with 0 obs, and average observations per species and site

source('C:/Users/user/OneDrive/PUBLICACIONES/BTU/montejo/functional_structure/scripts/import data.R')


# TRAITS
traits[traits$species=='a.his' & traits$forest=='Quercus pyrenaica',c('RDMC')] <- traits[traits$species=='a.his' & traits$forest=='Shrubland',c('RDMC')]
traits[traits$species=='c.sco' & traits$forest=='Q. pyrenaica-F. sylvatica',c('RDMC')] <- traits[traits$species=='c.sco' & traits$forest=='Quercus pyrenaica',c('RDMC')]
traits[traits$species=='q.pyr' & traits$forest=='Shrubland','LDMC'][4] <- traits %>% filter(species=='q.pyr') %>% dplyr::select(LDMC) %>% colMeans(na.rm=T)

# averages traits x species using a loop because aggregate and summarise miss some combinations
v_traits <- c("LDMC", "SLA", "SDMC", "RDMC", "SRA", "Rdi", "leaf_d13C", "leaf_CN", "leaf_N")
trait_means <- unique(traits[,c('species','forest')])
trait_means$LDMC <- NA
trait_means$SLA <- NA
trait_means$SDMC <- NA
trait_means$RDMC <- NA
trait_means$SRA <- NA
trait_means$Rdi <- NA
trait_means$leaf_d13C <- NA
trait_means$leaf_CN <- NA
trait_means$leaf_N <- NA

for (i in 1:nrow(trait_means)) {
  trait_means[i,v_traits] <- traits %>% filter(species==trait_means$species[i] & forest==trait_means$forest[i]) %>% dplyr::select(all_of(v_traits)) %>% colMeans(na.rm=T)
}
trait_means <- na.omit(trait_means) # remove malus


# ROOT
root_means <- trait_means[,1:4]
colnames(root_means) <- c("species", "forest", "root_N", "root_CN")

# averages traits x species using a loop because aggregate and summarise miss some combinations
v_traits <- c("root_N","root_CN")

for (i in 1:nrow(root_means)) {
  root_means[i,v_traits] <- root_traits %>% filter(species==root_means$species[i] & forest==root_means$forest[i]) %>% dplyr::select(all_of(v_traits)) %>% colMeans(na.rm=T)
}

root_means[root_means$species=='c.sco' & root_means$forest=='Quercus pyrenaica', c('root_N','root_CN')] <- root_means %>% subset(species=='c.sco') %>% dplyr::select(root_N,root_CN) %>% colMeans(na.rm=T)
root_means[root_means$species=='e.aus' & root_means$forest=='Quercus pyrenaica', c('root_N','root_CN')] <- root_means %>% subset(species=='e.arb') %>% dplyr::select(root_N,root_CN) %>% colMeans(na.rm=T)
root_means[root_means$species=='g.flo' & root_means$forest=='Quercus pyrenaica', c('root_N','root_CN')] <- root_means %>% subset(species=='g.flo') %>% dplyr::select(root_N,root_CN) %>% colMeans(na.rm=T)
root_means[root_means$species=='i.aqu' & root_means$forest=='Quercus pyrenaica', c('root_N','root_CN')] <- root_means %>% subset(species=='i.aqu') %>% dplyr::select(root_N,root_CN) %>% colMeans(na.rm=T)
root_means[root_means$species=='i.aqu' & root_means$forest=='Quercus petraea', c('root_N','root_CN')] <- root_means %>% subset(species=='i.aqu') %>% dplyr::select(root_N,root_CN) %>% colMeans(na.rm=T)
root_means[root_means$species=='j.com' & root_means$forest=='Quercus pyrenaica', c('root_N','root_CN')] <- root_means %>% subset(species=='j.com') %>% dplyr::select(root_N,root_CN) %>% colMeans(na.rm=T)
root_means[root_means$species=='p.avi' & root_means$forest=='Quercus pyrenaica', c('root_N','root_CN')] <- root_means %>% subset(species=='p.avi') %>% dplyr::select(root_N,root_CN) %>% colMeans(na.rm=T)
root_means[root_means$species=='s.auc' & root_means$forest=='Quercus petraea', c('root_N','root_CN')] <- root_means %>% subset(species=='s.auc') %>% dplyr::select(root_N,root_CN) %>% colMeans(na.rm=T)


# MERGE AND ADD HUBVAL

trait_means
root_means

# setdiff(paste(trait_means$species, trait_means$forest), paste(root_means$species, root_means$forest))
trait_means <- merge(trait_means, root_means, by=c('species','forest'))

trait_means$HubVal <- NA
v_traits <- colnames(HV_thick)[3]

for (i in 1:nrow(trait_means)) {
  sp <- trait_means$species[i]; fo <- trait_means$forest[i]
  if (paste(sp,fo) %in% paste(HV_thick$species,HV_thick$forest)) {
      trait_means[i,v_traits] <- HV_thick %>% filter(species%in%sp & forest%in%fo) %>% dplyr::select(all_of(v_traits)) %>% colMeans(na.rm=T)
  } else {
      trait_means[i,v_traits] <- HV_thick %>% filter(species==trait_means$species[i]) %>% dplyr::select(all_of(v_traits)) %>% colMeans(na.rm=T)
  }
}


# SAVE
write.table(trait_means, 'results/trait_means.txt')
