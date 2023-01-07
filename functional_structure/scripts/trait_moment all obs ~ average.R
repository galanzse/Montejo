
source('C:/Users/user/OneDrive/TESIS Y PUBLICACIONES/COTTBUS/montejo/functional_structure/scripts/import data.R')
source('C:/Users/user/OneDrive/TESIS Y PUBLICACIONES/COTTBUS/montejo/functional_structure/scripts/trait_moment function.R')

# we are going to compare trait moment indexes of our communities calculated using:
# 1/ all data and splitting species abundances by the number of observations
# 2/ trait averages


# trait averagerages
v_traits <- colnames(traits)[3:11]
traits[,v_traits] <- apply(traits[,v_traits], 2, scale)
averagerage_dist_weighted <- list() # abundances
traits <- traits %>% group_by(species, forest) %>% summarise(LDMC=mean(LDMC, na.rm=T),
                                                   SLA=mean(SLA, na.rm=T),
                                                   SDMC=mean(SDMC, na.rm=T),
                                                   RDMC=mean(RDMC, na.rm=T),
                                                   SRL=mean(SRL, na.rm=T),
                                                   TMDr=mean(TMDr, na.rm=T),
                                                   leaf_d13C=mean(leaf_d13C, na.rm=T),
                                                   leaf_C=mean(leaf_C, na.rm=T),
                                                   leaf_N=mean(leaf_N, na.rm=T))

sites <- read_excel("data2.xlsx",sheet = "sites") %>% as.data.frame()

for (t in 1:length(v_traits)) {
  sites$trait <- v_traits[t]
  sites$CWM <- NA
  sites$CWV <- NA
  sites$CWS <- NA
  sites$CWK <- NA
  sites_ab <- sites

  for (s in 1:nrow(sites)) {
    
    pl1 <- as.character(sites[s,'plot']) # plot
    ab1 <- abundances[pl1,]
    sp1 <- which(ab1>0) %>% names() # species names
    ab1 <- ab1[sp1] # abundaces

    # filter trait database
    df1 <- traits %>% filter(forest==sites[s,'forest'] & species%in%sp1) %>%
      dplyr::select(species,v_traits[t]) %>% na.omit()
    colnames(df1) <- c("species","trait")
    
    # remove species without trait data
    ab1 <- ab1[unique(df1$species)]
    
    # weight observations by missing data
    ab2 <- ab1 / sum(ab1)
    ab2 <- as.data.frame(ab2); ab2$species <- rownames(ab2)
    colnames(ab2) <- c('abundance','species')
    ab2 <- merge(df1, ab2, by='species')

    # run function
    sites_ab[s,c('CWM','CWV','CWS','CWK')] <- trait_moment(abundance=ab2$abundance, trait=ab2$trait)

  }
  
  # saverage results
  averagerage_dist_weighted[[v_traits[t]]] <- sites_ab

  # progress
  print(v_traits[t])
  
}

averagerage_dist_weighted <- do.call('rbind', averagerage_dist_weighted)
rownames(averagerage_dist_weighted) <- 1:nrow(averagerage_dist_weighted)
write.table(averagerage_dist_weighted, 'results/averagerage_dist_weighted.txt')


# comparison
alldata <- read.csv("results/moments_weighted.txt", sep="")
# rename indexes & merge
colnames(alldata)[4:7] <- c("CWM.all","CWV.all","CWS.all","CWK.all" )
colnames(averagerage_dist_weighted)[4:7] <- c("CWM.average","CWV.average","CWS.average","CWK.average" )
alldata <- merge(alldata, averagerage_dist_weighted, by=c('forest','plot','trait'))

alldata <- alldata %>% filter(CWK.average<200)

par(mfrow=c(2,2))
plot(CWM.all ~ CWM.average, data=alldata, main='CWM', xlab='Average', ylab='All observations')
abline(0,1,col='blue',lwd=2)
plot(CWV.all ~ CWV.average, data=alldata, main='CWV', xlab='Average', ylab='All observations')
abline(0,1,col='blue',lwd=2)
plot(CWS.all ~ CWS.average, data=alldata, main='CWS', xlab='Average', ylab='All observations')
abline(0,1,col='blue',lwd=2)
plot(CWK.all ~ CWK.average, data=alldata, main='CWK', xlab='Average', ylab='All observations')
abline(0,1,col='blue',lwd=2)

# colour by forest
ggplot(aes(x=CWM.all, y=CWM.average, color=forest), data=alldata) +
  geom_point() + theme_classic()
  # facet_wrap(~trait, scales='free_y') +
  # theme(legend.position='top')
ggplot(aes(x=CWV.all, y=CWV.average, color=forest), data=alldata) +
  geom_point() + theme_classic()
  # facet_wrap(~trait, scales='free_y') +
  # theme(legend.position='top')
ggplot(aes(x=CWS.all, y=CWS.average, color=forest), data=alldata) +
  geom_point() + theme_classic()
  # facet_wrap(~trait, scales='free_y') +
  # theme(legend.position='top')
ggplot(aes(x=CWK.all, y=CWK.average, color=forest), data=alldata) +
  geom_point() + theme_classic()
  # facet_wrap(~trait, scales='free_y') +
  # theme(legend.position='top')
