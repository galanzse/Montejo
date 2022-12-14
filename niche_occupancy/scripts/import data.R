

library(tidyverse)
library(readxl)


# abundances
abundances <- read_excel("data.xlsx", sheet = "abundancias")
str(abundances)
colnames(abundances)[colnames(abundances)=='bosque'] <- 'forest'
abundances$forest <- as.factor(abundances$forest)
levels(abundances$forest) <- c('Fagus sylvatica','Shrubland','Mixed forest 2','Quercus pyrenaica','Mixed forest 1','Quercus petraea')
abundances$forest <- factor(abundances$forest, levels=c("Fagus sylvatica","Quercus petraea","Mixed forest 1","Mixed forest 2","Quercus pyrenaica","Shrubland"))

# species vector
species <- colnames(abundances)[3:dim(abundances)[2]]

# presence/absence
presences <- abundances
presences[,species][presences[,species]>0] <- 1


# species richness plot
spp_richness <- presences
spp_richness$richness <- rowSums(spp_richness[,species])
spp_richness <- spp_richness[,c('forest','plot','richness')]
ggplot(aes(y=richness, x=forest, fill=forest), data=spp_richness) +
  geom_boxplot() + theme_classic() + ylim(2,16) +
  labs(y='Species richness (S)', x=NULL) +
  theme(legend.position = "none")


# traits
traits <- read_excel("data.xlsx", sheet = "traits_limpio")
str(traits)
colnames(traits)[colnames(traits)=='especie'] <- 'species'
colnames(traits)[colnames(traits)=='bosque'] <- 'forest'

traits$forest <- as.factor(traits$forest)
levels(traits$forest) <- c('Fagus sylvatica','Shrubland','Mixed forest 2','Quercus pyrenaica','Mixed forest 1','Quercus petraea')
traits$HubVal <- as.numeric(traits$HubVal)
traits$Thickness <- as.numeric(traits$Thickness)

# import species codes
spp_code <- read_excel("data.xlsx", sheet = "traits_3")[,1:2] %>% unique()
colnames(spp_code) <- c('code','species')
traits <- merge(spp_code, traits, by='species', all.y=T)

# NA explore
colSums(is.na(traits)) # miro si hay NAs

# check there is trait data for all species in abundance table
table(species %in% traits$code)
# delete malus
traits <- traits %>% subset(species != 'Malus sp.')

# impute average value of species and site for HubVal and Thickness
for (tr in c('Thickness','HubVal')) {
  for (sp in unique(traits$code)) {
    for (bq in unique(traits$forest)) {
      mn <- traits[traits$code==sp & traits$forest==bq,]
      if (dim(mn)[1] != 0) {
        traits[traits$code==sp & traits$forest==bq,tr][is.na(traits[traits$code==sp & traits$forest==bq ,tr])] <- mean(mn[,tr], na.rm=T)
      }
    }
  }
}

# traits to analyse
trait.ana <- c('LDMC', 'SLA', 'd13C', 'C.N', 'HubVal', 'SDMC', 'SRA', 'Rdi')


# clean dataset
traits <- traits[,c('species','code','forest',trait.ana)]


rm(mn, spp_code, bq, sp, tr)
