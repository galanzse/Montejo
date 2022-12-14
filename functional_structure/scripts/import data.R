
library(tidyverse)
library(readxl)

setwd('C:/Users/user/OneDrive/TESIS Y PUBLICACIONES/COTTBUS/montejo/functional_structure')

# create a matrix of species abundances
abundances <- read_excel("data2.xlsx",sheet = "abundances") %>% as.data.frame() # import data
rownames(abundances) <- abundances$species; abundances$species <- NULL # species as rownames
abundances <- as.matrix(abundances) %>% t() # sites as rows

# species vector
species <- colnames(abundances)

# create a presence/absence matrix
presences <- abundances
presences[,species][presences[,species]>0] <- 1


# traits
v_traits <- c('LDMC','SLA','Thickness','SDMC','HuberValue','RDMC','SRL','SRA','TMDr','AvgDiam','d13C','d15N','N','C','C.N')
traits <- read_excel("data2.xlsx", sheet = "traits")
traits <- traits %>% dplyr::select(species, forest, all_of(v_traits))
str(traits)

traits$forest <- as.factor(traits$forest)

# explore number of NAs per variable
colSums(is.na(traits))

# outliers
traits$SRL[traits$SRL>80] <- NA
traits$SLA[traits$SLA>45] <- NA
traits$Thickness[traits$Thickness>0.08] <- NA
traits$SDMC[traits$SDMC>0.7] <- NA
traits$HuberValue[traits$HuberValue>0.0015] <- NA
traits$RDMC[traits$RDMC<0.2] <- NA
traits$SRA[traits$SRA>42] <- NA
traits$AvgDiam[traits$AvgDiam>0.55] <- NA
traits$C.N[traits$C.N>55] <- NA

# pairs
pairs(traits[,v_traits], lower.panel=NULL)
