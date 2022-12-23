
library(tidyverse)
library(readxl)

setwd('C:/Users/user/OneDrive/TESIS Y PUBLICACIONES/COTTBUS/montejo/functional_structure')

# matrix of species abundances
abundances <- read_excel("data2.xlsx",sheet = "abundances") %>% as.data.frame() # import data
rownames(abundances) <- abundances$plot; abundances$plot <- NULL # species as rownames
abundances <- as.matrix(abundances) # sites as rows

# sites
sites <- read_excel("data2.xlsx",sheet = "sites") %>% as.data.frame()

# species vector
species <- colnames(abundances)

# create a presence/absence matrix
presences <- abundances
presences[,species][presences[,species]>0] <- 1

# traits
traits <- read_excel("data2.xlsx", sheet = "traits")
traits$forest <- as.factor(traits$forest)
str(traits)
root_traits <- read_excel("data2.xlsx", sheet = "root_traits")
root_traits$forest <- as.factor(root_traits$forest)
str(root_traits)

# explore number of NAs per variable
colSums(is.na(traits))
colSums(is.na(root_traits))

# outliers
traits$SRL[traits$SRL>80] <- NA
traits$SLA[traits$SLA>45] <- NA
traits$LDMC[traits$LDMC<0.15] <- NA
traits$LDMC[traits$LDMC>0.7] <- NA
traits$Thickness[traits$Thickness>0.08] <- NA
traits$SDMC[traits$SDMC>0.7] <- NA
traits$SDMC[traits$SDMC<0.15] <- NA
traits$HubVal[traits$HubVal>0.0015] <- NA
traits$RDMC[traits$RDMC<0.2] <- NA
traits$RDMC[traits$RDMC>0.9] <- NA
traits$SRA[traits$SRA>42] <- NA
traits$TMDr[traits$TMDr<0.2] <- NA
traits$Rdi[traits$Rdi>0.55] <- NA
traits$leaf_CN[traits$leaf_CN>50] <- NA
traits$leaf_d13C[traits$leaf_d13C<c(-35)] <- NA

root_traits$root_C[root_traits$root_C<20] <- NA
root_traits$root_C[root_traits$root_C>70] <- NA
root_traits$root_CN <- root_traits$root_C/root_traits$root_N

traits <- merge(traits, root_traits, by=c('species','forest'))
rm(root_traits)

# pairs
# pairs(traits[,3:20], lower.panel=NULL)

# par(mfrow=c(3,3))
# for (i in colnames(traits)[3:11]) { hist(traits[,i], 30, main=i) }
