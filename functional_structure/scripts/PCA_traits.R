
library(factoextra)
library(vegan)

source('scripts/import data.R')

# select traits
v_traits <- c('HubVal','LDMC','leaf_CN','leaf_d13C','Rdi','root_CN','SDMC','SLA','SRA')
traits <- traits %>% select(species, forest, v_traits)
str(traits)

# impute traits using species means at regional level
colSums(is.na(traits[,3:11]))
for (tr in v_traits) { for (r in 1:nrow(traits)) {
    if (is.na(traits[r,tr])) { traits[r,tr] <- traits %>% subset(species==traits$species[r]) %>% select(tr) %>% colMeans(na.rm=T) }
} }

# use congeneric species to impute the remaining NAs
traits$HubVal[traits$species=='c.pur'] <- traits$HubVal[traits$species=='c.sco'][1:length(traits$HubVal[traits$species=='c.pur'])]

# PCA
pca1 <- prcomp(traits[,v_traits], scale=T, center=T)
summary(pca1)
eigenvals(pca1)

# change PC2 so more positive means more acquisitive
pca1$rotation[,'PC2'] <- pca1$rotation[,'PC2'] * c(-1)
pca1$x[,'PC2'] <- pca1$x[,'PC2'] * c(-1)

fviz_pca_var(pca1, repel=T)

# merge data
pca1 <- cbind(traits[,1:2], pca1$x[,1:3])

# name functional axes
colnames(pca1)[3:5] <- c('Leaf_Economics','Root_Economics','HyArq')
