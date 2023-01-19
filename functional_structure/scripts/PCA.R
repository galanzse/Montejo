
library(tidyverse)
library(vegan)
library(factoextra)

# import data
moments_weighted <- read.csv("results/moments_weighted.txt", sep="")
# trait_dist_unweighted <- read.csv("results/trait_dist_unweighted.txt", sep="")

# PCA of weighted CWM
CWM_PCA <- moments_weighted[,1:4] %>% spread(key='trait', value='CWM') %>% na.omit()

# traits to use
v_traits <- colnames(CWM_PCA)[3:12]
# v_traits <- v_traits[!(v_traits%in%c(leaf_d15N','leaf_CN','root_CN','SRA','Rdi'))]

# there is something weird with Shrubland
# pairs(CWM_PCA[,3:15], lower.panel=NULL)

mipca <- prcomp(CWM_PCA[,v_traits], scale.=FALSE, center=FALSE) # already scaled traits
summary(mipca)
fviz_eig(mipca); eigenvals(mipca)
fviz_pca_var(mipca, col.var="contrib", repel=TRUE)
mipca$rotation

df_pca2 <- mipca$x %>% as.data.frame()
df_pca2$forest <- CWM_PCA$forest

ggplot(aes(x=PC1, y=PC2, colour=forest), data=df_pca2) +
  geom_point() +
  theme_classic()
