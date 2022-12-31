
library(tidyverse)

# import data
trait_dist_weighted <- read.csv("results/trait_dist_weighted.txt", sep="")
# trait_dist_unweighted <- read.csv("results/trait_dist_unweighted.txt", sep="")
# root_dist_weighted <- read.csv("results/root_dist_weighted.txt", sep="")
# root_dist_unweighted <- read.csv("results/root_dist_unweighted.txt", sep="")

# PCA of weighted CWM
df_pca <- trait_dist_weighted[,1:4] %>% spread(key='trait', value='CWM')

# there is something weird with Shrubland
pairs(df_pca[,v_traits], lower.panel=NULL)

mipca <- prcomp(df_pca[,v_traits], scale.=FALSE, center=FALSE) # already scaled traits
summary(mipca)
eigenvals(mipca)
fviz_pca_var(mipca, col.var="contrib", repel=TRUE)
mipca$rotation

df_pca2 <- mipca$x %>% as.data.frame()
df_pca2$forest <- df_pca$forest

ggplot(aes(x=PC1, y=PC2, colour=forest), data=df_pca2) +
  geom_point() +
  theme_classic()
