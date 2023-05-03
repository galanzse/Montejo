
# Optimization of bandwidth using the disjunt factor (Barros et al. 2016)

# The disjunct factor is the ratio between the size of the calculated hypervolume and the size of a hypervolume constructed from the same data with disjunct data points (i.e. no overlapping kernels). # Values > 0.9 indicate that the hypervolume has ‘holes’ and should be avoided by increasing the bandwidth value.

library(hypervolume)
source('scripts/import data.R')

# functional space
trait_pca <- prcomp(traits[,trait.ana], scale.=TRUE, center=TRUE)
df_trait_pca <- cbind(traits[,c("forest","code")], trait_pca$x[,which(eigenvals(trait_pca) > 1)])
df_trait_pca$forxspp <- paste(df_trait_pca$forest, df_trait_pca$code, sep=' ')


# estimate average bw using silverman-1d
bw_results <- vector()
for (i in 1:1000) {
  df_tr <- df_trait_pca %>% filter(forxspp == sample(df_trait_pca$forxspp, 1)) %>% dplyr::select(PC1,PC2,PC3)
  bw_results[i] <- estimate_bandwidth(df_tr, method="silverman-1d") %>% mean()
}
boxplot(bw_results, main=paste('mean bandwidth = ', round(mean(bw_results),2)))


# results
bw_results <- matrix(ncol=length(seq(0.1,1,0.05)), nrow=30)
colnames(bw_results) <- seq(0.1,1,0.05)

# starting bandwidth values of 0.1, increased in steps of 0.05, retain disjunct factor
for (c in 1:ncol(bw_results)) {
  bw <- as.numeric(colnames(bw_results)[c])
  
  for (r in 1:nrow(bw_results)) {
    df_tr <- df_trait_pca %>% filter(forxspp == sample(df_trait_pca$forxspp, 1)) %>% dplyr::select(PC1,PC2,PC3)
    hv1 <- hypervolume_box(data=df_tr, samples.per.point=1000, kde.bandwidth=bw, verbose=F)
    # formula extracted from hypervolume_box
    bw_results[r,c] <- hv1@Volume/(nrow(as.matrix(df_tr)) * prod(2 * rep(bw,3)))
  }
}

write.table(bw_results,"resultados/bw_results.txt")


# plot
bw_results <- bw_results %>% as.data.frame() %>% pivot_longer(1:19, names_to='bandwidth', values_to='disjunct_factor')
str(bw_results)
bw_results$bandwidth <- as.numeric(bw_results$bandwidth)

ggplot(aes(y=disjunct_factor, x=bandwidth), data=bw_results) +
  geom_point() + theme_classic() +
  geom_smooth(method='loess') +
  labs(x='Bandwidth', y='Disjunct factor') +
  geom_hline(yintercept=0.9, color='red', linetype="dashed") +
  geom_vline(xintercept=0.48, color='blue') + geom_vline(xintercept=0.55, color='green')
  
