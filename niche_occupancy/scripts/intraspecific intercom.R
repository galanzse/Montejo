
# following Benavides et al. 2019

source('scripts/abiotic.R')
library(hypervolume)
library(proxy)


# functional space
trait_pca <- prcomp(traits[,trait.ana], scale.=TRUE, center=TRUE)
df_trait_pca <- cbind(traits[,c("forest","code")], trait_pca$x[,which(eigenvals(trait_pca) > 1)])

# species with traits measured in 4 or more sites
table(rowSums(table(traits$code, traits$forest)) > 15)
spp_x4 <- which(rowSums(table(traits$code, traits$forest)) > 15) %>% names() # vector

# functional axes to consider
trait.ana <- c("PC1","PC2","PC3")

# observed differences
sppxfor <- df_trait_pca %>% filter(code %in% spp_x4) %>% dplyr::select(forest, code) %>% unique()
intr_inter_r <- sppxfor %>% dplyr::select(code) %>% unique()
intr_inter_r$diff.dist <- NA
intr_inter_r$diff.dist1 <- NA
intr_inter_r$diff.dist2 <- NA
intr_inter_r$diff.dist3 <- NA
intr_inter_r$diff.vol <- NA
intr_inter_r$diff.vol1 <- NA
intr_inter_r$diff.vol2 <- NA
intr_inter_r$diff.vol3 <- NA

# null model for each species
n.run <- 499
dist.null.m <- matrix(ncol=length(unique(sppxfor$code)), nrow=n.run) # niche shift
colnames(dist.null.m) <- unique(obs.centroids$code)
dist1.null.m <- dist.null.m
dist2.null.m <- dist.null.m
dist3.null.m <- dist.null.m

vol.null.m <- dist.null.m # niche breadth
vol1.null.m <- dist.null.m
vol2.null.m <- dist.null.m
vol3.null.m <- dist.null.m

for (r in 1:nrow(intr_inter_r)) {
  df_tm <- sppxfor %>% filter(code %in% intr_inter_r$code[r]) %>% unique()
  df_tm$vol <- NA
  df_tm$PC1 <- NA; df_tm$PC2 <- NA; df_tm$PC3 <- NA # shift
  df_tm$V1 <- NA; df_tm$V2 <- NA; df_tm$V3 <- NA # breadth
  
  # average observed distance
  for (f in 1:nrow(df_tm)) {
    df_hv <- df_trait_pca %>% filter(code==df_tm$code[f] & forest==df_tm$forest[f]) %>% dplyr::select(all_of(trait.ana))
    hn <- hypervolume_box(data=df_hv, samples.per.point=1000, kde.bandwidth=0.5, verbose=F)
    df_tm$vol[f] <- hn@Volume
    df_tm[f,c('PC1','PC2','PC3')] <- get_centroid(hn)
    df_tm[f,c('V1','V2','V3')] <- hypervolume_variable_importance(hn)
  }
  intr_inter_r$diff.dist[r] <- mean(dist(df_tm[,c('PC1','PC2','PC3')], method='euclidean'))
  intr_inter_r$diff.dist1[r] <- mean(dist(df_tm$PC1, method='euclidean'))
  intr_inter_r$diff.dist2[r] <- mean(dist(df_tm$PC2, method='euclidean'))
  intr_inter_r$diff.dist3[r] <- mean(dist(df_tm$PC3, method='euclidean'))
  
  intr_inter_r$diff.vol[r] <- mean(dist(df_tm$vol, method='euclidean'))
  intr_inter_r$diff.vol1[r] <- mean(dist(df_tm$V1, method='euclidean'))
  intr_inter_r$diff.vol2[r] <- mean(dist(df_tm$V2, method='euclidean'))
  intr_inter_r$diff.vol3[r] <- mean(dist(df_tm$V3, method='euclidean'))
  
  # null model
  for (n in 1:n.run) {
    # because the observed values are the average across 4-6 sites, I build the null model likewise
    df_hv <- df_trait_pca %>% filter(code==unique(df_tm$code)) %>% dplyr::select(all_of(trait.ana)) %>% slice_sample(n=5)
    h1 <- hypervolume_box(data=df_hv, samples.per.point=1000, kde.bandwidth=0.5, verbose=F)
    df_hv <- df_trait_pca %>% filter(code==unique(df_tm$code)) %>% dplyr::select(all_of(trait.ana)) %>% slice_sample(n=5)
    h2 <- hypervolume_box(data=df_hv, samples.per.point=1000, kde.bandwidth=0.5, verbose=F)
    df_hv <- df_trait_pca %>% filter(code==unique(df_tm$code)) %>% dplyr::select(all_of(trait.ana)) %>% slice_sample(n=5)
    h3 <- hypervolume_box(data=df_hv, samples.per.point=1000, kde.bandwidth=0.5, verbose=F)
    df_hv <- df_trait_pca %>% filter(code==unique(df_tm$code)) %>% dplyr::select(all_of(trait.ana)) %>% slice_sample(n=5)
    h4 <- hypervolume_box(data=df_hv, samples.per.point=1000, kde.bandwidth=0.5, verbose=F)
    df_hv <- df_trait_pca %>% filter(code==unique(df_tm$code)) %>% dplyr::select(all_of(trait.ana)) %>% slice_sample(n=5)
    h5 <- hypervolume_box(data=df_hv, samples.per.point=1000, kde.bandwidth=0.5, verbose=F)

    # shift
    t_dist <- rbind(get_centroid(h1),get_centroid(h2),get_centroid(h3), get_centroid(h4),get_centroid(h5))
    dist.null.m[n,unique(df_tm$code)] <- t_dist %>% dist(method="euclidean") %>% mean()
    dist1.null.m[n,unique(df_tm$code)] <- t_dist[,'PC1'] %>% dist(method="euclidean") %>% mean()
    dist2.null.m[n,unique(df_tm$code)] <- t_dist[,'PC2'] %>% dist(method="euclidean") %>% mean()
    dist3.null.m[n,unique(df_tm$code)] <- t_dist[,'PC3'] %>% dist(method="euclidean") %>% mean()

    # volume
    vol.null.m[n,unique(df_tm$code)] <- dist(c(h1@Volume, h2@Volume, h3@Volume, h4@Volume, h5@Volume), method="euclidean") %>% mean()
    v.imp <- rbind(hypervolume_variable_importance(h1), hypervolume_variable_importance(h2),
                  hypervolume_variable_importance(h3), hypervolume_variable_importance(h4), hypervolume_variable_importance(h5))
    vol1.null.m[n,unique(df_tm$code)] <- dist(v.imp[,'PC1'], method="euclidean") %>% mean()
    vol2.null.m[n,unique(df_tm$code)] <- dist(v.imp[,'PC2'], method="euclidean") %>% mean()
    vol3.null.m[n,unique(df_tm$code)] <- dist(v.imp[,'PC3'], method="euclidean") %>% mean()

  }
}

# SES
intr_inter_r$SES.dist <- (intr_inter_r$diff.dist - apply(dist.null.m, 2, mean)) / apply(dist.null.m, 2, sd)
intr_inter_r$SES.dist1 <- (intr_inter_r$diff.dist1 - apply(dist1.null.m, 2, mean)) / apply(dist1.null.m, 2, sd)
intr_inter_r$SES.dist2 <- (intr_inter_r$diff.dist2 - apply(dist2.null.m, 2, mean)) / apply(dist2.null.m, 2, sd)
intr_inter_r$SES.dist3 <- (intr_inter_r$diff.dist3 - apply(dist3.null.m, 2, mean)) / apply(dist3.null.m, 2, sd)
intr_inter_r$SES.vol <- (intr_inter_r$diff.vol - apply(vol.null.m, 2, mean)) / apply(vol.null.m, 2, sd)
intr_inter_r$SES.v1 <- (intr_inter_r$diff.vol1 - apply(vol1.null.m, 2, mean)) / apply(vol1.null.m, 2, sd)
intr_inter_r$SES.v2 <- (intr_inter_r$diff.vol2 - apply(vol2.null.m, 2, mean)) / apply(vol2.null.m, 2, sd)
intr_inter_r$SES.v3 <- (intr_inter_r$diff.vol3 - apply(vol3.null.m, 2, mean)) / apply(vol3.null.m, 2, sd)

# save
write.table(intr_inter_r, 'resultados/intr_inter_r.txt')
save(dist.null.m, dist1.null.m, dist2.null.m, dist3.null.m,
     vol.null.m, vol1.null.m, vol2.null.m, vol3.null.m, file="resultados/intras_inter_null.Rdata")


# table
intr_inter_l_sh <- intr_inter_r %>% pivot_longer(10:13, names_to='SES'); intr_inter_l_sh$var <- 'Shift'
intr_inter_l_br <- intr_inter_r %>% pivot_longer(14:17, names_to='SES'); intr_inter_l_br$var <- 'Breadth'

intr_inter_l_sh$SES <- as.factor(intr_inter_l_sh$SES)
levels(intr_inter_l_sh$SES) <- c('Total', 'Leaf economics', 'Root economics','Huber')

intr_inter_l_br$SES <- as.factor(intr_inter_l_br$SES)
levels(intr_inter_l_br$SES) <- c('Leaf economics', 'Root economics','Huber','Total')
intr_inter_l_br$SES  <- factor(intr_inter_l_br$SES, levels=c('Total','Leaf economics', 'Root economics','Huber'))

intr_inter_l <- rbind(intr_inter_l_sh[,c('SES','value','var')], intr_inter_l_br[,c('SES','value','var')])

# plot
ggplot(aes(y=value, x=SES), data=intr_inter_l) +
  geom_boxplot() + theme_classic() +
  ylim(-1, 3.5) + geom_hline(yintercept=0) +
  labs(y='Standardized Effect Size (SES)', x='') +
  geom_hline(yintercept=0) +
  facet_wrap(~var)


# tests
wilcox.test(intr_inter_r$SES.dist)
wilcox.test(intr_inter_r$SES.dist1)
wilcox.test(intr_inter_r$SES.dist2)
wilcox.test(intr_inter_r$SES.dist3)

wilcox.test(intr_inter_r$SES.vol)
wilcox.test(intr_inter_r$SES.v1)
wilcox.test(intr_inter_r$SES.v2)
wilcox.test(intr_inter_r$SES.v3)
