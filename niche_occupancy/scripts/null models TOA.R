
# This section of codes is to 
# 1. calculate functional niche occupancy metrics (T, O and A) each local community in each functional space (trait combination).
# 2. construct 1000 null communities for each local community and calculate the metrics for each null community.

source('scripts/import data.R')
source('scripts/TOA.R')

library(vegan)
library(hypervolume)
library(factoextra)
library(cooccur)


# defining 'restricted' regional pool ####
mat <- presences %>% as.data.frame()
rownames(mat) <- presences$plot
mat <- mat[,species] %>% as.matrix() %>% t()
cooc_mat <- cooccur(mat, spp_names=T)
plot(cooc_mat)
# matrix of pairs of observed cooccurrences
cooc_mat <- cooc_mat$results %>% dplyr::select(sp1_name, sp2_name, obs_cooccur)
# pairs of species which cooccur in atleast one plot
cooc_mat <- cooc_mat %>% filter(obs_cooccur > 0)
# retain pairs of observations which include the species of interest
cooc_mat <- cooc_mat %>% subset(sp1_name %in% c('f.syl','q.pet','q.pyr') | sp2_name %in% c('f.syl','q.pet','q.pyr'))

# pool
cooc_list <- list()

temp <- cooc_mat[,1:2] %>% dplyr::filter(sp1_name=='f.syl' | sp2_name=='f.syl') # fagus
cooc_list[['fagus']] <- c(temp$sp1_name, temp$sp2_name) %>% unique()

temp <- cooc_mat[,1:2] %>% dplyr::filter(sp1_name=='q.pyr' | sp2_name=='q.pyr') # pyrenaica
cooc_list[['pyrenaica']] <- c(temp$sp1_name, temp$sp2_name) %>% unique()

temp <- cooc_mat[,1:2] %>% dplyr::filter(sp1_name=='q.pet' | sp2_name=='q.pet') # oak
cooc_list[['oak']] <- c(temp$sp1_name, temp$sp2_name) %>% unique()

cooc_list[['shrubland']] <- species[species!='f.syl'] # shrubland

temp <- cooc_mat[,1:2] %>% dplyr::filter(sp1_name%in%c('f.syl','q.pet','q.pyr') | sp2_name%in%c('f.syl','q.pet','q.pyr')) # mixed
cooc_list[['mixed']] <- c(temp$sp1_name, temp$sp2_name) %>% unique()

temp <- cooc_mat[,1:2] %>% dplyr::filter(sp1_name%in%c('f.syl','q.pyr') | sp2_name%in%c('f.syl','q.pyr')) # pyrenaica.fagus
cooc_list[['pyrenaica.fagus']] <- c(temp$sp1_name, temp$sp2_name) %>% unique()


# PCA to compute hypervolumes ####
trait_pca <- prcomp(traits[,trait.ana], scale.=TRUE, center=TRUE)
summary(trait_pca)
eigenvals(trait_pca)
plot(trait_pca$x[,1:2])
fviz_pca_var(trait_pca, col.var="contrib", repel=TRUE)
trait_pca$rotation

# retain PCs with eigenvalues > 1
df_trait_pca <- cbind(traits[,c("forest","code")], trait_pca$x[,which(eigenvals(trait_pca) > 1)])
df_trait_pca$forxspp <- paste(df_trait_pca$forest, df_trait_pca$code, sep=' ')

# functional axes to analyse
trait.ana <- c("PC1","PC2","PC3")


# Li et al. 2017 [quantile threshold = 0.05, 1,000 MC samples x data point, fixed kernel band-width = 0.5 SD]
bw <- 0.5
sxp <- 1000
qt <- 0.05 # default (see Blonder et al 2018)

# global species pool with all datasets covering the selected traits to be analysed
df_trait_pca
# list of all local communities (plots)
plots <- presences$plot
# indicates the number of plots need to be analysed
print(length(plots))
# a matrix store the results
out_forest <- matrix(0,length(plots),11,dimnames=list(NULL,c("plot","S","T.ob","O.ob","A.ob","T.exp","O.exp","A.exp","T.sd","O.sd","A.sd")))
out_forest[,'plot'] <- plots
#
out_cocc <- out_forest
out_all <- out_forest

# N for null models
n.run <- 100
# save runs
forest.null <- list()
coocc.null <- list()
all.null <- list()

# for each local community (plot) construct "n.run" null communities
for(i in 1:nrow(out_forest)) {
  # species matrix
  myab <- abundances %>% subset(plot==out_forest[i,'plot'])
  myesp <- colnames(myab[,species])[colSums(myab[,species])>0]
  # species' names and trait values for each local community
  ob.com <- df_trait_pca %>% subset(forest==myab$forest & code %in% myesp) %>% dplyr::select(code, trait.ana)
  # species richness
  out_forest[i,"S"] <- length(unique(ob.com$code))
  # calculate T,O and A for observed communities
  out_forest[i,c('T.ob','O.ob','A.ob')] <- TOA(ob.com, bw=bw, sxp=sxp)
  
  
  # FIRST NULL MODEL: FOREST SPECIES AS REGIONAL POOL
  sp.region <- df_trait_pca %>% subset(forest==myab$forest) %>% dplyr::select(code) %>% unique() %>% deframe()
  null.b <- matrix(nrow=n.run, ncol=3)
  colnames(null.b) <- c("T.ob","O.ob","A.ob")
  # run null model
  for (j in 1:n.run) {
    # maintain species richness
    sp.id <- sp.region[sample(1:length(sp.region), out_forest[i,"S"], replace=F)]
    # select traits from the forest
    sam.com <- df_trait_pca %>% filter(forest==myab$forest) %>% subset(code %in% sp.id) %>% dplyr::select(code, trait.ana)
    # hypervolume
    null.b[j,c("T.ob","O.ob","A.ob")] <- TOA(sam.com, bw=bw, sxp=sxp)
  }
  # save raw results
  forest.null[[i]] <- null.b
  # calculate mean of T, O and A for the 1000 null communities
  out_forest[i,c('T.exp','O.exp','A.exp')] <- apply(null.b, 2, mean)
  # calculate standard deviation of T, O and A for the 1000 null communities
  out_forest[i,c('T.sd','O.sd','A.sd')] <- apply(null.b, 2, sd)
  
  
  # SECOND NULL MODEL: OBSERVED COOCCURRENCES AS REGIONAL POOL
  # sp.region <- cooc_list[[myab$forest]]
  # null.b <- matrix(nrow=n.run, ncol=3)
  # colnames(null.b) <- c("T.ob","O.ob","A.ob")
  # # run null model
  # for (j in 1:n.run) {
  #   # maintain species richness
  #   sp.id <- sp.region[sample(1:length(sp.region), out_forest[i,"S"], replace=F)]
  #   # randomly select traits from the dataset
  #   sam.com <- df_trait_pca %>% subset(code %in% sp.id) %>% dplyr::select(forest,code) %>% unique() %>% group_by(code) %>% slice_sample(n=1)
  #   sam.com$forxspp <- paste(sam.com$forest, sam.com$code, sep=' ')
  #   sam.com <- df_trait_pca %>% dplyr::filter(forxspp %in% sam.com$forxspp) %>% dplyr::select(code, trait.ana)
  #   # hypervolume
  #   null.b[j,c("T.ob","O.ob","A.ob")] <- TOA(sam.com, bw=bw, sxp=sxp)
  # }
  # # guardo los resultados brutos del modelo nulo
  # coocc.null[[i]] <- null.b
  # # calculate mean of T, O and A for the 1000 null communities
  # out_cocc[i,c('T.exp','O.exp','A.exp')] <- apply(null.b, 2, mean)
  # # calculate standard deviation of T, O and A for the 1000 null communities
  # out_cocc[i,c('T.sd','O.sd','A.sd')] <- apply(null.b, 2, sd)
  
  
  # THIRD NULL MODEL: OBSERVED COOCCURRENCES AS REGIONAL POOL
  sp.region <- species
  null.b <- matrix(nrow=n.run, ncol=3)
  colnames(null.b) <- c("T.ob","O.ob","A.ob")
  # run null model
  for (j in 1:n.run) {
    # maintain species richness
    sp.id <- sp.region[sample(1:length(sp.region), out_forest[i,"S"], replace=F)]
    # randomly select traits from the dataset
    sam.com <- df_trait_pca %>% subset(code %in% sp.id) %>% dplyr::select(forest,code) %>% unique() %>% group_by(code) %>% slice_sample(n=1)
    sam.com$forxspp <- paste(sam.com$forest, sam.com$code, sep=' ')
    sam.com <- df_trait_pca %>% dplyr::filter(forxspp %in% sam.com$forxspp) %>% dplyr::select(code, trait.ana)
    # hypervolume
    null.b[j,c("T.ob","O.ob","A.ob")] <- TOA(sam.com, bw=bw, sxp=sxp)
  }
  # save raw results
  all.null[[i]] <- null.b
  # calculate mean of T, O and A for the 1000 null communities
  out_all[i,c('T.exp','O.exp','A.exp')] <- apply(null.b, 2, mean)
  # calculate standard deviation of T, O and A for the 1000 null communities
  out_all[i,c('T.sd','O.sd','A.sd')] <- apply(null.b, 2, sd)

}

# complete observed values
out_cocc[,c('S','T.ob','O.ob','A.ob')] <- out_forest[,c('S','T.ob','O.ob','A.ob')]
out_all[,c('S','T.ob','O.ob','A.ob')] <- out_forest[,c('S','T.ob','O.ob','A.ob')]

# add forest + null model info
out_forest <- merge(presences[,c('forest','plot')], out_forest, by='plot'); out_forest$null <- 'Forest'
out_cocc <- merge(presences[,c('forest','plot')], out_cocc, by='plot'); out_cocc$null <- 'Cooccurrences'
out_all <- merge(presences[,c('forest','plot')], out_all, by='plot'); out_all$null <- 'All'


# Effect Sizes ####
out_forest$T.SES <- (out_forest$T.ob - out_forest$T.exp) / out_forest$T.sd
out_forest$O.SES <- (out_forest$O.ob - out_forest$O.exp) / out_forest$O.sd
out_forest$A.SES <- (out_forest$A.ob - out_forest$A.exp) / out_forest$A.sd

out_cocc$T.SES <- (out_cocc$T.ob - out_cocc$T.exp) / out_cocc$T.sd
out_cocc$O.SES <- (out_cocc$O.ob - out_cocc$O.exp) / out_cocc$O.sd
out_cocc$A.SES <- (out_cocc$A.ob - out_cocc$A.exp) / out_cocc$A.sd

out_all$T.SES <- (out_all$T.ob - out_all$T.exp) / out_all$T.sd
out_all$O.SES <- (out_all$O.ob - out_all$O.exp) / out_all$O.sd
out_all$A.SES <- (out_all$A.ob - out_all$A.exp) / out_all$A.sd


# save results ####
write.table(out_forest,"resultados/out_forest.txt")
save(forest.null, file="resultados/forest.null.Rdata")
write.table(out_cocc,"resultados/out_cocc.txt")
save(coocc.null, file="resultados/coocc.null.Rdata")
write.table(out_all,"resultados/out_all.txt")
save(all.null, file="resultados/all.null.Rdata")
