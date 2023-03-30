
setwd('C:/Users/user/OneDrive/PUBLICACIONES/BTU/montejo/functional_structure')
source('scripts/import data.R')

# variance partitining considering species abundances
# de la Riva et a. 2016 (Oikos); modified from Lepš et al. (2011)

# template
temp <- read.csv("results/moments_weighted.txt", sep="")[,1:4] %>% spread(key=trait, value=CWM)
temp[,3:14] <- NA
head(temp)
  
v_traits <- colnames(temp)[3:14]

# calculate five types of CWM (parameters) using:

  # 1) specific (forest averages): itv + turnover
site_means <- read.csv("results/trait_means.txt", sep="")
specific <- temp
  # 2) fixed (region averages): turnover
montejo_means <- site_means %>% dplyr::select(-forest) %>% group_by(species) %>%
  summarise(HubVal=mean(HubVal, na.rm=T),
            LDMC=mean(LDMC, na.rm=T),
            leaf_CN=mean(leaf_CN, na.rm=T),
            leaf_d13C=mean(leaf_d13C, na.rm=T),
            leaf_N=mean(leaf_N, na.rm=T),
            Rdi=mean(Rdi, na.rm=T),
            RDMC=mean(RDMC, na.rm=T),
            root_CN=mean(root_CN, na.rm=T),
            root_N=mean(root_N, na.rm=T),
            SDMC=mean(SDMC, na.rm=T),
            SLA=mean(SLA, na.rm=T),
            SRA=mean(SRA, na.rm=T))

fixed <- temp
  # 3) unweighted (presence data): abundance
unweighted <- temp

for (i in 1:nrow(specific)) {
  
  # species
  p1 <- abundances[as.character(fixed$plot[i]),] %>% as.data.frame()
  colnames(p1) <- 'abundance'; p1$species <- rownames(p1); rownames(p1) <- NULL
  p1 <- p1 %>% filter(abundance>0)

    for (t in colnames(specific)[3:14]) {
    
      # specific
      t1 <- site_means %>% filter(species %in% p1$species & forest == specific$forest[i]) %>% dplyr::select(species, t)
      colnames(t1) <- c('species','trait')
      t1 <- merge(t1, p1, by = intersect(names(t1), names(p1)))
      specific[i,t] <- sum(t1$trait * t1$abundance / sum(t1$abundance))
      
      # fixed
      t1 <- montejo_means %>% filter(species %in% p1$species) %>% dplyr::select(species, t)
      colnames(t1) <- c('species','trait')
      t1 <- merge(t1, p1, by = intersect(names(t1), names(p1)))
      fixed[i,t] <- sum(t1$trait * t1$abundance / sum(t1$abundance))
      
      # unweighted
      unweighted[i,t] <- mean(t1$trait)
      
  }
}

  # 4) intraspecific parameter (specific - fixed): itv
intraspecific <- temp
for (t in colnames(intraspecific)[3:14]) {  intraspecific[,t] <- specific[,t] - fixed[,t] }

  # 5) species_abundance parameter (fixed – unweighted): abundance
sp_abundance <- temp
for (t in colnames(sp_abundance)[3:14]) {  sp_abundance[,t] <- unweighted[,t] - fixed[,t] }

# save
var_partitioning <- list()
var_partitioning[['specific']] <- specific
var_partitioning[['fixed']] <- fixed
var_partitioning[['unweighted']] <- unweighted
var_partitioning[['intraspecific']] <- intraspecific
var_partitioning[['sp_abundance']] <- sp_abundance

# save
save(var_partitioning, file="results/var_partitioning.RData")


