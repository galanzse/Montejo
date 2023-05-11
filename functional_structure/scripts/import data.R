
library(tidyverse)
library(readxl)
library(GGally)

setwd('C:/Users/user/OneDrive/PUBLICACIONES/BTU/montejo/functional_structure')

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
presences[presences>0] <- 1


# traits
traits <- read_excel("data2.xlsx", sheet = "traits")
traits$indiv <- NULL
traits$forest <- as.factor(traits$forest)
str(traits)

# boxplot(scale(traits[,3:15]), range=2)
traits$LDMC[traits$LDMC %in% boxplot.stats(traits$LDMC, coef=2)$out] <- NA
traits$SLA[traits$SLA %in% boxplot.stats(traits$SLA, coef=2)$out] <- NA
traits$SDMC[traits$SDMC %in% boxplot.stats(traits$SDMC, coef=2)$out] <- NA
traits$RDMC[traits$RDMC %in% boxplot.stats(traits$RDMC, coef=2)$out] <- NA
traits$SRL[traits$SRL %in% boxplot.stats(traits$SRL, coef=2)$out] <- NA
traits$SRA[traits$SRA %in% boxplot.stats(traits$SRA, coef=2)$out] <- NA
traits$TMDr[traits$TMDr %in% boxplot.stats(traits$TMDr, coef=2)$out] <- NA
traits$Rdi[traits$Rdi %in% boxplot.stats(traits$Rdi, coef=2)$out] <- NA
traits$leaf_d13C[traits$leaf_d13C %in% boxplot.stats(traits$leaf_d13C, coef=2)$out] <- NA
traits$leaf_d15N[traits$leaf_d15N %in% boxplot.stats(traits$leaf_d15N, coef=2)$out] <- NA
traits$leaf_C[traits$leaf_C %in% boxplot.stats(traits$leaf_C, coef=2.2)$out] <- NA #
traits$leaf_N[traits$leaf_N %in% boxplot.stats(traits$leaf_N, coef=2.2)$out] <- NA #
traits$leaf_CN[traits$leaf_CN %in% boxplot.stats(traits$leaf_CN, coef=2)$out] <- NA


# hubervalue, measured at the regional level
HV_thick <- read_excel("data2.xlsx", sheet = "HubVal_Thickness") %>% group_by(species) %>%
  summarise(HubVal = mean(HubVal, na.rm=T))
str(HV_thick)

# boxplot(HV_thick$HubVal, range=2)
HV_thick$HubVal[HV_thick$HubVal %in% boxplot.stats(HV_thick$HubVal, coef=2)$out] <- NA


# roots, measured polling samples per site
root_traits <- read_excel("data2.xlsx", sheet = "root_traits") %>% group_by(species, forest) %>%
  summarise(root_C = mean(root_C, na.rm=T),
            root_N = mean(root_N, na.rm=T),
            root_CN = mean(root_CN, na.rm=T) )
root_traits$forest <- as.factor(root_traits$forest)
str(root_traits)

# boxplot(scale(root_traits[,3:5]), range=2)
root_traits$root_C[root_traits$root_C %in% boxplot.stats(root_traits$root_C, coef=2)$out] <- NA
root_traits$root_N[root_traits$root_N %in% boxplot.stats(root_traits$root_N, coef=2)$out] <- NA
root_traits$root_CN[root_traits$root_CN %in% boxplot.stats(root_traits$root_CN, coef=2)$out] <- NA


# explore number of NAs per variable
colSums(is.na(traits))
colSums(is.na(HV_thick))
colSums(is.na(root_traits))


# merge tables
traits2 <- merge(traits, root_traits, by=c('species','forest'), all.x=T)
traits2 <- merge(traits2, HV_thick, by='species', all.x=T)
traits <- traits2
rm(traits2, root_traits, HV_thick)

# trait correlation: ggpairs with species means
# ggpairs(traits[,3:19])

# remove malus
traits <- traits %>% subset(species!='malus')

# explore the cumulative abundance we can consider for assemble without imputation
for (i in 1:nrow(sites)) {
  
  # select plot
  p1 <- abundances[as.character(sites$plot[i]),]
  p1 <- p1[p1>0]
  # compute relative abundances
  p1 <- p1/sum(p1)
  
  # species with trait data available
  sp <- traits %>% filter(species %in% names(p1) & forest == sites$forest[i]) %>%
    dplyr::select(species) %>% unique %>% deframe()
  
  # calculate cumulative cover for species with trait data, and nº of missing species
  sites$cum_cover[i] <- sum(p1[sp])
  sites$sp_missing[i] <- length(p1) - length(sp)
}

cumulative_cover <- sites
write.table(cumulative_cover, 'results/cumulative_cover.txt')

# impute q.pyr en shrubland
q.pyr <- traits %>% subset(species=='q.pyr') %>% sample_n(5)
q.pyr$forest <- 'Shrubland'
traits <- rbind(traits, q.pyr)
rm(q.pyr)

sites$cum_cover <- NULL
sites$sp_missing <- NULL
