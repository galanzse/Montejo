
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
presences[presences>0] <- 1

# traits
traits <- read_excel("data2.xlsx", sheet = "traits")
traits$forest <- as.factor(traits$forest)
str(traits)

# hubervalue & thickness
HV_thick <- read_excel("data2.xlsx", sheet = "HubVal_Thickness")
HV_thick$forest <- as.factor(HV_thick$forest)
str(HV_thick)

# roots, use average per site because most samples were averaged in the lab
root_traits <- read_excel("data2.xlsx", sheet = "root_traits") %>% group_by(species, forest) %>%
  summarise(root_C = mean(root_C, na.rm=T),
            root_N = mean(root_N, na.rm=T),
            root_CN = mean(root_CN, na.rm=T) )
root_traits$forest <- as.factor(root_traits$forest)
str(root_traits)

# explore number of NAs per variable
colSums(is.na(traits))
colSums(is.na(HV_thick))
colSums(is.na(root_traits))

# remove outliers
traits$LDMC[traits$LDMC %in% boxplot.stats(traits$LDMC)$out] <- NA
traits$SLA[traits$SLA %in% boxplot.stats(traits$SLA)$out] <- NA
traits$SDMC[traits$SDMC %in% boxplot.stats(traits$SDMC)$out] <- NA
traits$RDMC[traits$RDMC %in% boxplot.stats(traits$RDMC)$out] <- NA
traits$SRL[traits$SRL %in% boxplot.stats(traits$SRL)$out] <- NA
traits$SRA[traits$SRA %in% boxplot.stats(traits$SRA)$out] <- NA
traits$TMDr[traits$TMDr %in% boxplot.stats(traits$TMDr)$out] <- NA
traits$Rdi[traits$Rdi %in% boxplot.stats(traits$Rdi)$out] <- NA
traits$leaf_d13C[traits$leaf_d13C %in% boxplot.stats(traits$leaf_d13C)$out] <- NA
traits$leaf_d15N[traits$leaf_d15N %in% boxplot.stats(traits$leaf_d15N)$out] <- NA
traits$leaf_C[traits$leaf_C %in% boxplot.stats(traits$leaf_C)$out] <- NA
traits$leaf_N[traits$leaf_N > 4.5] <- NA # less sensitive threshold
traits$leaf_CN[traits$leaf_CN %in% boxplot.stats(traits$leaf_CN)$out] <- NA
HV_thick$HubVal[HV_thick$HubVal > 0.0015] <- NA  # less sensitive threshold
root_traits$root_C[root_traits$root_C %in% boxplot.stats(root_traits$root_C)$out] <- NA
root_traits$root_N[root_traits$root_N < 0.7 | root_traits$root_N > 1.6] <- NA
root_traits$root_CN[root_traits$root_CN %in% boxplot.stats(root_traits$root_CN)$out] <- NA

# pairs
t1 <- merge(aggregate(.~species+forest, traits, mean),
            aggregate(.~species+forest, root_traits, mean), by=c('species','forest'), all=T)
# pairs(t1[,3:18], upper.panel=NULL)
# pairs(HV_thick[,3:4], upper.panel=NULL)

# select: LDMC, SLA, SDMC, Hv, RDMC, SRA, Rdiam, d13C, leaf_CN y y root_CN
traits <- traits %>% dplyr::select(species, forest, LDMC, SLA, SDMC, RDMC, SRA, Rdi, leaf_d13C, leaf_CN) %>% as.data.frame()
root_traits <- root_traits %>% dplyr::select(species, forest, root_CN) %>% as.data.frame() %>% na.omit()
HV_thick <- HV_thick %>% dplyr::select(species, forest, HubVal) %>% as.data.frame() %>% na.omit()

# par(mfrow=c(4,4))
# for (i in colnames(traits)[3:11]) { hist(deframe(traits[,i]), 30, main=i) }
# for (i in colnames(root_traits)[3:4]) { hist(deframe(root_traits[,i]), 30, main=i) }
# for (i in colnames(HV_thick)[3:4]) { hist(deframe(HV_thick[,i]), 30, main=i) }


# explore the cumulative abundance we can consider for assemble without imputation
sites$cum_cover <- NA
for (i in 1:nrow(sites)) {
  
  # select plot
  p1 <- abundances[as.character(sites$plot[i]),]
  p1 <- p1[p1>0]
  # compute relative abundances
  p1 <- p1/sum(p1)
  
  # species with trait data available
  sp <- traits %>% filter(species %in% names(p1) & forest == sites$forest[i]) %>%
    dplyr::select(species) %>% unique %>% deframe()
  
  # calculate cumulative cover for species with trait data
  sites$cum_cover[i] <- sum(p1[sp])
  
}

cumulative_cover <- sites
write.table(cumulative_cover, 'results/cumulative_cover.txt')
