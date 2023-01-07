
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
traits$leaf_N[traits$leaf_N %in% boxplot.stats(traits$leaf_N)$out] <- NA
traits$leaf_CN[traits$leaf_CN %in% boxplot.stats(traits$leaf_CN)$out] <- NA
HV_thick$HubVal[HV_thick$HubVal > 0.0015] <- NA
root_traits$root_C[root_traits$root_C %in% boxplot.stats(root_traits$root_C)$out] <- NA
root_traits$root_N[root_traits$root_N %in% boxplot.stats(root_traits$root_N)$out] <- NA
root_traits$root_CN[root_traits$root_CN %in% boxplot.stats(root_traits$root_CN)$out] <- NA

# pairs
t1 <- merge(aggregate(.~species+forest, traits, mean),
            aggregate(.~species+forest, root_traits, mean), by=c('species','forest'), all=T)
# pairs(t1[,3:18], upper.panel=NULL)
# pairs(HV_thick[,3:4], upper.panel=NULL)

# remove correlated variables
traits <- traits %>% dplyr::select(-c(SRA,Rdi,leaf_CN,leaf_d15N))
root_traits <- root_traits %>% dplyr::select(-root_CN)
HV_thick

# par(mfrow=c(4,4))
# for (i in colnames(traits)[3:11]) { hist(deframe(traits[,i]), 30, main=i) }
# for (i in colnames(root_traits)[3:4]) { hist(deframe(root_traits[,i]), 30, main=i) }
# for (i in colnames(HV_thick)[3:4]) { hist(deframe(HV_thick[,i]), 30, main=i) }
