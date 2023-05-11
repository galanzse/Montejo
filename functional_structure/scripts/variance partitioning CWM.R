
# variance partitioning considering species abundances
library(cati)

source('scripts/PCA_traits.R')


# 1/ Prepare the tables
temp <- read.csv("results/lg_moments_weighted.txt", sep="") %>% subset(moment=='CWM') %>%
  select(-moment) %>%
  pivot_wider(names_from='trait', values_from='value')

temp[,3:5] <- NA
head(temp)
  
v_traits <- colnames(temp)[3:5]

# calculate five types of CWM (parameters) using:

  # 1) specific (forest averages): itv + turnover
specific <- temp
site_means <- pca1 %>% group_by(species, forest) %>%
  summarise(Leaf_Economics=mean(Leaf_Economics),
            Root_Economics=mean(Root_Economics),
            HyArq=mean(HyArq))

  # 2) fixed (region averages): turnover
fixed <- temp
montejo_means <- pca1 %>% group_by(species) %>%
  summarise(Leaf_Economics=mean(Leaf_Economics),
            Root_Economics=mean(Root_Economics),
            HyArq=mean(HyArq))

  # 3) unweighted (presence data): abundance
unweighted <- temp

for (i in 1:nrow(temp)) {
  
  # species
  p1 <- abundances[as.character(temp$plot[i]),] %>% as.data.frame()
  colnames(p1) <- 'abundance'; p1$species <- rownames(p1); rownames(p1) <- NULL
  p1 <- p1 %>% filter(abundance>0)

    for (t in c("Leaf_Economics","Root_Economics","HyArq")) {
    
      # specific
      t1 <- site_means %>% filter(species %in% p1$species & forest == specific$forest[i]) %>% dplyr::select(species, t)
      colnames(t1) <- c('species','trait')
      t1 <- merge(t1, p1, by='species')
      specific[i,t] <- sum(t1$trait * t1$abundance / sum(t1$abundance))
      
      # fixed
      t1 <- montejo_means %>% filter(species %in% p1$species) %>% dplyr::select(species, t)
      colnames(t1) <- c('species','trait')
      t1 <- merge(t1, p1, by='species')
      fixed[i,t] <- sum(t1$trait * t1$abundance / sum(t1$abundance))
      
      # unweighted
      t1 <- montejo_means %>% filter(species %in% p1$species) %>% dplyr::select(species, t)
      colnames(t1) <- c('species','trait')
      t1 <- merge(t1, p1, by='species')
      unweighted[i,t] <- mean(t1$trait)
      
  }
}

  # 4) intraspecific parameter (specific - fixed): itv
intraspecific <- temp
intraspecific[,'Leaf_Economics'] <- specific[,'Leaf_Economics'] - fixed[,'Leaf_Economics']
intraspecific[,'Root_Economics'] <- specific[,'Root_Economics'] - fixed[,'Root_Economics']
intraspecific[,'HyArq'] <- specific[,'HyArq'] - fixed[,'HyArq']

  # 5) species_abundance parameter (fixed – unweighted): abundance
sp_abundance <- temp
sp_abundance[,'Leaf_Economics'] <- fixed[,'Leaf_Economics'] - unweighted[,'Leaf_Economics']
sp_abundance[,'Root_Economics'] <- fixed[,'Root_Economics'] - unweighted[,'HyArq']
sp_abundance[,'HyArq'] <- fixed[,'HyArq'] - unweighted[,'HyArq']

# save
var_partitioning <- list()
var_partitioning[['specific']] <- specific
var_partitioning[['fixed']] <- fixed
var_partitioning[['unweighted']] <- unweighted
var_partitioning[['intraspecific']] <- intraspecific
var_partitioning[['sp_abundance']] <- sp_abundance

# save
save(var_partitioning, file="results/var_partitioning.RData")


# 2/ check the package 'cati'
specific_avg <- var_partitioning[['specific']][,c('plot','forest','Leaf_Economics')]; colnames(specific_avg)[3] <- 'specific_avg'
const_avg <- var_partitioning[['fixed']][,c('plot','forest','Leaf_Economics')]; colnames(const_avg)[3] <- 'const_avg'
d1 <- merge(specific_avg, const_avg, by=c('plot','forest'))
traitflex.anova(~forest, specific_avg, const_avg, data=d1)

specific_avg <- var_partitioning[['specific']][,c('plot','forest','Root_Economics')]; colnames(specific_avg)[3] <- 'specific_avg'
const_avg <- var_partitioning[['fixed']][,c('plot','forest','Root_Economics')]; colnames(const_avg)[3] <- 'const_avg'
d1 <- merge(specific_avg, const_avg, by=c('plot','forest'))
traitflex.anova(~forest, specific_avg, const_avg, data=d1)

specific_avg <- var_partitioning[['specific']][,c('plot','forest','HyArq')]; colnames(specific_avg)[3] <- 'specific_avg'
const_avg <- var_partitioning[['fixed']][,c('plot','forest','HyArq')]; colnames(const_avg)[3] <- 'const_avg'
d1 <- merge(specific_avg, const_avg, by=c('plot','forest'))
traitflex.anova(~forest, specific_avg, const_avg, data=d1)


# 3/ de la Riva et a. 2016 (Oikos) modified from Lepš et al. (2011)
SS_part <- matrix(nrow=7, ncol=3) %>% as.data.frame() # sum of squares
colnames(SS_part) <- c('Leaf_Economics','Root_Economics','HyArq')
rownames(SS_part) <- c("specific","fixed","unweighted","intraspecific","sp_abundance","covSSI","covSSII")

for (p in c('specific','fixed','unweighted','intraspecific','sp_abundance')) {
  t1 <- var_partitioning[[p]]
  for (t in c('Leaf_Economics','Root_Economics','HyArq')) {
  a1 <- lm(deframe(t1[,t]) ~ deframe(t1[,'forest'])) %>% anova()
  SS_part[p,t] <- a1$`Sum Sq`[1]
  }
}

# covariance
SS_part['covSSII',] <- SS_part['fixed',] - SS_part['unweighted',] - SS_part['sp_abundance',]
SS_part['covSSI',] <- SS_part['specific',] - SS_part['fixed',] - SS_part['intraspecific',] - SS_part['covSSII',]

# recalculate specific
SS_part['specific',] <- colSums(SS_part[-which(rownames(SS_part)=='specific'),])

# barplot
SS_part$Leaf_Economics <- SS_part$Leaf_Economics / SS_part$Leaf_Economics[1] * 100
SS_part$Root_Economics <- SS_part$Root_Economics / SS_part$Root_Economics[1] * 100
SS_part$HyArq <- SS_part$HyArq / SS_part$HyArq[1] * 100

# with covariance
SS_part$parameter <- rownames(SS_part)
SS_part <- SS_part %>% subset(!(parameter%in%c('specific','fixed'))) %>% pivot_longer(1:3, names_to='trait')
SS_part <- SS_part %>% subset(trait!='HyArq')

SS_part$trait <- as.factor(SS_part$trait)
levels(SS_part$trait)
SS_part$trait <- factor(SS_part$trait, levels=c("Root_Economics","Leaf_Economics"))

SS_part$parameter[SS_part$parameter=='unweighted'] <- 'occurrence'
SS_part$parameter[SS_part$parameter=='sp_abundance'] <- 'abundance'

SS_part$parameter <- as.factor(SS_part$parameter)
levels(SS_part$parameter)
SS_part$parameter <- factor(SS_part$parameter, levels=c("occurrence","abundance","intraspecific","covSSI","covSSII"))

ggplot(aes(fill=parameter, y=trait, x=value), data=SS_part) + 
  geom_bar(position="stack", stat="identity",colour="black") +
  xlab('% variance explained') + ylab(NULL) +
  theme(legend.position='top', legend.title=element_blank(), legend.box='vertical') +
  scale_fill_grey(start=0.2, end=1)

# without covariance
SS_part <- SS_part %>% subset(!(parameter%in%c('covSSI','covSSII')))

ggplot(aes(fill=parameter, y=trait, x=value), data=SS_part) + 
  geom_bar(position="stack", stat="identity",colour="black") +
  xlab('% variance explained') + ylab(NULL) +
  theme(legend.position='top', legend.title=element_blank(), legend.box='vertical') +
  scale_fill_grey(start=0.2, end=1)
