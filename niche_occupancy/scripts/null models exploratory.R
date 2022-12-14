

library(tidyverse)
library(ggpubr)
source('scripts/import data.R')


# import data
out_forest <- read.table("resultados/out_forest.txt")
out_cocc <- read.table("resultados/out_cocc.txt")
out_all <- read.table("resultados/out_all.txt")


# observed values
g1 <- ggplot(aes(y=T.ob, x=S), data=out_forest) +
  geom_point() + theme_classic() +
  labs(y='Total Functional Volume (T)', x='Species Richness (S)') +
  geom_smooth(method='lm') +
  scale_x_continuous(breaks=seq(0,12,2))
g2 <- ggplot(aes(y=O.ob, x=S), data=out_forest) +
  geom_point() + theme_classic() + ylim(0,15) +
  labs(y='Functional Overlap (O)', x='Species Richness (S)') +
  geom_smooth(method='lm') +
  scale_x_continuous(breaks=seq(0,12,2))
g3 <- ggplot(aes(y=A.ob, x=S), data=out_forest) +
  geom_point() + theme_classic() +
  labs(y='Average Functional Volume (A)', x='Species Richness (S)') +
  geom_smooth(method='lm') +
  scale_x_continuous(breaks=seq(0,12,2))

ggarrange(g1, g2, g3, ncol=3, labels=c('(a)','(b)','(c)'), hjust=0.05)

lm(A.ob ~ S, data=out_forest) %>% anova()


# niche differentiation
# out_forest$ND.ob <- out_forest$T.ob - out_forest$O.ob
# ggplot(aes(y=ND.ob, x=S), data=out_forest) +
#   geom_point() + theme_classic() +
#   labs(y='Niche Differentiation (ND)', x='Species Richness (S)') +
#   geom_smooth(method='lm') +
#   scale_x_continuous(breaks=seq(0,12,2))


# observed values
out_merged <- out_forest
out_merged$forest <- as.factor(out_merged$forest) # organize factor variables
out_merged$forest <- factor(out_merged$forest, levels=c("Fagus sylvatica","Quercus petraea","Mixed forest 1","Mixed forest 2","Quercus pyrenaica","Shrubland"))

g1 <- ggplot(aes(y=T.ob, x=forest, fill=forest), data=out_merged) +
  geom_boxplot() + theme_classic() +
  labs(y='T observed', x=NULL) +
  theme(legend.position="none", legend.title=element_blank())
g2 <- ggplot(aes(y=O.ob, x=forest, fill=forest), data=out_merged) +
  geom_boxplot() + theme_classic() +
  labs(y='O observed', x=NULL) +
  theme(legend.position="none", legend.title=element_blank())
g3 <- ggplot(aes(y=A.ob, x=forest, fill=forest), data=out_merged) +
  geom_boxplot() + theme_classic() +
  labs(y='A observed', x=NULL) +
  theme(legend.position="none", legend.title=element_blank())

ggarrange(g1, g2, g3, hjust=0.03,
          labels=c('(a)','(b)','(c)'), nrow=3)


# expected values
out_merged <- rbind(out_forest, out_all) # merge
SES_long <- out_merged %>% pivot_longer(13:15, names_to='SES', values_to='value') # long format
str(SES_long)

SES_long$forest <- as.factor(SES_long$forest) # organize factor variables
SES_long$forest <- factor(SES_long$forest, levels=c("Fagus sylvatica","Quercus petraea","Mixed forest 1","Mixed forest 2","Quercus pyrenaica","Shrubland"))

SES_long$null <- as.factor(SES_long$null)
SES_long$null <- factor(SES_long$null, levels=c("Forest","All"))
levels(SES_long$null) <- c('Community pool', 'Regional pool')

SES_long$SES <- as.factor(SES_long$SES)
levels(SES_long$SES) <- c('A','O','T')
SES_long$SES <- factor(SES_long$SES, levels=c("T","O","A"))


ggplot(aes(y=value, x=forest, fill=null), data=SES_long) +
  geom_boxplot() + theme_classic() +
  geom_hline(yintercept=0) +
  labs(y='Standardized Effect Size (SES)', x=NULL) +
  facet_grid(rows=vars(SES)) +
  theme(legend.position="top", legend.title=element_blank())
  

# statistical differences
out2 <- read.table("TOA_plot", header = TRUE, sep="", dec=".")
out_wil <- merge(abiotico,out2, by="plot") 

out_sp<-out_wil %>%  filter(bosque=="robledal")
a<-wilcox.test(out_sp$A.ob,out_sp$A.exp, paired=T)
a2<-t.test(out_sp$A.ob,out_sp$A.exp, paired=T)

