
# Import and explore abiotic data

library(readxl)
library(tidyverse)
library(ggforce)
library(vegan)

# traits and abundances
source('scripts/import data.R')

# import abiotic data
abiotic <- read_excel("data.xlsx", sheet = "abiotic")
colnames(abiotic)[colnames(abiotic)=='bosque'] <- 'forest'
str(abiotic)
abiotic$forest <- as.factor(abiotic$forest)
levels(abiotic$forest) <- c('Fagus sylvatica','Shrubland','Mixed forest 2','Quercus pyrenaica','Mixed forest 1','Quercus petraea')
abiotic$forest <- factor(abiotic$forest, levels=c("Fagus sylvatica","Quercus petraea","Mixed forest 1","Mixed forest 2","Quercus pyrenaica","Shrubland"))

#  correlation among variables
pairs(abiotic[,3:13], lower.panel=NULL, cex.labels=3)

# PCA
abiot.pca <- prcomp(abiotic[,3:10], scale.=TRUE, center=TRUE)
summary(abiot.pca)
eigenvals(abiot.pca)
fviz_pca_var(abiot.pca, col.var="contrib", repel=TRUE)

# dataframe
abiotic <- cbind(abiotic, abiot.pca$x[,1:3])


# community function ~ abiotic data
out_all <- read.csv("resultados/out_all.txt", sep="")
commxabiotic <- merge(abiotic, out_all, by=c('plot','forest')) # merge


ggplot(aes(y=T.ob, x=PC1, color=forest), data=commxabiotic) + geom_point()
ggplot(aes(y=O.ob, x=PC1, color=forest), data=commxabiotic) + geom_point()
ggplot(aes(y=A.ob, x=PC1, color=forest), data=commxabiotic) + geom_point()

ggplot(aes(y=T.SES, x=PC1, color=forest), data=commxabiotic) + geom_point()
ggplot(aes(y=O.SES, x=PC1, color=forest), data=commxabiotic) + geom_point()
ggplot(aes(y=A.SES, x=PC1, color=forest), data=commxabiotic) + geom_point()
