
# adapt unction of Li et al

# Function to calculate for the metrics (T, O and A) for a given community
# x is dataframe of species' names and their trait values for a given community
# bw is the bandwidth
TOA <- function(x, bw, sxp) {
  n <- ncol(x)
  # list the species in community x
  sp <- table(x$code)
  # a vector store the T, O and A
  out <- numeric()
  # a vector to store the functional volumes of each species in community x
  vol <- numeric()
  # the trait values of the first species
  x1 <- droplevels(x[x$code==names(sp)[1], trait.ana])
  # functional volume of the first species
  hv1 <- hypervolume_box(data=x1, samples.per.point=sxp, kde.bandwidth=bw, verbose=F)
  vol[1] <- hv1@Volume
  # union is calculated in a sequential way
  for(j in 2:length(sp)) {
    # the trait values of each following species in community x
    x2 <- droplevels(x[x$code==names(sp)[j], trait.ana])
    # functional volumes of the following species in community x
    hv2 <- hypervolume_box(data=x2, samples.per.point=sxp, kde.bandwidth=bw, verbose=F)
    vol[j] <- hv2@Volume
    # The union of the union of functional volumes for previous species and functional volume of next species 
    hv1 <- hypervolume_set(hv1, hv2, verbose=F, check.memory=F)@HVList[[4]]
  }
  # total functional volume (T)
  out[1]<-hv1@Volume
  # functional overlap (O)
  out[2]<-sum(vol)-out[1]
  # average functional volume (A)
  out[3]<-mean(vol)
  return(out)
}
