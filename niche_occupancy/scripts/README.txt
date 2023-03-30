#### Scripts for "Living at the edge: the functional niche occupation of woody plant communities in the Submediterranean ecotone" ####

Scripts should be run in this order:

1/ import data.R

Import trait and presence data, NA imputation

2/ TOA.R

Import and adapt function from Li et al. 2018

3/ bw_optimisation.R

Optimise bandwidth for our dataset

4/ null models TOA.R

Compute observed values of T, O nad A; and run two null models: using the species present at a forest type, or all species in Montejo.

5/ null models exploratory.R

Plots and SES of null models.

6/ abiotic.R

Import and explore abiotic data.

7/ intraspecific intercom.R

Estimate the relative importance of niche shifts and expansion following Benvides et al. 2019.


REFERENCES
Li Y, Shipley B, Price JN et al. (2018) Habitat filtering determines the functional niche occupancy of plant communities worldwide. Journal of Ecology 106: 1001–1009. https://doi.org/10.1111/1365-2745.12802 

Benavides R, Scherer-Lorenzen M, Valladares F (2019) The functional trait space of tree species is influenced by the species richness of the canopy and the type of forest. Oikos 128: 1435–1445. https://doi.org/10.1111/oik.06348 