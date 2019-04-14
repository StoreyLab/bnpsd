# Model and Simulate Admixed Populations with `bnpsd`

The `bnpsd` ("Balding-Nichols Pritchard-Stephens-Donnelly") R package facilitates construction of admixed population structures and simulation of allele frequencies and genotypes from the BN-PSD admixture model.
This model combines the Balding-Nichols (BN) allele frequency model for the intermediate subpopulations with the Pritchard-Stephens-Donnelly (PSD) model of individual-specific admixture proportions.
This model enables the simulation of complex population structures, ideal for illustrating challenges in kinship coefficient and FST estimation.
Note that simulated loci are drawn independently (in linkage equilibrium).

## Installation

The stable version of the package is now on CRAN and can be installed using
```R
install.packages("bnpsd")
```

The current development version can be installed from the GitHub repository using `devtools`:
```R
install.packages("devtools") # if needed
library(devtools)
install_github('StoreyLab/bnpsd', build_opts = c())
```

You can see the package vignette, which has more detailed documentation, by typing this into your R session:
```R
vignette('bnpsd')
```

## Example

This is a quick overview of the main `bnpsd` functions.

Define the population structure (in this case for 1D admixture scenario).
```R
library(bnpsd)
# dimensions of data/model
m <- 10 # number of loci
n <- 5 # number of individuals
k <- 2 # number of intermediate subpops

# define population structure
inbr_subpops <- c(0.1, 0.3) # FST values for k=2 subpopulations
sigma <- 1 # dispersion parameter of intermediate subpops
# admixture proportions from 1D geography
admix_proportions <- admix_prop_1d_linear(n, k, sigma)

# get pop structure parameters of the admixed individuals
coancestry <- coanc_admix(admix_proportions, inbr_subpops) # the coancestry matrix
Fst <- fst(admix_proportions, inbr_subpops) # FST of admixed individuals
```

Draw random allele frequencies and genotypes from this population structure.
```R
# draw all random allele freqs and genotypes
out <- rbnpsd(admix_proportions, inbr_subpops, m)
X <- out$X # genotypes
P <- out$P # IAFs (individual-specific AFs)
B <- out$B # intermediate AFs
pAnc <- out$Pa # ancestral AFs

# OR... draw each vector or matrix separately
# provided for additional flexibility
pAnc <- rpanc(m) # "anc"estral AFs
B <- rpint(pAnc, inbr_subpops) # "int"ermediate AFs
P <- rpiaf(B, admix_proportions) # "IAF"s (individual-specific AFs)
X <- rgeno(P) # "geno"types
```

## Citations

Ochoa, Alejandro, and John D. Storey. 2016a. "FST And Kinship for Arbitrary Population Structures I: Generalized Definitions." bioRxiv [doi:10.1101/083915](http://doi.org/10.1101/083915).

Ochoa, Alejandro, and John D. Storey. 2016b. "FST And Kinship for Arbitrary Population Structures II: Method of Moments Estimators." bioRxiv [doi:10.1101/083923](http://doi.org/10.1101/083923).
