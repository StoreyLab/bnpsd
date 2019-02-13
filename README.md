Model and Simulate Admixed Populations with bnpsd
===

The `bnpsd` ("Balding-Nichols Pritchard-Stephens-Donnelly") R package facilitates construction of admixed population structures and simulation of allele frequencies and genotypes from the BN-PSD admixture model.
This model combines the Balding-Nichols (BN) allele frequency model for the intermediate subpopulations with the Pritchard-Stephens-Donnelly (PSD) model of individual-specific admixture proportions.
This model enables the simulation of complex population structures, ideal for illustrating challenges in kinship coefficient and FST estimation.
Note that simulated loci are drawn independently (in linkage equilibrium).

Installation
===

The stable version of the package is now on CRAN and can be installed using
```R
install.packages("bnpsd")
```

The current development version can be installed from the GitHub repository using `devtools`:
```R
install.packages("devtools") # if needed
library(devtools)
devtools::install_github('StoreyLab/bnpsd', build_opts=c())
```

You can see the package vignette, which has additional documentation, by typing this into your R session:
``` r
vignette('bnpsd')
```

Synopsis of commands
===

This is a quick overview of the main `bnpsd` functions.

Define the population structure (in this case for 1D admixture scenario).
```R
library(bnpsd)
# dimensions of data/model
m <- 10 # number of loci
n <- 5 # number of individuals
k <- 2 # number of intermediate subpops

# define population structure
F <- c(0.1, 0.3) # FST values for k=2 subpopulations
sigma <- 1 # dispersion parameter of intermediate subpops
Q <- q1d(n, k, sigma) # admixture proportions from 1D geography

# get pop structure parameters of the admixed individuals
Theta <- coanc(Q,F) # the coancestry matrix
Fst <- fst(Q,F) # Fst
```

Draw random allele frequencies and genotypes from this population structure.
```R
# draw all random allele freqs and genotypes
out <- rbnpsd(Q, F, m)
X <- out$X # genotypes
P <- out$P # IAFs (individual-specific AFs)
B <- out$B # intermediate AFs
pAnc <- out$Pa # ancestral AFs

# OR... draw each vector or matrix separately
# provided for additional flexibility
pAnc <- rpanc(m) # "anc"estral AFs
B <- rpint(pAnc, F) # "int"ermediate AFs
P <- rpiaf(B, Q) # "IAF"s (individual-specific AFs)
X <- rgeno(P) # "geno"types
```


More details
===

Please see the bnpsd vignette in R for a description of the key parameters and more detailed examples.

Citations
===

Ochoa, Alejandro, and John D. Storey. 2016a. "FST And Kinship for Arbitrary Population Structures I: Generalized Definitions." bioRxiv [doi:10.1101/083915](http://doi.org/10.1101/083915).

Ochoa, Alejandro, and John D. Storey. 2016b. "FST And Kinship for Arbitrary Population Structures II: Method of Moments Estimators." bioRxiv [doi:10.1101/083923](http://doi.org/10.1101/083923).
