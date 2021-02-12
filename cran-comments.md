## Test environments
* local R:             x86_64-redhat-linux-gnu R 4.0.3 (2020-10-10)
* local R-devel:       x86_64-redhat-linux-gnu R Under development (unstable) (2021-02-09 r79976)
* local R-devel-noLD:  x86_64-redhat-linux-gnu R Under development (unstable) (2021-02-09 r79976)

TODO

* win-builder devel:   x86_64-w64-mingw32      R Under development (unstable) (2020-01-03 r77630)
* win-builder release: x86_64-w64-mingw32      R 3.6.2 (2019-12-12)
* rhub:       (debian) x86_64-pc-linux-gnu     R Under development (unstable) (2020-01-03 r77629)

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs. 

## Downstream dependencies
Tested `popkin` and `ggmix` (only downstream dependencies) and found no errors or warnings related to `popkin`.

* The latest `ggmix` (from GitHub) failed a test unrelated to `popkin`, which I also found reported on CRAN:
  * Failure (test-KKT.R:32:3): Check predict and coef methods with multiple s values
    all(abs(kkt)[-1] < 0.02) is not TRUE
  * https://www.stats.ox.ac.uk/pub/bdr/M1mac/ggmix.out
