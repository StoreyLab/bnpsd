# 2017-09-11 - bnpsd 1.0.0.9000

* Public release!

# 2018-01-15 - bnpsd 1.0.1

* Minor non-code changes for first CRAN submission.

# 2018-02-01 - bnpsd 1.0.1.9000

* README.md now contains instructions for installing from CRAN as well as from GitHub.

# 2019-02-06 - bnpsd 1.0.2.9000

* Added option `noFixed` to function `rbnpsd` to redraw loci that were drawn fixed for a single allele.
These loci are not polymorphic so they would normally not be considered in analyses.

* Added function `fixed_loci` to test for fixed loci within rbnpsd.

# 2019-02-07 - bnpsd 1.0.3.9000

* Added function `coanc_to_kinship` to easily obtain kinship matrices from coancestry matrices.

# 2019-02-11 - bnpsd 1.0.4

* Second CRAN submission.

# 2019-02-13 - bnpsd 1.0.4.9000

* Converted the vignette from PDF to HTML

# 2019-04-11 - bnpsd 1.0.5.9000

* `qis` now returns a numeric admixture proportions matrix (used to be logical).
* `q1d` and `q1dc` now handle `sigma = 0` special case.
* `q1d` and `q1dc` now provide more informative out-of-bounds messages when `sigma` is missing (and `s` is provided)
* `sigma` root finding in `q1d` and `q1dc` (when `s` is provided) is now more robust, explicitly tested at boundaries (min `s > 0` achieved at `sigma = 0` and max `s = 1` achieved at `sigma = Inf`).
  * Removed arguments `interval` and `tol` from both `q1d` and `q1dc` (users would never need to set them now that procedure is more robust).
* Updated coding style, renamed some internal functions and variables.

# 2019-04-16 - bnpsd 1.1.0.9000

* Renamed most functions for clarity:
  * `coanc` -> `coanc_admix`
  * `q1d` -> `admix_prop_1d_linear`
  * `q1dc` -> `admix_prop_1d_circular`
  * `qis` -> `admix_prop_indep_subpops`
  * `rpanc` -> `draw_p_anc`
  * `rpint` -> `draw_p_subpops`
  * `rpiaf` -> `make_p_ind_admix`
  * `rgeno` -> `draw_genotypes_admix`
  * `rbnpsd` -> `draw_all_admix`
  * `fst` -> `fst_admix` (no deprecated version available in this case, to eliminate conflict with `popkin::fst`)
  * Functions with old names remain for now as deprecated functions (to be removed in the future).
* Renamed several recurrent argument names for clarity:
  * `Q` -> `admix_proportions`
  * `F` -> `coanc_subpops` (if general matrix is accepted), `inbr_subpops` (vector or scalar versions required)
  * `s` -> `bias_coeff`
  * `w` -> `weights`
  * `Theta` -> `coancestry`
  * `m` -> `m_loci`
  * `n` -> `n_ind`
  * `k` -> `k_subpops`
  * `pAnc` -> `p_anc`
  * `B` -> `p_subpops`
  * `P` -> `p_ind`
  * Deprecated functions still accept old argument names.
* Fixed a `sigma = 0` bug in `admix_prop_1d_circular`.
* Changed default values for `draw_all_admix` (compared to deprecated `rbnpsd`, which retains old defaults):
  * `require_polymorphic_loci` (old `noFixed`) is now `TRUE` by default.
  * `want_p_ind` and `want_p_subpops` (old `wantP` and `wantB`) are now `FALSE` by default.
  * Names (following above conventions) and order of items in return list changed.
* `draw_p_subpops` now admits scalar inputs `p_anc` and `inbr_subpops`, while number of loci and number of subpopulations can be provided as additional options.
* Added more input checks to functions, informative error messages.
* Updated documentation, particularly on whether intermediate subpopulation coancestries are accepted generally (`coanc_subpops`) or if the diagonal matrix case is required (specified as vector or scalar `inbr_subpops`).

# 2019-05-15 - bnpsd 1.1.1

* Third CRAN submission.
* Added ORCIDs to authors.
* Corrected doc typos.
* Adjusted layout of subpopulations and individuals (default limits) for circular 1D geography (`admix_prop_1d_circular`) to prevent overlapping individuals on the edges, and to better agree visually with the linear version (`admix_prop_1d_linear`).

# 2019-06-05 - bnpsd 1.1.2

* Non-code changes:
  * Edited .Rbuildignore to stop ignoring README; also removed non-existent files from list
  * Removed unused .travis.yml and bnpsd.Rproj files

# 2019-08-13 - bnpsd 1.1.2.9000

* Improved memory efficiency of default `draw_genotypes_admix`
  * Old approach was by default very memory-hungry (created IAF matrix whole when admixture proportions were provided).
    The `low_mem` option could be set but filled slowly by locus only.
  * New approach is always low-memory (so the explicit option was removed).
    This was made faster by filling by individual when there are fewer individuals than loci, or filling by locus otherwise, therefore always vectorizing as much as possible.
	Test showed this was always as fast as the original full memory approach, so the latter was removed as an option.
* `draw_all_admix` is also now automatically low-memory whenever `want_p_ind = FALSE`, and the explicit `low_mem` option has also been removed.
* Updated documentation to use RMarkdown
* Other code tidying

# 2019-09-06 - bnpsd 1.1.3.9000

* Added option `beta` in function `draw_p_anc` to trigger a symmetric Beta distribution for the ancestral allele frequencies, with the desired shape parameter.
  The `beta` option can also be set on the wrapper function `draw_all_admix`.
  This option allows simulation of a distribution heavier on rare variants (when `beta` is much smaller than 1), more similar to real human data.

# 2019-12-17 - bnpsd 1.2.0

* Fourth CRAN submission.
* Removed deprecated function names: `q1dc`, `q1d`, `qis`, `coanc`, `rbnpsd`, `rgeno`, `rpanc`, `rpint`, `rpiaf`.
* Moved logo to `man/figures/`
* Minor Roxygen-related updates.

# 2020-01-08 - bnpsd 1.2.1

* Fourth CRAN submission, second attempt.
* Fixed a rare bug in `bias_coeff_admix_fit`, which caused it to die if the desired bias coefficient was an extreme value (particularly `1`).
  The error message was: `f() values at end points not of opposite sign`.
  The actual bug was not observed in the regular R build, but rather in a limited precision setting where R was configured with `--disable-long-double`.

# 2020-01-08 - bnpsd 1.2.1.9000

* Added option `p_anc` to function `draw_all_admix`, to specify desired ancestral allele frequencies instead of having the code generate it randomly (default).
* Added details for documentation of function `draw_p_subpops.R`, clarifying that input `p_anc` can be scalar.

# 2021-01-21 - bnpsd 1.2.2.9000

* Function `draw_all_admix`: when option `p_anc` is provided as scalar and `want_p_anc = TRUE`, now the return value is always a vector (in this case the input scalar value repeated `m_loci` times).  The previous behavior was to return `p_anc` as scalar if that was the input, which could be problematic for downstream applications.

# 2021-02-11 - bnpsd 1.2.3

* 5th CRAN submission
* Functions `admix_prop_1d_linear` and `admix_prop_1d_circular` had these changes:
  - The optional parameters `bias_coeff`, `coanc_subpops` and `fst` now have default values (of `NA`, `NULL`, and `NA`, respectively) instead of missing, and these "missing" values can be passed to get the same behavior as if they hadn't been passed at all.
  - Their documentation has been clarified.
  - Improved internal code to handle edge case `bias_coeff = 1` (to fix an issue only observed on Apple M1).
* Function `admix_prop_indep_subpops`: default value for the optional parameter `subpops` is now made more clear in arguments definition.
* Simplified documentation (most functions) by clarifying language, using markdown roxygen, and replacing all LaTeX equations with simpler code equations.
* Updated paper citations in `DESCRIPTION`, `README.md` and the vignette, to point to the published method in PLoS Genetics.

# 2021-02-16 - bnpsd 1.2.3.9000

* Documentation updates:
  - Fixed links to functions, in many cases these were broken because of incompatible mixed Rd and markdown syntax (now markdown is used more fully).

# 2021-03-24 - bnpsd 1.3.0.9000

- Added support for intermediate subpopulations related by a tree
  - New function `draw_p_subpops_tree` is the tree version of `draw_p_subpops`.
  - New function `coanc_tree` calculates the true coancestry matrix corresponding to the subpopulations related by a tree.
  - Function `draw_all_admix` has new argument `tree_subpops` that can be used in place of `inbr_subpops` (to simulated subpopulation allele frequencies using `draw_p_subpops_tree` instead of `draw_p_subpops`).
  - Note: These other functions work for trees (without change) because they accept arbitrary coancestry matrices (param `coanc_subpops`) as input, so they work if they are passed the matrix that `coanc_tree` returns: `coanc_admix`, `fst_admix`, `admix_prop_1d_linear`, `admix_prop_1d_circular`.
- Functions `admix_prop_1d_linear` and `admix_prop_1d_circular`, when `sigma` is missing (and therefore fit to a desired `coanc_subpops`, `fst`, and `bias_coeff`), now additionally return multiplicative `factor` used to rescale `coanc_subpops`.
