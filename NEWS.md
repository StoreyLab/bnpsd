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
