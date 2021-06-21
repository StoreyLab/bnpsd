# bnpsd 1.0.0.9000 (2017-09-11)

* Public release!

# bnpsd 1.0.1 (2018-01-15)

* Minor non-code changes for first CRAN submission.

# bnpsd 1.0.1.9000 (2018-02-01)

* README.md now contains instructions for installing from CRAN as well as from GitHub.

# bnpsd 1.0.2.9000 (2019-02-06)

* Added option `noFixed` to function `rbnpsd` to redraw loci that were drawn fixed for a single allele.
These loci are not polymorphic so they would normally not be considered in analyses.

* Added function `fixed_loci` to test for fixed loci within rbnpsd.

# bnpsd 1.0.3.9000 (2019-02-07)

* Added function `coanc_to_kinship` to easily obtain kinship matrices from coancestry matrices.

# bnpsd 1.0.4 (2019-02-11)

* Second CRAN submission.

# bnpsd 1.0.4.9000 (2019-02-13)

* Converted the vignette from PDF to HTML

# bnpsd 1.0.5.9000 (2019-04-11)

* `qis` now returns a numeric admixture proportions matrix (used to be logical).
* `q1d` and `q1dc` now handle `sigma = 0` special case.
* `q1d` and `q1dc` now provide more informative out-of-bounds messages when `sigma` is missing (and `s` is provided)
* `sigma` root finding in `q1d` and `q1dc` (when `s` is provided) is now more robust, explicitly tested at boundaries (min `s > 0` achieved at `sigma = 0` and max `s = 1` achieved at `sigma = Inf`).
  * Removed arguments `interval` and `tol` from both `q1d` and `q1dc` (users would never need to set them now that procedure is more robust).
* Updated coding style, renamed some internal functions and variables.

# bnpsd 1.1.0.9000 (2019-04-16)

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

# bnpsd 1.1.1 (2019-05-15)

* Third CRAN submission.
* Added ORCIDs to authors.
* Corrected doc typos.
* Adjusted layout of subpopulations and individuals (default limits) for circular 1D geography (`admix_prop_1d_circular`) to prevent overlapping individuals on the edges, and to better agree visually with the linear version (`admix_prop_1d_linear`).

# bnpsd 1.1.2 (2019-06-05)

* Non-code changes:
  * Edited .Rbuildignore to stop ignoring README; also removed non-existent files from list
  * Removed unused .travis.yml and bnpsd.Rproj files

# bnpsd 1.1.2.9000 (2019-08-13)

* Improved memory efficiency of default `draw_genotypes_admix`
  * Old approach was by default very memory-hungry (created IAF matrix whole when admixture proportions were provided).
    The `low_mem` option could be set but filled slowly by locus only.
  * New approach is always low-memory (so the explicit option was removed).
    This was made faster by filling by individual when there are fewer individuals than loci, or filling by locus otherwise, therefore always vectorizing as much as possible.
	Test showed this was always as fast as the original full memory approach, so the latter was removed as an option.
* `draw_all_admix` is also now automatically low-memory whenever `want_p_ind = FALSE`, and the explicit `low_mem` option has also been removed.
* Updated documentation to use RMarkdown
* Other code tidying

# bnpsd 1.1.3.9000 (2019-09-06)

* Added option `beta` in function `draw_p_anc` to trigger a symmetric Beta distribution for the ancestral allele frequencies, with the desired shape parameter.
  The `beta` option can also be set on the wrapper function `draw_all_admix`.
  This option allows simulation of a distribution heavier on rare variants (when `beta` is much smaller than 1), more similar to real human data.

# bnpsd 1.2.0 (2019-12-17)

* Fourth CRAN submission.
* Removed deprecated function names: `q1dc`, `q1d`, `qis`, `coanc`, `rbnpsd`, `rgeno`, `rpanc`, `rpint`, `rpiaf`.
* Moved logo to `man/figures/`
* Minor Roxygen-related updates.

# bnpsd 1.2.1 (2020-01-08)

* Fourth CRAN submission, second attempt.
* Fixed a rare bug in `bias_coeff_admix_fit`, which caused it to die if the desired bias coefficient was an extreme value (particularly `1`).
  The error message was: `f() values at end points not of opposite sign`.
  The actual bug was not observed in the regular R build, but rather in a limited precision setting where R was configured with `--disable-long-double`.

# bnpsd 1.2.1.9000 (2020-01-08)

* Added option `p_anc` to function `draw_all_admix`, to specify desired ancestral allele frequencies instead of having the code generate it randomly (default).
* Added details for documentation of function `draw_p_subpops.R`, clarifying that input `p_anc` can be scalar.

# bnpsd 1.2.2.9000 (2021-01-21)

* Function `draw_all_admix`: when option `p_anc` is provided as scalar and `want_p_anc = TRUE`, now the return value is always a vector (in this case the input scalar value repeated `m_loci` times).  The previous behavior was to return `p_anc` as scalar if that was the input, which could be problematic for downstream applications.

# bnpsd 1.2.3 (2021-02-11)

* 5th CRAN submission
* Functions `admix_prop_1d_linear` and `admix_prop_1d_circular` had these changes:
  - The optional parameters `bias_coeff`, `coanc_subpops` and `fst` now have default values (of `NA`, `NULL`, and `NA`, respectively) instead of missing, and these "missing" values can be passed to get the same behavior as if they hadn't been passed at all.
  - Their documentation has been clarified.
  - Improved internal code to handle edge case `bias_coeff = 1` (to fix an issue only observed on Apple M1).
* Function `admix_prop_indep_subpops`: default value for the optional parameter `subpops` is now made more clear in arguments definition.
* Simplified documentation (most functions) by clarifying language, using markdown roxygen, and replacing all LaTeX equations with simpler code equations.
* Updated paper citations in `DESCRIPTION`, `README.md` and the vignette, to point to the published method in PLoS Genetics.

# bnpsd 1.2.3.9000 (2021-02-16)

* Documentation updates:
  - Fixed links to functions, in many cases these were broken because of incompatible mixed Rd and markdown syntax (now markdown is used more fully).

# bnpsd 1.3.0.9000 (2021-03-24)

- Added support for intermediate subpopulations related by a tree
  - New function `draw_p_subpops_tree` is the tree version of `draw_p_subpops`.
  - New function `coanc_tree` calculates the true coancestry matrix corresponding to the subpopulations related by a tree.
  - Function `draw_all_admix` has new argument `tree_subpops` that can be used in place of `inbr_subpops` (to simulated subpopulation allele frequencies using `draw_p_subpops_tree` instead of `draw_p_subpops`).
  - Note: These other functions work for trees (without change) because they accept arbitrary coancestry matrices (param `coanc_subpops`) as input, so they work if they are passed the matrix that `coanc_tree` returns: `coanc_admix`, `fst_admix`, `admix_prop_1d_linear`, `admix_prop_1d_circular`.
- Functions `admix_prop_1d_linear` and `admix_prop_1d_circular`, when `sigma` is missing (and therefore fit to a desired `coanc_subpops`, `fst`, and `bias_coeff`), now additionally return multiplicative `factor` used to rescale `coanc_subpops`.

# bnpsd 1.3.1.9000 (2021-04-17)

It's Fangorn Forest around here with all the tree updates!

- Added these functions:
  - `fit_tree` for fitting trees to coancestry matrices!
  - `scale_tree` to easily scale coancestry trees and check for out-of-bounds values.
  - `tree_additive` for calculating "additive" edges for probabilistic edge coancestry trees, and also the reverse function .
    - This already existed as an internal, unexported function used mainly by `coanc_tree`, but now it's renamed, exported, and well documented!
- Added support of `$root.edge` to tree `phylo` objects passed to these functions:
  - `coanc_tree`: edge is a shared covariance value affecting all subpopulations.
  - `draw_all_admix` and `draw_p_subpops_tree`: if root edge is present, functions warn that it will be ignored.
- Functions `admix_prop_1d_linear` and `admix_prop_1d_circular`: debugged an edge case where `sigma` is small but not zero and numerically-calculated densities all come out to zero in a given row of the `admix_proportions` matrix (for `admix_prop_1d_circular` infinite values also arise), which used to lead to NAs upon row normalization; now for those rows, the closest ancestry (by coordinate distance) gets assigned the full admixture fraction (just as for independent subpopulations/`sigma = 0`).

# bnpsd 1.3.2.9000 (2021-04-22)

- Updated various functions to transfer names between inputs and outputs as it makes sense
  - Functions `admix_prop_1d_linear`, `admix_prop_1d_circular` now copy names from the input `coanc_subpops` (vector and matrix versions, only required when fitting `bias_coeff`) to the columns of the output `admix_proportions` matrix.
  - Function `draw_genotypes_admix` now copies row and column names from input matrix `p_ind` (or rownames from `p_ind` and column names from the rownames of `admix_proportions` when the latter is provided) to output genotype matrix
  - Function `draw_p_subpops` now copies names from `p_anc` to rows, names from `inbr_subpops` to columns, when present and of the right dimensions.
  - Function `draw_p_subpops_tree` now copies names from `p_anc` to rows.  Names from `tree_subpops` were already copied to columns before.
  - All other functions already transfered names as desired/appropriate.  Tests were added for these functions to ensure that this is so.
- Updated various functions to stop if there are paired names for two objects that are both non-NULL and disagree, as this suggests that the data is misaligned or incompatible.
  - Functions `coanc_admix` and `fst_admix` stop if the column names of `admix_proportions` and the names of `coanc_subpops` disagree.
  - Function `draw_all_admix` stops if the column names of `admix_proportions` and the names of either `inbr_subpops` or `tree_subpops` disagree.
  - Function `draw_genotypes_admix`, when `admix_proportions` is passed, stops if the column names of `admix_proportions` and `p_ind` disagree.
  - Function `make_p_ind_admix` stops if the column names of `admix_proportions` and `p_subpops` disagree.
- Function `tree_additive` now has option `force`, which when `TRUE` simply proceeds without stopping if additive edges were already present (in `tree$edge.length.add`, which is ignored and overwritten).

# bnpsd 1.3.3.9000 (2021-04-29)

New functions and bug fixes dealing with reordering tree edges and tips.

- Added function `tree_reindex_tips` for ensuring that tip order agrees in both the internal labels vector and the edge matrix.
  Such lack of agreement is generally possible (technically the tree is the same for arbitrary orders of edges in the edge matrix).
  However, such a disagreement causes visual disagreement in plots (for example, trees are plotted in the order of the edge matrix, versus coancestry matrices are ordered as in the tip labels vector instead), which can now be fixed in general.
- Added function `tree_reorder` for reordering tree edges and tips to agree as much as possible with a desired tip order.
  The heuristic finds the exact solution if it exists, otherwise returns a reasonable order close to the desired order.
  Tip order in labels and edge matrix agree (via `tree_reindex_tips`).
- Function `fit_tree` now outputs trees with tip order that better agrees with the input data, and tip order in labels vector and edge matrix now agree (via `tree_reorder`).
- Several functions now work with trees whose edges are arbitrarily ordered, particularly when they do not move out from the root (i.e. reverse postorder):
  - Function `tree_additive`.
    Before this bug fix, some trees could trigger the error message "Error: Node index 6 was not assigned coancestry from root! (unexpected)", where "6" could be other numbers.
  - Function `draw_p_subpops_tree`.
	Before this bug fix, some trees could trigger the error message "Error: The root node index in `tree_subpops$edge` (9) does not match `k_subpops + 1` (6) where `k_subpops` is the number of tips!  Is the `tree_subpops` object malformed?", where "9" and "6" could be other numbers.  Other possible error messages contain "Parent node index 6 has not yet been processed ..." or "Child" instead of "Parent", where "6" could be other numbers.
  - Internal functions used by `fit_tree` had related fixes, but overall `fit_tree` appears to have had no bugs because users cannot provide trees, and the tree-building algorithm does not produce scrambled edges that would have caused problems.

# bnpsd 1.3.4.9000 (2021-05-12)

- Functions `fixed_loci` and `draw_all_admix` have a new parameter `maf_min` that, when greater than zero, allows for treating rare variants as fixed.
  In `draw_all_admix`, this now allows for simulating loci with frequency-based ascertainment bias.

# bnpsd 1.3.5.9000 (2021-05-14)

- Fixed a rare bug in `draw_all_admix` that could cause a "stack overflow" error.
  The function used to call itself recursively if `require_polymorphic_loci = TRUE`, and in cases where there are very rare allele frequencies or high `maf_min` the number of recursions could be so large that it triggered this error.
  Now the function has a `while` loop, and does not recurse more than one level at the time; there is no limit to the number of iterations and no errors occur inherently due to large numbers of iterations.

# bnpsd 1.3.6.9000 (2021-06-02)

- Function `fit_tree` internally simplified to use `stats::hclust`, which also results in a small runtime gain.
  The new code (when `method = "mcquitty"`, which is default) gives the same answers as before (in other words, the original algorithm was a special case of hierarchical clustering).
  - New option `method` is passed to `hclust`.
    Although all `hclust` methods are allowed, for this application the only ones that make sense are "mcquitty" (WPGMA) and "average" (UPGMA).
	In internal evaluations, both algorithms had similar accuracy and runtime, but only "mcquitty" exactly recapitulates the original algorithm.

# bnpsd 1.3.7.9000 (2021-06-04)

- Updated citations in `inst/CITATION` (missed last time I updated them in other locations).

# bnpsd 1.3.8.9000 (2021-06-21)

- Added function `undiff_af` for creating "undifferentiated" allele frequency distributions based on real data but with a lower variance (more concentrated around 0.5) according to a given FST, useful for simulating data trying to match real data.
- Added `LICENSE.md`.
- Reformatted this `NEWS.md` slightly to improve its automatic parsing.
