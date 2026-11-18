# SelectSim (development version)

* New three-state mode (`selectX(..., clonality = TRUE)`): alteration matrices
  encode 0 = absent, 1 = subclonal, 2 = clonal, and SelectSim additionally tests
  co-clonality (`CC`), co-subclonality (`SS`), directional partial mutual
  exclusivity (`CS`, `SC`) and discordant clonality (`DISC`). Results are
  returned in `result$clonality`; `result$result` (co-occurrence / mutual
  exclusivity) is unchanged.
* New clonality null model (`clonality_null_model()`) preserving per-gene clonal
  and subclonal counts and correcting for the per-sample clonal propensity
  (`sample_clonal_propensity()`); `clonality.null = "conditional"` (default,
  tests clonality given co-occurrence) or `"joint"`.
* New `clonality.table()`, `clonality_patterns()` and `maf2clonal_tmb()`.
* `maf2gam()` gains `ccf.col` / `clonal.threshold` and `clonality.col` /
  `clonal.labels` to build three-state matrices.
* `new.AL.general()` binarises each alteration-type matrix (no change for 0/1 input).
* Uncertain clonality: states in (1, 2) encode a clonal probability (1 + p) that
  is propagated to the clonality tests (expected counts; the null re-assigns the
  sorted per-gene probabilities along the residual ranking, identical to the
  hard-call null for p in {0, 1}); NA marks present alterations whose clonality
  is excluded from the clonality tests. `maf2gam()` and `maf2clonal_tmb()` gain
  `clonal.prob.col`. The clonality FDR uses standardised effect sizes (`zES`,
  `zFDR`), and the sample clonal propensity is only estimated from independent
  (non-GAM) clonal/subclonal counts.

# SelectSim 0.1.6

* Removed `reshape2` dependency; replaced with base R equivalents (`tapply`, `as.data.frame.table`).
* Fixed `theme_ridges()` partial argument match (`center` → `center_axis_labels`).
* Fixed parallel log file leak in `null_model_parallel` using `tempfile()` and `on.exit()`.
* Refactored `selectX()` to eliminate duplicate verbose/non-verbose code paths.
* Added `utils::globalVariables()` declarations to suppress R CMD check NOTEs for ggplot2 aesthetics.
* Dropped `SystemRequirements: C++14` from DESCRIPTION.
* Added GitHub community health files: `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, issue templates, PR template, and CI workflow.
* Expanded test coverage with `test-gam_utils.R`.
* Repository migrated: development moved from personal branches to [CSOgroup/SelectSim](https://github.com/CSOgroup/SelectSim); `dev` branch merged into `main` for initial public release.
* Bumped version to 0.1.6 across all source files
* Updated citation (`inst/CITATION`, `README`) to the published Nature Genetics article: Iyer A, Petrovic M, Sesia D, Nanni L, Mina M, Ciriello G (2026). Evolving patterns of co-mutations from tumor initiation to metastatic progression. *Nature Genetics*. DOI: 10.1038/s41588-026-02661-4

# SelectSim 0.0.1.3

* Added a `NEWS.md` file to track changes to the package.
* Renamed the package to SelectSim by creating a new Git repository.
* Created the CSOgroup repository and moved the code there.
* Created the GitHub website.

# SelectSim 0.0.1.4

* Added Mijan to the author list.
* Removed the C/C++ code dependency to avoid installation difficulties across systems.
  * Moved to `Matrix` package functions and removed the `Rcpp` functions and compiled code.
* Updated the website and vignette accordingly. 

# SelectSim 0.0.1.5

* Fixed bug of outlier functions and added C/C++ code back.
* More description in functions