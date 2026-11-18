###
# Author  : Arvind Iyer, Miljan Petrovic
# Project : SelectSim
# Desc    : Main file which to run the SelecSim algoritm via calling selectX function to create alteration object with background model and funtion to generate the table.
# Version : 0.1.6
# Notes:
# - Better Error message and running text
# - Edge case: When sample size in less than 2 there is error in computation (need to fix a number to do this analysis)
# - parallel::makeCluster(2, setup_strategy = "sequential") a possible fix to remove the error of connection problem (https://github.com/rstudio/rstudio/issues/6692)
# - Try to fix the parallel processing issuses
###

###
#' SelectX main function from SelectSim to create alteration object with background model
#'
#' @description
#' `selectX()` takes a list object which consist of genome alteration matrix and tumor mutation burden.
#'
#' @import foreach
#' @import doParallel
#' @importFrom  stats as.dist
#' @importFrom  stats density
#' @importFrom  stats ecdf
#' @importFrom  stats runif
#' @importFrom  stats setNames
#' @param M a list data object which consist of gams and tmbs.
#' @param sample.class sample covariates as named list.
#' @param alteration.class alteration covariates as named list.
#' @param n.cores no of cores.
#' @param min.freq number of samples for features to be at least mutated in.
#' @param n.permut number of simulations.
#' @param lambda lambda parameter.
#' @param tau tau (fold change) parameter.
#' @param folder folder path to store the results.
#' @param save.object store the SelectX object.
#' @param verbose print the time and each steps.
#' @param estimate_pairwise Compute pairwise p-value
#' @param maxFDR FDR value
#' @param seed a random seed
#' @param clonality Logical; if \code{TRUE} the matrices in \code{M$M} are
#'   three-state matrices (0 = absent, 1 = subclonal, 2 = clonal) and, in addition
#'   to co-occurrence / mutual exclusivity, SelectSim tests co-clonality,
#'   co-subclonality and partial mutual exclusivity (discordant clonality); see
#'   \code{\link{clonality_patterns}}. The TMB data frames may include
#'   \code{clonal} and \code{subclonal} columns with clonal and subclonal counts
#'   of mutations not included in \code{M} (e.g. non-driver mutations), used to
#'   correct for the sample clonal propensity; without them no correction is applied.
#' @param clonality.null \code{"conditional"} (default): clonal labels are
#'   re-assigned on the observed alterations, so clonality patterns are tested
#'   given co-occurrence. \code{"joint"}: labels are assigned on the simulated
#'   alterations of the binary null model.
#' @param min.overlap Minimum number of samples with both alterations present
#'   for a pair to be tested for clonality patterns.
#' @return result  a SelectSim object with background model and other info along with result table.
#'   With \code{clonality = TRUE} the list also contains \code{clonality}, a
#'   long-format table with one row per pair and clonality pattern (see
#'   \code{\link{clonality.table}}).
#'
#' @examples
#' \donttest{
#' data(luad_run_data, package = "SelectSim")
#' result <- selectX(
#'   M                = luad_run_data$M,
#'   sample.class     = luad_run_data$sample.class,
#'   alteration.class = luad_run_data$alteration.class,
#'   n.cores          = 1,
#'   min.freq         = 10,
#'   n.permut         = 10,
#'   verbose          = FALSE
#' )
#' head(result$result)
#' }
#'
#' @export
selectX <- function(M,
                    sample.class,
                    alteration.class,
                    n.cores = 1,
                    min.freq = 10,
                    n.permut = 1000,
                    lambda = 0.3,
                    tau = 1,
                    save.object = FALSE,
                    folder = "./",
                    verbose = TRUE,
                    estimate_pairwise = FALSE,
                    maxFDR = 0.25,
                    seed = 42,
                    clonality = FALSE,
                    clonality.null = c("conditional", "joint"),
                    min.overlap = 3) {
  clonality.null <- match.arg(clonality.null)
  set.seed(seed)
  use_tictoc <- verbose && requireNamespace("tictoc", quietly = TRUE)
  if (verbose) {
    print(paste("#### Creating SelectSim object ####"))
    if (use_tictoc) tictoc::tic("Total time taken:")
  }

  if (verbose) {
    print(paste("Step1-> Parsing and Filtering GAM..."))
    if (use_tictoc) tictoc::tic("Time:")
  }
  al <- new.AL.general(M, alteration.class, sample.class, min.freq, verbose, clonality = clonality)
  if (verbose) {
    print(paste("-> Alteration Landscape object created"))
    if (use_tictoc) tictoc::toc()
  }

  if (verbose) {
    print(paste("Step2-> Generating Template object..."))
    if (use_tictoc) tictoc::tic("Time:")
  }
  temp.data <- template.obj.gen(al)
  if (verbose) {
    print(paste("-> Template object created"))
    if (use_tictoc) tictoc::toc()
  }

  if (verbose) {
    print(paste("Step3-> Generating sample weight matrix..."))
    if (use_tictoc) tictoc::tic("Time:")
  }
  W <- generateW_block(al, lambda, tau)
  if (verbose) {
    print(paste("-> Weight Matrix created"))
    if (use_tictoc) tictoc::toc()
  }

  if (verbose) {
    print(paste("Step4-> Generating null model..."))
    if (use_tictoc) tictoc::tic("Time:")
  }
  sim <- null_model_parallel(al, temp.data$temp_mat, W$W, n.cores, n.permut)
  obj <- list("al" = al, "W" = W, "T" = temp.data, "null" = sim, "nSim" = n.permut)
  if (verbose) print(paste("-> Removing the outliers matrix from null model..."))
  outliers <- retrieveOutliers(obj = obj, nSim = n.permut)
  if (length(which(outliers)) == 0) {
    if (verbose) print(paste("Removed null-matrix:", length(which(outliers)), sep = " "))
    obj$nSim <- n.permut
  } else {
    if (verbose) print(paste(" Removed null-matrix:", length(which(outliers)), sep = " "))
    obj$null <- obj$null[which(!outliers)]
    if (verbose) print(paste(" Updated the null-model and nSim variables..."))
    obj$nSim <- length(which(!outliers))
  }
  if (verbose) {
    if (use_tictoc) tictoc::toc()
    print(paste("-> Null model generated"))
    print(paste("### SelectSim object created ###"))
  }

  if (verbose) {
    print(paste("#### Computing EDs on the dataset ####"))
    if (use_tictoc) tictoc::tic("Time:")
  }
  al <- obj$al
  als <- al.stats(obj)
  alp <- al.pairwise.alteration.stats(obj, als, do.blocks = FALSE)
  als[["alteration.pairwise"]] <- alp
  als[["alteration.pairwise"]]$sample.blocks <- NULL
  for (i in names(als$sample.blocks)) als$sample.blocks[[i]][["alteration.pairwise"]] <- alp$sample.blocks[[i]]
  obs.co <- as.matrix(als$alteration.pairwise$overlap)
  wobs.co <- as.matrix(als$alteration.pairwise$w_overlap)
  robs.co <- r.am.pairwise.alteration.overlap(null = obj$null, n.permut = obj$nSim, n.cores = 1)
  wrobs.co <- w.r.am.pairwise.alteration.overlap(null = obj$null, W = obj$W$W, n.permut = obj$nSim, n.cores = 1)
  selectX_result <- interaction.table(al, als, obs.co, wobs.co, robs.co, wrobs.co,
    null = obj$null, maxFDR = maxFDR, n.cores = 1,
    estimate_pairwise = estimate_pairwise, n.permut = obj$nSim
  )
  obj$robs.co <- robs.co
  obj$wrobs.co <- wrobs.co
  if (verbose) {
    if (use_tictoc) tictoc::toc()
    print(paste("#### EDs computed ####"))
  }

  clonality_result <- NULL
  if (clonality) {
    if (verbose) {
      print(paste("#### Testing clonality patterns (", clonality.null, " null) ####", sep = ""))
      if (use_tictoc) tictoc::tic("Time:")
    }
    obj$clonal.null <- clonality_null_model(al,
      null = obj$null, mode = clonality.null,
      n.permut = obj$nSim, seed = seed
    )
    clonality_result <- clonality.table(al,
      W = obj$W$W, clonal.null = obj$clonal.null,
      maxFDR = maxFDR, min.overlap = min.overlap
    )
    if (verbose) {
      if (use_tictoc) tictoc::toc()
      print(paste("#### Clonality patterns computed ####"))
    }
  }

  if (use_tictoc) tictoc::toc()

  if (save.object) {
    saveRDS(obj, file = paste(folder, "selectsim_object.rds", sep = ""))
    saveRDS(selectX_result, file = paste(folder, "selectsim_results.rds", sep = ""))
    if (clonality) saveRDS(clonality_result, file = paste(folder, "selectsim_clonality_results.rds", sep = ""))
  }
  out <- list(obj = obj, result = selectX_result)
  if (clonality) out$clonality <- clonality_result
  return(out)
}

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

## usethis namespace: start
#' @useDynLib SelectSim, .registration = TRUE
## usethis namespace: end
NULL
