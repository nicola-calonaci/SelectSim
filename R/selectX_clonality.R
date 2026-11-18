###
# Author  : SelectSim contributors
# Project : SelectSim
# Desc    : Three-state (absent / subclonal / clonal) extension of SelectSim.
#           Tests co-clonality, co-subclonality and partial mutual exclusivity
#           (discordant clonality) between pairs of alterations.
# Version : 0.1.6.9000
#
# State encoding used throughout:
#   0          absent
#   1          subclonal
#   2          clonal
#   (1, 2)     present, with probability of being clonal p = state - 1
#              (propagation of uncertainty in the clonality call)
#   NA         present, clonality excluded from the clonality tests
#              (e.g. ambiguous calls); still counted as present for CO / ME
#
# Null model (see ?clonality_null_model):
#   Presence of alterations is either kept as observed ("conditional", default)
#   or taken from the binary SelectSim null ("joint"). Given presence, the
#   samples where a gene is present are ranked by q_gs - U_gs with
#   q_gs = phi_g * c_s / mean_block(c) (phi_g gene clonal fraction, c_s sample
#   clonal propensity), and the gene's observed clonal probabilities, sorted in
#   decreasing order, are assigned along this ranking. With hard calls (p in
#   {0, 1}) this is exactly "label the k_g^C top-ranked samples clonal",
#   mirroring the residual-ranking scheme of the binary null; in general it
#   preserves the multiset of clonal probabilities (hence the expected clonal
#   and subclonal counts) of every gene.
###

#' Clonality pattern codes
#'
#' Codes of the pairwise clonality patterns tested by SelectSim when
#' \code{clonality = TRUE}:
#' \describe{
#'   \item{CC}{co-clonality: both alterations clonal}
#'   \item{SS}{co-subclonality: both alterations subclonal}
#'   \item{CS}{SFE_1 clonal and SFE_2 subclonal (directional partial mutual exclusivity)}
#'   \item{SC}{SFE_1 subclonal and SFE_2 clonal (directional partial mutual exclusivity)}
#'   \item{DISC}{discordant clonality, CS + SC (partial mutual exclusivity)}
#' }
#' Total mutual exclusivity (one altered, the other not) and co-occurrence are
#' the binary SelectSim patterns reported in \code{result$result}.
#'
#' @return Named character vector with pattern descriptions.
#'
#' @examples
#' clonality_patterns()
#'
#' @export
clonality_patterns <- function() {
  c(
    CC = "co-clonal",
    SS = "co-subclonal",
    CS = "SFE_1 clonal, SFE_2 subclonal",
    SC = "SFE_1 subclonal, SFE_2 clonal",
    DISC = "discordant clonality (partial ME)"
  )
}

#' Validate a three-state alteration matrix
#'
#' @param m Matrix with values 0, values in [1, 2] or NA (see file header).
#' @param name Name used in error messages.
#' @return Invisibly \code{TRUE}; stops with an error otherwise.
#' @noRd
check_state_matrix <- function(m, name = "matrix") {
  v <- as.vector(as.matrix(m))
  v <- v[!is.na(v)]
  if (!all(v == 0 | (v >= 1 & v <= 2))) {
    stop(paste0(
      "Alteration matrix '", name, "' contains values outside {0} U [1, 2]. ",
      "Encode alterations as 0 = absent, 1 = subclonal, 2 = clonal, ",
      "1 + p for a clonal probability p, or NA for present with clonality excluded."
    ))
  }
  invisible(TRUE)
}

#' Merge two three-state matrices (e.g. two alteration types of the same gene)
#'
#' The gene is clonal in a sample if any of its alterations is clonal: with
#' clonal probabilities \eqn{p_a, p_b} the merged probability is
#' \eqn{1 - (1 - p_a)(1 - p_b)} (the maximum, for hard calls). If one of the
#' alterations has excluded clonality (NA), the merged entry is NA unless the
#' other alteration is certainly clonal.
#' @noRd
merge_states <- function(a, b) {
  absent_a <- !is.na(a) & a == 0
  absent_b <- !is.na(b) & b == 0
  out <- a
  out[absent_a] <- b[absent_a]
  both <- !absent_a & !absent_b
  comb <- 2 - (2 - a) * (2 - b) # 1 + (1 - (1 - p_a) * (1 - p_b))
  certain <- (!is.na(a) & a == 2) | (!is.na(b) & b == 2)
  comb[certain] <- 2
  out[both] <- comb[both]
  out
}

#' Estimate the per-sample clonal propensity
#'
#' @description
#' The clonal propensity \eqn{c_s} of a sample is the fraction of its mutations
#' called clonal. It captures sample-level technical and biological effects
#' (purity, sequencing depth, intra-tumour heterogeneity) that change how many
#' subclonal mutations are detectable. It is estimated from the \code{clonal}
#' and \code{subclonal} columns of the TMB data frames, which should count
#' mutations \emph{other than} the alterations in the GAM (e.g. all non-driver
#' mutations, see \code{\link{maf2clonal_tmb}}): estimating it from the GAM
#' itself leaks the pairwise signal into the null and inflates false
#' discoveries. Estimates are shrunk toward the cohort clonal fraction
#' \eqn{\bar p}: \eqn{c_s = (n^C_s + w \bar p) / (n_s + w)}. Without these
#' columns no sample correction is applied (constant propensity).
#'
#' @param state Three-state matrix (features x samples).
#' @param tmb Named list of TMB data frames (as in \code{al$tmb}), optional.
#' @param prior.weight Shrinkage pseudo-count \eqn{w} (in number of mutations).
#' @return Named numeric vector of clonal propensities, with attribute
#'   \code{source} (\code{"tmb"} or \code{"none"}).
#'
#' @examples
#' state <- matrix(c(2, 1, 0, 2, 2, 1), nrow = 2,
#'                 dimnames = list(c("geneA", "geneB"), c("s1", "s2", "s3")))
#' tmb <- list(mut = data.frame(sample = c("s1", "s2", "s3"), mutation = c(10, 20, 5),
#'                              clonal = c(8, 5, 4), subclonal = c(2, 15, 1)))
#' sample_clonal_propensity(state, tmb)
#'
#' @export
sample_clonal_propensity <- function(state, tmb = NULL, prior.weight = 10) {
  samples <- colnames(state)
  use_tmb <- !is.null(tmb) && length(tmb) > 0 &&
    all(sapply(tmb, function(t) all(c("clonal", "subclonal") %in% colnames(t))))
  if (!use_tmb) {
    known <- !is.na(state) & state > 0
    tot <- sum(known)
    prop <- rep(if (tot > 0) sum(state[known] - 1) / tot else 0.5, length(samples))
    names(prop) <- samples
    attr(prop, "source") <- "none"
    return(prop)
  }
  n_c <- rep(0, length(samples))
  n_s <- rep(0, length(samples))
  names(n_c) <- names(n_s) <- samples
  for (t in tmb) {
    rownames(t) <- t$sample
    n_c <- n_c + t[samples, "clonal"]
    n_s <- n_s + t[samples, "subclonal"]
  }
  tot <- n_c + n_s
  p_bar <- if (sum(tot) > 0) sum(n_c) / sum(tot) else 0.5
  prop <- (n_c + prior.weight * p_bar) / (tot + prior.weight)
  names(prop) <- samples
  attr(prop, "source") <- "tmb"
  return(prop)
}

#' Generate the clonality null model
#'
#' @description
#' Re-assigns clonality to alterations while preserving, for every gene, the
#' observed clonal probabilities of its alterations. For each gene \eqn{g} and
#' sample \eqn{s} where the gene is present (with non-excluded clonality), a
#' score \eqn{q_{gs} - U_{gs}} is computed with
#' \eqn{q_{gs} = \min(1, \phi_g c_s / \bar{c}_b)}, \eqn{\phi_g} the gene clonal
#' fraction, \eqn{c_s} the sample clonal propensity and \eqn{\bar{c}_b} its mean
#' within the sample block. The gene's clonal probabilities, sorted in
#' decreasing order, are then assigned to the samples in decreasing order of
#' score. With hard calls this labels the \eqn{k_g^C} highest scoring samples
#' clonal. Entries with excluded clonality (NA) keep their position in the
#' \code{"conditional"} null and are placed at random among the present samples
#' in the \code{"joint"} null.
#'
#' @param al Alteration landscape object created with \code{clonality = TRUE}.
#' @param null Binary null model (list of simulated matrices); required when
#'   \code{mode = "joint"}.
#' @param mode \code{"conditional"} (default) keeps the observed presence matrix
#'   so clonality patterns are tested given co-occurrence; \code{"joint"}
#'   assigns clonality on the simulated presence matrices of the binary null.
#' @param n.permut Number of simulations (ignored for \code{"joint"}, where one
#'   labelling is generated per binary null matrix).
#' @param seed Random seed.
#' @return List of clonal-probability matrices (features x samples; logical
#'   when all calls are hard), with attribute \code{known}: the matrix of
#'   present entries with non-excluded clonality (\code{"conditional"}) or a
#'   list of such matrices, one per simulation (\code{"joint"}).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' st <- matrix(sample(0:2, 20 * 60, replace = TRUE, prob = c(.6, .2, .2)), 20, 60,
#'              dimnames = list(paste0("g", 1:20), paste0("s", 1:60)))
#' tmb <- data.frame(sample = colnames(st), mutation = colSums(st > 0))
#' rownames(tmb) <- tmb$sample
#' al <- new.AL.general(list(M = list(mut = st), tmb = list(mut = tmb)),
#'                      min.freq = 2, clonality = TRUE)
#' cn <- clonality_null_model(al, n.permut = 5)
#' length(cn)
#' }
#'
#' @export
clonality_null_model <- function(al,
                                 null = NULL,
                                 mode = c("conditional", "joint"),
                                 n.permut = 1000,
                                 seed = 42) {
  mode <- match.arg(mode)
  if (is.null(al$clonality)) stop("al was not created with clonality = TRUE.")
  set.seed(seed)
  state <- al$clonality$state
  genes <- rownames(state)
  samples <- colnames(state)
  known_obs <- !is.na(state) & state > 0
  n_excluded <- rowSums(is.na(state))
  P_obs <- ifelse(known_obs, state - 1, 0)
  hard <- all(P_obs[known_obs] %in% c(0, 1))
  n_known <- rowSums(known_obs)
  phi <- ifelse(n_known > 0, rowSums(P_obs) / n_known, 0)
  p_sorted <- lapply(seq_along(genes), function(g) sort(P_obs[g, known_obs[g, ]], decreasing = TRUE))

  # sample propensity normalised within sample blocks
  prop <- al$clonality$sample.propensity[samples]
  blocks <- get.blocks(al)
  rel <- prop
  for (b in names(blocks$sample.blocks)) {
    sb <- names(blocks$sample.blocks[[b]])
    rel[sb] <- prop[sb] / mean(prop[sb])
  }
  Q <- tcrossprod(phi, rel)
  Q[Q > 1] <- 1
  dimnames(Q) <- list(genes, samples)

  label_once <- function(K) {
    # K: logical matrix of present entries with non-excluded clonality
    R <- Q - matrix(runif(length(Q)), nrow(Q), ncol(Q))
    C <- matrix(0, nrow(Q), ncol(Q), dimnames = dimnames(Q))
    for (g in seq_len(nrow(Q))) {
      pg <- p_sorted[[g]]
      if (length(pg) == 0 || pg[1] == 0) next
      cand <- which(K[g, ])
      if (length(cand) == 0) next
      n <- min(length(pg), length(cand))
      top <- cand[order(R[g, cand], decreasing = TRUE)[seq_len(n)]]
      C[g, top] <- pg[seq_len(n)]
    }
    if (hard) C > 0 else C
  }

  if (mode == "conditional") {
    out <- lapply(seq_len(n.permut), function(i) label_once(known_obs))
    attr(out, "known") <- known_obs
  } else {
    if (is.null(null)) stop("mode = 'joint' requires the binary null model.")
    known <- lapply(null, function(m) {
      m <- as.matrix(m)
      if (!is.null(colnames(m))) m <- m[genes, samples, drop = FALSE]
      K <- m > 0
      for (g in which(n_excluded > 0)) {
        pres <- which(K[g, ])
        K[g, pres[sample.int(length(pres), min(n_excluded[g], length(pres)))]] <- FALSE
      }
      K
    })
    out <- lapply(known, label_once)
    attr(out, "known") <- known
  }
  attr(out, "mode") <- mode
  return(out)
}

#' Indices of the unordered feature pairs
#'
#' Returns the pair indices in the same order used by
#' \code{interaction.table} (lower triangle, column-major): \code{a} is the
#' index of SFE_1 and \code{b} the index of SFE_2.
#' @param n Number of features.
#' @return Two-column integer matrix.
#' @noRd
pair_indices <- function(n) {
  idx <- which(lower.tri(matrix(0, n, n)), arr.ind = TRUE)
  # lower.tri gives (row = i, col = j) with i > j in column-major order
  cbind(a = idx[, "col"], b = idx[, "row"])
}

#' Compute pairwise clonality-pattern overlaps
#'
#' @param C Clonal indicator or probability matrix (features x samples).
#' @param S Subclonal indicator or probability matrix (features x samples).
#' @param W Weight matrix (features x samples) or \code{NULL} for counts.
#' @param ab Pair index matrix from \code{pair_indices}.
#' @return Named list of numeric vectors (one value per pair): CC, SS, CS, SC, DISC.
#' @noRd
clonality_pair_overlaps <- function(C, S, W = NULL, ab) {
  C <- C * 1
  S <- S * 1
  WC <- if (is.null(W)) C else W * C
  WS <- if (is.null(W)) S else W * S
  O_cc <- tcrossprod(WC, C)
  O_ss <- tcrossprod(WS, S)
  O_cs <- tcrossprod(WC, S) # [i, j]: i clonal, j subclonal
  a <- ab[, 1]
  b <- ab[, 2]
  cs <- O_cs[cbind(a, b)]
  sc <- O_cs[cbind(b, a)]
  list(CC = O_cc[cbind(a, b)], SS = O_ss[cbind(a, b)], CS = cs, SC = sc, DISC = cs + sc)
}

#' Build the clonality interaction table
#'
#' @description
#' For every pair of alterations and every clonality pattern (see
#' \code{\link{clonality_patterns}}) compares the observed TMB-weighted count of
#' samples showing the pattern with its expectation under the clonality null
#' model. The effect size \code{wES} is the observed minus expected weighted
#' count (as for co-occurrence); it is standardised by the null standard
#' deviation of the pair (\code{zES}) so that pairs with different numbers of
#' co-mutated samples are comparable, and the FDR is estimated by comparing the
#' observed \code{|zES|} with the pooled null \code{|zES|} of all tested pairs
#' (\code{\link{estimateFDR2}}); \code{FDR} flags \code{zFDR <= maxFDR}. In
#' addition, \code{p_value} is the two-sided permutation p-value of each pair;
#' its resolution is limited by \code{n.permut} (the smallest attainable value
#' is \code{2 / (n.permut + 1)}), so a multiple-testing adjustment of it is
#' only informative with many permutations.
#' With clonal probabilities (states in (1, 2)) all counts are expected counts,
#' e.g. CC = \eqn{\sum_s W_s p_{as} p_{bs}}. Samples where either alteration has
#' excluded clonality (NA) do not contribute; \code{overlap_known} is the number
#' of co-mutated samples where both clonalities are available, and equals
#' CC + SS + CS + SC (unweighted). Pairs with \code{overlap_known} below
#' \code{min.overlap}, or whose pattern count cannot vary under the null, are
#' not tested (\code{NA}).
#'
#' @param al Alteration landscape object created with \code{clonality = TRUE}.
#' @param W Weight matrix (features x samples).
#' @param clonal.null Output of \code{\link{clonality_null_model}}.
#' @param maxFDR FDR cutoff.
#' @param min.overlap Minimum number of co-mutated samples with available
#'   clonality for a pair to be tested (pairs below it get \code{NA}).
#' @return Long-format data frame, one row per pair and pattern.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' st <- matrix(sample(0:2, 20 * 60, replace = TRUE, prob = c(.6, .2, .2)), 20, 60,
#'              dimnames = list(paste0("g", 1:20), paste0("s", 1:60)))
#' tmb <- data.frame(sample = colnames(st), mutation = colSums(st > 0))
#' rownames(tmb) <- tmb$sample
#' res <- selectX(list(M = list(mut = st), tmb = list(mut = tmb)),
#'                sample.class = NULL, alteration.class = NULL,
#'                min.freq = 2, n.permut = 20, verbose = FALSE, clonality = TRUE)
#' head(res$clonality)
#' }
#'
#' @export
clonality.table <- function(al,
                            W,
                            clonal.null,
                            maxFDR = 0.25,
                            min.overlap = 3) {
  state <- al$clonality$state
  genes <- rownames(state)
  samples <- colnames(state)
  W <- as.matrix(W)[genes, samples, drop = FALSE]
  n <- length(genes)
  ab <- pair_indices(n)
  nP <- nrow(ab)
  nSim <- length(clonal.null)
  mode <- attr(clonal.null, "mode")
  known <- attr(clonal.null, "known")

  present_obs <- is.na(state) | state > 0
  known_obs <- !is.na(state) & state > 0
  C_obs <- ifelse(known_obs, state - 1, 0)
  S_obs <- known_obs - C_obs
  obs_w <- clonality_pair_overlaps(C_obs, S_obs, W, ab)
  obs_n <- clonality_pair_overlaps(C_obs, S_obs, NULL, ab)
  overlap <- tcrossprod(present_obs * 1)[cbind(ab[, 1], ab[, 2])]
  overlap_known <- tcrossprod(known_obs * 1)[cbind(ab[, 1], ab[, 2])]

  pats <- names(clonality_patterns())
  null_w <- lapply(pats, function(p) matrix(0, nP, nSim))
  names(null_w) <- pats
  null_n_sum <- lapply(pats, function(p) rep(0, nP))
  names(null_n_sum) <- pats
  for (k in seq_len(nSim)) {
    Ck <- clonal.null[[k]] * 1
    Kk <- if (mode == "joint") known[[k]] else known
    Sk <- Kk - Ck
    ow <- clonality_pair_overlaps(Ck, Sk, W, ab)
    on <- clonality_pair_overlaps(Ck, Sk, NULL, ab)
    for (p in pats) {
      null_w[[p]][, k] <- ow[[p]]
      null_n_sum[[p]] <- null_n_sum[[p]] + on[[p]]
    }
  }

  k_c <- rowSums(C_obs)
  k_s <- rowSums(S_obs)
  k_x <- rowSums(is.na(state))
  a <- ab[, 1]
  b <- ab[, 2]
  tabs <- lapply(pats, function(p) {
    mean_null <- rowMeans(null_w[[p]])
    sd_null <- apply(null_w[[p]], 1, stats::sd)
    tested <- overlap_known >= min.overlap & sd_null > 1e-12
    wES <- (obs_w[[p]] - mean_null) * sin(pi / 4)
    zES <- ifelse(tested, (obs_w[[p]] - mean_null) / sd_null, NA_real_)
    res <- data.frame(
      SFE_1 = genes[a],
      SFE_2 = genes[b],
      name = paste(genes[a], genes[b], sep = " - "),
      pattern = p,
      description = unname(clonality_patterns()[p]),
      clonal_1 = unname(k_c[a]),
      subclonal_1 = unname(k_s[a]),
      clonal_2 = unname(k_c[b]),
      subclonal_2 = unname(k_s[b]),
      excluded_1 = unname(k_x[a]),
      excluded_2 = unname(k_x[b]),
      overlap = overlap,
      overlap_known = overlap_known,
      count = obs_n[[p]],
      r_count = null_n_sum[[p]] / nSim,
      w_count = obs_w[[p]],
      w_r_count = mean_null,
      wES = wES,
      zES = zES,
      zFDR = NA_real_,
      p_value = NA_real_,
      stringsAsFactors = FALSE
    )
    if (sum(tested) > 0) {
      z_null <- (null_w[[p]][tested, , drop = FALSE] - mean_null[tested]) / sd_null[tested]
      res$zFDR[tested] <- estimateFDR2(
        obs = abs(zES[tested]),
        exp = abs(as.vector(z_null)),
        nSim = nSim,
        maxFDR = maxFDR
      )
      # per-pair two-sided permutation p-value
      nw <- null_w[[p]][tested, , drop = FALSE]
      ob <- obs_w[[p]][tested]
      eps <- 1e-9
      p_lo <- (1 + rowSums(nw <= ob + eps)) / (nSim + 1)
      p_hi <- (1 + rowSums(nw >= ob - eps)) / (nSim + 1)
      res$p_value[tested] <- pmin(1, 2 * pmin(p_lo, p_hi))
    }
    res
  })
  results <- do.call(rbind, tabs)
  results$direction <- ifelse(is.na(results$zES), NA,
    ifelse(results$zES > 0, "enriched", ifelse(results$zES < 0, "depleted", "none"))
  )
  results$FDR <- !is.na(results$zFDR) & results$zFDR <= maxFDR
  results <- results[order(is.na(results$zES), -abs(results$zES)), ]
  rownames(results) <- paste(results$name, results$pattern, sep = " | ")
  return(results)
}
