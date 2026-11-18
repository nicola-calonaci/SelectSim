make_clonal_data <- function(seed = 7, ng = 30, ns = 300) {
  set.seed(seed)
  st <- matrix(0, ng, ns, dimnames = list(paste0("g", 1:ng), paste0("s", 1:ns)))
  for (g in 1:ng) {
    p <- which(runif(ns) < 0.15)
    st[g, p] <- ifelse(runif(length(p)) < 0.5, 2, 1)
  }
  sh <- sample(ns, 60); st["g1", sh] <- 2; st["g2", sh] <- 2 # co-clonal
  sh <- sample(ns, 60); st["g3", sh] <- 1; st["g4", sh] <- 1 # co-subclonal
  sh <- sample(ns, 60); st["g5", sh] <- 2; st["g6", sh] <- 1 # g5 clonal / g6 subclonal
  o <- setdiff(1:ns, which(st["g1", ] > 0)); st["g1", sample(o, 40)] <- 1
  o <- setdiff(1:ns, which(st["g2", ] > 0)); st["g2", sample(o, 40)] <- 1
  o <- setdiff(1:ns, which(st["g3", ] > 0)); st["g3", sample(o, 40)] <- 2
  o <- setdiff(1:ns, which(st["g4", ] > 0)); st["g4", sample(o, 40)] <- 2
  tmb <- data.frame(sample = colnames(st), mutation = colSums(st > 0) * 10 + rpois(ns, 50))
  tmb$clonal <- round(tmb$mutation * 0.6)
  tmb$subclonal <- tmb$mutation - tmb$clonal
  rownames(tmb) <- tmb$sample
  list(M = list(mut = st), tmb = list(mut = tmb))
}

test_that("new.AL.general validates three-state input", {
  M <- make_clonal_data()
  bad <- M
  bad$M$mut[1, 1] <- 3
  expect_error(new.AL.general(bad, min.freq = 5, clonality = TRUE), "outside")
  bad$M$mut[1, 1] <- 0.5
  expect_error(new.AL.general(bad, min.freq = 5, clonality = TRUE), "outside")
  al <- new.AL.general(M, min.freq = 5, clonality = TRUE)
  expect_true(all(al$am$full %in% c(0, 1)))
  expect_equal(dim(al$clonality$state), dim(al$am$full))
  expect_equal(attr(al$clonality$sample.propensity, "source"), "tmb")
  M$tmb$mut$clonal = NULL
  al0 <- new.AL.general(M, min.freq = 5, clonality = TRUE)
  expect_equal(attr(al0$clonality$sample.propensity, "source"), "none")
  expect_length(unique(al0$clonality$sample.propensity), 1)
})

test_that("alteration types are merged by highest state", {
  m1 <- matrix(c(1, 0, 2, 0), 2, dimnames = list(c("a", "b"), c("s1", "s2")))
  m2 <- matrix(c(2, 1, 1, 0), 2, dimnames = list(c("a", "b"), c("s1", "s2")))
  tmb <- data.frame(sample = c("s1", "s2"), mutation = c(5, 5), row.names = c("s1", "s2"))
  al <- new.AL.general(list(M = list(x = m1, y = m2), tmb = list(x = tmb, y = tmb)),
    min.freq = 0, clonality = TRUE
  )
  expect_equal(unname(al$clonality$state["a", ]), c(2, 2))
  expect_equal(al$am$x["a", "s2"], 1)
})

test_that("conditional clonality null preserves per-gene clonal counts", {
  M <- make_clonal_data()
  al <- new.AL.general(M, min.freq = 5, clonality = TRUE)
  cn <- clonality_null_model(al, n.permut = 20)
  kc <- rowSums(al$clonality$state == 2)
  pres <- al$clonality$state > 0
  for (C in cn) {
    expect_equal(rowSums(C), kc)
    expect_true(all(!C | pres))
  }
})

test_that("selectX detects planted clonality patterns and keeps binary results", {
  M <- make_clonal_data()
  r <- selectX(M,
    sample.class = NULL, alteration.class = NULL, min.freq = 5,
    n.permut = 100, verbose = FALSE, clonality = TRUE
  )
  ct <- r$clonality
  expect_setequal(unique(ct$pattern), names(clonality_patterns()))
  expect_true(all(ct$p_value[!is.na(ct$p_value)] <= 1))
  expect_lte(ct["g1 - g2 | CC", "p_value"], 2 / (r$obj$nSim + 1) + 1e-9) # smallest attainable
  sig <- ct[ct$FDR & ct$direction == "enriched", ]
  expect_true("g1 - g2 | CC" %in% rownames(sig))
  expect_true("g3 - g4 | SS" %in% rownames(sig))
  expect_true("g5 - g6 | CS" %in% rownames(sig))
  expect_false("g5 - g6 | SC" %in% rownames(sig))
  # the four disjoint patterns sum to the co-occurrence count
  x <- ct[ct$pattern %in% c("CC", "SS", "CS", "SC"), ]
  s <- tapply(x$count, x$name, sum)
  ov <- tapply(x$overlap, x$name, `[`, 1)
  expect_equal(as.vector(s), as.vector(ov[names(s)]))

  r0 <- selectX(M,
    sample.class = NULL, alteration.class = NULL, min.freq = 5,
    n.permut = 100, verbose = FALSE
  )
  expect_identical(r0$result, r$result)
  expect_null(r0$clonality)
})

test_that("joint clonality null runs", {
  M <- make_clonal_data()
  r <- selectX(M,
    sample.class = NULL, alteration.class = NULL, min.freq = 5,
    n.permut = 50, verbose = FALSE, clonality = TRUE, clonality.null = "joint"
  )
  expect_equal(attr(r$obj$clonal.null, "mode"), "joint")
  expect_true("g1 - g2 | CC" %in% rownames(r$clonality[r$clonality$FDR, ]))
})

test_that("maf2gam and maf2clonal_tmb build three-state inputs", {
  maf <- data.frame(
    Tumor_Sample_Barcode = c("s1", "s1", "s2", "s3"),
    Hugo_Symbol = c("A", "A", "B", "A"), HGVSp_Short = "x",
    ccf = c(0.3, 1, 0.5, NA)
  )
  g <- maf2gam(maf, ccf.col = "ccf")
  expect_equal(g["s1", "A"], 2)
  expect_equal(g["s2", "B"], 1)
  expect_equal(g["s2", "A"], 0)
  expect_true(is.na(g["s3", "A"]))
  lab <- maf2gam(transform(maf, cl = c("subclonal", "subclonal", "clonal", "clonal")),
    clonality.col = "cl"
  )
  expect_equal(lab["s1", "A"], 1)
  t <- maf2clonal_tmb(maf, samples = c("s1", "s2", "s4"), ccf.col = "ccf")
  expect_equal(t$clonal, c(1, 0, 0))
  expect_equal(t$subclonal, c(1, 1, 0))
})

test_that("states are merged as 'clonal if any alteration is clonal'", {
  a <- matrix(c(0, 1.5, 2, NA, 1.2, NA), 1)
  b <- matrix(c(1.4, 1.5, 1, 1.3, 0, 2), 1)
  expect_equal(as.vector(merge_states(a, b)), c(1.4, 1.75, 2, NA, 1.2, 2))
  expect_error(check_state_matrix(matrix(0.5)), "outside")
  expect_silent(check_state_matrix(matrix(c(0, 1, 1.3, 2, NA))))
})

test_that("clonal probabilities and excluded entries are propagated", {
  M <- make_clonal_data()
  st <- M$M$mut
  set.seed(3)
  pres <- which(st > 0)
  # uncertain calls: clonal prob 0.8 / 0.2 instead of 1 / 0 for a third of entries
  soft <- sample(pres, length(pres) %/% 3)
  st[soft] <- ifelse(st[soft] == 2, 1.8, 1.2)
  # excluded clonality for 5% of entries
  excl <- sample(setdiff(pres, soft), length(pres) %/% 20)
  st[excl] <- NA
  M$M$mut <- st
  al <- new.AL.general(M, min.freq = 5, clonality = TRUE)
  expect_equal(sum(is.na(al$clonality$state)), sum(is.na(st[rownames(al$am$full), ])))
  expect_true(all(al$am$full[is.na(al$clonality$state)] == 1)) # excluded = still present

  cn <- clonality_null_model(al, n.permut = 10)
  known <- attr(cn, "known")
  p_obs <- ifelse(known, al$clonality$state - 1, 0)
  for (C in cn) {
    expect_equal(rowSums(C), rowSums(p_obs)) # expected clonal count per gene preserved
    expect_true(all(C[!known] == 0))
  }

  r <- selectX(M,
    sample.class = NULL, alteration.class = NULL, min.freq = 5,
    n.permut = 100, verbose = FALSE, clonality = TRUE
  )
  ct <- r$clonality
  expect_true(all(ct$overlap_known <= ct$overlap))
  x <- ct[ct$pattern %in% c("CC", "SS", "CS", "SC"), ]
  s <- tapply(x$count, x$name, sum)
  ov <- tapply(x$overlap_known, x$name, `[`, 1)
  expect_equal(as.vector(s), as.vector(ov[names(s)]))
  expect_true("g1 - g2 | CC" %in% rownames(ct[ct$FDR & ct$direction == "enriched", ]))

  rj <- selectX(M,
    sample.class = NULL, alteration.class = NULL, min.freq = 5,
    n.permut = 30, verbose = FALSE, clonality = TRUE, clonality.null = "joint"
  )
  kj <- attr(rj$obj$clonal.null, "known")
  expect_equal(rowSums(kj[[1]]), rowSums(known))
})

test_that("maf2gam and maf2clonal_tmb accept clonal probabilities", {
  maf <- data.frame(
    Tumor_Sample_Barcode = c("s1", "s1", "s2", "s3", "s3"),
    Hugo_Symbol = c("A", "A", "B", "A", "A"), HGVSp_Short = "x",
    p = c(0.5, 0.5, 0.2, NA, 1)
  )
  g <- maf2gam(maf, clonal.prob.col = "p")
  expect_equal(g["s1", "A"], 1.75)
  expect_equal(g["s2", "B"], 1.2)
  expect_equal(g["s3", "A"], 2)
  t <- maf2clonal_tmb(maf, samples = c("s1", "s2", "s3"), clonal.prob.col = "p")
  expect_equal(t$clonal, c(1, 0.2, 1))
  expect_equal(t$subclonal, c(1, 0.8, 0))
  expect_equal(t$mutation, c(2, 1, 2))
})
