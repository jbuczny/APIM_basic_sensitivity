
## title: "APIM meta-analysis (study-level)"
## author: "Anonymized"
## date: "15-10-2026"

# ============================================================
#  APIMeM Power at a FIXED Number of Dyads -- TWO PARALLEL MEDIATORS
#  Actor-Partner Interdependence Mediation Model
#  (distinguishable dyads; Ledermann, Macho & Kenny, 2011)
#  Simulation-based (Monte Carlo) power analysis using lavaan
# ============================================================
#
#  Variables (member 1 = e.g. men, member 2 = e.g. women)
#    X1, X2 = predictors
#    M1, M2 = mediator M (e.g. stress)     measured for both members
#    W1, W2 = mediator W (e.g. support)    measured for both members
#    Y1, Y2 = outcomes
#
#  Paths ("A" = actor, "P" = partner; M or W = which mediator)
#    a paths X -> M   aMA1: X1->M1  aMA2: X2->M2  aMP1: X2->M1  aMP2: X1->M2
#    a paths X -> W   aWA1: X1->W1  aWA2: X2->W2  aWP1: X2->W1  aWP2: X1->W2
#    b paths M -> Y   bMA1: M1->Y1  bMA2: M2->Y2  bMP1: M2->Y1  bMP2: M1->Y2
#    b paths W -> Y   bWA1: W1->Y1  bWA2: W2->Y2  bWP1: W2->Y1  bWP2: W1->Y2
#    c paths X -> Y   cA1:  X1->Y1  cA2:  X2->Y2  cP1:  X2->Y1  cP2:  X1->Y2
#    (index = member whose mediator/outcome is being predicted)
#
#  Simple indirect effects, per mediator (shown for M; W is analogous)
#  (type = a-path type, then b-path type; index = outcome member)
#    M.AA1 = aMA1*bMA1  X1->M1->Y1     M.AA2 = aMA2*bMA2  X2->M2->Y2
#    M.PP1 = aMP2*bMP1  X1->M2->Y1     M.PP2 = aMP1*bMP2  X2->M1->Y2
#    M.PA1 = aMP1*bMA1  X2->M1->Y1     M.PA2 = aMP2*bMA2  X1->M2->Y2
#    M.AP1 = aMA2*bMP1  X2->M2->Y1     M.AP2 = aMA1*bMP2  X1->M1->Y2
#  Combined effects, for each X->Y relation
#  (TA1: X1->Y1, TP1: X2->Y1, TA2: X2->Y2, TP2: X1->Y2)
#    M.TA1 = M.AA1 + M.PP1          total indirect via M
#    W.TA1 = W.AA1 + W.PP1          total indirect via W
#    TA1   = M.TA1 + W.TA1          total indirect via both
#    D.TA1 = M.TA1 - W.TA1          contrast: M vs W
#
#  Two layers of "Monte Carlo" in this script:
#   OUTER: the power simulation. Many data sets are drawn from the
#          population model; power = proportion of data sets in
#          which an effect is detected.
#   INNER: the test of each indirect effect within one data set,
#          a Monte Carlo confidence interval (Preacher & Selig,
#          2012): draw path estimates from their joint sampling
#          distribution, form the products and sums, and check
#          whether the 95% percentile interval excludes 0.
#
#  Outputs, for EACH data set listed in `datasets`
#   (A) Achieved power for your hypothesized set of betas: all
#       paths, simple indirect, total indirect, and contrasts.
#   (B) Power curves: each focal simple indirect effect is swept
#       across a grid of betas (others at hypothesized values).
#   (C) The power-curve plot, exported to a PNG file.
#  With more than one data set, side-by-side comparison tables are
#  printed at the end.
#
#  Each data set is analysed ON ITS OWN, as if it were the only
#  sample. Pooling several data sets into one analysis (data-set
#  intercepts, paths estimated within data sets) needs a multilevel
#  model and is outside the scope of this script.
#
#  Runtime at the defaults (16 focal effects x 11 betas x 1000
#  reps) is several minutes PER DATA SET; see the notes on
#  reducing it.
# ------------------------------------------------------------

# install.packages(c("lavaan", "future.apply"))   # run once if needed
library(lavaan)
library(future.apply)

# ---- 1. Settings --------------------------------------------
#  Data sets ---------------------------------------------------
#  One full analysis per entry. The names label the output and the
#  PNG files, so keep them short (letters, digits, _ or -). Add or
#  remove entries freely, e.g. c(D1 = 130) or c(A = 90, B = 130, C = 174).
datasets <- c(D1 = 130, D2 = 174)    # number of dyads in each data set

#  Plot export -------------------------------------------------
export_png <- TRUE           # one PNG per data set
out_dir    <- "power_plots"  # created if needed; relative to getwd()
png_width  <- 13             # inches
png_height <- 5.5            # inches
png_res    <- 300            # dots per inch
show_plots <- TRUE           # also draw each plot in the Plots pane

#  Hypothesized standardized betas
paths_hyp <- c(
  aMA1 = 0.30, aMA2 = 0.30, aMP1 = 0.15, aMP2 = 0.15,   # X -> M
  aWA1 = 0.20, aWA2 = 0.20, aWP1 = 0.10, aWP2 = 0.10,   # X -> W
  bMA1 = 0.30, bMA2 = 0.30, bMP1 = 0.15, bMP2 = 0.15,   # M -> Y
  bWA1 = 0.20, bWA2 = 0.20, bWP1 = 0.10, bWP2 = 0.10,   # W -> Y
  cA1  = 0.10, cA2  = 0.10, cP1  = 0.05, cP2  = 0.05)   # X -> Y

rxx    <- 0.00   # correlation between partners' predictors (X1, X2)
rmm    <- 0.30   # residual cor. M1-M2 (partners, mediator M)
rww    <- 0.30   # residual cor. W1-W2 (partners, mediator W)
rmw_in <- 0.30   # residual cor. M1-W1 and M2-W2 (same person)
rmw_x  <- 0.10   # residual cor. M1-W2 and M2-W1 (cross person)
ryy    <- 0.30   # residual cor. Y1-Y2

#  Sweep for output (B)
beta_grid <- seq(0, 0.50, by = 0.05)
focal_ind <- c(paste0("M.", c("AA1", "AA2", "PP1", "PP2", "AP1", "AP2", "PA1", "PA2")),
               paste0("W.", c("AA1", "AA2", "PP1", "PP2", "AP1", "AP2", "PA1", "PA2")))
sweep_component <- "both"
#   "both": a = b = beta   (X axis = beta of each component path;
#                           implied indirect effect = beta^2)
#   "a"   : only the a path is swept, b stays at hypothesized value
#   "b"   : only the b path is swept, a stays at hypothesized value

#  Test for simple indirect effects
ind_test <- "mc"   # "mc"   : Monte Carlo CI (Preacher & Selig, 2012)
                   # "joint": joint significance (a and b both p < alpha)
                   # Totals and contrasts always use the MC CI.
mc_draws <- 2000   # draws for the inner MC CI (fresh in every rep)

#  Parallel processing -----------------------------------------
parallel  <- TRUE     # FALSE -> plan(sequential)
n_workers <- NULL     # NULL = parallelly::availableCores() - 1
set_plan  <- TRUE     # FALSE: leave a plan() you set yourself alone
                      #        (e.g. in an Rmd setup chunk)

parallel_over <- "reps"
#  Both settings give results IDENTICAL to the sequential script,
#  for any number of workers. That is a property of the seeding
#  scheme already used here: replication i always starts from
#  set.seed(seed + i), so no replication and no scenario depends on
#  another's RNG state, and the work can be split any way at all.
#   "reps" (default): replications are spread over workers. Output
#       (A) is parallelised too, and load balancing is finer.
#   "scenarios": one scenario per worker. Slightly lower dispatch
#       overhead, but (A) is a single scenario and runs on one core.
#       Prefer this if you have many scenarios and few reps.

n_reps <- 1000     # outer replications per scenario
alpha  <- 0.05
target <- 0.80
seed   <- 2024

#  Estimation engine: "fast" = closed-form ML, identical to lavaan
#  for this saturated model (checked below); "lavaan" = fit the SEM
#  every replication (much slower). If you change apim_model, use
#  "lavaan".
engine <- "fast"

# ---- 2. Definitions -----------------------------------------
stopifnot(parallel_over %in% c("reps", "scenarios"))

stopifnot(is.numeric(datasets), length(datasets) >= 1,
          all(datasets == round(datasets)), all(datasets >= 10))
if (is.null(names(datasets))) names(datasets) <- rep("", length(datasets))
unnamed <- !nzchar(names(datasets))
names(datasets)[unnamed] <- paste0("N", datasets[unnamed])
if (anyDuplicated(names(datasets)))
  stop("Every entry of `datasets` needs a unique name.")

#  Parallel back end. availableCores() is used rather than
#  detectCores(): it respects cgroup/container limits, HPC
#  scheduler allocations and R's own options, so it will not
#  over-subscribe the machine. multicore (forking) is cheaper to
#  start but is unavailable on Windows and inside RStudio;
#  supportsMulticore() falls back to multisession there.
if (set_plan) {
  n_workers <- if (is.null(n_workers))
    max(1L, parallelly::availableCores() - 1L) else as.integer(n_workers)
  if (!parallel || n_workers < 2L)          plan(sequential)
  else if (parallelly::supportsMulticore()) plan(multicore,    workers = n_workers)
  else                                      plan(multisession, workers = n_workers)
}
plan_name <- function() {
  cl <- class(plan())
  hit <- intersect(c("sequential", "multicore", "multisession", "cluster"), cl)
  if (length(hit)) hit[1] else cl[1]
}
cat(sprintf("Workers: %d (%s), parallel over %s\n",
            nbrOfWorkers(), plan_name(), parallel_over))

path_names <- names(paths_hyp)
ab_names   <- grep("^[ab]", path_names, value = TRUE)

base_ind <- data.frame(
  tag   = c("AA1", "PP1", "PA1", "AP1", "AA2", "PP2", "PA2", "AP2"),
  a     = c("A1",  "P2",  "P1",  "A2",  "A2",  "P1",  "P2",  "A1"),
  b     = c("A1",  "P1",  "A1",  "P1",  "A2",  "P2",  "A2",  "P2"),
  route = c("X1->#1->Y1", "X1->#2->Y1", "X2->#1->Y1", "X2->#2->Y1",
            "X2->#2->Y2", "X2->#1->Y2", "X1->#2->Y2", "X1->#1->Y2"),
  stringsAsFactors = FALSE)
ind_def <- do.call(rbind, lapply(c("M", "W"), function(m) data.frame(
  name  = paste0(m, ".", base_ind$tag),
  a     = paste0("a", m, base_ind$a),
  b     = paste0("b", m, base_ind$b),
  route = gsub("#", m, base_ind$route, fixed = TRUE),
  med   = m, stringsAsFactors = FALSE)))

#  Combined effects as weights over the 16 simple indirect effects
rel_def <- list(TA1 = c("AA1", "PP1"), TP1 = c("PA1", "AP1"),
                TA2 = c("AA2", "PP2"), TP2 = c("PA2", "AP2"))
rel_route <- c(TA1 = "X1->Y1", TP1 = "X2->Y1", TA2 = "X2->Y2", TP2 = "X1->Y2")
wvec <- function(m_w, w_w, tags) {
  setNames(ifelse(ind_def$med == "M", m_w, w_w) *
           (sub("^.\\.", "", ind_def$name) %in% tags), ind_def$name)
}
comb_def <- c(
  setNames(lapply(rel_def, function(t) wvec(1, 0, t)),  paste0("M.", names(rel_def))),
  setNames(lapply(rel_def, function(t) wvec(0, 1, t)),  paste0("W.", names(rel_def))),
  setNames(lapply(rel_def, function(t) wvec(1, 1, t)),  names(rel_def)),
  setNames(lapply(rel_def, function(t) wvec(1, -1, t)), paste0("D.", names(rel_def))))
Wmat <- do.call(cbind, comb_def)
comb_type <- rep(c("total indirect via M", "total indirect via W",
                   "total indirect", "contrast M - W"), each = 4)

eff_names <- c(path_names, ind_def$name, names(comb_def))

# ---- 3. lavaan model: distinguishable-dyad APIMeM -----------
#  Saturated (df = 0): all direct paths, all mediator residual
#  covariances, correlated outcome residuals.
apim_model <- '
  M1 ~ aMA1*X1 + aMP1*X2
  M2 ~ aMA2*X2 + aMP2*X1
  W1 ~ aWA1*X1 + aWP1*X2
  W2 ~ aWA2*X2 + aWP2*X1
  Y1 ~ cA1*X1 + cP1*X2 + bMA1*M1 + bMP1*M2 + bWA1*W1 + bWP1*W2
  Y2 ~ cA2*X2 + cP2*X1 + bMA2*M2 + bMP2*M1 + bWA2*W2 + bWP2*W1
  M1 ~~ M2 + W1 + W2
  M2 ~~ W1 + W2
  W1 ~~ W2
  Y1 ~~ Y2
  X1 ~~ X2
'

# ---- 4. Population model and data generation ----------------
#  All variables have variance 1, so the betas are standardized.
#  Residual variances are solved from the implied covariances.
#  Mediator order: M1, M2, W1, W2.
Rmed <- matrix(c(1,      rmm,    rmw_in, rmw_x,
                 rmm,    1,      rmw_x,  rmw_in,
                 rmw_in, rmw_x,  1,      rww,
                 rmw_x,  rmw_in, rww,    1), 4, byrow = TRUE)
if (any(eigen(Rmed, only.values = TRUE)$values <= 0))
  stop("Mediator residual correlations are not positive definite.")

setup_pop <- function(pr) {
  A  <- rbind(c(pr[["aMA1"]], pr[["aMP1"]]),          # M1 ~ X1 + X2
              c(pr[["aMP2"]], pr[["aMA2"]]),          # M2
              c(pr[["aWA1"]], pr[["aWP1"]]),          # W1
              c(pr[["aWP2"]], pr[["aWA2"]]))          # W2
  Cm <- rbind(c(pr[["cA1"]], pr[["cP1"]], pr[["bMA1"]], pr[["bMP1"]],
                pr[["bWA1"]], pr[["bWP1"]]),          # Y1 ~ X1 X2 M1 M2 W1 W2
              c(pr[["cP2"]], pr[["cA2"]], pr[["bMP2"]], pr[["bMA2"]],
                pr[["bWP2"]], pr[["bWA2"]]))          # Y2
  Sx  <- matrix(c(1, rxx, rxx, 1), 2)
  vM  <- 1 - diag(A %*% Sx %*% t(A))
  if (any(vM <= 0)) return(list(ok = FALSE))
  PsiM <- outer(sqrt(vM), sqrt(vM)) * Rmed
  SXM  <- Sx %*% t(A)
  SP   <- rbind(cbind(Sx, SXM), cbind(t(SXM), A %*% Sx %*% t(A) + PsiM))
  vY   <- 1 - diag(Cm %*% SP %*% t(Cm))
  if (any(vY <= 0)) return(list(ok = FALSE))
  PsiY <- outer(sqrt(vY), sqrt(vY)) * matrix(c(1, ryy, ryy, 1), 2)
  list(ok = TRUE, A = A, Cm = Cm,
       LX = chol(Sx), LM = chol(PsiM), LY = chol(PsiY))
}

sim_data <- function(n, S) {
  X <- matrix(rnorm(2 * n), n) %*% S$LX
  M <- X %*% t(S$A) + matrix(rnorm(4 * n), n) %*% S$LM   # M1 M2 W1 W2
  Y <- cbind(X, M) %*% t(S$Cm) + matrix(rnorm(2 * n), n) %*% S$LY
  list(X = X, M = M, Y = Y)
}

# ---- 5. Estimation engines: return estimates + ACOV ---------
fit_lavaan <- function(d) {
  df <- data.frame(X1 = d$X[, 1], X2 = d$X[, 2],
                   M1 = d$M[, 1], M2 = d$M[, 2], W1 = d$M[, 3], W2 = d$M[, 4],
                   Y1 = d$Y[, 1], Y2 = d$Y[, 2])
  fit <- try(sem(apim_model, data = df), silent = TRUE)
  if (inherits(fit, "try-error") || !lavInspect(fit, "converged"))
    return(NULL)
  list(est = coef(fit)[path_names],
       V   = unclass(vcov(fit))[path_names, path_names])
}

#  Closed form: each block (mediator equations, outcome equations)
#  is a system with identical regressors, so ML = OLS per equation,
#  residual covariance = SSE/N, and ACOV = Sigma_res (x) (P'P)^-1.
#  The two blocks are independent (recursive, saturated model).
ols_block <- function(P, Y) {
  n  <- nrow(P)
  Pc <- P - rep(colMeans(P), each = n)
  Yc <- Y - rep(colMeans(Y), each = n)
  Pi <- solve(crossprod(Pc))
  B  <- Pi %*% crossprod(Pc, Yc)
  list(B = B, V = kronecker(crossprod(Yc - Pc %*% B) / n, Pi))
}
nm_fast <- c("aMA1", "aMP1", "aMP2", "aMA2",              # vec(B), mediator block
             "aWA1", "aWP1", "aWP2", "aWA2",
             "cA1", "cP1", "bMA1", "bMP1", "bWA1", "bWP1",  # vec(B), outcome block
             "cP2", "cA2", "bMP2", "bMA2", "bWP2", "bWA2")
fit_fast <- function(d) {
  m <- ols_block(d$X, d$M)
  y <- ols_block(cbind(d$X, d$M), d$Y)
  est <- setNames(c(m$B, y$B), nm_fast)
  V <- matrix(0, 20, 20, dimnames = list(nm_fast, nm_fast))
  V[1:8, 1:8] <- m$V; V[9:20, 9:20] <- y$V
  list(est = est[path_names], V = V[path_names, path_names])
}

fit_fun <- switch(engine, lavaan = fit_lavaan, fast = fit_fast,
                  stop("engine must be \"lavaan\" or \"fast\""))

S_hyp <- setup_pop(paths_hyp)
if (!S_hyp$ok) stop("Hypothesized betas imply residual variance <= 0.")

if (engine == "fast") {                  # safety check vs lavaan
  set.seed(seed)
  d_chk <- sim_data(datasets[[1]], S_hyp)
  f1 <- fit_fast(d_chk); f2 <- fit_lavaan(d_chk)
  dev <- max(abs(f1$est - f2$est), max(abs(f1$V - f2$V)) / max(abs(f2$V)))
  if (!is.finite(dev) || dev > 1e-4)
    stop("Fast engine does not match lavaan (max rel. diff = ",
         signif(dev, 3), "). Set engine <- \"lavaan\".")
  cat(sprintf("Fast engine verified against lavaan (N = %d, max rel. diff = %.1e)\n",
              datasets[[1]], dev))
}

# ---- 6. One replication and power for one scenario ----------
#  Inner MC draws are fresh in every replication, so their
#  Monte Carlo error averages out over the outer replications.
n_ab <- length(ab_names)
ia <- match(ind_def$a, ab_names)         # column positions of a and b
ib <- match(ind_def$b, ab_names)

mc_sig <- function(Dm) {                 # percentile CI excludes 0?
  pp <- colMeans(Dm > 0)
  pp > 1 - alpha / 2 | pp < alpha / 2
}

one_rep <- function(S, n) {
  f <- fit_fun(sim_data(n, S))
  if (is.null(f)) return(rep(NA_real_, length(eff_names)))
  p <- 2 * pnorm(-abs(f$est / sqrt(diag(f$V))))
  #  draws ~ MVN(est, V): [1, Z] %*% rbind(est, chol(V))
  Z1  <- cbind(1, matrix(rnorm(mc_draws * n_ab), mc_draws))
  D   <- Z1 %*% rbind(unname(f$est[ab_names]),
                      chol(f$V[ab_names, ab_names]))
  IND <- D[, ia] * D[, ib]
  sig_ind <- if (ind_test == "joint") {
    p[ind_def$a] < alpha & p[ind_def$b] < alpha
  } else mc_sig(IND)
  as.numeric(c(p < alpha, sig_ind, mc_sig(IND %*% Wmat)))
}

run_power <- function(pr, n) {
  S <- setup_pop(pr)
  if (!S$ok) stop("Implied residual variance <= 0.")
  #  Replication i always starts from seed + i: the same data sets
  #  are used in every scenario (common random numbers), whatever
  #  the engine or number of inner draws. Because the seed is set
  #  INSIDE the loop, replications are independent of each other's
  #  RNG state, so future.seed = NULL is correct here: no seeding by
  #  future is wanted, and splitting the loop over workers cannot
  #  change a single number.
  rep_apply <- if (parallel_over == "reps")
    function(X, FUN, FV) future_vapply(X, FUN, FV, future.seed = NULL,
                                       future.packages = "lavaan") else vapply
  hits <- rep_apply(seq_len(n_reps), function(i) {
    set.seed(seed + i); one_rep(S, n)
  }, numeric(length(eff_names)))
  n_ok  <- rowSums(!is.na(hits))
  power <- rowMeans(hits, na.rm = TRUE)
  names(power) <- names(n_ok) <- eff_names
  list(power = power, mcse = sqrt(power * (1 - power) / n_ok))
}

# ---- 7. Scenarios for output (B) ----------------------------
#  The grid depends only on the betas, not on N, so it is built
#  once and reused for every data set.
comp_to_set <- function(k) {
  switch(sweep_component,
         both = c(ind_def$a[k], ind_def$b[k]),
         a    = ind_def$a[k],
         b    = ind_def$b[k])
}
scen <- expand.grid(beta = beta_grid, ind = focal_ind,
                    stringsAsFactors = FALSE)
scen_pars <- lapply(seq_len(nrow(scen)), function(k) {
  pr <- paths_hyp
  pr[comp_to_set(match(scen$ind[k], ind_def$name))] <- scen$beta[k]
  pr
})
valid <- vapply(scen_pars, function(pr) setup_pop(pr)$ok, logical(1))
if (any(!valid)) {
  warning(sum(!valid), " scenario(s) dropped: implied residual variance <= 0")
  scen <- scen[valid, ]; scen_pars <- scen_pars[valid]
}

#  Progress messages from inside a worker do not reach the console,
#  so they are only printed when the scenarios run sequentially.
scen_apply <- if (parallel_over == "scenarios")
  function(X, FUN) future_lapply(X, FUN, future.seed = NULL,
                                 future.packages = "lavaan") else lapply

# ---- 8. Analysis of one data set: (A), (B), MDE -------------
#  Minimal detectable effect (linear interpolation)
mde <- function(beta, power, target) {
  o <- order(beta); beta <- beta[o]; power <- power[o]
  i <- which(power >= target)[1]
  if (is.na(i)) return(NA_real_)
  if (i == 1)   return(beta[1])
  beta[i - 1] + (target - power[i - 1]) *
    (beta[i] - beta[i - 1]) / (power[i] - power[i - 1])
}

ind_val  <- setNames(paths_hyp[ind_def$a] * paths_hyp[ind_def$b], ind_def$name)
comb_val <- drop(ind_val %*% Wmat)

analyse_dataset <- function(n, label) {
  cat(sprintf("\n%s\n  DATA SET %s: N = %d dyads\n%s\n",
              strrep("#", 62), label, n, strrep("#", 62)))

  #  (A) achieved power at the hypothesized betas
  cat(sprintf("\nRunning (A): N = %d dyads, %d reps (engine: %s)...\n",
              n, n_reps, engine))
  res_A <- run_power(paths_hyp, n)
  tab_A <- data.frame(
    dataset = label, n = n,
    effect  = eff_names,
    type    = c(rep("path", length(path_names)),
                rep("simple indirect", nrow(ind_def)), comb_type),
    route   = c(rep("", length(path_names)), ind_def$route,
                rep(rel_route, 4)),
    value   = c(paths_hyp, ind_val, comb_val),
    power   = res_A$power,
    mcse    = res_A$mcse, row.names = NULL)
  cat(sprintf("\n=== (A) Achieved power, %s: N = %d dyads ===\n", label, n))
  print(format(tab_A[, -(1:2)], digits = 3), row.names = FALSE)
  #  Note: a contrast whose true value is 0 (equal M and W effects)
  #  shows its Type I error rate, not power.

  #  (B) power across the effect-size grid
  cat(sprintf("\nRunning (B): %d scenarios x %d reps...\n", nrow(scen), n_reps))
  res_B <- do.call(rbind, scen_apply(seq_len(nrow(scen)), function(k) {
    r <- run_power(scen_pars[[k]], n)
    if (interactive() && parallel_over != "scenarios") {
      cat(sprintf("\r  %s: scenario %d / %d", label, k, nrow(scen)))
      flush.console()
    }
    data.frame(dataset = label, n = n, ind = scen$ind[k], beta = scen$beta[k],
               power = r$power[[scen$ind[k]]], mcse = r$mcse[[scen$ind[k]]])
  }))
  cat("\n")

  cat(sprintf("\n=== (B) Power by beta, %s: N = %d (sweep: %s) ===\n",
              label, n, sweep_component))
  for (m in c("M", "W")) {
    sub_B <- res_B[startsWith(res_B$ind, paste0(m, ".")), c("ind", "beta", "power")]
    if (nrow(sub_B) == 0) next
    cat(sprintf("-- via mediator %s --\n", m))
    print(reshape(sub_B, idvar = "beta", timevar = "ind", direction = "wide"),
          row.names = FALSE, digits = 3)
  }

  #  minimal detectable effect
  v <- vapply(focal_ind, function(ind) {
    s <- res_B[res_B$ind == ind, ]; mde(s$beta, s$power, target)
  }, numeric(1))
  mde_tab <- data.frame(dataset = label, n = n, ind = focal_ind,
                        route = ind_def$route[match(focal_ind, ind_def$name)],
                        mde_beta = unname(v),
                        mde_indirect = if (sweep_component == "both") unname(v^2)
                                       else NA_real_, row.names = NULL)
  cat(sprintf("\n=== Smallest beta with power >= %.2f, %s: N = %d ===\n",
              target, label, n))
  for (k in seq_along(focal_ind)) {
    cat(sprintf("  %-6s %-11s : %s\n", focal_ind[k], mde_tab$route[k],
        if (is.na(v[k])) sprintf("> %.2f (max tested)", max(beta_grid)) else
        if (sweep_component == "both")
          sprintf("beta = %.3f  (indirect = %.3f)", v[k], v[k]^2) else
          sprintf("beta = %.3f", v[k])))
  }

  list(tab_A = tab_A, res_B = res_B, mde = mde_tab, power_A = res_A$power)
}

# ---- 9. Plot: one panel per mediator + a legend panel -------
#  colour = indirect-effect type, line type = outcome member
type_col <- c(AA = "#1f5f99", PP = "#7d3c98", AP = "#b03a2e", PA = "#d68910")
col_of <- function(ind) type_col[substr(ind, 3, 4)]
lty_of <- function(ind) ifelse(substr(ind, 5, 5) == "1", 1, 2)

xlab_txt <- switch(sweep_component,
  both = "Beta of each component path (a = b = beta)",
  a    = "Beta of the a path (b at hypothesized value)",
  b    = "Beta of the b path (a at hypothesized value)")

meds <- unique(substr(focal_ind, 1, 1))
mlab <- function(k) paste0(meds, k, collapse = "/")      # e.g. "M1/W1"
legend_items <- list(
  legend = c(sprintf("AA  X1->%s->Y1", mlab(1)), sprintf("PP  X1->%s->Y1", mlab(2)),
             sprintf("AP  X2->%s->Y1", mlab(2)), sprintf("PA  X2->%s->Y1", mlab(1)),
             "solid: outcome Y1", "dashed: outcome Y2",
             "hypothesized betas (A)"),
  col = c(type_col, "grey30", "grey30", "grey30"),
  lty = c(1, 1, 1, 1, 1, 2, NA), lwd = c(2, 2, 2, 2, 1, 1, 2),
  pch = c(NA, NA, NA, NA, NA, NA, 23), pt.bg = "white",
  bty = "n", y.intersp = 1.3)

#  Error bars. A point with power of exactly 0 or 1 has mcse = 0,
#  so its bar has zero length; arrows() warns about those and skips
#  them. They are filtered out here rather than silenced, so that
#  any OTHER warning the plot raises still reaches you.
ci_bars <- function(x, power, mcse, col) {
  lo <- pmax(0, power - 1.96 * mcse)
  hi <- pmin(1, power + 1.96 * mcse)
  i  <- which(is.finite(lo) & is.finite(hi) & hi > lo)
  if (length(i))
    arrows(x[i], lo[i], x[i], hi[i],
           length = 0.02, angle = 90, code = 3, col = col)
}

#  The legend gets a column of its own (layout), instead of being
#  drawn into the outer margin of the last panel. It can no longer
#  run off the edge of the device, and its text is shrunk until it
#  fits that column, so it displays the same in the Plots pane and
#  in the PNG.
plot_power <- function(res_B, power_A, n, label) {
  op <- par(no.readonly = TRUE); on.exit(par(op))
  layout(matrix(seq_len(length(meds) + 1), nrow = 1),
         widths = c(rep(1, length(meds)), 0.6))
  #  Text size follows the device width: full size from 11 inches up
  #  (the PNG default is 13), smaller in a narrow Plots pane so that
  #  axis labels are not clipped. layout() alone would fix it at 0.66.
  par(oma = c(2.5, 0, 2.5, 0), mar = c(5, 4.5, 5.5, 1),
      cex = min(1, max(0.55, par("din")[1] / 11)))

  for (m in meds) {
    inds <- focal_ind[startsWith(focal_ind, paste0(m, "."))]
    plot(NA, xlim = range(res_B$beta), ylim = c(0, 1), las = 1,
         xlab = xlab_txt, ylab = "Power")
    title(paste("via mediator", m), line = if (sweep_component == "both") 3.3 else 1.5)
    if (sweep_component == "both") {
      at <- pretty(res_B$beta)
      axis(3, at = at, labels = format(round(at^2, 3)), cex.axis = 0.8)
      mtext("Implied indirect effect (a x b)", side = 3, line = 2, cex = 0.75)
    }
    abline(h = target, lty = 2, col = "grey50")
    abline(h = alpha,  lty = 3, col = "grey70")
    for (ind in inds) {
      s <- res_B[res_B$ind == ind, ]; s <- s[order(s$beta), ]
      ci_bars(s$beta, s$power, s$mcse, col_of(ind))
      lines(s$beta, s$power, col = col_of(ind), lty = lty_of(ind), lwd = 2)
      points(s$beta, s$power, col = col_of(ind), pch = 16, cex = 0.8)
      #  (A) marker where the hypothesized values sit on the swept axis
      #  (in "both" mode only when hypothesized a == b)
      k  <- match(ind, ind_def$name)
      xa <- paths_hyp[[ind_def$a[k]]]; xb <- paths_hyp[[ind_def$b[k]]]
      x  <- switch(sweep_component, both = if (xa == xb) xa else NA, a = xa, b = xb)
      if (!is.na(x))
        points(x, power_A[[ind]], pch = 23, cex = 1.6, lwd = 2,
               col = col_of(ind), bg = "white")
    }
  }

  #  legend panel: largest text size that still fits the column
  par(mar = c(0, 0.5, 0, 0.5))
  plot.new()
  usr <- par("usr"); cex_leg <- 1
  repeat {
    r <- do.call(legend, c(list(x = "center", cex = cex_leg, plot = FALSE),
                           legend_items))$rect
    if ((r$w <= 0.97 * diff(usr[1:2]) && r$h <= 0.97 * diff(usr[3:4])) ||
        cex_leg <= 0.4) break
    cex_leg <- cex_leg - 0.05
  }
  do.call(legend, c(list(x = "center", cex = cex_leg), legend_items))

  #  titles in the outer margins (footnote shrunk to the device width)
  mtext(sprintf("APIMeM power, two parallel mediators, %s: N = %d dyads", label, n),
        outer = TRUE, side = 3, line = 0.8, font = 2, cex = 1.1)
  foot <- sprintf(paste0("test: %s | rxx = %.2f, rmm = %.2f, rww = %.2f, ",
                         "rmw = %.2f/%.2f, ryy = %.2f | %d reps | ",
                         "other paths at hypothesized values"),
                  ifelse(ind_test == "mc", "Monte Carlo CI", "joint significance"),
                  rxx, rmm, rww, rmw_in, rmw_x, ryy, n_reps)
  w1 <- strwidth(foot, units = "inches", cex = 1 / par("cex"))   # at mtext cex = 1
  mtext(foot, outer = TRUE, side = 1, line = 1,
        cex = min(0.75, 0.97 * par("din")[1] / w1))
}

safe_name <- function(s) gsub("[^A-Za-z0-9_-]+", "_", s)

# ---- 10. Run every data set ---------------------------------
if (export_png) dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

results <- setNames(vector("list", length(datasets)), names(datasets))
for (j in seq_along(datasets)) {
  lab <- names(datasets)[j]; n <- datasets[[j]]
  res <- analyse_dataset(n, lab)
  results[[lab]] <- res

  if (show_plots)
    tryCatch(plot_power(res$res_B, res$power_A, n, lab),
             error = function(e) message(
               "On-screen plot for ", lab, " skipped (", conditionMessage(e),
               "). Enlarge the Plots pane and rerun plot_power(); ",
               "the PNG export is not affected."))

  if (export_png) {
    f <- file.path(out_dir, sprintf("apimem2_power_%s_N%d.png", safe_name(lab), n))
    png(f, width = png_width, height = png_height, units = "in", res = png_res)
    tryCatch(plot_power(res$res_B, res$power_A, n, lab), finally = dev.off())
    cat(sprintf("\n  plot saved: %s\n", normalizePath(f)))
  }
}

#  All data sets in one table each (column `dataset` tells them apart)
all_A   <- do.call(rbind, c(lapply(results, `[[`, "tab_A"), make.row.names = FALSE))
all_B   <- do.call(rbind, c(lapply(results, `[[`, "res_B"), make.row.names = FALSE))
all_mde <- do.call(rbind, c(lapply(results, `[[`, "mde"),   make.row.names = FALSE))

#  Side-by-side comparison across data sets
col_lab <- sprintf("%s_N%d", names(datasets), datasets)
power_by_dataset <- data.frame(
  results[[1]]$tab_A[, c("effect", "type", "route", "value")],
  setNames(as.data.frame(do.call(cbind, lapply(results, `[[`, "power_A"))), col_lab),
  check.names = FALSE, row.names = NULL)
mde_by_dataset <- data.frame(
  ind   = focal_ind,
  route = ind_def$route[match(focal_ind, ind_def$name)],
  setNames(as.data.frame(do.call(cbind, lapply(results, function(r) r$mde$mde_beta))),
           col_lab),
  check.names = FALSE, row.names = NULL)

if (length(datasets) > 1) {
  cat(sprintf("\n%s\n  COMPARISON ACROSS DATA SETS\n%s\n",
              strrep("#", 62), strrep("#", 62)))
  cat("\n=== (A) Achieved power at the hypothesized betas ===\n")
  print(format(power_by_dataset, digits = 3), row.names = FALSE)
  cat(sprintf("\n=== Smallest beta with power >= %.2f (sweep: %s) ===\n",
              target, sweep_component))
  print(format(mde_by_dataset, digits = 3), row.names = FALSE)
  cat("  (NA = target not reached within beta_grid)\n")
}

write.csv(all_A,            "apimem2_power_A.csv",            row.names = FALSE)
write.csv(all_B,            "apimem2_power_B.csv",            row.names = FALSE)
write.csv(power_by_dataset, "apimem2_power_by_dataset.csv",   row.names = FALSE)

# ============================================================
#  NOTES
#  - Data sets: each entry of `datasets` is analysed on its own,
#    exactly as the single-N version would analyse it. The seeds are
#    the same for every data set (replication i starts from
#    seed + i), so rerunning gives identical numbers, but the data
#    sets are not nested: the N = 130 samples are not the first 130
#    dyads of the N = 174 samples.
#  - Plots: the PNGs are written to out_dir, one per data set, named
#    apimem2_power_<name>_N<n>.png; the full path is printed. To redraw
#    one in the Plots pane later:
#      plot_power(results$D1$res_B, results$D1$power_A, 130, "D1")
#    If the pane is very small R may refuse ("figure margins too
#    large"); the script catches that and the PNG is still written.
#  - Power for a simple indirect effect is roughly the product of
#    the power for its a and b components. Partner-partner (PP)
#    effects multiply two typically small paths.
#  - Correlated mediators (rmw_in, rmw_x) make the b paths of M
#    and W collinear in the Y equations, which inflates their SEs
#    and lowers power for the specific indirect effects. Choose
#    these values carefully; they matter more here than in the
#    one-mediator model.
#  - Contrasts (D.*) test whether the indirect effect via M differs
#    from that via W. Their power depends on the difference in
#    hypothesized effects and is usually low.
#  - At beta = 0 in "both" mode, a = b = 0 and the indirect test
#    is far below alpha (a known property of product tests).
#  - To shorten the runtime: sweep only the effects you need
#    (e.g. focal_ind <- c("M.AA1", "M.PP1", "W.AA1", "W.PP1")),
#    use a coarser beta_grid, or use n_reps = 500 while exploring.
#    Runtime grows linearly with the number of data sets.
#  - For the actual analysis, a percentile bootstrap in lavaan
#    (se = "bootstrap") gives power close to the MC CI used here.
#  - Parallelism pays off in proportion to the cost of one
#    replication. With engine = "fast" a replication is a closed-form
#    fit plus the inner Monte Carlo draws; with engine = "lavaan" it
#    is a full SEM fit and the speed-up is close to linear.
#  - plan(sequential) after the run releases multisession workers.
#  - To confirm nothing changed, run once with parallel <- FALSE and
#    once with parallel <- TRUE and compare all_B: the two are
#    identical, because each replication re-seeds itself.
#  - Monte Carlo SE of a power estimate = sqrt(p(1-p)/reps).
# ============================================================
