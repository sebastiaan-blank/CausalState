# Shared policy handling preserves factor labels, declared levels, and ordering.
g_shift_values <- function(A_obs, spec, t) {
  A_shift <- A_obs
  if (is.null(spec) || !length(spec)) return(A_shift)
  a_names <- names(A_obs)
  n <- nrow(A_obs)
  values <- function(x, a) {
    if (!length(x) || !length(x) %in% c(1L, n))
      stop(sprintf("policy_spec_fun must return one value or %d values at t=%d", n, t),
           call. = FALSE)
    prototype <- A_obs[[a]]
    if (is.factor(prototype)) {
      labels <- rep_len(as.character(x), n)
      if (anyNA(labels))
        stop(sprintf("policy_spec_fun returned missing value(s) for '%s' at t=%d", a, t),
             call. = FALSE)
      if (any(!labels %in% levels(prototype)))
        stop(sprintf("policy_spec_fun returned values outside factor levels for '%s' at t=%d", a, t),
             call. = FALSE)
      return(factor(labels, levels = levels(prototype), ordered = is.ordered(prototype)))
    }
    if ((!is.numeric(x) && !is.logical(x)) || is.factor(x))
      stop(sprintf("policy_spec_fun must return numeric values for numeric treatment '%s' at t=%d",
                   a, t), call. = FALSE)
    rep_len(as.numeric(x), n)
  }
  if (is.data.frame(spec)) {
    for (a in intersect(names(spec), a_names)) A_shift[[a]] <- values(spec[[a]], a)
  } else {
    for (a in intersect(names(spec), a_names)) {
      s <- spec[[a]]
      if (is.atomic(s)) {
        A_shift[[a]] <- values(s, a)
      } else if (is.list(s)) {
        m <- if (!is.null(s$idx) && length(s$idx) == n) as.logical(s$idx) else rep(TRUE, n)
        m[is.na(m)] <- FALSE
        if (!any(m)) next
        val <- if (!is.null(s$value)) s$value else if (!is.null(s$mean)) s$mean else s$setto
        A_shift[[a]][m] <- values(val, a)[m]
      }
    }
  }
  for (a in a_names) {
    invalid <- if (is.factor(A_shift[[a]])) anyNA(A_shift[[a]]) else any(!is.finite(A_shift[[a]]))
    if (invalid)
      stop(sprintf("policy_spec_fun returned non-finite or missing value(s) for '%s' at t=%d", a, t),
           call. = FALSE)
  }
  A_shift
}

# Map each possible source dose to its policy destination for each history.
# Policies must depend only on the current treatment and pre-treatment history.
ordinal_policy_map <- function(block, t, a_name, levels, policy_spec_fun) {
  n <- nrow(block)
  destinations <- matrix(NA_integer_, n, length(levels))
  for (j in seq_along(levels)) {
    candidate <- data.table::copy(block)
    value <- factor(rep(levels[j], n), levels = levels, ordered = is.ordered(block[[a_name]]))
    data.table::set(candidate, j = a_name, value = value)
    obs <- as.data.frame(candidate[, a_name, with = FALSE])
    shifted <- g_shift_values(obs, policy_spec_fun(candidate, t, a_name), t)[[a_name]]
    destinations[, j] <- match(as.character(shifted), levels)
    if (anyNA(destinations[, j]))
      stop(sprintf(paste0("Ordinal propensity policy returned values outside factor levels ",
                          "at t=%d when '%s'=%s. Return a deterministic factor level."),
                   t, a_name, levels[j]), call. = FALSE)
  }
  destinations
}

ordinal_propensity_fit <- function(A, X, newX, levels, sl_g, inner_v,
                                   ids, cluster, seed, fold, t, bounds, sl_workers) {
  # A contains declared category indices, never numeric conversions of labels.
  counts <- tabulate(A, nbins = length(levels))
  if (any(counts == 0L))
    stop(sprintf("Ordinal propensity at t=%d, fold=%d: no training observations for dose(s) %s.",
                 t, fold, paste(levels[counts == 0L], collapse = ", ")), call. = FALSE)

  # Continuation probabilities q_j = P(A > a_j | A >= a_j, H).
  # Multiplication gives a coherent probability mass function without requiring
  # common slopes across thresholds or monotonicity repairs to fitted CDFs.
  p <- matrix(0, nrow(newX), length(levels))
  remaining <- rep(1, nrow(newX))
  sl_rows <- vector("list", length(levels) - 1L)
  for (j in seq_len(length(levels) - 1L)) {
    risk <- which(A >= j)
    fit <- sl_block_fit(
      Y = as.integer(A[risk] > j), X = X[risk, , drop = FALSE],
      family = stats::binomial(), sl_lib = sl_g, v_inner = inner_v,
      id_inner = if (is.null(cluster)) ids[risk] else NULL,
      cluster_inner = if (!is.null(cluster)) cluster[risk] else NULL,
      seed_inner = est_seed_for(seed, fold = fold, t = t, component = "g_prop", phase = "fit") + j - 1L,
      seed_cv_rows = est_seed_for(seed, fold = fold, t = t, component = "g_prop", phase = "cv_rows") + j - 1L,
      seed_sl_fit = est_seed_for(seed, fold = fold, t = t, component = "g_prop", phase = "sl_fit") + j - 1L,
      method = SuperLearner::method.NNLS, sl_workers = sl_workers
    )
    q <- as.numeric(sl_predict(fit$fit, newdata = newX))
    if (length(q) != nrow(newX) || any(!is.finite(q)))
      stop(sprintf("Non-finite ordinal propensity predictions at t=%d, fold=%d, dose=%s",
                   t, fold, levels[j]), call. = FALSE)
    q <- pmin(1 - bounds, pmax(bounds, q))
    p[, j] <- remaining * (1 - q)
    remaining <- remaining * q
    if (!is.null(fit$fit$coef)) sl_rows[[j]] <- data.table::data.table(
      t = as.integer(t), outer_fold = as.integer(fold), threshold = levels[j],
      learner = names(fit$fit$coef), weight = as.numeric(fit$fit$coef),
      cvRisk = as.numeric(fit$fit$cvRisk)
    )
  }
  p[, length(levels)] <- remaining
  list(prob = p, sl_tab = data.table::rbindlist(sl_rows, use.names = TRUE, fill = TRUE))
}

ordinal_policy_ratio <- function(prob, destinations, observed, bounds) {
  # g*(a | h) = sum_b I{d(b,h) = a} g(b | h), including every preimage.
  g_star <- rowSums(prob * (destinations == observed))
  g_obs <- prob[cbind(seq_len(nrow(prob)), observed)]
  r <- g_star / pmax(g_obs, bounds)
  if (any(!is.finite(r)))
    stop("Ordinal propensity produced non-finite density ratios; increase bounds.", call. = FALSE)
  r
}
