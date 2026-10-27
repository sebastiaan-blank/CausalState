# Factor columns declare categorical treatment; their levels define a fixed schema.
factor_treatment_levels <- function(data, a_names) {
  out <- list()
  for (a in a_names) {
    x <- data[[a]]
    if (is.factor(x)) {
      lev <- levels(x)
      if (length(lev) < 2L || anyNA(lev))
        stop("Treatment factor '", a, "' must declare at least two non-missing levels.",
             call. = FALSE)
      out[[a]] <- lev
    } else if (!is.numeric(x) && !is.logical(x)) {
      stop("Treatment '", a, "' must be numeric or a factor. Declare categorical treatment as a factor.",
           call. = FALSE)
    }
  }
  out
}

factor_treatment_ordered <- function(data, treatment_levels) {
  vapply(names(treatment_levels), function(a) is.ordered(data[[a]]), logical(1L))
}

resolve_treatment_levels <- function(data, a_names, weight_object = NULL) {
  inferred <- factor_treatment_levels(data, a_names)
  if (is.list(weight_object) && !is.data.frame(weight_object)) {
    stored <- weight_object$settings$treatment_levels
    ordered <- weight_object$settings$treatment_ordered
    if (!is.null(stored) && !identical(stored, inferred))
      stop("Treatment factor levels differ from weight_object. Use the same factors and refit fit_g().",
           call. = FALSE)
    if (!is.null(ordered) &&
        !identical(ordered, factor_treatment_ordered(data, inferred)))
      stop("Treatment factor ordering differs from weight_object. Refit fit_g() after changing ordered status.",
           call. = FALSE)
  }
  inferred
}

encode_treatment_design <- function(X, treatment_levels) {
  X <- as.data.frame(X)
  if (!length(treatment_levels)) return(X)
  out <- list()
  for (nm in names(X)) {
    source <- if (nm %in% names(treatment_levels)) nm else {
      matches <- names(treatment_levels)[vapply(names(treatment_levels), function(a) {
        startsWith(nm, a) &&
          grepl("^(_lag|_t)[0-9]+$", substring(nm, nchar(a) + 1L))
      }, logical(1L))]
      if (length(matches) > 1L)
        stop("Ambiguous treatment history column: ", nm, call. = FALSE)
      if (length(matches)) matches else NULL
    }
    lev <- if (length(source)) treatment_levels[[source]] else NULL
    if (is.null(lev)) {
      out[[nm]] <- X[[nm]]
      next
    }
    x <- X[[nm]]
    if (!is.factor(x) || !identical(levels(x), lev))
      stop("Categorical treatment predictor '", nm, "' must retain its declared factor levels.",
           call. = FALSE)
    cn <- paste0(nm, "__level_", seq_along(lev)[-1L])
    if (any(cn %in% names(X)))
      stop("Treatment indicator names conflict with existing predictor columns for '", nm, "'.",
           call. = FALSE)
    out[cn] <- lapply(lev[-1L], function(value) as.numeric(as.character(x) == value))
  }
  if (anyDuplicated(names(out)))
    stop("Treatment encoding produced duplicated predictor names.", call. = FALSE)
  as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
}
