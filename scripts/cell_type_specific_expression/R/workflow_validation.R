validate_boolean_flag <- function(value, label) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("%s must be true or false", label), call. = FALSE)
  }
  value
}

validate_proportion_mode <- function(precomputed_proportions_defined) {
  precomputed_proportions_defined <- validate_boolean_flag(
    precomputed_proportions_defined,
    "precomputed_proportions_defined"
  )
  if (precomputed_proportions_defined) {
    list(selected_mode = "precomputed", estimate_proportions = FALSE)
  } else {
    list(selected_mode = "hspe", estimate_proportions = TRUE)
  }
}
