#' @include import.R
NULL

#' @noRd
supported_solver_list <- function() {
  x <- c("LPSYMPHONY", "RSYMPHONY", "LPSOLVE", "GUROBI", "RGLPK")
  return(x)
}

#' @noRd
supported_item_selection_method_list <- function() {
  x <- c("MFI", "MPWI", "EB", "FB", "GFI", "FIXED")
  return(x)
}

#' @noRd
supported_refresh_policy_list <- function() {
  x <- c(
    "ALWAYS", "POSITION",
    "INTERVAL", "THRESHOLD", "INTERVAL-THRESHOLD",
    "STIMULUS", "SET", "PASSAGE"
  )
  return(x)
}

#' @noRd
supported_exposure_control_method_list <- function() {
  x <- c(
    "NONE", "ELIGIBILITY",
    "BIGM", "BIGM-BAYESIAN",
    "ALPHA-STRATIFICATION",
    "PROGRESSIVE-RESTRICTED",
    "MULTIPLE-OBJECTIVE",
    "HYBRID-SE",
    "HYBRID-ES"
  )
  return(x)
}

#' @noRd
supported_theta_estimation_method_list <- function() {
  x <- c(
    "EAP",
    "MLE", "MLEF",
    "EB", "FB"
  )
  return(x)
}

#' @noRd
supported_prior_distribution_list <- function() {
  x <- c("NORMAL", "UNIFORM")
  return(x)
}

#' @noRd
accepts_from_list <- function(x) {
  x <- paste(x, collapse = ", ")
  return(x)
}
