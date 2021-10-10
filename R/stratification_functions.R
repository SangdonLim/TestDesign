#' @include shadow_functions.R
NULL

#' @noRd
applyStratificationConstraintsToXdata <- function(xdata, idx_stratum, o, constraints) {

  idx_exclude <- constraints@item_attrib@data$STRATUM != idx_stratum

  # don't exclude administered items; otherwise it leads to incompatible constraints
  idx_exclude[na.omit(o@administered_item_index)] <- FALSE

  # augment constraints
  xmat_augment <- matrix(0, 1, constraints@nv)
  xmat_augment[, idx_exclude] <- 1
  xdir_augment <- "=="
  xrhs_augment <- 0

  xdata_augmented = list(
    xmat = rbind(xmat_augment, xdata$xmat),
    xdir =     c(xdir_augment, xdata$xdir),
    xrhs =     c(xrhs_augment, xdata$xrhs))

  return(xdata_augmented)

}

#' @noRd
getStratumForCurrentPosition <- function(position, constants, constraints) {

  n_strata <- constants$n_strata
  stratum_size <- constants$test_length / n_strata
  stratum_begins_from <- seq(1, constants$test_length, stratum_size)
  idx_stratum <- sum(position >= stratum_begins_from)

  return(idx_stratum)

}
