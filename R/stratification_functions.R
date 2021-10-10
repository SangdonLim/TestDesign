#' @include shadow_functions.R
NULL

#' @noRd
dbind <- function(...) {

  x <- list(...)

  n_rows <- lapply(x, function(xx) dim(xx)[1])
  n_cols <- lapply(x, function(xx) dim(xx)[2])
  n_rows_total <- sum(unlist(n_rows))
  n_cols_total <- sum(unlist(n_cols))

  o <- matrix(0, n_rows_total, n_cols_total)

  leftpad_row <- 0
  leftpad_col <- 0

  for (i in 1:length(x)) {
    o[
      leftpad_row + (1:n_rows[[i]]),
      leftpad_col + (1:n_cols[[i]])
    ] <- x[[i]]
    leftpad_row <- leftpad_row + n_rows[[i]]
    leftpad_col <- leftpad_col + n_cols[[i]]
  }

  return(o)

}

#' Stratify item pool
#'
#' \code{\link{stratifyItemPool}} is function to perform item pool stratification, using the linear programming method described in Chang & van der Linden (2003).
#'
#' @param cfg a \code{\linkS4class{config_Static}} object. Used for solver configuration.
#' @param item_pool the item pool to be stratified.
#' @param target_a target a-parameter values to use in stratification. This determines the number of produced strata.
#' @param target_b target b-parameter values to use in stratification. This affects the evenness of b-distributions between strata. This does not affect the number of produced strata.
#'
#' @return \code{\link{stratifyItemPool}} returns a vector containing stratum indices for each item.
#'
#' @examples
#'
#' @template alpha-ref
#'
#' @export
setGeneric(
  name = "stratifyItemPool",
  def = function(cfg, item_pool, target_a, target_b) {
    standardGeneric("stratifyItemPool")
  }
)

#' @export
setMethod(
  f = "stratifyItemPool",
  signature = "config_Static",
  definition = function(cfg, item_pool, target_a, target_b) {

    bin_targets <-
      expand.grid(
        a = target_a,
        b = target_b
      )
    n_bins <- dim(bin_targets)[1]
    ni <- item_pool@ni

    pool_a <- item_pool@ipar[, 1]
    pool_b <- item_pool@ipar[, 2]

    scaling_factor <- computeScalingFactor(pool_a, pool_b)

    d_ib <- computeDistanceToTargetItemParameters(
      pool_a, pool_b, bin_targets, scaling_factor
    )

    # bin assignment constraint: an item must be assigned to no more than one bin
    mat_ba <- matrix(0, ni, length(d_ib))
    for (i in 1:ni) {
      mat_ba[i, seq(i, length(d_ib), ni)] <- 1
    }
    dir_ba <- rep("==", ni)
    rhs_ba <- rep(1   , ni)

    # bin size constraint:
    bin_size <- ni / n_bins
    if (floor(bin_size) == bin_size) {
      mat_bs <- matrix(0, n_bins, length(d_ib))
      for (b in 1:n_bins) {
        mat_bs[b, ((b - 1) * ni) + 1:ni] <- 1
      }
      dir_bs <- rep("=="    , n_bins)
      rhs_bs <- rep(bin_size, n_bins)
    }
    if (floor(bin_size) != bin_size) {
      warning(sprintf(
        "bin size %1.2f (%s / %s) is not an integer, using the two nearest integer values instead",
        bin_size,
        ni, n_bins
      ))
      mat_bs <- matrix(0, n_bins * 2, length(d_ib))
      for (b in 1:n_bins) {
        mat_bs[((b - 1) * 2) + 1, ((b - 1) * ni) + 1:ni] <- 1
        mat_bs[((b - 1) * 2) + 2, ((b - 1) * ni) + 1:ni] <- 1
      }
      dir_bs <- rep(c(">=", "<="), n_bins)
      rhs_bs <- rep(c(floor(bin_size), ceiling(bin_size)), n_bins)
    }

    mat <- rbind(mat_ba, mat_bs)
    dir <-     c(dir_ba, dir_bs)
    rhs <-     c(rhs_ba, rhs_bs)

    types <- rep("B", length(d_ib))

    solution <- runMIP(
      toupper(cfg@MIP$solver),
      d_ib, mat, dir, rhs,
      maximize = FALSE,
      types = types,
      verbosity = cfg@MIP$verbosity,
      time_limit = cfg@MIP$time_limit,
      gap_limit_abs = cfg@MIP$gap_limit_abs,
      gap_limit = cfg@MIP$gap_limit
    )

    solution_per_bin <- splitSolutionToBins(solution$solution, n_bins, ni)
    items_per_stratum <- aggregateBinsToStrata(solution_per_bin, bin_targets)

    idx_stratum <- numeric(item_pool@ni)

    for (idx_a in 1:length(target_a)) {
      idx_stratum[items_per_stratum[[idx_a]]] <- idx_a
    }

    return(idx_stratum)

  }
)

#' @noRd
computeScalingFactor <- function(pool_a, pool_b) {
  scaling_factor <-
    (max(pool_b) - min(pool_b)) /
    (max(pool_a) - min(pool_a))
  return(scaling_factor)
}

#' @noRd
computeDistanceToTargetItemParameters <- function(pool_a, pool_b, bin_targets, scaling_factor) {

  n_bins <- dim(bin_targets)[1]

  d_ib <- list()
  for (b in 1:n_bins) {
    d_ib[[b]] <-
      sqrt(
        (((pool_b - bin_targets$b[b]) ** 2)) +
        (((pool_a - bin_targets$a[b]) ** 2) * (scaling_factor ** 2))
      )
  }
  d_ib <- unlist(d_ib)

  return(d_ib)

}

#' @noRd
splitSolutionToBins <- function(solution, n_bins, ni_per_bin) {
  pool_size <- length(solution) / n_bins
  o <- split(solution, ceiling(seq_along(solution) / pool_size))
  o <- lapply(o, function(x) x[1:ni_per_bin])
  return(o)
}

#' @noRd
aggregateBinsToStrata <- function(solution_per_bin, bin_targets) {

  target_a <- unique(bin_targets$a)

  o <- list()

  for (idx_a in 1:length(target_a)) {
    bin_idx <- which(bin_targets$a == target_a[idx_a])
    x <- solution_per_bin[bin_idx]
    x <- do.call(rbind, x)
    x <- apply(x, 2, sum)
    o[[idx_a]] <- list()
    o[[idx_a]] <- which(x == 1)
  }

  return(o)

}

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
