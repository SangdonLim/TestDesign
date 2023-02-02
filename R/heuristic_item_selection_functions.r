#' @include static_functions.R
NULL

#' @noRd
computeWeightedDeviation <- function(v, lb, ub, w) {
  d <-
    ((lb - v) * (v < lb)) +
    ((v - ub) * (v > ub))
  wd <- w * d
  return(sum(wd))
}

#' @noRd
selectItemUsingWeightedDeviation <- function(
  info, position, o, constants, interim_theta, constraints,
  previous_selection, groupings_record
) {

  administered_items <- o@administered_item_index[0:(position - 1)]

  candidate_tests <- lapply(
    1:constants$ni,
    function(candidate_item) {
      c(administered_items, candidate_item)
    }
  )

  # filter out items in completed domains --------------------------------------
  if (constants$group_by_domain) {

    idx_to_filter_out <- constraints@item_attrib@data$DOMAIN %in% groupings_record$completed_domains
    candidate_tests[idx_to_filter_out] <- NA

  }

  # if previous selection was the last item in the domain ----------------------
  if (constants$group_by_domain) {

    if (!is.null(previous_selection$domain_selected)) {

      if (
        previous_selection$domain_selected %in%
        groupings_record$completed_domains
      ) {
        previous_selection$domain_selected <- NULL
      }

    }

  }

  # if domain is selected, filter out items in other domains -------------------
  if (!is.null(previous_selection$domain_selected)) {
    n_domain_items <- sum(
      o@administered_domain_index[0:(position - 1)] ==
      previous_selection$domain_selected
    )
    # if domain item limit has not been reached yet ----------------------------
    if (n_domain_items < 15) {
      idx_to_filter_out <- constraints@item_attrib@data$DOMAIN != previous_selection$domain_selected
      candidate_tests[idx_to_filter_out] <- NA
    }
  }

  # filter out administered items ----------------------------------------------
  candidate_tests[administered_items] <- NA

  wd_c <- lapply(
    candidate_tests,
    function(candidate_test) {
      if (all(is.na(candidate_test))) {
        return(Inf)
      }
      v <- constants$h$value[, candidate_test, drop = FALSE]
      v <- apply(v, 1, sum, simplify = TRUE)
      wd <- computeWeightedDeviation(
        v,
        constants$h$lb,
        constants$h$ub,
        constants$h$weight
      )
      return(wd)
    }
  )
  wd_c <- unlist(wd_c)

  # add information from administered items -----------------------------------
  if (length(administered_items) > 0) {
    info_from_administered_items <-
      Reduce("+", info[administered_items])
    info <- lapply(
      info, function(x) {
        info_from_administered_items + x
      }
    )
  }

  phi <- interim_theta$prior_par[[o@simulee_id]]$prior_sigma
  inv_phi <- solve(phi)

  w_ti <- lapply(info, function(x) x + inv_phi)
  det_w_ti <- lapply(w_ti, det)
  det_w_ti <- unlist(det_w_ti)

  info_lb     <- constants$heuristic_information_lb
  info_weight <- constants$heuristic_information_weight

  wd_ti <- info_weight *
    (info_lb - det_w_ti) * (det_w_ti < info_lb)

  wd <- wd_c + wd_ti

  o <- list()
  o$item_selected <- which.min(wd)

  if (constants$group_by_domain) {
    o$domain_selected <- constraints@item_attrib@data$DOMAIN[o$item_selected]
  }

  return(o)

}

#' @noRd
selectItemUsingMaximumPriorityIndex <- function(
  info, position, o, constants, interim_theta, constraints,
  previous_selection, groupings_record
) {

  administered_items <- o@administered_item_index[0:(position - 1)]

  current_value <- apply(
    constants$h$value[, administered_items, drop = FALSE],
    1, sum
  )

  all_lb_reached <- all(current_value >= constants$h$lb)

  if (!all_lb_reached) {
    # phase one: target lb
    current_value <-
      (constants$h$lb - current_value) / constants$h$lb
  }
  if (all_lb_reached) {
    # phase two: target ub
    current_value <-
      (constants$h$ub - current_value) / constants$h$ub
  }

  wd_c <- constants$h$weight * current_value
  x_ci <- constants$h$value > 0
  wd_ci <- wd_c[, 1] ** x_ci
  wd_i <- apply(wd_ci, 2, prod)

  # add information from administered items ------------------------------------
  if (length(administered_items) > 0) {
    info_from_administered_items <-
    Reduce("+", info[administered_items])
    info <- lapply(
      info, function(x) {
        info_from_administered_items + x
      }
    )
  }

  phi <- interim_theta$prior_par[[o@simulee_id]]$prior_sigma
  inv_phi <- solve(phi)

  w_ti <- lapply(info, function(x) x + inv_phi)
  det_w_ti <- lapply(w_ti, det)
  det_w_ti <- unlist(det_w_ti)

  wd <- det_w_ti * wd_i

  # filter out items in completed domains --------------------------------------
  if (constants$group_by_domain) {

    idx_to_filter_out <- constraints@item_attrib@data$DOMAIN %in% groupings_record$completed_domains
    wd[idx_to_filter_out] <- -Inf

  }

  # if previous selection was the last item in the domain ----------------------
  if (constants$group_by_domain) {

    if (!is.null(previous_selection$domain_selected)) {

      if (
        previous_selection$domain_selected %in%
        groupings_record$completed_domains
      ) {
        previous_selection$domain_selected <- NULL
      }

    }

  }

  # if domain is selected, filter out items in other domains -------------------
  if (!is.null(previous_selection$domain_selected)) {
    n_domain_items <- sum(
      o@administered_domain_index[0:(position - 1)] ==
      previous_selection$domain_selected
    )
    # if domain item limit has not been reached yet ----------------------------
    if (n_domain_items < 15) {
      idx_to_filter_out <- constraints@item_attrib@data$DOMAIN != previous_selection$domain_selected
      wd[idx_to_filter_out] <- -Inf
    }
  }

  # filter out administered items ----------------------------------------------
  wd[administered_items] <- -Inf

  o <- list()
  o$item_selected <- which.max(wd)

  if (constants$group_by_domain) {
    o$domain_selected <- constraints@item_attrib@data$DOMAIN[o$item_selected]
  }

  return(o)

}
