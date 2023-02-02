#' @include shadow_functions.R
NULL

#' @noRd
getConstants <- function(constraints, config, arg_data, true_theta, max_info) {

  o <- list()
  o$ni <- constraints@ni
  o$ns <- constraints@ns
  o$nv <- constraints@nv
  o$nd <- constraints@pool@nd

  o$theta_q <- config@theta_grid
  if (inherits(o$theta_q, "numeric")) {
    o$theta_q <- matrix(o$theta_q, , 1)
  }
  o$nq      <- nrow(o$theta_q)
  o$min_q   <- min(o$theta_q)
  o$max_q   <- max(o$theta_q)

  if (!is.null(arg_data)) {
    o$nj <- nrow(arg_data)
  }
  if (!is.null(true_theta)) {
    o$nj <- nrow(true_theta)
  }
  if (is.null(o$nj)) {
    stop("either 'data' or 'true_theta' must be supplied")
  }

  content_balancing_method <- toupper(config@content_balancing$method)
  if (content_balancing_method %in% c("STA", "SHADOW", "SHADOWTEST", "SHADOW TEST")) {
    if (is.null(constraints)) {
      stop(sprintf("config@content_balancing: 'constraints' must be supplied when $method is '%s'", content_balancing_method))
    }
    o$use_shadowtest    <- TRUE
    o$group_by_stimulus <- constraints@set_based
    o$group_by_domain   <- constraints@group_by_domain
    o$test_length       <- constraints@test_length
    o$min_ni            <- constraints@test_length
    o$max_ni            <- constraints@test_length
    o$max_se            <- NULL
  }
  if (content_balancing_method %in% "HEURISTIC") {
    o$use_shadowtest    <- FALSE
    o$group_by_stimulus <- FALSE
    o$group_by_domain   <- constraints@group_by_domain
    o$test_length       <- constraints@test_length
    o$min_ni            <- constraints@test_length
    o$max_ni            <- constraints@test_length
    o$max_se            <- NULL
  }
  if (content_balancing_method %in% "NONE") {
    o$use_shadowtest    <- FALSE
    o$group_by_stimulus <- FALSE
    o$group_by_domain   <- FALSE
    o$test_length       <- NULL
    o$min_ni            <- config@stopping_criterion$min_ni
    o$max_ni            <- config@stopping_criterion$max_ni
    o$max_se            <- config@stopping_criterion$se_threshold
  }

  # parse constraints for heuristic methods ------------------------------------
  if (content_balancing_method %in% "HEURISTIC") {

    constraints_data <- constraints@constraints

    h <- list()

    for (idx_c in 1:nrow(constraints_data)) {

      if (!constraints_data$TYPE[idx_c] %in% c("NUMBER", "SUM")) {
        stop("constraints on other than the number of items are not supported on heuristic methods")
      }

      if (constraints_data$TYPE[idx_c] == "NUMBER") {
        if (constraints_data$CONDITION[idx_c] == "") {
          # skip test length constraint
          next
        }
        v <- with(constraints@item_attrib@data, eval(parse(text = constraints_data$CONDITION[idx_c])))
      }
      if (constraints_data$TYPE[idx_c] == "SUM") {
        v <- with(constraints@item_attrib@data, eval(parse(text = constraints_data$CONDITION[idx_c])))
      }

      h[[idx_c]] <- list()
      h[[idx_c]]$value  <- v
      h[[idx_c]]$p      <- mean(v)
      h[[idx_c]]$lb     <- constraints_data$LB[idx_c]
      h[[idx_c]]$ub     <- constraints_data$UB[idx_c]
      h[[idx_c]]$lp     <- constraints_data$LB[idx_c] / o$test_length
      h[[idx_c]]$up     <- constraints_data$UB[idx_c] / o$test_length
      h[[idx_c]]$mp     <- mean(c(h[[idx_c]]$lp, h[[idx_c]]$up))
      h[[idx_c]]$weight <- constraints_data$WEIGHT[idx_c]

    }

    hh <- list()
    hh$value  <- do.call(rbind, lapply(h, function(x) x$value))
    hh$value_binary <- hh$value > 0
    hh$p      <- do.call(rbind, lapply(h, function(x) x$p))
    hh$lb     <- do.call(rbind, lapply(h, function(x) x$lb))
    hh$ub     <- do.call(rbind, lapply(h, function(x) x$ub))
    hh$lp     <- do.call(rbind, lapply(h, function(x) x$lp))
    hh$up     <- do.call(rbind, lapply(h, function(x) x$up))
    hh$mp     <- do.call(rbind, lapply(h, function(x) x$mp))
    hh$ld     <- hh$lp - hh$mp
    hh$ud     <- hh$up - hh$mp
    hh$weight <- do.call(rbind, lapply(h, function(x) x$weight))
    hh$nc     <- nrow(hh$value)
    hh$k      <- 2

    o$h <- hh

  }

  # parse options for heuristic methods ----------------------------------------
  o$heuristic_information_lb     <- config@content_balancing$heuristic_information_lb
  o$heuristic_information_weight <- config@content_balancing$heuristic_information_weight

  o$exclude_method <- toupper(config@exclude_policy$method)
  o$exclude_M      <- config@exclude_policy$M

  o$use_hand_scored <- !is.null(config@interim_theta$hand_scored_attribute)
  if (o$use_hand_scored) {
    if (length(config@interim_theta$hand_scored_attribute) != 1) {
      stop(sprintf(
        "config@interim_theta$hand_scored_attribute: too many values (expecting one item attribute name)"
      ))
    }
  }
  if (o$use_hand_scored) {
    if (!config@interim_theta$hand_scored_attribute %in% names(constraints@item_attrib)) {
      stop(sprintf(
        "column not found in item attribute table: '%s'",
        config@interim_theta$hand_scored_attribute
      ))
    }
  }
  if (o$use_hand_scored) {
    o$item_is_hand_scored <-
      constraints@item_attrib@data[[config@interim_theta$hand_scored_attribute]]
  }
  if (o$use_hand_scored) {
    # normalize
    if (all(sort(unique(o$item_is_hand_scored)) == c("N", "Y"))) {
      o$item_is_hand_scored <- o$item_is_hand_scored == "Y"
    }
    if (all(sort(unique(o$item_is_hand_scored)) == c("n", "y"))) {
      o$item_is_hand_scored <- o$item_is_hand_scored == "y"
    }
    if (all(sort(unique(o$item_is_hand_scored)) == c("NO", "YES"))) {
      o$item_is_hand_scored <- o$item_is_hand_scored == "YES"
    }
    if (all(sort(unique(o$item_is_hand_scored)) == c("no", "yes"))) {
      o$item_is_hand_scored <- o$item_is_hand_scored == "yes"
    }
    if (all(sort(unique(o$item_is_hand_scored)) == c("0", "1"))) {
      o$item_is_hand_scored <- o$item_is_hand_scored == "1"
    }
    if (!all(sort(unique(o$item_is_hand_scored)) == c(FALSE, TRUE))) {
      stop(sprintf(
        "hand_scored_attribute '%s': {%s} could not be normalized to {FALSE, TRUE} (expecting logical values in item attribute '%s')",
        config@interim_theta$hand_scored_attribute,
        paste0(sort(unique(o$item_is_hand_scored)), collapse = ", "),
        config@interim_theta$hand_scored_attribute
      ))
    }
  }

  o$exposure_control_method <- toupper(config@exposure_control$method)
  if (o$exposure_control_method %in% c("ELIGIBILITY", "BIGM", "BIGM-BAYESIAN")) {
    o$use_eligibility_control <- TRUE
  } else {
    o$use_eligibility_control <- FALSE
  }
  o$exposure_M <- config@exposure_control$M
  if (is.null(o$exposure_M)) {
    o$exposure_M <- max_info + 1
  }

  o$max_exposure_rate   <- config@exposure_control$max_exposure_rate
  o$fading_factor       <- config@exposure_control$fading_factor
  o$acceleration_factor <- config@exposure_control$acceleration_factor
  o$n_segment           <- config@exposure_control$n_segment
  o$segment_cut         <- config@exposure_control$segment_cut
  o$cut_lower           <- o$segment_cut[(1:o$n_segment)]
  o$cut_upper           <- o$segment_cut[(1:o$n_segment) + 1]

  if (!length(o$max_exposure_rate) %in% c(1, o$n_segment)) {
    stop("length(max_exposure_rate) must be 1 or n_segment")
  }

  return(o)

}

#' @noRd
sanitizeModel <- function(model) {
  model[which(model == "item_1PL")] <- 1
  model[which(model == "item_2PL")] <- 2
  model[which(model == "item_3PL")] <- 3
  model[which(model == "item_PC")]  <- 4
  model[which(model == "item_GPC")] <- 5
  model[which(model == "item_GR")]  <- 6
  model[which(model == "item_M2PL")] <- 102
  model[which(model == "item_M3PL")] <- 103
  model[which(model == "item_MGPC")] <- 105
  model[which(model == "item_MGR")]  <- 106
  model <- as.numeric(model)
  return(model)
}

#' @noRd
getInfoFixedTheta <- function(item_selection, constants, item_pool, model) {

  nj <- constants$nj
  o <- list()

  if (!is.null(item_selection$fixed_theta)) {
    if (length(item_selection$fixed_theta) == 1) {
      o$info_fixed_theta <- lapply(seq_len(nj), function(j) calc_info(item_selection$fixed_theta, item_pool@ipar, item_pool@NCAT, model))
      o$select_at_fixed_theta <- TRUE
    }
    if (length(item_selection$fixed_theta) == nj) {
      o$info_fixed_theta <- lapply(seq_len(nj), function(j) calc_info(item_selection$fixed_theta[j], item_pool@ipar, item_pool@NCAT, model))
      o$select_at_fixed_theta <- TRUE
    }
    if (is.null(o$info_fixed_theta)) {
      stop("config@item_selection: length($fixed_theta) must be either 1 or nj")
    }
  } else {
    o$select_at_fixed_theta <- FALSE
  }

  return(o)

}

#' @noRd
computeInfoAtCurrentTheta <- function(
  item_selection,
  j,
  current_theta,
  item_pool,
  model_code,
  info_fixed_theta,
  info_grid,
  prob_grid,
  item_parameter_sample
) {

  item_method <- toupper(item_selection$method)
  info_type   <- toupper(item_selection$info_type)

  if (item_method == "FIXED") {
    info <- info_fixed_theta[[j]]
    return(info)
  }
  if (item_method == "MFI") {
    if (item_pool@nd == 1) {
      info <- calc_info(
        current_theta$theta,
        item_pool@ipar,
        item_pool@NCAT,
        model_code
      )
      return(info)
    }
  }
  if (item_method == "GFI") {
    if (item_pool@nd == 1) {
      info <- calc_info(
        current_theta$theta,
        item_pool@ipar,
        item_pool@NCAT,
        model_code
      )
      return(info)
    }
  }
  if (item_method == "DIRINFO-45") {
    if (item_pool@nd > 1) {
      alpha_vec <- a_to_alpha(rep(1, item_pool@nd))
      info <- calc_thisdirinfo(
        current_theta$theta,
        item_pool@ipar,
        item_pool@nd,
        item_pool@NCAT,
        model_code,
        alpha_vec
      )
      return(info)
    }
  }
  if (item_method == "DIRINFO-THISANGLE") {
    if (item_pool@nd > 1) {
      info <- calc_thisdirinfo(
        current_theta$theta,
        item_pool@ipar,
        item_pool@nd,
        item_pool@NCAT,
        model_code,
        item_selection$alpha_vec
      )
      return(info)
    }
  }
  if (item_method == "MKL") {
    info <- calcKL(
      item_pool,
      matrix(current_theta$theta, 1, ),
      item_selection$KL_width,
      item_selection$KL_quadrature_unit,
      item_selection$KL_use_ellipse
    )
    return(info)
  }
  if (item_method == "MMI") {
    info <- calcMI(
      prob_grid,
      current_theta$posterior
    )
    return(info)
  }
  if (item_method == "MPWI") {
    info <- as.vector(matrix(current_theta$posterior, nrow = 1) %*% info_grid)
    return(info)
  }
  if (item_method == "EB") {
    info <- calc_info_EB(
      matrix(current_theta$posterior_sample),
      item_pool@ipar, item_pool@NCAT, model_code)[, 1]
    return(info)
  }
  if (item_method == "FB" & info_type == "FISHER") {
    info <- calc_info_FB(
      matrix(current_theta$posterior_sample),
      item_parameter_sample, item_pool@NCAT, model_code)[, 1]
    return(info)
  }
  if (item_method == "FB" & info_type %in% c("MI", "MUTUAL")) {
    info <- calc_MI_FB(
      matrix(current_theta$posterior_sample),
      item_parameter_sample, item_pool@NCAT, model_code)[, 1]
    return(info)
  }
  if (item_method == "D-OPTIMAL") {
    info <- calc_info_matrix(
      current_theta$theta,
      item_pool@ipar,
      item_pool@nd,
      item_pool@NCAT,
      model_code
    )
    return(info)
  }
  if (item_method == "WEIGHTED-DEVIATION") {
    info <- calc_info_matrix(
      current_theta$theta,
      item_pool@ipar,
      item_pool@nd,
      item_pool@NCAT,
      model_code
    )
    return(info)
  }
  if (item_method == "MAXIMUM-PRIORITY-INDEX") {
    info <- calc_info_matrix(
      current_theta$theta,
      item_pool@ipar,
      item_pool@nd,
      item_pool@NCAT,
      model_code
    )
    return(info)
  }
}

#' @noRd
initializeCompletedGroupingsRecord <- function() {
  o <- list()
  o$completed_stimulus_index <- NULL
  o$completed_stimulus_size  <- NULL
  o$completed_domains <- NULL
  return(o)
}

#' @noRd
selectItem <- function(info, position, o) {

  info[o@administered_item_index[0:(position - 1)]] <- -1

  info_index    <- order(info, decreasing = TRUE)
  item_selected <- info_index[1]

  if (item_selected %in% o@administered_item_index[0:(position - 1)]) {
    stop(sprintf("item %i has been already administered", item_selected))
  }

  return(item_selected)

}
