#' @include shadow_functions.R
NULL

#' @noRd
assembleShadowTest <- function(
  j, position, o,
  eligible_flag,
  exclude_index,
  stimulus_record,
  info,
  config,
  constants,
  constraints
) {

  administered_stimulus_index <- na.omit(unique(o@administered_stimulus_index))

  xdata         <- getXdataOfAdministered(constants, position, o, stimulus_record, constraints)
  xdata_exclude <- getXdataOfExcludedEntry(constants, exclude_index[[j]])
  xdata         <- combineXdata(xdata, xdata_exclude)

  if (constants$use_eligibility_control) {

    # Get eligible items in the current theta segment
    current_segment <- o@theta_segment_index[position]
    eligible_flag_in_current_theta_segment <- getEligibleFlagInSegment(eligible_flag, current_segment, constants)
    eligible_flag_in_current_theta_segment <- flagAdministeredAsEligible(eligible_flag_in_current_theta_segment, o, position, constants)
    # Allow administered items to be selected, because otherwise constraints will contradict each other

  }

  if (constants$use_eligibility_control && constants$exposure_control_method %in% c("ELIGIBILITY")) {

    xdata_elg  <- applyEligibilityConstraintsToXdata(xdata, eligible_flag_in_current_theta_segment, constants, constraints)
    shadowtest <- runAssembly(config, constraints, xdata = xdata_elg, objective = info)
    is_optimal <- isShadowtestOptimal(shadowtest)

    if (is_optimal) {
      shadowtest$feasible <- TRUE
      return(shadowtest)
    }

    # If not optimal, retry without xmat
    shadowtest <- runAssembly(config, constraints, xdata = xdata, objective = info)
    shadowtest$feasible <- FALSE
    return(shadowtest)

  }

  if (constants$use_eligibility_control && constants$exposure_control_method %in% c("BIGM", "BIGM-BAYESIAN")) {

    # Do Big-M based exposure control: penalize item info
    info <- applyEligibilityConstraintsToInfo(
      info, eligible_flag_in_current_theta_segment, config, constants
    )

    shadowtest <- runAssembly(config, constraints, xdata = xdata, objective = info)
    shadowtest$feasible <- TRUE
    return(shadowtest)

  }

  if (constants$use_eligibility_control && constants$exposure_control_method %in% c("ALPHA-STRATIFICATION")) {

    # Do alpha-stratification
    idx_stratum   <- getStratumForCurrentPosition(position, constants, constraints)
    xdata_str     <- applyStratificationConstraintsToXdata(xdata, idx_stratum, o, constraints)

    shadowtest <- runAssembly(config, constraints, xdata = xdata_str, objective = info)
    is_optimal <- isShadowtestOptimal(shadowtest)

    if (is_optimal) {
      shadowtest$feasible <- TRUE
      return(shadowtest)
    }

    # If not optimal, retry without str
    shadowtest <- runAssembly(config, constraints, xdata = xdata, objective = info)
    shadowtest$feasible <- FALSE
    return(shadowtest)

  }

  if (constants$use_eligibility_control && constants$exposure_control_method %in% c("PROGRESSIVE-RESTRICTED")) {

    # Progressive component
    weight <- (position - 1) / constants$min_ni
    # numerator is the number of previous items; weight is 0 for the first item
    # TODO: prevent doing this when min_ni != max_ni
    info_unused <- info
    info_unused[o@administered_item_index, ] <- 0
    max_info_of_unused_items <- max(info_unused)
    # unused items are items excluding items administered to this examinee
    random_component <- runif(constants$ni, 0, max_info_of_unused_items)
    weighted_info <- (weight * info) + ((1 - weight) * random_component)

    # Restricted component
    xdata_elg  <- applyEligibilityConstraintsToXdata(xdata, eligible_flag_in_current_theta_segment, constants, constraints)

    # Main assembly
    shadowtest <- runAssembly(config, constraints, xdata = xdata_elg, objective = weighted_info)
    is_optimal <- isShadowtestOptimal(shadowtest)

    if (is_optimal) {
      shadowtest$feasible <- TRUE
      return(shadowtest)
    }

    # If not optimal, retry without xmat
    shadowtest <- runAssembly(config, constraints, xdata = xdata, objective = weighted_info)
    shadowtest$feasible <- FALSE
    return(shadowtest)

  }

  if (!constants$use_eligibility_control) {

    shadowtest <- runAssembly(config, constraints, xdata = xdata, objective = info)
    shadowtest$feasible <- TRUE
    return(shadowtest)

  }

}

#' @noRd
isShadowtestOptimal <- function(shadowtest) {
  return(isOptimal(shadowtest$status, shadowtest$solver))
}
