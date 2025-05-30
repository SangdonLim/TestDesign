#' @include static_class.R
NULL

#' @rdname createShadowTestConfig
#' @export
setClass("config_Shadow",
  slots = c(
    item_selection     = "list",
    content_balancing  = "list",
    MIP                = "list",
    MCMC               = "list",
    exclude_policy     = "list",
    refresh_policy     = "list",
    exposure_control   = "list",
    overlap_control    = "list",
    stopping_criterion = "list",
    interim_theta      = "list",
    final_theta        = "list",
    theta_grid         = "numeric"
  ),
  prototype = list(
    item_selection = list(
      method                    = "MFI",
      info_type                 = "FISHER",
      initial_theta             = NULL,
      fixed_theta               = NULL,
      target_value              = NULL
    ),
    content_balancing = list(
      method                    = "STA"
    ),
    MIP = list(
      solver                    = "LPSOLVE",
      verbosity                 = -2,
      time_limit                = 60,
      gap_limit                 = .05,
      gap_limit_abs             = .05,
      obj_tol                   = .05,
      retry                     = 5
    ),
    MCMC = list(
      burn_in                   = 100,
      post_burn_in              = 500,
      thin                      = 1,
      jump_factor               = 2.4
    ),
    exclude_policy = list(
      method = "HARD",
      M = 1000
    ),
    refresh_policy = list(
      method                    = "ALWAYS",
      interval                  = 1,
      threshold                 = 0.1,
      position                  = 1
    ),
    exposure_control = list(
      method                    = "ELIGIBILITY",
      M                         = NULL,
      max_exposure_rate         = 0.25,
      acceleration_factor       = 1.0,
      n_segment                 = 7,
      first_segment             = NULL,
      segment_cut               = c(-Inf, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, Inf),
      initial_eligibility_stats = NULL,
      fading_factor             = 0.999,
      diagnostic_stats          = FALSE
    ),
    overlap_control = list(
      method                    = "NONE",
      M                         = 100,
      max_overlap_rate          = 0.20
    ),
    stopping_criterion = list(
      method                    = "FIXED",
      test_length               = NULL,
      min_ni                    = NULL,
      max_ni                    = NULL,
      se_threshold              = NULL
    ),
    interim_theta = list(
      method                    = "EAP",
      shrinkage_correction      = FALSE,
      prior_dist                = "NORMAL",
      prior_par                 = c(0, 1),
      bound_ML                  = c(-4, 4),
      truncate_ML               = FALSE,
      max_iter                  = 50,
      crit                      = 0.001,
      max_change                = 1.0,
      use_step_size             = FALSE,
      step_size                 = 0.5,
      do_Fisher                 = TRUE,
      fence_slope               = 5,
      fence_difficulty          = c(-5, 5),
      hand_scored_attribute     = NULL
    ),
    final_theta = list(
      method                    = "EAP",
      shrinkage_correction      = FALSE,
      prior_dist                = "NORMAL",
      prior_par                 = c(0, 1),
      bound_ML                  = c(-4, 4),
      truncate_ML               = FALSE,
      max_iter                  = 50,
      crit                      = 0.001,
      max_change                = 1.0,
      use_step_size             = FALSE,
      step_size                 = 0.5,
      do_Fisher                 = TRUE,
      fence_slope               = 5,
      fence_difficulty          = c(-5, 5)
    ),
    theta_grid                  = seq(-4, 4, .1)
  ),
  validity = function(object) {

    err <- NULL
    if (!object@MIP$solver %in% c("RSYMPHONY", "HIGHS", "GUROBI", "LPSOLVE", "RGLPK")) {
      # only capitalized names are valid values;
      # the rest of the package assumes this is capitalized
      msg <- sprintf("config@MIP: unrecognized $solver '%s' (accepts RSYMPHONY, HIGHS, GUROBI, LPSOLVE, or RGLPK)", object@MIP$solver)
      err <- c(err, msg)
    }

    for (solver_name in c("Rsymphony", "highs", "gurobi", "Rglpk")) {
      if (object@MIP$solver == toupper(solver_name)) {
        if (!requireNamespace(solver_name, quietly = TRUE)) {
          msg <- sprintf("config@MIP: could not find the specified solver package '%s'", solver_name)
          err <- c(err, msg)
        }
      }
    }

    if (object@MIP$solver == "HIGHS") {
      warning(sprintf(
        "HIGHS solver is unstable as of 1.10.0.1 (installed version: %s) - it sometimes returns 0.5 as a solution for binary decision variables",
        packageVersion("highs")
      ))
    }

    if (!toupper(object@item_selection$method) %in% c("MFI", "MPWI", "EB", "FB", "GFI", "FIXED", "RANDOM")) {
      msg <- sprintf("config@item_selection: unrecognized $method '%s' (accepts MFI, MPWI, EB, FB, GFI, FIXED or RANDOM)", object@item_selection$method)
      err <- c(err, msg)
    }
    if (toupper(object@item_selection$method) %in% c("FIXED")) {
      if (is.null(object@item_selection$fixed_theta)) {
        msg <- sprintf("config@item_selection: $method 'FIXED' requires $fixed_theta")
        err <- c(err, msg)
      }
    }

    if (!object@content_balancing$method %in% c("NONE", "STA")) {
      msg <- sprintf("config@content_balancing: unrecognized $method '%s' (accepts NONE, or STA)", object@content_balancing$method)
      err <- c(err, msg)
    }
    if (!object@exclude_policy$method %in%
      c("HARD", "SOFT")) {
      msg <- sprintf("config@exclude_policy: unrecognized $method '%s' (accepts HARD, or SOFT)", object@exclude_policy$method)
      err <- c(err, msg)
    }
    if (object@exclude_policy$method == "SOFT") {
      if (!is.numeric(object@exclude_policy$M)) {
        if (!is.null(object@exclude_policy$M)) {
          msg <- sprintf("$method 'SOFT' requires $M to be a positive value")
          err <- c(err, msg)
        }
      }
      if (is.numeric(object@exclude_policy$M)) {
        if (object@exclude_policy$M < 0) {
          msg <- sprintf("$method 'SOFT' requires $M to be a positive value")
          err <- c(err, msg)
        }
      }
    }
    if (!object@refresh_policy$method %in%
      c("ALWAYS", "NONE", "POSITION", "INTERVAL", "THRESHOLD", "INTERVAL-THRESHOLD", "STIMULUS", "SET", "PASSAGE")) {
      msg <- sprintf("config@refresh_policy: unrecognized $method '%s'", object@refresh_policy$method)
      err <- c(err, msg)
    }
    if (!object@exposure_control$method %in% c("NONE", "ELIGIBILITY", "BIGM", "BIGM-BAYESIAN")) {
      msg <- sprintf("config@exposure_control: unrecognized $method '%s' (accepts NONE, ELIGIBILITY, BIGM, or BIGM-BAYESIAN)", object@exposure_control$method)
      err <- c(err, msg)
    }
    if (!object@overlap_control$method %in% c("NONE", "ELIGIBILITY", "BIGM", "BIGM-BAYESIAN")) {
      msg <- sprintf("config@overlap_control: unrecognized $method '%s' (accepts NONE, BIGM, or BIGM-BAYESIAN)", object@overlap_control$method)
      err <- c(err, msg)
    }
    if (object@exposure_control$method %in% c("BIGM", "BIGM-BAYESIAN")) {
      if (!is.numeric(object@exposure_control$M)) {
        if (!is.null(object@exposure_control$M)) {
          msg <- sprintf("$method 'BIGM', 'BIGM-BAYESIAM' requires $M to be a positive value")
          err <- c(err, msg)
        }
      }
      if (is.numeric(object@exposure_control$M)) {
        if (object@exposure_control$M < 0) {
          msg <- sprintf("$method 'BIGM', 'BIGM-BAYESIAM' requires $M to be a positive value")
          err <- c(err, msg)
        }
      }
    }
    if (object@overlap_control$method %in% c("BIGM", "BIGM-BAYESIAN")) {
      if (!is.numeric(object@overlap_control$M)) {
        if (!is.null(object@overlap_control$M)) {
          msg <- sprintf("$method 'BIGM', 'BIGM-BAYESIAM' requires $M to be a positive value")
          err <- c(err, msg)
        }
      }
      if (is.numeric(object@overlap_control$M)) {
        if (object@overlap_control$M < 0) {
          msg <- sprintf("$method 'BIGM', 'BIGM-BAYESIAM' requires $M to be a positive value")
          err <- c(err, msg)
        }
      }
    }
    if (object@overlap_control$method == "ELIGIBILITY") {
      if (object@exposure_control$method != "ELIGIBILITY") {
        msg <- sprintf("config@overlap_control: $method 'ELIGIBILITY' also requires @exposure_control$method 'ELIGIBILITY'")
        err <- c(err, msg)
      }
    }
    if (toupper(object@item_selection$method) %in% c("GFI") &
      object@exposure_control$method %in% c("ELIGIBILITY")) {
      msg <- sprintf("config@exposure_control: $method 'ELIGIBILITY' is not compatible with @item_selection$method 'GFI'")
      err <- c(err, msg)
    }
    if (!length(object@exposure_control$max_exposure_rate) %in% c(1, object@exposure_control$n_segment)) {
      msg <- "config@exposure_control: length($max_exposure_rate) must be 1 or $n_segment"
      err <- c(err, msg)
    }
    if (object@exposure_control$n_segment != length(object@exposure_control$segment_cut) - 1) {
      msg <- "config@exposure_control: $n_segment must be equal to length($segment_cut) - 1"
      err <- c(err, msg)
    }
    if (!length(object@exposure_control$max_exposure_rate) %in% c(1, object@exposure_control$n_segment)) {
      msg <- sprintf("config@exposure_control: unexpected length($max_exposure_rate) %s (must be 1 or $n_segment)", length(object@exposure_control$max_exposure_rate))
      err <- c(err, msg)
    }
    if (!is.null(object@exposure_control$first_segment)) {
      if (any(
        !(object@exposure_control$first_segment %in% 1:object@exposure_control$n_segment)
      )) {
        msg <- sprintf(
          "config@exposure_control: out-of-range values in $first_segment %s (must be within 1:n_segment; n_segment is %s)",
          setdiff(object@exposure_control$first_segment, 1:object@exposure_control$n_segment),
          object@exposure_control$n_segment
        )
        err <- c(err, msg)
      }
    }
    if (!object@stopping_criterion$method %in% c("FIXED")) {
      msg <- sprintf("config@stopping_criterion: unrecognized $method '%s'", object@stopping_criterion$method)
      err <- c(err, msg)
    }
    if (!object@interim_theta$method %in% c("EAP", "MLE", "MLEF", "EB", "FB", "CARRYOVER")) {
      msg <- sprintf("config@interim_theta: unrecognized $method '%s' (accepts EAP, MLE, MLEF, EB, FB, or CARRYOVER)", object@interim_theta$method)
      err <- c(err, msg)
    }
    if (!object@interim_theta$prior_dist %in% c("NORMAL", "UNIFORM")) {
      msg <- sprintf("config@interim_theta: unrecognized $prior_dist '%s' (accepts NORMAL or UNIFORM)", object@interim_theta$prior_dist)
      err <- c(err, msg)
    }
    if (!object@final_theta$method %in% c("EAP", "MLE", "MLEF", "EB", "FB", "CARRYOVER")) {
      msg <- sprintf("config@final_theta: unrecognized $method '%s' (accepts EAP, MLE, MLEF, EB, FB, or CARRYOVER)", object@final_theta$method)
      err <- c(err, msg)
    }
    if (toupper(object@final_theta$method) == "EAP") {
      if (!toupper(object@final_theta$prior_dist) %in% c("NORMAL", "UNIFORM")) {
        msg <- sprintf("config@final_theta: unrecognized $prior_dist '%s' (when $method is EAP, accepts NORMAL or UNIFORM)", object@final_theta$prior_dist)
        err <- c(err, msg)
      }
    }
    if ((object@exposure_control$method == c("BIGM-BAYESIAN")) &&
      (!object@interim_theta$method %in% c("EB", "FB"))) {
      err <- c(err, "config@exposure_control: $method 'BIGM-BAYESIAN' requires interim_theta$method to be EB or FB")
    }
    if ((object@overlap_control$method == c("BIGM-BAYESIAN")) &&
        (!object@interim_theta$method %in% c("EB", "FB"))) {
      err <- c(err, "config@overlap_control: $method 'BIGM-BAYESIAN' requires interim_theta$method to be EB or FB")
    }
    if (length(err) == 0) {
      return(TRUE)
    } else {
      return(err)
    }
  }
)

#' Create a config_Shadow object
#'
#' \code{\link{createShadowTestConfig}} is a config function for creating a \code{\linkS4class{config_Shadow}} object for shadowtest assembly.
#' Default values are used for any unspecified parameters/slots.
#'
#' @param item_selection a named list containing item selection criteria.
#' \itemize{
#'   \item{\code{method}} the type of selection criteria. Accepts \code{MFI, MPWI, FB, EB, GFI}. (default = \code{MFI})
#'   \item{\code{info_type}} the type of information. Accepts \code{FISHER}. (default = \code{FISHER})
#'   \item{\code{initial_theta}} (optional) initial theta values to use.
#'   \item{\code{fixed_theta}} (optional) fixed theta values to use throughout all item positions.
#'   \item{\code{target_value}} (optional) the target value to use for \code{method = 'GFI'}.
#' }
#' @param content_balancing a named list containing content balancing options.
#' \itemize{
#'   \item{\code{method}} the type of balancing method. Accepts \code{NONE, STA}. (default = \code{STA})
#' }
#' @param MIP a named list containing solver options.
#' \itemize{
#'   \item{\code{solver}} the type of solver. Accepts \code{Rsymphony, highs, gurobi, lpSolve, Rglpk}. (default = \code{LPSOLVE})
#'   \item{\code{verbosity}} verbosity level of the solver. (default = \code{-2})
#'   \item{\code{time_limit}} time limit in seconds. Used in solvers \code{Rsymphony, gurobi, Rglpk}. (default = \code{60})
#'   \item{\code{gap_limit}} search termination criterion. Gap limit in relative scale passed onto the solver. Used in solver \code{gurobi}. (default = \code{.05})
#'   \item{\code{gap_limit_abs}} search termination criterion. Gap limit in absolute scale passed onto the solver. Used in solvers \code{Rsymphony}. (default = \code{0.05})
#'   \item{\code{obj_tol}} search termination criterion. The lower bound to use on the minimax deviation variable. Used when \code{item_selection$method} is \code{GFI}, and ignored otherwise. (default = \code{0.05})
#'   \item{\code{retry}} number of times to retry running the solver if the solver returns no solution. Some solvers incorrectly return no solution even when a solution exists. This is the number of attempts to verify that the problem is indeed infeasible in such cases. Set to \code{0} to not retry. (default = \code{5})
#' }
#' @param MCMC a named list containing Markov-chain Monte Carlo configurations for obtaining posterior samples.
#' \itemize{
#'   \item{\code{burn_in}} the number of chains from the start to discard. (default = \code{100})
#'   \item{\code{post_burn_in}} the number of chains to use after discarding the first \code{burn_in} chains. (default = \code{500})
#'   \item{\code{thin}} thinning interval to apply. \code{1} represents no thinning. (default = \code{1})
#'   \item{\code{jump_factor}} the jump (scaling) factor for the proposal distribution. \code{1} represents no jumping. (default = \code{2.4})
#' }
#' @param exclude_policy a named list containing the exclude policy for use with the \code{exclude} argument in \code{\link{Shadow}}.
#' \itemize{
#'   \item{\code{method}} the type of policy. Accepts \code{HARD, SOFT}. (default = \code{HARD})
#'   \item{\code{M}} the Big M penalty to use on item information. Used in the \code{SOFT} method.
#' }
#' @param refresh_policy a named list containing the refresh policy for when to obtain a new shadowtest.
#' \itemize{
#'   \item{\code{method}} the type of policy. Accepts \code{ALWAYS, POSITION, INTERVAL, THRESHOLD, INTERVAL-THRESHOLD, STIMULUS, SET, PASSAGE}. (default = \code{ALWAYS})
#'   \item{\code{interval}} used in methods \code{INTERVAL, INTERVAL-THRESHOLD}. Set to 1 to refresh at each position, 2 to refresh at every two positions, and so on. (default = \code{1})
#'   \item{\code{threshold}} used in methods \code{THRESHOLD, INTERVAL-THRESHOLD}. The absolute change in between interim theta estimates to trigger the refresh. (default = \code{0.1})
#'   \item{\code{position}} used in methods \code{POSITION}. Item positions to trigger the refresh. (default = \code{1})
#' }
#' @param exposure_control a named list containing exposure control settings.
#' \itemize{
#'   \item{\code{method}} the type of exposure control method. Accepts \code{NONE, ELIGIBILITY, BIGM, BIGM-BAYESIAN}. (default = \code{ELIGIBILITY})
#'   \item{\code{M}} used in methods \code{BIGM, BIGM-BAYESIAN}. the Big M penalty to use on item information.
#'   \item{\code{max_exposure_rate}} target exposure rates for each segment. (default = \code{rep(0.25, 7)})
#'   \item{\code{acceleration_factor}} the acceleration factor to apply. (default = \code{1})
#'   \item{\code{n_segment}} the number of theta segments to use. (default = \code{7})
#'   \item{\code{first_segment}} (optional) the theta segment assumed at the beginning of test for all participants.
#'   \item{\code{segment_cut}} theta segment cuts. (default = \code{c(-Inf, seq(-2.5, 2.5, 1), Inf)})
#'   \item{\code{initial_eligibility_stats}} (optional) initial eligibility statistics to use.
#'   \item{\code{fading_factor}} the fading factor to apply. (default = \code{.999})
#'   \item{\code{diagnostic_stats}} set to \code{TRUE} to generate segment-wise diagnostic statistics. (default = \code{FALSE})
#' }
#' @param overlap_control a named list containing overlap control settings.
#' \itemize{
#'   \item{\code{method}} the type of overlap control method. Accepts \code{NONE, ELIGIBILITY, BIGM, BIGM-BAYESIAN}. (default = \code{NONE})
#'   \item{\code{M}} used in methods \code{BIGM, BIGM-BAYESIAN}. the Big M penalty to use on item information.
#'   \item{\code{max_overlap_rate}} target overlap rate. (default = \code{0.20})
#' }
#' @param stopping_criterion a named list containing stopping criterion.
#' \itemize{
#'   \item{\code{method}} the type of stopping criterion. Accepts \code{FIXED}. (default = \code{FIXED})
#'   \item{\code{test_length}} test length.
#'   \item{\code{min_ni}} the maximum number of items to administer.
#'   \item{\code{max_ni}} the minimum number of items to administer.
#'   \item{\code{se_threshold}} standard error threshold. Item administration is stopped when theta estimate standard error becomes lower than this value.
#' }
#' @param interim_theta a named list containing interim theta estimation options.
#' \itemize{
#'   \item{\code{method}} the type of estimation. Accepts \code{EAP, MLE, MLEF, EB, FB, CARRYOVER}. (default = \code{EAP})
#'   \item{\code{shrinkage_correction}} set \code{TRUE} to apply shrinkage correction. Used when \code{method} is \code{EAP}. (default = \code{FALSE})
#'   \item{\code{prior_dist}} the type of prior distribution. Accepts \code{NORMAL, UNIFORM}. (default = \code{NORMAL})
#'   \item{\code{prior_par}} distribution parameters for \code{prior_dist}. (default = \code{c(0, 1)})
#'   \item{\code{bound_ML}} theta bound in \code{c(lower_bound, upper_bound)} format. Used when \code{method} is \code{MLE}. (default = \code{-4, 4})
#'   \item{\code{truncate_ML}} set \code{TRUE} to truncate ML estimate within \code{bound_ML}. (default = \code{FALSE})
#'   \item{\code{max_iter}} maximum number of Newton-Raphson iterations. Used when \code{method} is \code{MLE}. (default = \code{50})
#'   \item{\code{crit}} convergence criterion. Used when \code{method} is \code{MLE}. (default = \code{1e-03})
#'   \item{\code{max_change}} maximum change in ML estimates between iterations. Changes exceeding this value is clipped to this value. Used when \code{method} is \code{MLE}. (default = \code{1.0})
#'   \item{\code{use_step_size}} set \code{TRUE} to use \code{step_size}. Used when \code{method} is \code{MLE} or \code{MLEF}. (default = \code{FALSE})
#'   \item{\code{step_size}} upper bound to impose on the absolute change in initial theta and estimated theta. Absolute changes exceeding this value will be capped to \code{step_size}. Used when \code{method} is \code{MLE} or \code{MLEF}. (default = \code{0.5})
#'   \item{\code{do_Fisher}} set \code{TRUE} to use Fisher's method of scoring. Used when \code{method} is \code{MLE}. (default = \code{TRUE})
#'   \item{\code{fence_slope}} slope parameter to use for \code{method = 'MLEF'}. This must have two values in total, for the lower and upper bound item respectively. Use one value to use the same value for both bounds. (default = \code{5})
#'   \item{\code{fence_difficulty}} difficulty parameters to use for \code{method = 'MLEF'}. This must have two values in total, for the lower and upper bound item respectively. (default = \code{c(-5, 5)})
#'   \item{\code{hand_scored_attribute}} (optional) the item attribute name for whether each item is hand-scored or not. The attribute should have \code{TRUE} (hand-scored) and \code{FALSE} (machine-scored) values. If a hand-scored item is administered to an examinee, the previous interim theta (or the starting theta if this occurs for the first item) is reused without updating the estimate.
#' }
#' @param final_theta a named list containing final theta estimation options.
#' \itemize{
#'   \item{\code{method}} the type of estimation. Accepts \code{EAP, MLE, MLEF, EB, FB, CARRYOVER}. (default = \code{EAP})
#'   \item{\code{shrinkage_correction}} set \code{TRUE} to apply shrinkage correction. Used when \code{method} is \code{EAP}. (default = \code{FALSE})
#'   \item{\code{prior_dist}} the type of prior distribution. Accepts \code{NORMAL, UNIFORM}. (default = \code{NORMAL})
#'   \item{\code{prior_par}} distribution parameters for \code{prior_dist}. (default = \code{c(0, 1)})
#'   \item{\code{bound_ML}} theta bound in \code{c(lower_bound, upper_bound)} format. Used when \code{method} is \code{MLE}. (default = \code{-4, 4})
#'   \item{\code{truncate_ML}} set \code{TRUE} to truncate ML estimate within \code{bound_ML}. (default = \code{FALSE})
#'   \item{\code{max_iter}} maximum number of Newton-Raphson iterations. Used when \code{method} is \code{MLE}. (default = \code{50})
#'   \item{\code{crit}} convergence criterion. Used when \code{method} is \code{MLE}. (default = \code{1e-03})
#'   \item{\code{max_change}} maximum change in ML estimates between iterations. Changes exceeding this value is clipped to this value. Used when \code{method} is \code{MLE}. (default = \code{1.0})
#'   \item{\code{use_step_size}} set \code{TRUE} to use \code{step_size}. Used when \code{method} is \code{MLE} or \code{MLEF}.  (default = \code{FALSE})
#'   \item{\code{step_size}} upper bound to impose on the absolute change in initial theta and estimated theta. Absolute changes exceeding this value will be capped to \code{step_size}. Used when \code{method} is \code{MLE} or \code{MLEF}. (default = \code{0.5})
#'   \item{\code{do_Fisher}} set \code{TRUE} to use Fisher's method of scoring. Used when \code{method} is \code{MLE}. (default = \code{TRUE})
#'   \item{\code{fence_slope}} slope parameter to use for \code{method = 'MLEF'}. This must have two values in total, for the lower and upper bound item respectively. Use one value to use the same value for both bounds. (default = \code{5})
#'   \item{\code{fence_difficulty}} difficulty parameters to use for \code{method = 'MLEF'}. This must have two values in total, for the lower and upper bound item respectively. (default = \code{c(-5, 5)})
#' }
#' @param theta_grid the theta grid to use as quadrature points.
#'
#' @examples
#' cfg1 <- createShadowTestConfig(refresh_policy = list(
#'   method = "STIMULUS"
#' ))
#' cfg2 <- createShadowTestConfig(refresh_policy = list(
#'   method = "POSITION",
#'   position = c(1, 5, 9)
#' ))
#' @rdname createShadowTestConfig
#' @export
createShadowTestConfig <- function(
  item_selection = NULL, content_balancing = NULL, MIP = NULL, MCMC = NULL,
  exclude_policy = NULL, refresh_policy = NULL, exposure_control = NULL, overlap_control = NULL,
  stopping_criterion = NULL, interim_theta = NULL, final_theta = NULL, theta_grid = seq(-4, 4, .1)) {
  cfg <- new("config_Shadow")

  arg_names <- c(
    "item_selection", "content_balancing", "MIP", "MCMC",
    "exclude_policy", "refresh_policy", "exposure_control", "overlap_control",
    "stopping_criterion", "interim_theta", "final_theta"
  )
  obj_names <- c()
  for (arg in arg_names) {
    if (!is.null(eval(parse(text = arg)))) {
      accepted_slots <-
        eval(parse(text = sprintf("names(cfg@%s)", arg)))
      supplied_slots <-
        eval(parse(text = sprintf("names(%s)", arg)))
      leftovers <- setdiff(supplied_slots, accepted_slots)
      if (length(leftovers) > 0) {
        for (x in leftovers) {
          stop(
            sprintf(
              "cfg@%s: slot '%s' is unused",
              arg, x
            )
          )
        }
      }
      for (entry in accepted_slots) {
        entry_l <- paste0("cfg@", arg, "$", entry)
        entry_r <- paste0(arg, "$", entry)
        tmp <- eval(parse(text = entry_r))
        if (!is.null(tmp)) {
          eval(parse(text = paste0(entry_l, " <- ", entry_r)))
        }
      }
    }
  }

  # ensure the solver name is capitalized here;
  # the rest of the package assumes this is already done
  cfg@MIP$solver <- toupper(cfg@MIP$solver)

  if (!is.null(theta_grid)) {
    cfg@theta_grid <- theta_grid
  }
  if (length(cfg@exposure_control$max_exposure_rate) == 1) {
    cfg@exposure_control$max_exposure_rate <- rep(
      cfg@exposure_control$max_exposure_rate,
      cfg@exposure_control$n_segment
    )
  }
  v <- validObject(cfg)
  if (v) {
    return(cfg)
  }
}

#' Class 'output_Shadow_all': a set of adaptive assembly solutions
#'
#' \code{\linkS4class{output_Shadow_all}} is an S4 class for representing a set of adaptive assembly solutions.
#'
#' \describe{
#'   \item{\emph{notations}}{\itemize{
#'     \item{\emph{ni} denotes the number of items in the \code{\linkS4class{item_pool}} object.}
#'     \item{\emph{ns} denotes the number of stimuli.}
#'     \item{\emph{nj} denotes the number of participants.}
#'   }}
#' }
#'
#' @slot call the function call used for obtaining this object.
#' @slot output a length-*nj* list of \code{\linkS4class{output_Shadow}} objects, containing the assembly results for each participant.
#' @slot final_theta_est a length-*nj* vector containing final theta estimates for each participant.
#' @slot final_se_est a length-*nj* vector standard errors of the final theta estimates for each participant.
#' @slot exposure_rate a matrix containing item-level exposure rates of all items in the pool. Also contains stimulus-level exposure rates if the assembly was set-based.
#' @slot usage_matrix a *nj* by (*ni* + *ns*) matrix representing whether the item/stimulus was administered to each participant. Stimuli representations are appended to the right side of the matrix.
#' @slot cumulative_usage_matrix a *nj* by (*ni* + *ns*) matrix representing the number of times the item/stimulus was administered to each participant over multiple administrations.
#' @slot true_segment_count a length-*nj* vector containing the how many examinees are now in their segment based on the true theta. This will tend to increase. This can be reproduced with true theta values alone.
#' @slot est_segment_count a length-*nj* vector containing the how many examinees are now in their segment based on the estimated theta. This will tend to increase. This can be reproduced with estimated theta values alone.
#' @slot eligibility_stats exposure record for diagnostics.
#' @slot check_eligibility_stats detailed segment-wise exposure record for diagnostics. available when \code{config_Shadow@exposure_control$diagnostic_stats} is \code{TRUE}.
#' @slot no_fading_eligibility_stats detailed segment-wise exposure record without fading for diagnostics. available when \code{config_Shadow@exposure_control$diagnostic_stats} is \code{TRUE}.
#' @slot freq_infeasible a table representing the number of times the assembly was initially infeasible.
#' @slot pool the \code{\linkS4class{item_pool}} used in the assembly.
#' @slot config the \code{\linkS4class{config_Shadow}} used in the assembly.
#' @slot constraints the \code{\linkS4class{constraints}} used in the assembly.
#' @slot true_theta the \code{true_theta} argument used in the assembly.
#' @slot data the \code{data} argument used in the assembly.
#' @slot prior the \code{prior} argument used in the assembly.
#' @slot prior_par the \code{prior_par} argument used in the assembly.
#' @slot adaptivity a list of adaptivity indices.
#' @slot simulation_constants a list containing simulation constants parsed from input.
#'
#' @export
setClass("output_Shadow_all",
  slots = c(
    call                        = "call",
    output                      = "list_or_null",
    final_theta_est             = "matrix_or_numeric_or_null",
    final_se_est                = "matrix_or_numeric_or_null",
    exposure_rate               = "matrix_or_null",
    usage_matrix                = "matrix_or_null",
    cumulative_usage_matrix     = "matrix_or_null",
    true_segment_count          = "numeric_or_null",
    est_segment_count           = "numeric_or_null",
    eligibility_stats           = "list_or_null",
    check_eligibility_stats     = "list_or_null",
    no_fading_eligibility_stats = "list_or_null",
    freq_infeasible             = "table",
    pool                        = "item_pool",
    config                      = "config_Shadow",
    constraints                 = "constraints",
    data                        = "matrix_or_null",
    true_theta                  = "matrix_or_numeric_or_null",
    prior                       = "matrix_or_numeric_or_null",
    prior_par                   = "matrix_or_numeric_or_null",
    adaptivity                  = "list_or_null",
    simulation_constants        = "list_or_null"
  ),
  prototype = list(
    output                      = NULL,
    final_theta_est             = NULL,
    final_se_est                = NULL,
    exposure_rate               = NULL,
    usage_matrix                = NULL,
    cumulative_usage_matrix     = NULL,
    true_segment_count          = NULL,
    est_segment_count           = NULL,
    eligibility_stats           = NULL,
    check_eligibility_stats     = NULL,
    no_fading_eligibility_stats = NULL,
    freq_infeasible             = new("table"),
    pool                        = new("item_pool"),
    config                      = new("config_Shadow"),
    constraints                 = new("constraints"),
    data                        = NULL,
    true_theta                  = NULL,
    prior                       = NULL,
    prior_par                   = NULL,
    adaptivity                  = NULL,
    simulation_constants        = NULL
  ),
  validity = function(object) {
    return(TRUE)
  }
)

#' Class 'output_Shadow': adaptive assembly solution for one simulee
#'
#' \code{\linkS4class{output_Shadow}} is an S4 class for representing the adaptive assembly solution for one simulee.
#'
#' @slot simulee_id the numeric ID of the simulee.
#' @slot true_theta the true theta of the simulee, if was specified.
#' @slot true_theta_segment the segment number of the true theta.
#' @slot final_theta_est final theta estimate.
#' @slot final_se_est the standard error of \code{final_theta_est}.
#' @slot administered_item_index item IDs administered at each position.
#' @slot administered_item_resp item responses from the simulee at each position.
#' @slot administered_item_ncat the number of categories of each administered item.
#' @slot administered_stimulus_index stimulus IDs administered at each position.
#' @slot shadow_test_refreshed \code{TRUE} indicates the shadowtest was refreshed for the position.
#' @slot shadow_test_feasible \code{TRUE} indicates the MIP was feasible with all constraints.
#' @slot solve_time elapsed time in running the solver at each position.
#' @slot initial_theta_est initial theta estimate.
#' @slot interim_theta_est interim theta estimates at each position.
#' @slot interim_se_est the standard error of the interim estimate at each position.
#' @slot theta_segment_index segment numbers of interim theta estimates.
#' @slot prior prior distribution, if was specified.
#' @slot prior_par prior parameters, if were specified.
#' @slot posterior the posterior distribution after completing test.
#' @slot posterior_sample posterior samples of interim theta before the estimation of final theta. \code{mean(posterior_sample) == interim_theta_est[test_length]} holds.
#' @slot likelihood the likelihood distribution after completing test.
#' @slot shadow_test the list containing the item IDs within the shadowtest used in each position.
#' @slot max_cat_pool the maximum number of response categories the item pool had.
#' @slot ni_pool the total number of items the item pool had.
#' @slot ns_pool the total number of stimuli the item pool had.
#' @slot test_length_constraints the test length constraint used in assembly.
#' @slot set_based whether the item pool was set-based.
#' @slot item_index_by_stimulus the list of items by each stimulus the item pool had.
#'
#' @export
setClass("output_Shadow",
  slots = c(
    simulee_id                  = "numeric",
    true_theta                  = "matrix_or_numeric_or_null",
    true_theta_segment          = "numeric_or_null",
    final_theta_est             = "matrix_or_numeric",
    final_se_est                = "matrix_or_numeric",
    administered_item_index     = "numeric",
    administered_item_resp      = "numeric",
    administered_item_ncat      = "numeric",
    administered_stimulus_index = "numeric",
    shadow_test_refreshed       = "logical",
    shadow_test_feasible        = "logical",
    solve_time                  = "numeric",
    initial_theta_est           = "list_or_matrix_or_numeric",
    interim_theta_est           = "matrix_or_numeric",
    interim_se_est              = "matrix_or_numeric",
    theta_segment_index         = "numeric",
    prior                       = "numeric",
    prior_par                   = "numeric",
    posterior                   = "numeric",
    posterior_sample            = "numeric",
    likelihood                  = "numeric",
    shadow_test                 = "list",
    max_cat_pool                = "numeric",
    ni_pool                     = "numeric",
    ns_pool                     = "numeric",
    test_length_constraints     = "numeric",
    set_based                   = "logical",
    item_index_by_stimulus      = "list_or_null"
  ),
  prototype = list(
    simulee_id                  = numeric(0),
    true_theta                  = NULL,
    true_theta_segment          = NULL,
    final_theta_est             = numeric(0),
    final_se_est                = numeric(0),
    administered_item_index     = numeric(0),
    administered_item_resp      = numeric(0),
    administered_item_ncat      = numeric(0),
    administered_stimulus_index = numeric(0),
    shadow_test_refreshed       = logical(0),
    shadow_test_feasible        = logical(0),
    solve_time                  = numeric(0),
    initial_theta_est           = list(0),
    interim_theta_est           = numeric(0),
    interim_se_est              = numeric(0),
    theta_segment_index         = numeric(0),
    prior                       = numeric(0),
    prior_par                   = numeric(0),
    posterior                   = numeric(0),
    posterior_sample            = numeric(0),
    likelihood                  = numeric(0),
    shadow_test                 = list(0),
    max_cat_pool                = numeric(0),
    ni_pool                     = numeric(0),
    ns_pool                     = numeric(0),
    test_length_constraints     = numeric(0),
    set_based                   = logical(0),
    item_index_by_stimulus      = NULL
  ),
  validity = function(object) {
    return(TRUE)
  }
)

#' An S4 class for representing the exposure rate plot
#'
#' @noRd
setClass("exposure_rate_plot",
  slots = c(
    plot = "recordedplot_or_null",
    item_exposure_rate               = "numeric_or_null",
    item_exposure_rate_segment       = "list_or_null",
    item_exposure_rate_segment_final = "list_or_null",
    stim_exposure_rate               = "numeric_or_null",
    stim_exposure_rate_segment       = "list_or_null",
    stim_exposure_rate_segment_final = "list_or_null",
    segment_rate_table               = "dataframe_or_null",
    n_segment                        = "numeric_or_null",
    segment_n                        = "numeric_or_null",
    segment_cut                      = "numeric_or_null",
    segment_label                    = "character_or_null"
  ),
  prototype = list(
    plot = NULL,
    item_exposure_rate               = NULL,
    item_exposure_rate_segment       = NULL,
    item_exposure_rate_segment_final = NULL,
    stim_exposure_rate               = NULL,
    stim_exposure_rate_segment       = NULL,
    stim_exposure_rate_segment_final = NULL,
    segment_rate_table               = NULL,
    n_segment                        = NULL,
    segment_n                        = NULL,
    segment_cut                      = NULL,
    segment_label                    = NULL
  ),
  validity = function(object) {
    return(TRUE)
  }
)
