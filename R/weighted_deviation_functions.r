#' @include static_functions.R
NULL

#' @noRd
selectItemUsingWeightedDeviation <- function(info, position, o, constraints) {

  # runs too slow, needs optimization
  # maybe generate the attribute contribution vector for each item in advance

  ni <- constraints@ni

  if (position == 1) {
    administered_test <- c()
    administered_test_attributes <- getSolutionAttributes(constraints, administered_test, FALSE)
  }
  if (position > 1) {
    administered_test <- o@administered_item_index[1:(position - 1)]
    administered_test_attributes <- getSolutionAttributes(constraints, administered_test, FALSE)
  }

  wd <- lapply(
    1:ni,
    function(x) {

      if (position == 1) {
        candidate_test <- x
      }
      if (position > 1) {
        if (x %in% administered_test) {
          return(Inf)
        }
        candidate_test <- c(administered_test, x)
      }

      candidate_item_attributes <- getSolutionAttributes(constraints, x, FALSE)
      candidate_test_attributes <-
        administered_test_attributes$solution +
        candidate_item_attributes$solution

      d_lowerbound_violation <-
        (administered_test_attributes$LB - candidate_test_attributes) *
        (administered_test_attributes$LB > candidate_test_attributes)
      d_upperbound_violation <-
        (candidate_test_attributes - administered_test_attributes$UB) *
        (candidate_test_attributes > administered_test_attributes$UB)
      candidate_test_info <- sum(info[candidate_test, ])
      d_info_violation <- 100000 - candidate_test_info

      wd <-
        sum(d_lowerbound_violation, na.rm = TRUE) +
        sum(d_upperbound_violation, na.rm = TRUE) +
        d_info_violation

      # TODO: add weights

      return(wd)
      
    }
  )

  wd <- unlist(wd)
  selected_item <- which.min(wd)
  return(selected_item)

}
