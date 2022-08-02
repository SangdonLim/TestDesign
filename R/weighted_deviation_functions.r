#' @include static_functions.R
NULL

#' @noRd
selectItemUsingWeightedDeviation <- function(info, position, o, attribute_contribution, constraints) {

  ni <- constraints@ni

  if (position == 1) {
    administered_test <- c()
    administered_test_attributes <- attribute_contribution[[1]] * 0
  }
  if (position > 1) {
    administered_test <- o@administered_item_index[1:(position - 1)]
    administered_test_attributes <- do.call(rbind, attribute_contribution[administered_test])
    administered_test_attributes <- apply(administered_test_attributes, 2, sum)
  }

  wd <- lapply(
    1:ni,
    function(x) {

      if (x %in% administered_test) {
        return(Inf)
      }

      candidate_test <- c(administered_test, x)

      candidate_item_attributes <- attribute_contribution[[x]]
      candidate_test_attributes <-
        administered_test_attributes +
        candidate_item_attributes

      d_lowerbound_violation <-
        (constraints@constraints$LB - candidate_test_attributes) *
        (constraints@constraints$LB > candidate_test_attributes)
      d_upperbound_violation <-
        (candidate_test_attributes - constraints@constraints$UB) *
        (candidate_test_attributes > constraints@constraints$UB)
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
