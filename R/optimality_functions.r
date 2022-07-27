#' @include static_functions.R
NULL

#' @noRd
selectItemUsingDOptimality <- function(info, position, o, constants) {

  if (position == 1) {
    info_matrix_of_administered_items <- matrix(0, constants$nd, constants$nd)
  }

  if (position > 1) {
    info_matrix_of_administered_items <- info[o@administered_item_index[1:(position - 1)]]
    info_matrix_of_administered_items <- Reduce("+", info_matrix_of_administered_items)
  }

  determinant_of_candidate_items <-
    lapply(1:constants$ni, function(x) {
      if (x %in% o@administered_item_index[1:(position - 1)]) {
        return(-Inf)
      }
      info_matrix_of_candidate_test <- info_matrix_of_administered_items + info[[x]]
      d <- det(info_matrix_of_candidate_test)
      return(d)
  })

  determinant_of_candidate_items <- unlist(determinant_of_candidate_items)
  selected_item <- which.max(determinant_of_candidate_items)

  return(selected_item)

}
