
# - `gridsearch` ------------------------------------------------------------- #
# - function: matrix smoothing ----------------------------------------------- #

### note: this function
###       (a) takes an n x m matrix as its input,
###       (b) averages each cell by using the values from its neighboring cells.

# FUNCTION CALL -------------------------------------------------------------- #

#' Smooth Matrix
#'
#' This function takes an `n x m` matrix as its input and averages each cell by using the values from its neighboring cells.
#'
#' @param mat An `n x m` matrix.
#' @param n The radius term for the neighbors.
#'
#' @return A matrix.
#' @noRd
#' @keywords internal
#'

smoothMatrix <- function(mat, n = 1) {

  # --- get dimensions
  nrow <- nrow(mat)
  ncol <- ncol(mat)

  # --- an empty matrix for results
  result <- matrix(NA, nrow = nrow, ncol = ncol)

  # --- loop through each cell
  for (i in 1:nrow) {
    for (j in 1:ncol) {

      # define the neighbor range, with edge cases handled
      row_range <- max(1, i - n):min(nrow, i + n)
      col_range <- max(1, j - n):min(ncol, j + n)

      # extract values of the cell and its neighbors
      neighborhood <- mat[row_range, col_range]

      # Calculate the average
      result[i, j] <- mean(neighborhood, na.rm = TRUE)
    }
  }
  return(result)
}

# ---------------------------------------------------------------------------- #
