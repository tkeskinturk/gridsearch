
# - `gridsearch` ------------------------------------------------------------- #
# - function: varying slopes ------------------------------------------------- #

### note: this function
###       (a) estimates bounded predictions at t = 0 and t = 1 at person level,
###       (b) calculates model-predicted change scores.

# FUNCTION CALL -------------------------------------------------------------- #

#' Calculate Varying Slopes
#'
#' This function estimates individual-level logistic regression models by regressing binary outcomes on a normalized time variable.
#' @param data A dataframe.
#' @keywords internal
#' @param pname The unit identifier.
#' @param tname The time identifier.
#' @param yname The outcome identifier.
#'
#' @return A data frame.
#' @noRd
#' @keywords internal
#'

varyingSlopes <-

  function(data,
           pname,
           tname,
           yname) {

    ## data file
    df <- data |>
      dplyr::select(# get person identifier
        pid = {
          {
            pname
          }
        }, # get time identifier
        t = {
          {
            tname
          }
        }, # get outcome identifier
        y  = {
          {
            yname
          }
        })

    ## get unique `t`
    if (df |>
      dplyr::summarize(pattern = list(.data$t), .by = 'pid') |>
      dplyr::count(.data$pattern) |>
      nrow() != 1) stop(
        "Error: The dataset is not balanced.")
    X <- stats::model.matrix(~ t,
                            data = df |>
                              dplyr::filter(.data$pid == 1))

    ## get unique `y`
    Y <- df |>
      dplyr::summarize(pattern = list(.data$y), .by = 'pid') |>
      dplyr::count(.data$pattern)

    ## function for probs
    bounded_change <- function(X, Y) {
      fit <- fastglm::fastglmPure(
        X,
        Y,
        family = stats::binomial(),
        method = 1)
      prob1 <- stats::plogis(fit$coefficients[1])
      prob2 <- stats::plogis(fit$coefficients[1] + fit$coefficients[2])
      return( round(prob2 - prob1, 3) )
    }

    ## pattern fits
    Y <-
      suppressWarnings(Y |>
                         dplyr::mutate(estimate = purrr::map_dbl(
                           .x = .data$pattern, .f = ~ bounded_change(X = X, Y = .)
                         )) |>
                         dplyr::select(-.data$n))

    ## match by patterns
    estimates_coef <- df |>
      dplyr::summarize(pattern = list(.data$y), .by = 'pid') |>
      dplyr::left_join(Y, by = "pattern") |>
      dplyr::select(.data$pid, .data$estimate)

    return(estimates_coef)

  }

# ---------------------------------------------------------------------------- #
