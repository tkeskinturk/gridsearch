
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
#' @param pname The unit identifier.
#' @param tname The time identifier.
#' @param yname The outcome identifier.
#'
#' @return A data frame.
#' @export
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
        }) |>
      dplyr::mutate(t = scales::rescale(t)) |>
      tidyr::drop_na()

    ## get unique model matrix
    X = stats::model.matrix( ~ t, data = df |>
                               dplyr::filter(pid == 1) |> dplyr::distinct(pid, t))

    ## function for probs
    bounded_change <- function(d) {
      fit <- fastglm::fastglmPure(
        X,
        d$y,
        family = stats::binomial(),
        method = 1)
      prob1 <- stats::plogis(fit$coefficients[1])
      prob2 <- stats::plogis(fit$coefficients[1] + fit$coefficients[2])
      return( round(prob2 - prob1, 3) )
    }

    # extract coefficient per pid
    estimates_coef <-
      suppressWarnings(
        data.table::as.data.table(df)[, estimate := bounded_change(.SD), by = pid]
        ) |>
      tibble::as_tibble() |>
      dplyr::distinct(pid, estimate)

    return(estimates_coef)

  }

# ---------------------------------------------------------------------------- #
