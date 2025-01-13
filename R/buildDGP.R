
# - `gridsearch` ------------------------------------------------------------- #
# - function: build DGP ------------------------------------------------------ #

### note: this function
###       (a) generates a panel data from the specified parameter values,
###       (b) tabulates the response patterns across time.

# FUNCTION CALL -------------------------------------------------------------- #

#' Build Panel Data from a DGP
#'
#' This function builds a simulated dataset and tabulates the simulated trajectories using the parameters supplied in the function call.
#' @param n The number of individuals in the DGP.
#' @param t The number of time periods.
#' @param rate The percent of individuals changing across the panel.
#' @param balance_dir The direction of change, where 0 codes negative change, 1 codes positive change, and all values in-between codes the percentage of changers changing in the positive direction.
#' @param balance_res The marginal distribution of the outcome (effectively regulating the percent distribution of 0s and 1s).
#' @param strength The strength of change in the latent variable.
#' @param reliable The reliability score of the outcome measurement.
#' @param export If `export = TRUE`, the function exports the simulated dataset.
#' @param status if `status = TRUE`, the function exports a column for changer status.
#'
#' @return A dataframe or a list.
#'
#' @export
#'

buildDGP <-

  function(

    # ----------------------------------------------------- #
    # FUNCTION ARGUMENTS AND WARNINGS                       #
    # ----------------------------------------------------- #

    n = 1000,
    t = 3,
    rate = .5,
    balance_dir = .5,
    balance_res = .5,
    strength = .5,
    reliable = 1,
    export = FALSE,
    status = FALSE) {

    if (n <= 0) stop("Error: We need at least some people. Check your `n` call.")
    if (t <= 1) stop("Error: We need at least 2 time periods to generate panel data. Try again.")
    if (rate < 0 | rate > 1) stop("Error: Rate of change must range between 0 and 1.")
    if (balance_dir <  0 | balance_dir >  1) stop("Error: Directionality balance must be a number between 0 and 1.")
    if (balance_res <= 0 | balance_res >= 1) stop("Error: Outcome balance must be between 0 and 1 (not inclusive).")
    if (strength < 0) stop("Error: Change strength can't be negative (remember: it's 'strength')")
    if (reliable < 0) stop("Error: Reliability can't be negative. Try again.")

    # ----------------------------------------------------- #
    # GENERATE DATA                                         #
    # ----------------------------------------------------- #

    ## latent variable
    mlv <- matrix(nrow = n, ncol = t) # latent variable
    mt <- matrix(nrow = n, ncol = t) # store survey answer at time t
    mlv[,] <- stats::rnorm(n) # generate latent variable

    ## add changer status
    changers <- # actors who change
      sample(1:n, size = rate * n, replace = FALSE)
    change_wave <- # in which wave
      sample(2:t, size = rate * n, replace = TRUE)
    change_dir <- # direction of change
      sample(
        c(-1, 1),
        size = rate * n,
        replace = TRUE,
        prob = c((1 - balance_dir), balance_dir)
      )
    change_value <- # amount of change
      change_dir * strength
    j <- paste0("V", 1:t) # columns to join on

    ## iterate over time periods
    for (i in 2:t) {
      mlv[changers, i] <-
        mlv[changers, i - 1] + ifelse(change_wave == i, change_value, 0)
    }

    ## get the new error
    new_error <- stats::sd(mlv)

    ## reliability
    mlv <- mlv * sqrt(reliable) + stats::rnorm(n * t, 0, new_error) * sqrt(1 - reliable)

    ## binarize the latent variable
    thres <- stats::quantile(mlv, prob = 1 - balance_res)
    mt <- ifelse(mlv <= thres, 0, 1)
    mt <- data.table::data.table(mt) # for data.table
    data <- mt # for listing later
    mt <- mt[, .N, by=j]
    mt[, "N" := .SD[["N"]] / sum(.SD[["N"]])]
    mt <- as.data.frame(mt)

    ## generate the changer column
    changer_col <- rep(0, n)
    if (status) {
      changer_col[changers] <- change_dir # assign change direction
      data$changer <- changer_col # add changer column
    }

    if (export == TRUE) {

      return(
        list(patterns = mt, data = data)
      )

    }

    else {

      return(
        mt
      )

    }

  }

# ---------------------------------------------------------------------------- #
