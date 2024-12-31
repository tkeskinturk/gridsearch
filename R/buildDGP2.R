
# - `gridsearch` ------------------------------------------------------------- #
# - function: build DGP ------------------------------------------------------ #

### note: this function
###       (a) generates panel data from specified parameters values,
###       (b) estimates user-specified models of change,
###       (c) stores (1) patterns, (2) model results OR (3) data, if specified.

# FUNCTION CALL -------------------------------------------------------------- #

#' Build Panel Data from a DGP
#'
#' This function (1) builds a simulated dataset using the parameters supplied in the function call, (2) generates distinct outcome patterns by concatenating the outcomes across time periods for each individual and provides counts for each pattern, (3) calculates slope coefficients of change at the individual level.
#' @param n The number of units in the DGP.
#' @param t The number of time periods.
#' @param rate The percent of units changing across the panel.
#' @param balance_dir The direction of change, where 0 codes negative change, 1 codes positive change, and all values in-between codes the percentage of changers changing in the positive direction.
#' @param balance_res The marginal distribution of the outcome (effectively regulating the percent distribution of 0s and 1s).
#' @param strength The strength of change in the latent variable.
#' @param reliable The reliability score of the outcome measurement.
#' @param export If `export = TRUE`, the function exports the simulated dataset.
#' @param patterns If `patterns = TRUE`, the function exports the concatenated panel patterns.
#' @param slopes If `slopes = TRUE`, the function exports slope coefficients.
#'
#' @return A dataframe or a list.
#'
#' @details # Examples

#' @export
#'

buildDGP2 <-

  function(

    # ----------------------------------------------------- #
    # FUNCTION ARGUMENTS AND WARNINGS                       #
    # ----------------------------------------------------- #

    ## the number of Ns
    n = 1000,
    ## the number of Ts
    t = 3,
    ## the rate of change in the population
    rate = .5,
    ## % balance in changes
    balance_dir = .5,
    ## % balance in `y` = 1
    balance_res = 0.5,
    ## strength of change (realized in latent scores)
    strength = .5,
    ## reliability for latent scores
    reliable = 1) {

    if (n <= 0) stop("Error: We need at least some people. Check your `n` call.")
    if (t <= 1) stop("Error: We need at least 2 time periods to generate panel data. Try again.")
    if (rate < 0) stop("Error: Rate of change must range between 0 and 1.")
    if (balance_dir <  0 | balance_dir >  1) stop("Error: Directionality balance must be a number between 0 and 1.")
    if (balance_res <= 0 | balance_res >= 1) stop("Error: Outcome balance must be between 0 and 1 (not inclusive).")
    if (strength < 0) stop("Error: Change strength can't be negative (remember: it's 'strength')")
    if (reliable < 0) stop("Error: Reliability can't be negative. Try again.")

    # ----------------------------------------------------- #
    # BUILDING BLOCKS                                       #
    # ----------------------------------------------------- #

    mlv <- matrix(nrow = n, ncol = t) # latent variable
    mt <- matrix(nrow = n, ncol = t) # store survey answer at time t
    mlv[,] <- rnorm(n) # generate latent variable


    changers <- sample(1:n, size = rate*n, replace = FALSE) # actors who change
    change_wave <- sample(2:t, size = rate*n, replace = TRUE) # in which wave
    change_dir <- sample(c(-1,1), size = rate*n, replace = TRUE,
                         prob = c((1 - balance_dir), balance_dir)) # direction of change
    change_value <- change_dir*strength # amount of change
    j <- paste0("V",1:t) # columns to join on


    # Iterate over time periods
    for(i in 2:t) {

      mlv[changers, i] <- mlv[changers, i - 1] + ifelse(change_wave == i, change_value, 0)

    }

    new_error <- sd(mlv)

    # Reliability
    mlv <- mlv * sqrt(reliable) + rnorm(n*t, 0, new_error) * sqrt(1 - reliable)

    thres <- quantile(mlv, prob = 1 - balance_res)
    mt <- ifelse(mlv <= thres, 0, 1)
    mt <- data.table::data.table(mt) # for data.table
    mt <- mt[, .N, by=j]
    mt[,N := N / sum(N)]
    mt <- as.data.frame(mt)



    ## spit out the data
      return(mt)
    }

# ---------------------------------------------------------------------------- #
