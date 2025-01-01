
# - `gridsearch` ------------------------------------------------------------- #
# - function: varying slopes ------------------------------------------------- #

### note: this function
###       (a) estimates bounded predictions at t = 0, t = 1 for trajectories,
###       (b) calculates model-predicted change scores.

# FUNCTION CALL -------------------------------------------------------------- #

#' Calculate Varying Slopes
#'
#' This function estimates trajectory-level logistic regression models by regressing binary outcomes on a normalized time variable.
#' @param t The number of time periods.
#' @noRd
#' @keywords internal
#' @return A data frame.
#'

varyingSlopes <-

  function(t) {

    ## possible trajectories
    tx <- seq(0, 1, length.out = t)
    patterns <- do.call(expand.grid, rep(list(c(0, 1)), t))

    ## an empty vector
    estimate <- c()

    ## loop over possible time trajectories
    for (i in 1:nrow(patterns)) {
      fit <- suppressWarnings(
        stats::glm(as.numeric(patterns[i, ]) ~ tx, family = stats::binomial()))
      prob1 <- stats::plogis(fit$coefficients[1])
      prob2 <- stats::plogis(fit$coefficients[1] + fit$coefficients[2])
      estimate[i] <- round(prob2 - prob1, 3)
    }

    ## generate a patterns dataframe
    patterns <- patterns |>
      cbind("estimate" = estimate)
    colnames(patterns) <- c(paste0("V", 1:t), "estimate")

    return(patterns)
  }

# ---------------------------------------------------------------------------- #
