
# - `gridsearch` ------------------------------------------------------------- #
# - function: grid search ---------------------------------------------------- #

### note: this function
###       (a) takes a panel dataframe and processes it,
###       (b) performs an ABC grid search via tabulated response patterns/slopes,
###       (c) returns accepted samples for IC, PC, BAL & REL from the search.

# FUNCTION CALL -------------------------------------------------------------- #

#' Grid Search for Adjudicating DGPs
#'
#' This function uses an Approximate Bayesian Computation (ABC) algorithm to calculate the extent to which a dataframe can be approximated by a known DGP. It implements simulations across plausible values---rate of change, strength of change, directionality, and reliability---and provides a list of accepted samples from this procedure.
#' @param data A panel dataframe in the long format.
#' @param yname The outcome identifier.
#' @param tname The time identifier.
#' @param pname The unit identifier.
#' @param n_samples The number of samples from priors.
#' @param ic_min Minimum strength of change parameter.
#' @param ic_max Maximum strength of change parameter.
#' @param pc_min Minimum rate of change parameter.
#' @param pc_max Maximum rate of change parameter.
#' @param bal_min Minimum balance parameter.
#' @param bal_max Maximum balance parameter.
#' @param rel_min Minimum reliability parameter.
#' @param rel_max Maximum reliability parameter.
#' @param fix Fix any parameter? (values: "ic_sample," "pc_sample," etc.).
#' @param fix_at At what value to fix it on?
#' @param verbose Whether to see detailed messages.
#' @return A data frame.
#'
#' @export
#'

gridSearch <-

  function(

    # ----------------------------------------------------- #
    # FUNCTION ARGUMENTS AND WARNINGS                       #
    # ----------------------------------------------------- #

    data,
    yname,
    tname,
    pname,
    n_samples = 1000,
    ic_min = 0,
    ic_max = 2,
    pc_min = 0,
    pc_max = 1,
    bal_min = 0,
    bal_max = 1,
    rel_min = 0,
    rel_max = 1,
    fix = "none",
    fix_at = 1,
    verbose = TRUE
    ) {

    # ----------------------------------------------------- #
    # DATA PROCESSING                                       #
    # ----------------------------------------------------- #

    data <- data.table::as.data.table(data)

    # use the column names provided in the function arguments
    data$y <- data[[yname]] # access the column for `yname`
    data$t <- data[[tname]] # access the column for `tname`
    data$p <- data[[pname]] # access the column for `pname`

    # subset the relevant columns
    data <- data[, .SD, .SDcols = c("y", "t", "p")]

    # check if the dataframe supplied is balanced
    data <- stats::na.omit(data) ## missing
    is_balanced <- all(table(data$p) == (table(data$p))[1])
    if (!is_balanced) {
      stop("Error: The panel data is not balanced.")
    }
    b <- mean(data$y) ## the distribution of observed ys
    p <- table(data$t)[[1]] ## the number of unique individuals
    t <- table(data$p)[[1]] ## the number of time periods

    # reshape the data from long to wide format
    data <- data.table::dcast(data, p ~ t, value.var = "y")

    # rename the columns for easy handling
    j <- paste0("V", 1:t) ## columns to join on
    data.table::setnames(data, old = as.character(1:t), new = paste0("V", 1:t))
    data <- data[, .N, by=j]
    data[, "N" := .SD[["N"]] / sum(.SD[["N"]])]

    # ----------------------------------------------------- #
    # PREPARATIONS                                          #
    # ----------------------------------------------------- #

    ## preps
    samples <- list()

    ## slopes
      ref <- varyingSlopes(t)
      data <- collapse::join(data,
                             ref,
                             on = j,
                             how = "left",
                             verbose = FALSE)
      data <- stats::aggregate(N ~ estimate, data, sum)

    # ----------------------------------------------------- #
    # ABC                                                   #
    # ----------------------------------------------------- #

    for (i in 1:n_samples) {

      ### --- sample from prior
      ic_sample <- stats::runif(n = 1,
                                min = ic_min,
                                max = ic_max)
      pc_sample <- stats::runif(n = 1,
                                min = pc_min,
                                max = pc_max)
      bal_sample <- stats::runif(n = 1,
                                 min = bal_min,
                                 max = bal_max)
      rel_sample <- stats::runif(n = 1,
                                 min = rel_min,
                                 max = rel_max)
      assign(fix, fix_at)

      ### --- simulate the data
      res <- buildDGP(
        n = p,
        t = t,
        rate = pc_sample,
        balance_dir = bal_sample,
        balance_res = b,
        strength = ic_sample,
        reliable = rel_sample
      )

      ### --- join slopes
      res <- collapse::join(res,
                            ref,
                            on = j,
                            how = "left",
                            verbose = FALSE)
      res <- stats::aggregate(N ~ estimate, res, sum)

      ### --- ks test
      ks <- suppressWarnings(stats::ks.test(
          x = rep(res$estimate, res$N * 1000),
          y = rep(data$estimate, data$N * 1000)
        ))
        error <- ks$statistic


      ### --- accept/reject
        samples[[i]] <-
          data.table::data.table(ic_sample, pc_sample, bal_sample, rel_sample, error)


      if (i %% 1000 == 0 && verbose == TRUE) {
        print(paste0(round(i / n_samples * 100, 1), "%"))
      }

    }

    samples <- data.table::rbindlist(samples)

    return(samples)

  }

# ---------------------------------------------------------------------------- #
