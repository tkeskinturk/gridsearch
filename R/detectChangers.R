
# - `gridsearch` ------------------------------------------------------------- #
# - function: detectChangers ------------------------------------------------- #

### note: this function
###       (a) receives parameter arguments for generating a simulation data,
###       (b) calculates the extent to which we recover "changers."

# FUNCTION CALL -------------------------------------------------------------- #

#' Detect the Accuracy of Change Classification
#'
#' This function builds a simulated dataset and classifies individuals based on their change scores to see whether we can recover the changers from the true DGP.
#' @param n The number of individuals in the DGP.
#' @param t The number of time periods.
#' @param rate The percent of individuals changing across the panel.
#' @param balance_dir The direction of change, where 0 means all negative change, 1 means all positive change, and values in-between means the percentage of changers changing positively.
#' @param balance_res The marginal distribution of the outcome (the percent distribution of 0s and 1).
#' @param strength The strength of change in the latent variable.
#' @param reliable The reliability score of the outcome measurement.
#' @param n_samples The number of simulated datasets.
#' @param verbose See progress messages in the console.
#' @return A dataframe or a list.
#'
#' @export
#'

detectChangers <-

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
    n_samples = 1000,
    verbose = TRUE) {

    if (n <= 0) stop("Error: We need at least some people. Check your `n` call.")
    if (t <= 1) stop("Error: We need at least 2 time periods to generate panel data. Try again.")
    if (rate < 0 | rate > 1) stop("Error: Rate of change must range between 0 and 1.")
    if (balance_dir <  0 | balance_dir >  1) stop("Error: Directionality balance must range between 0 and 1.")
    if (balance_res <= 0 | balance_res >= 1) stop("Error: Outcome balance must be between 0 and 1 (not inclusive).")
    if (strength < 0) stop("Error: Change strength can't be negative (remember: it's 'strength')")
    if (reliable < 0) stop("Error: Reliability can't be negative. Try again.")

    # ----------------------------------------------------- #
    # GENERATE DATA                                         #
    # ----------------------------------------------------- #

    d <- buildDGP(
      n = n,
      t = t,
      rate = rate,
      balance_dir = balance_dir,
      balance_res = balance_res,
      strength = strength,
      reliable = reliable,
      export = TRUE,
      status = TRUE
    )$data

    # ----------------------------------------------------- #
    # DATASETS                                              #
    # ----------------------------------------------------- #

    sim <- vector()
    ref <- varyingSlopes(t)
    j <- paste0("V", 1:t)

    for (i in 1:n_samples) {

      # --- generate a data
      data <- buildDGP(
        n = n,
        t = t,
        rate = rate,
        balance_dir = balance_dir,
        balance_res = balance_res,
        strength = strength,
        reliable = reliable,
        export = TRUE,
        status = TRUE
      )$data

      # --- add slopes and order
      data <- collapse::join(data,
                             ref,
                             on = j,
                             how = "left",
                             verbose = FALSE)
      data <- data[order(-data$estimate), ]

      # --- assignment

      ## changer counts
      changer_n <- n * rate
      changer_u <- changer_n * balance_dir
      changer_d <- changer_n * (1 - balance_dir)

      ## assign values
      if (balance_dir != 0 & balance_dir != 1) {
        data$class <- 0
        data$class[1:changer_u] <- 1
        data$class[(nrow(data) - changer_d + 1):nrow(data)] <- -1
      }

      if (balance_dir == 1) {
        data$class <- 0
        data$class[1:changer_u] <- 1
      }

      if (balance_dir == 0) {
        data$class <- 0
        data$class[(nrow(data) - changer_d + 1):nrow(data)] <- -1
      }

      # --- classification
      list1 <- c(-1, 0, 1)
      list2 <- c(-1, 0, 1)
      table <- expand.grid(changer = list1, class = list2)
      counts <- as.data.frame(table(data$changer, data$class))
      colnames(counts) <- c("changer", "class", "Freq")
      table <- merge(table, counts, by = c("changer", "class"), all.x = TRUE)
      table$Freq[is.na(table$Freq)] <- 0

      # --- calculate accuracy
      acc1 <-
        table$Freq[table$changer == -1 & table$class == -1] /
        sum(table$Freq[table$changer == -1])
      acc2 <-
        table$Freq[table$changer == 1 & table$class == 1] /
        sum(table$Freq[table$changer == 1])

      acc <- mean(c(acc1, acc2), na.rm = TRUE)

      sim[[i]] <- acc

      if (i %% 100 == 0 & verbose == TRUE) {
        print(paste0(round(i / n_samples * 100, 1), "%"))
      }

      rm(data)

    }

    return(sim)


  }

# ---------------------------------------------------------------------------- #
