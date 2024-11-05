
# - `gridsearch` ------------------------------------------------------------- #
# - function: build DGP ------------------------------------------------------ #

### note: this function
###       (a) generates panel data from specified parameters values,
###       (b) estimates user-specified models of change,
###       (c) stores (1) patterns, (2) model results OR (3) data, if specified.

# FUNCTION CALL -------------------------------------------------------------- #

#' Build DGP
#'
#' This function (1) builds simulated datasets using the parameter inputs in the function call, (2) generates panel patterns by concatenating the outcomes across time periods for each individual, and (3) calculates slope coefficients of change at the individual level.
#' @param n The number of units in the DGP.
#' @param t The number of time periods.
#' @param rate The percent of units changing across the panel period.
#' @param balance_dir The direction of change, where 0 codes negative change, 1 codes positive change, and all values in-between codes the percentage of changers changing in the positive direction.
#' @param balance_res The marginal distribution of the outcome (effectively indexing the 0-1 cut for the latent variable).
#' @param strength The strength of change in the latent variable.
#' @param reliable The reliability score of the outcome measurement.
#' @param export If `export = TRUE`, the function exports the simulated dataset.
#' @param patterns If `patterns = TRUE`, the function exports the concatenated panel patterns.
#' @param slopes If `slopes = TRUE`, the function exports slope coefficients.
#'
#' @return A list.
#' @export
#'

buildDGP <-

  function(

    # ----------------------------------------------------- #
    # FUNCTION ARGUMENTS AND WARNINGS                       #
    # ----------------------------------------------------- #

    ## the number of Ns
    n = 1000,
    ## the number of Ts
    t = 3,
    ## the rate of change in the population
    rate = 0.25,
    ## % balance in changes
    balance_dir = 0.5,
    ## % balance in `y` = 1
    balance_res = 0.5,
    ## strength of change (realized in latent scores)
    strength = 1,
    ## reliability for latent scores
    reliable = 0.8,
    ## return options for components
    export = TRUE, patterns = FALSE, slopes = FALSE) {

    if (n <= 0) stop("Error: We need at least some people. Check your `n` call.")
    if (t <= 1) stop("Error: We need at least 2 time periods to generate panel data. Try again.")
    if (rate < 0) stop("Error: Rate of change must range between 0 and 1.")
    if (balance_dir <  0 | balance_dir >  1) stop("Error: Directionality balance must be a number between 0 and 1.")
    if (balance_res <= 0 | balance_res >= 1) stop("Error: Outcome balance must be between 0 and 1 (not inclusive).")
    if (strength < 0) stop("Error: Change strength can't be negative (remember: it's 'strength')")
    if (reliable < 0) stop("Error: Reliability can't be negative. Try again.")
    if (export == FALSE & patterns == FALSE & slopes == FALSE) stop("Error: You said no data, no patterns, or no model. What do you want?")

    # ----------------------------------------------------- #
    # BUILDING BLOCKS                                       #
    # ----------------------------------------------------- #

    # --- part 1: generate actors
    u <-
      tibble::tibble(
        ## actors
        pid = c(1:n),
        ## central tendency
        u = stats::rnorm(n = n, mean = 0, sd = 1)
      )

    # --- part 2: generate the longitudinal data setup

    ## 2.1: build the window
    span_window <- seq(
      ## first t
      from = 0,
      ## final t
      to =   1,
      ## the number of ts
      length.out = t)

    ## 2.2: generate pid x time grid
    data <- tidyr::expand_grid(u, t = span_window)

    # --- part 3: assign change status

    ## 3.1: assign change status
    if (rate == 0) {
      data <- data |>
        dplyr::mutate(changer = 0) |> # placeholder for few clutter
        dplyr::mutate(timing = 0)
    } else {
      data <- data |>
        dplyr::left_join(tibble::tibble(
          pid = 1:n,
          changer = stats::rbinom(n = n, size = 1, prob = rate)
        ), by = "pid")
    }

    ## 3.2: assign change timing
    rate_window <- span_window[-1]
    if (rate != 0) {
      data <- data |>
        dplyr::left_join(
          data |>
            dplyr::filter(.data$changer == 1) |>
            dplyr::distinct(.data$pid, .data$changer) |>
            dplyr::select(-"changer") |>
            dplyr::mutate(timing = sample(
              rate_window, size = dplyr::n(), replace = TRUE
            )),
          by = "pid"
        ) |>
        dplyr::mutate(timing = ifelse(is.na(.data$timing) == TRUE, 0, .data$timing))
    }

    # ----------------------------------------------------- #
    # RESPONSE CONSTRUCTION                                 #
    # ----------------------------------------------------- #

    # --- part 4: generate the observed scores

    ## 4.1: add change scores
    data <- data |>
      dplyr::left_join(
        u |> dplyr::select("pid") |>
          dplyr::mutate(upper = stats::rbinom(
            n = n, size = 1, prob = balance_dir
          )) |>
          dplyr::mutate(upper = ifelse(.data$upper == 1, 1, -1)),
        by = "pid"
      ) |>
      dplyr::mutate(change = dplyr::case_when(
        .data$changer == 0 ~ 0,
        .data$changer == 1 & .data$t >= .data$timing ~ 1,
        .default = 0)
      ) |>
      dplyr::mutate(u = .data$u + strength * .data$change * .data$upper) |>
      dplyr::select(-"upper", -"change")

    ## 4.2: extract the new error variance
    new_error <- stats::sd(data$u)

    ## 4.3: generate the realized y scores
    data <- data |>
      dplyr::mutate(y =
                      ## true scores
                      (.data$u * sqrt(reliable))
                    +
                      ## error
                      stats::rnorm(n = n, mean = 0, sd = new_error)
                    *
                      sqrt((1 - reliable)))

    ## 4.4: binarize observed scores
    data <- data |>
      dplyr::mutate(y_obs =
                      ifelse(.data$y <= stats::quantile(.data$y, prob = 1 - balance_res), 0, 1))

    # --- part 5: organize and spit out if necessary

    ## 5.1: organize the data
    data <- data |>
      dplyr::select("pid", "changer", "timing", "t", "y_true" = "u", "y_obs")

    ## 5.2: spit out the data
    if (export   == TRUE   &
        patterns == FALSE  &
        slopes   == FALSE) {
      return(data)
    }

    # ----------------------------------------------------- #
    # PATTERNS                                              #
    # ----------------------------------------------------- #

    data_patterns <-
      ## theoretical patterns
      tibble::tibble(
        patterns =
          expand.grid(replicate(t, 0:1, simplify = F)) |>
          tibble::as_tibble() |>
          tidyr::unite(col = "position", sep = "") |>
          dplyr::pull("position") |>
          sort()
      ) |>
      ## empirical patterns
      dplyr::left_join(
        data |>
          dplyr::select("pid", "t", "y_obs") |>
          tidyr::pivot_wider(names_from = "t", values_from = "y_obs") |>
          janitor::clean_names() |>
          tidyr::unite("patterns", tidyr::starts_with("x"), sep = "") |>
          dplyr::summarize(sim_counts = dplyr::n(), .by = "patterns") |>
          dplyr::arrange(.data$patterns),
        by = "patterns"
      ) |>
      ## fill out the NA terms
      dplyr::mutate(sim_counts = ifelse(is.na(.data$sim_counts) == TRUE, 0, .data$sim_counts))

    if (export   == TRUE   &
        patterns == TRUE   &
        slopes   == FALSE) {
      return(list(data = data,
                  patterns = data_patterns))
    }

    if (export   == FALSE  &
        patterns == TRUE   &
        slopes   == FALSE) {
      return(data_patterns)
    }

    # ----------------------------------------------------- #
    # SLOPES                                                #
    # ----------------------------------------------------- #

    estimates_coef <- gridsearch::varyingSlopes(data,
                                                yname = "y_obs",
                                                tname = "t",
                                                pname = "pid")

    # --- part M3: export
    if (export   == FALSE  &
        patterns == FALSE  &
        slopes   == TRUE)  {
      return(estimates_coef)
    }

    if (export   == FALSE  &
        patterns == TRUE   &
        slopes   == TRUE)  {
      return(list(patterns = data_patterns,
                  slopes = estimates_coef))
    }

    if (export   == TRUE   &
        patterns == FALSE  &
        slopes   == TRUE)  {
      return(list(data = data,
                  slopes = estimates_coef))
    }

    if (export   == TRUE   &
        patterns == TRUE   &
        slopes   == TRUE)  {
      return(list(data = data,
                  patterns = data_patterns,
                  slopes = estimates_coef))
    }

  }

# ---------------------------------------------------------------------------- #
