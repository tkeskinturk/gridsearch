
# - `gridsearch` ------------------------------------------------------------- #
# - function: simulate changers ---------------------------------------------- #

### note: this function
###       (a) takes a set of parameters from a specific DGP and simulates data,
###       (b) calculates accuracy scores for detecting changers.

# FUNCTION CALL -------------------------------------------------------------- #

#' Simulate Changers
#'
#' This function gets several DGP parameters and calculates kappa, sensitivity, and specificity scores for classifying individual changers in simulated data frames.
#' @param n The number of units in the DGP.
#' @param t The number of time periods.
#' @param rate The percent of units changing across the panel period.
#' @param balance_dir The direction of change, where 0 codes negative change, 1 codes positive change, and all values in-between codes the percentage of changers changing in the positive direction.
#' @param balance_res The marginal distribution of the outcome (effectively indexing the 0-1 cut for the latent variable).
#' @param strength The strength of change in the latent variable.
#' @param reliable The reliability score of the outcome measurement.
#' @param nrep The number of simulation runs.
#' @param seed Set seed.
#'
#' @return A data frame.
#' @export
#'

simulateChangers <-

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
    ## the number of simulations and set seed
    nrep = 1000, seed = 11235

  ) {

    if(n <= 0)
      return(stop(
        "Error: We need at least some people. Check your `n` call."))
    if(t <= 1)
      return(stop(
        "Error: We need at least 2 time periods to generate panel data. Try again."))
    if(rate < 0)
      return(stop(
        "Error: Rate of change must range between 0 and 1."))
    if(balance_dir <  0 | balance_dir >  1)
      return(stop(
        "Error: Directionality balance must be a number between 0 and 1."))
    if(balance_res <= 0 | balance_res >= 1)
      return(stop(
        "Error: Outcome balance must be between 0 and 1 (not inclusive)."))
    if(strength < 0)
      return(stop(
        "Error: Change strength can't be negative (remember: it's 'strength')"))
    if(reliable < 0)
      return(stop(
        "Error: Reliability can't be negative. Try again."))

    future::plan(future::multisession)
    set.seed(seed)

    # ----------------------------------------------------- #
    # BUILDING BLOCKS                                       #
    # ----------------------------------------------------- #

    # --- parameter space

    data <-
      tidyr::expand_grid(
        p_n = n,
        p_t = t,
        p_rate = rate,
        p_dir = balance_dir,
        p_res = balance_res,
        p_strength = strength,
        p_rel = reliable,
        sims = 1:nrep
      )

    # --- simulate panel datasets from the parameters
    data <- data |>
      dplyr::mutate(
        data =
          furrr::future_pmap(
            ## mapping list
            .l = list(p_n, p_t, p_strength, p_rate, p_dir, p_res, p_rel),
            ## refer to list based on index
            .f = purrr::possibly(
              ~ gridsearch::buildDGP(
                # varying parameters
                n = ..1,
                t = ..2,
                strength = ..3,
                rate = ..4,
                balance_dir = ..5,
                balance_res = ..6,
                reliable = ..7,
                export = TRUE,
                patterns = FALSE,
                slopes = TRUE
              )
            ),
            ## for reproducibility
            .options = furrr::furrr_options(seed = TRUE),
            .progress = TRUE
          )
      )

    # ----------------------------------------------------- #
    # CALCULATE CONFUSION                                   #
    # ----------------------------------------------------- #

    data <- data |>
      tidyr::unnest_wider(col = data) |>

      # step 1: classify
      dplyr::mutate(
        correspondence = purrr::map2(
          .x = data,
          .y = slopes,
          .f =
            ~ dplyr::left_join(.x |>
                                 dplyr::distinct(pid, changer), .y, by = "pid") |>
            dplyr::mutate(
              predicted =
                dplyr::case_when(
                  (estimate > 0) == TRUE
                  &
                    (abs(estimate - 0) >= abs(estimate - strength)) == TRUE ~
                    1,
                  (estimate < 0) == TRUE
                  &
                    (abs(estimate - 0) >= abs(estimate + strength)) == TRUE  ~
                    1,
                  TRUE ~ 0
                )
            ) |>
            dplyr::select(pid, changer, predicted)
        )
      ) |>

      # step 2: confuse
      dplyr::mutate(confusion = purrr::map(
        .x = correspondence,
        .f = ~  tibble::tibble(changer =
                                 c(0, 0, 1, 1),
                               predicted =
                                 c(0, 1, 0, 1)) |>
          dplyr::left_join(
            . |>
              dplyr::count(changer, predicted),
            by = c("changer", "predicted")
          ) |>
          dplyr::mutate(n = ifelse(is.na(n) == TRUE, 0, n))
      )) |>

      # step 3: organize
      dplyr::mutate(confused =
                      purrr::map(.x = confusion, .f = ~ confusionMatrix(.))) |>
      dplyr::select(sims, confused) |> tidyr::unnest_wider(confused)

    return(data)

  }

# ---------------------------------------------------------------------------- #
