
# - `gridsearch` ------------------------------------------------------------- #
# - function: grid search ---------------------------------------------------- #

### note: this function
###       (a) takes a panel dataframe (with `pid`, `t`, and `y`),
###       (b) iteratively prune to recover best DGPs,
###       (c) calculates error by measuring pattern distance from the DGPs.

# FUNCTION CALL -------------------------------------------------------------- #

#' Grid Search for Adjudicating DGPs
#'
#' This function uses a three-step grid search algorithm to calculate the extent to which a dataframe can be approximated by a known DGP. It implements fake-data simulations from a grid of plausible values---rate of change ranging from 0% to 100%, strength of change ranging from 0 SD to 2 SD, and 5 directional cases (everyone changing negatively, everyone changing positively, half changing negatively and half changing positively, 25% changing negatively and 75% changing positively, and 75% changing negatively and 25% changing positively), calculates the distribution of response patterns or slopes, and provides an error term that summarizes the distance of the DGP from the observed values.
#' @param data A dataframe.
#' @param yname The outcome identifier.
#' @param tname The time identifier.
#' @param pname The unit identifier.
#' @param pattern If `pattern = contingency`, the calculations are based on the use of a contingency table of panel patterns; if `pattern = slopes`, function uses individual slope coefficients.
#' @param steps The number of simulated datasets.
#' @param reliability The assumed reliability score of the outcome measurement.
#' @param seed Seed for reproducibility.
#' @param workers The number of workers for parallelization (note that parallelization is highly recommended for reasonable duration).
#'
#' @return A data frame.
#'
#' @details # Examples
#'
#' **Basic [gridSearch()] call:**
#' ```{r, comment = "#>", collapse = TRUE}
#' set.seed(1871)
#' data <- buildDGP(n = 50,
#'                  t = 3,
#'                  rate = 0.5,
#'                  balance_dir = 1,
#'                  balance_res = 0.5,
#'                  strength = 1,
#'                  reliable = 0.9,
#'                  export = TRUE,
#'                  patterns = FALSE,
#'                  slopes = FALSE)
#' gridSearch(data = data,
#'            yname = 'y_obs',
#'            tname = 't',
#'            pname = 'pid',
#'            pattern = 'contingency',
#'            steps = 1,
#'            reliability = 0.9,
#'            seed = 11234,
#'            workers = 10)
#' ```
#'
#' @export
#'

gridSearch <-

  function(

    # ----------------------------------------------------- #
    # FUNCTION ARGUMENTS AND WARNINGS                       #
    # ----------------------------------------------------- #

    # dataframe
    data,
    # outcome variable
    yname,
    # time variable
    tname,
    # person identifier
    pname,
    # pattern match
    pattern = c("contingency", "slopes"),
    # simulation runs
    steps = 1,
    # reliability
    reliability = 0.9,
    # seed
    seed = 11235,
    # the number of workers
    workers = 8

  ) {

    # --- dataframe

    df <- data |>
      dplyr::select(
        # get person identifier
        pid = dplyr::all_of(pname),
        # get time identifier
        t = dplyr::all_of(tname),
        # get outcome identifier
        y  = dplyr::all_of(yname)) |>
      dplyr::mutate(t = scales::rescale(.data$t)) |>
      tidyr::drop_na()

    # --- messages!

    message("This function drops missing values. Tread carefully.")

    if(
      sort(unique(df$y))[1] == 0 &
      sort(unique(df$y))[2] == 1 &
      length(unique(df$y)) == 2) {
    } else{
      stop("Error: The outcome needs to be binary and coded as 0 and 1.")
    }
    if(length(unique(df$t)) <= 1) {
      stop("Error: There needs to be at least two time periods.")
    }

    pattern = match.arg(pattern)

    # ----------------------------------------------------- #
    # BUILDING BLOCKS                                       #
    # ----------------------------------------------------- #

    # --- part 1: parallelize

    set.seed(seed)
    future::plan(future::multisession, workers = workers)

    # --- part 2: prepare the reference

    if (pattern == "contingency") {

      rf <- df |>
        tidyr::pivot_wider(names_from = "t", values_from = "y") |>
        tidyr::unite("patterns", -tidyr::starts_with("pid"), sep = "") |>
        dplyr::summarize(reference = dplyr::n(), .by = "patterns") |>
        dplyr::arrange(.data$patterns) ## contingency of 0s and 1s across waves

    }

    if (pattern == "slopes") {

      rf <- varyingSlopes(df,
                          pname = "pid",
                          tname = "t",
                          yname = "y")

    }

    # --- part 3: fake-data grid

    fd <- tidyr::expand_grid(
      p_n =
        dplyr::distinct(df, .data$pid) |> nrow(),
      p_t =
        dplyr::distinct(df, .data$t)   |> nrow(),
      p_rate =
        round(seq(
          from = 0,
          to   = 1,
          length.out = 41
        ), 3),
      p_strength =
        round(seq(
          from = 0.05,
          to   = 2,
          length.out = 40
        ), 3),
      p_dir = c(0, 0.25, 0.50, 0.75, 1),
      p_res = mean(df$y),
      sim = 1:steps
    )

    # --- part 4: function for simulations

    parallel_grid <- function(dat) {

      if (pattern == "contingency") {

        dat <- dat |>
          dplyr::mutate(
            data =
              furrr::future_pmap(
                ## mapping list
                .l = list(.data$p_n,
                          .data$p_t,
                          .data$p_strength,
                          .data$p_rate,
                          .data$p_dir,
                          .data$p_res),
                ## refer to list based on index
                .f = ~ buildDGP(
                  # varying parameters
                  n = ..1,
                  t = ..2,
                  strength = ..3,
                  rate = ..4,
                  balance_dir = ..5,
                  balance_res = ..6,
                  reliable = reliability,
                  export = FALSE,
                  patterns = TRUE,
                  slopes = FALSE
                ),
                ## for reproducibility
                .options = furrr::furrr_options(seed = TRUE),
                .progress = TRUE
              )
          ) |>
          ## get the distance
          tidyr::unnest(data) |>
          dplyr::left_join(rf, by = "patterns") |>
          dplyr::mutate(reference = ifelse(is.na(.data$reference) == TRUE,
                                           0,
                                           .data$reference)) |>
          dplyr::mutate(deviation = abs(.data$n - .data$reference)) |>
          dplyr::summarize(
            deviation_sum = sum(.data$deviation),
            .by = c("p_rate", "p_strength", "p_dir", "sim")
          )

      }

      if (pattern == "slopes") {

        dat <- dat |>
          dplyr::mutate(
            data =
              furrr::future_pmap_dbl(
                ## mapping list
                .l = list(.data$p_n,
                          .data$p_t,
                          .data$p_strength,
                          .data$p_rate,
                          .data$p_dir,
                          .data$p_res),
                ## refer to list based on index
                .f = ~ buildDGP(
                  # varying parameters
                  n = ..1,
                  t = ..2,
                  strength = ..3,
                  rate = ..4,
                  balance_dir = ..5,
                  balance_res = ..6,
                  reliable = reliability,
                  export = FALSE,
                  patterns = FALSE,
                  slopes = TRUE) |>
                  dplyr::summarize(ks = stats::ks.test(.data$estimate,
                                                       rf$estimate,
                                                       exact = TRUE)$statistic) |>
                  dplyr::pull(ks),
                ## for reproducibility
                .options = furrr::furrr_options(seed = TRUE),
                .progress = TRUE
              )
          ) |>
          ## get the distance
          dplyr::summarize(
            deviation_sum = sum(.data$data),
            .by = c("p_rate", "p_strength", "p_dir", "sim")
          )
      }

      return(dat)

    }

    # ----------------------------------------------------- #
    # GRID SEARCH ALGORITHM                                 #
    # ----------------------------------------------------- #

    cat("\n We now start the calculations. There will be 3 steps.")

    # ----------------------------------------------------- #

    # --- STEP 1

    message("\n Step 1 for grid search...")

    fd <- parallel_grid(fd)

    # ----------------------------------------------------- #

    # --- STEP 2

    message("\n Step 2 for grid search...")

    ## empty dataframe
    fd_smoothed <-
      tibble::tibble(
        p_rate = double(),
        p_strength = double(),
        p_dir = double(),
        error_smoothed = double(),
        sim = double()
      )

    ## loop across simulations
    for (sims in c(1:steps)) {
      ## loop across grids
      for (dir in c(0, 0.25, 0.5, 0.75, 1)) {
        ## construct a matrix
        mat <- fd |>
          dplyr::filter(.data$sim == sims) |>
          dplyr::filter(.data$p_dir == dir) |>
          dplyr::select(-.data$sim, -.data$p_dir) |>
          tidyr::pivot_wider(names_from = "p_strength", values_from = "deviation_sum") |>
          textshape::column_to_rownames('p_rate') |>
          as.matrix()

        ## apply neighborhood smoothing
        mat <- round(smoothMatrix(mat, n = 2), 2)

        ## add strength values
        colnames(mat) <-
          round(seq(
            from = 0.05,
            to   = 2,
            length.out = 40
          ), 3)

        ## reshape and organize the data
        mat <- mat |>
          tibble::as_tibble(.data) |>
          dplyr::mutate(p_rate = round(seq(
            from = 0,
            to   = 1,
            length.out = 41
          ), 3)) |>
          dplyr::relocate(.data$p_rate) |>
          tidyr::pivot_longer(
            cols = -.data$p_rate,
            names_to = "p_strength",
            values_to = "error_smoothed"
          ) |>
          dplyr::mutate(p_strength = as.numeric(.data$p_strength),
                        p_rate     = as.numeric(.data$p_rate)) |>
          dplyr::mutate(p_dir = dir) |>
          dplyr::mutate(sim = sims) |>
          dplyr::relocate(.data$p_dir, .before = "error_smoothed")

        ## bind the new loop step
        fd_smoothed <- fd_smoothed |> dplyr::bind_rows(mat)

      }
    }

    ## join
    fd <- fd |>
      dplyr::left_join(fd_smoothed, by = c("p_rate", "p_strength", "p_dir", "sim")) |>
      dplyr::select(.data$p_rate,
                    .data$p_strength,
                    .data$p_dir,
                    error = .data$error_smoothed)
    rm(fd_smoothed, mat) # clean-up

    # ----------------------------------------------------- #

    # --- STEP 3

    message("\n Step 3 for grid search...")

    fd <- fd |>
      dplyr::mutate(p_dir = as.factor(.data$p_dir)) # for interactions

    ## model
    m <- mgcv::gam(
      error ~ p_dir + te(p_rate, p_strength, by = p_dir), data = fd)

    # predict
    fd$pred <- as.vector(mgcv::predict.gam(m))
    fd <- fd |>
      dplyr::summarize(error = mean(.data$pred),
                       .by = c("p_rate", "p_strength", "p_dir"))

    # ----------------------------------------------------- #
    # EXPORT                                                #
    # ----------------------------------------------------- #

    fd <- fd |>
      ## cleanup
      dplyr::mutate(p_dir = factor(
        .data$p_dir,
        levels = c(0, 0.25, 0.5, 0.75, 1),
        labels = c(
          "100% Down",
          "75% Down-25% Up",
          "50% Down-50% Up",
          "25% Down-75% Up",
          "100% Up"
        )
      )) |>
      ## selection
      dplyr::select(
        rate = "p_rate",
        strength = "p_strength",
        direction = "p_dir",
        "error"
      )

    if (pattern == "contingency") {

      fd <- fd |>
        dplyr::mutate(pattern = "contingency") |>
        dplyr::mutate(n = df |> dplyr::distinct(.data$pid) |> nrow())

      return(fd)

    }

    if (pattern == "slopes") {

      fd <- fd |>
        dplyr::mutate(pattern = "slopes") |>
        dplyr::mutate(n = df |> dplyr::distinct(.data$pid) |> nrow())

      return(fd)

    }

  }

# ---------------------------------------------------------------------------- #
