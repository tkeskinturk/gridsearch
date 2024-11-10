
# - `gridsearch` ------------------------------------------------------------- #
# - function: grid search ---------------------------------------------------- #

### note: this function
###       (a) takes a panel dataframe (with `pid`, `t`, and `y`),
###       (b) iteratively prune to recover best DGPs,
###       (c) calculates error by measuring pattern distance from the DGPs.

# FUNCTION CALL -------------------------------------------------------------- #

#' Grid Search for Adjudicating DGPs
#'
#' This function uses a three-step grid search algorithm to calculate the extent to which a dataframe can be approximated by a known DGP. It implements fake-data simulations from a grid of plausible values---rate of change ranging from 0% to 100%, strength of change ranging from 0 SD to 2 SD, and 5 directional cases (everyone changing negatively, everyone changing positively, half changing negatively and half changing positively, 25% changing negatively and 75% changing positively, and 75% changing negatively and 25% changing positive), calculates the distribution of response patterns or slopes, and provides an error term that summarizes the distance of the DGP from the observed values.
#' @param data A dataframe.
#' @param yname The outcome identifier.
#' @param tname The time identifier.
#' @param pname The unit identifier.
#' @param pattern If `pattern = contingency`, the calculations are based on the use of a contingency table of panel patterns; if `pattern = slopes`, function uses individual slope coefficients.
#' @param step1 The number of simulated datasets in Step 1.
#' @param step2 The number of simulated datasets in Step 2.
#' @param step3 The number of simulated datasets in Step 3.
#' @param reliability The assumed reliability score of the outcome measurement.
#' @param seed Seed for reproducibility.
#' @param workers The number of workers for parallelization (note that parallelization is highly recommended for reasonable duration).
#'
#' @return A data frame.
#'
#' **Basic [gridSearch()] call:**
#' ```{r, comment = "#>", collapse = TRUE}
#' set.seed(11235)
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
#'            step1 = 1,
#'            step2 = 2,
#'            step3 = 3,
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
    # step 1 runs
    step1 = 1,
    # step 2 runs
    step2 = 30,
    # step 3 runs
    step3 = 500,
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
    if(step1 > step2 | step2 > step3 | step1 > step3) {
      stop("Error: Each `step` should have higher runs than the previous ones.")
    }

    # ----------------------------------------------------- #
    # BUILDING BLOCKS                                       #
    # ----------------------------------------------------- #

    # --- part 1: initialize

    ## session
    future::plan(future::multisession, workers = workers)
    set.seed(seed)
    pattern = match.arg(pattern)

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
        seq(
          from = 0,
          to   = 1,
          length.out = 21),
      p_strength =
        seq(
          from = 0.1,
          to   = 2,
          length.out = 20),
      p_dir = c(0, 0.25, 0.50, 0.75, 1),
      p_res = mean(df$y)) |>
      ## round it up
      dplyr::mutate(p_rate = round(.data$p_rate, 2),
                    p_strength = round(.data$p_strength, 1))

    # --- part 4: function for parallelization

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
                  patterns = FALSE,
                  slopes = TRUE) |>
                  dplyr::summarize(ks = stats::ks.test(.data$estimate,
                                                       rf$estimate,
                                                       exact = TRUE)$statistic),
                ## for reproducibility
                .options = furrr::furrr_options(seed = TRUE),
                .progress = TRUE
              )
          )
      }

      return(dat)

    }

    # ----------------------------------------------------- #
    # GRID SEARCH ALGORITHM: CONTINGENCY                    #
    # ----------------------------------------------------- #

    if (pattern == "contingency") {

      cat("\n We now start the calculations. There will be 3 steps.")

      # ----------------------------------------------------- #

      # --- STEP 1

      message("\n Step 1 for grid search...")

      fd_1 <-
        tidyr::expand_grid(fd, sim = c(1:step1))
      fd_1 <- parallel_grid(fd_1)
      fd_1 <- fd_1 |>
        dplyr::summarize(
          error_1 = mean(.data$deviation_sum),
          .by = c("p_rate", "p_strength", "p_dir")
        )

      # --- STEP 1, SMOOTHING
      fd_smoothed <-
        tibble::tibble(
          ## empty vectors!
          p_rate = double(),
          p_strength = double(),
          p_dir = double(),
          error_1 = double()
        )

      ## loop across grids
      for (dir in c(0, 0.25, 0.5, 0.75, 1)) {

        ## get matrix representation
        mat <- fd_1 |>
          dplyr::filter(.data$p_dir == dir) |>
          dplyr::select(-.data$p_dir) |>
          tidyr::pivot_wider(names_from = "p_strength", values_from = "error_1") |>
          textshape::column_to_rownames('p_rate') |>
          as.matrix()

        ## smooth the matrix with neighboring cells and round
        mat <- round(smoothMatrix(mat), 2)

        ## reshape to the initial tidy version
        colnames(mat) <- round(seq(
          from = 0.1,
          to = 2,
          length.out = 20
        ), 2)
        mat <- mat |>
          tibble::as_tibble(.data) |>
          dplyr::mutate(p_rate = round(seq(
            from = 0,
            to   = 1,
            length.out = 21
          ), 2)) |>
          dplyr::relocate(.data$p_rate) |>
          tidyr::pivot_longer(cols = -.data$p_rate,
                              names_to = "p_strength",
                              values_to = "error_1") |>
          dplyr::mutate(p_strength = as.numeric(.data$p_strength),
                        p_rate     = as.numeric(.data$p_rate)) |>
          dplyr::mutate(p_dir = dir) |>
          dplyr::relocate(.data$p_dir, .before = "error_1")
        fd_smoothed <- fd_smoothed |> dplyr::bind_rows(mat)

      }

      # --- STEP 1, MAPPING
      fd <- fd |>
        dplyr::left_join(fd_smoothed, by = c("p_rate", "p_strength", "p_dir"))
      rm(fd_smoothed, mat) # clean-up

      # ----------------------------------------------------- #

      # --- STEP 2

      message("\n Step 2 for grid search...")

      fd_2 <-
        tidyr::expand_grid(fd |>
                             dplyr::arrange(.data$error_1) |>
                             dplyr::slice(1:1050),
                           sim = c(1:(step2 - step1)))
      fd_2 <- parallel_grid(fd_2)
      fd_2 <- fd_2 |>
        dplyr::summarize(
          error_2 = mean(.data$deviation_sum),
          .by = c("p_rate", "p_strength", "p_dir")
        )

      # --- STEP 2, MAPPING
      fd <- fd |>
        dplyr::left_join(fd_2, by = c("p_rate", "p_strength", "p_dir")) |>
        ## update the error term for group 2
        dplyr::mutate(error_2 = dplyr::case_when(
          is.na(.data$error_1) ~ NA_real_,
          .default =
            .data$error_1 * (step1 / step2)
          +
            .data$error_2 * (step2 - step1) / step2
        ))

      # ----------------------------------------------------- #

      # --- STEP 3

      message("\n Step 3 for grid search...")

      fd_3 <-
        tidyr::expand_grid(fd |>
                             dplyr::arrange(.data$error_2) |>
                             dplyr::slice(1:420),
                           sim = c(1:(step3 - step2)))
      fd_3 <- parallel_grid(fd_3)
      fd_3 <- fd_3 |>
        dplyr::summarize(
          error_3 = mean(.data$deviation_sum),
          .by = c("p_rate", "p_strength", "p_dir")
        )

      # --- STEP 3, MAPPING
      fd <- fd |>
        dplyr::left_join(fd_3, by = c("p_rate", "p_strength", "p_dir")) |>
        ## update the error term for group 2
        dplyr::mutate(error_3 = dplyr::case_when(
          is.na(.data$error_2) ~ NA_real_,
          .default =
            .data$error_2 * (step2 / step3)
          +
            .data$error_3 * (step3 - step2) / step3
        ))

    }

    # ----------------------------------------------------- #
    # GRID SEARCH ALGORITHM: SLOPES                         #
    # ----------------------------------------------------- #

    if (pattern == "slopes") {

      cat("\n We now start the calculations. There will be 3 steps.")

      # ----------------------------------------------------- #

      # --- STEP 1

      message("\n Step 1 for grid search...")

      fd_1 <-
        tidyr::expand_grid(fd, sim = c(1:step1))
      fd_1 <- parallel_grid(fd_1)
      fd_1 <- fd_1 |>
        tidyr::unnest(cols = data) |>
        dplyr::summarize(ks_stat1 = mean(.data$ks),
                         .by = c("p_rate", "p_strength", "p_dir"))

      # --- STEP 1, SMOOTHING
      fd_smoothed <-
        tibble::tibble(
          ## empty vectors!
          p_rate = double(),
          p_strength = double(),
          p_dir = double(),
          ks_stat1 = double()
        )

      ## loop across grids
      for (dir in c(0, 0.25, 0.5, 0.75, 1)) {

        ## get matrix representation
        mat <- fd_1 |>
          dplyr::filter(.data$p_dir == dir) |>
          dplyr::select(-.data$p_dir) |>
          tidyr::pivot_wider(names_from = "p_strength", values_from = "ks_stat1") |>
          textshape::column_to_rownames('p_rate') |>
          as.matrix()

        ## smooth the matrix with neighboring cells and round
        mat <- round(smoothMatrix(mat), 2)

        ## reshape to the initial tidy version
        colnames(mat) <- round(seq(
          from = 0.1,
          to = 2,
          length.out = 20
        ), 2)
        mat <- mat |>
          tibble::as_tibble(.data) |>
          dplyr::mutate(p_rate = round(seq(
            from = 0,
            to   = 1,
            length.out = 21
          ), 2)) |>
          dplyr::relocate(.data$p_rate) |>
          tidyr::pivot_longer(cols = -.data$p_rate,
                              names_to = "p_strength",
                              values_to = "ks_stat1") |>
          dplyr::mutate(p_strength = as.numeric(.data$p_strength),
                        p_rate     = as.numeric(.data$p_rate)) |>
          dplyr::mutate(p_dir = dir) |>
          dplyr::relocate(.data$p_dir, .before = "ks_stat1")
        fd_smoothed <- fd_smoothed |> dplyr::bind_rows(mat)

      }

      # --- STEP 1, MAPPING
      fd <- fd |>
        dplyr::left_join(fd_smoothed, by = c("p_rate", "p_strength", "p_dir"))
      rm(fd_smoothed, mat) # clean-up

      # ----------------------------------------------------- #

      # --- STEP 2

      message("\n Step 2 for grid search...")

      fd_2 <-
        tidyr::expand_grid(fd |>
                             dplyr::arrange(.data$ks_stat1) |>
                             dplyr::slice(1:1050),
                           sim = c(1:(step2 - step1)))
      fd_2 <- parallel_grid(fd_2)
      fd_2 <- fd_2 |>
        tidyr::unnest(cols = data) |>
        dplyr::summarize(ks_stat2 = mean(.data$ks),
                         .by = c("p_rate", "p_strength", "p_dir"))

      # --- STEP 2, MAPPING
      fd <- fd |>
        dplyr::left_join(fd_2, by = c("p_rate", "p_strength", "p_dir")) |>
        ## update the error term for group 2
        dplyr::mutate(ks_stat2 = dplyr::case_when(
          is.na(.data$ks_stat1) ~ NA_real_,
          .default =
            .data$ks_stat1 * (step1 / step2)
          +
            .data$ks_stat2 * (step2 - step1) / step2
        ))

      # ----------------------------------------------------- #

      # --- STEP 3

      message("\n Step 3 for grid search...")

      fd_3 <-
        tidyr::expand_grid(fd |>
                             dplyr::arrange(.data$ks_stat2) |>
                             dplyr::slice(1:420),
                           sim = c(1:(step3 - step2)))
      fd_3 <- parallel_grid(fd_3)
      fd_3 <- fd_3 |>
        tidyr::unnest(cols = data) |>
        dplyr::summarize(ks_stat3 = mean(.data$ks),
                         .by = c("p_rate", "p_strength", "p_dir"))

      # --- STEP 3, MAPPING
      fd <- fd |>
        dplyr::left_join(fd_3, by = c("p_rate", "p_strength", "p_dir")) |>
        ## update the error term for group 2
        dplyr::mutate(ks_stat3 = dplyr::case_when(
          is.na(.data$ks_stat2) ~ NA_real_,
          .default =
            .data$ks_stat2 * (step2 / step3)
          +
            .data$ks_stat3 * (step3 - step2) / step3
        ))

    }

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
      ))

    if (pattern == "contingency") {

      fd <- fd |>
        dplyr::rename(error = "error_3") |>
        dplyr::mutate(error = ifelse(is.na(.data$error),
                                     .data$error_2,
                                     .data$error)) |>
        dplyr::mutate(error = ifelse(is.na(.data$error),
                                     .data$error_1,
                                     .data$error)) |>
        dplyr::select(
          rate = "p_rate",
          strength = "p_strength",
          direction = "p_dir",
          "error"
        ) |>
        dplyr::mutate(pattern = "contingency") |>
        dplyr::mutate(n = df |> dplyr::distinct(.data$pid) |> nrow())

      return(fd)

    }

    if (pattern == "slopes") {

      fd <- fd |>
        dplyr::rename(ks_stat = "ks_stat3") |>
        dplyr::mutate(ks_stat = ifelse(is.na(.data$ks_stat),
                                       .data$ks_stat2,
                                       .data$ks_stat)) |>
        dplyr::mutate(ks_stat = ifelse(is.na(.data$ks_stat),
                                       .data$ks_stat1,
                                       .data$ks_stat)) |>
        dplyr::select(
          rate = "p_rate",
          strength = "p_strength",
          direction = "p_dir",
          error = "ks_stat"
        ) |>
        dplyr::mutate(pattern = "slopes") |>
        dplyr::mutate(n = df |> dplyr::distinct(.data$pid) |> nrow())

      return(fd)

    }

  }

# ---------------------------------------------------------------------------- #
