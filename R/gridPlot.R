
# - `gridsearch` ------------------------------------------------------------- #
# - function: grid plot ------------------------------------------------------ #

### note: this function
###       (a) takes the resulting data frame from `gridSearch` function,
###       (b) plots a heatmap that shows the distribution of model fits.

# PART 1: FUNCTION ----------------------------------------------------------- #

#' Plot DGP Grids
#'
#' This function plots the distribution of error scores resulting from `gridSearch`.
#' @param data A data frame generated by `gridSearch`.
#' @param cut1 Cut error values for DGPs: Group 1.
#' @param cut2 Cut error values for DGPs: Group 2.
#' @param cut3 Cut error values for DGPs: Group 3.
#'
#' @return A `ggplot2` object.
#' @export
#'

gridPlot <-

  function(data,
           cut1 = 0.05,
           cut2 = 0.10,
           cut3 = 0.20) {

    lvls <- paste("Group", 1:4)

    if (data |>
        dplyr::select(dplyr::any_of(c(
          "rate", "strength", "direction", "error", "pattern", "n"
        ))) |>
        ncol() != 6)
      stop("Error: One or more columns from `gridSearch` is missing.")

    ## organize data
    data <- data |>
      dplyr::mutate(
        group = dplyr::case_when(
          error <= stats::quantile(.data$error, cut1, na.rm = TRUE) ~
            "Group 1",
          error <= stats::quantile(.data$error, cut2, na.rm = TRUE) ~
            "Group 2",
          error <= stats::quantile(.data$error, cut3, na.rm = TRUE) ~
            "Group 3",
          TRUE ~ "Group 4"
        )
      ) |>
      dplyr::mutate(group = factor(.data$group, levels = lvls))

    ## check the raw values
    if (unique(data$pattern) == "contingency" & data |>
        dplyr::mutate(error = (.data$error / .data$n) / 2) |>
        dplyr::pull(.data$error) |>
        min() > 0.10) {
      message(
        "\n Warning: DGPs are not very close to the observed values.\n Tread carefully."
      )
    }
    if (unique(data$pattern) == "slopes" & data |>
        dplyr::pull(.data$error) |>
        min() > 0.10) {
      message(
        "\n Warning: DGPs are not very close to the observed values.\n Tread carefully."
      )
    }

    ## facet design
    design <-
      c(
      "
      AABBCC
      #DDEE#
      ")

    ## plot
    plot <- data |>
      ggplot2::ggplot(ggplot2::aes(x = .data$rate, y = .data$strength)) +
      ggplot2::scale_x_continuous(labels = scales::label_percent()) +
      ggplot2::scale_y_continuous(labels = scales::pretty_breaks()) +
      ggplot2::labs(x = "Rate of Change", y = "Strength of Change") +
      hrbrthemes::theme_ipsum_rc(grid = "none") +
      ggplot2::geom_tile(ggplot2::aes(fill = .data$group),
                         color = "black",
                         size = 0.01) +
      ggplot2::scale_fill_manual(
        values = c("#000000", "#777777", "#eeeeee", "transparent"), drop = FALSE
        ) +
      ggplot2::theme(legend.position = "none")

    ## design
    plot + ggh4x::facet_manual(~direction, design = design)

  }

# ---------------------------------------------------------------------------- #
