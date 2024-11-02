
# data table
.datatable.aware = TRUE

# global variables
## usethis namespace: start
#' @importFrom rlang :=
#' @importFrom rlang .data
## usethis namespace: end

globalVariables(unique(
  c(
    # buildDGP:
    "change",
    "changer",
    "pid",
    "position",
    "sim_counts",
    "timing",
    "upper",
    "y",
    "y_obs",
    # gridPlot:
    "error",
    "group",
    "rate",
    "strength",
    # gridSearch:
    "deviation_sum",
    "error",
    "error_1",
    "error_2",
    "error_3",
    "ks",
    "ks_stat",
    "ks_stat1",
    "ks_stat2",
    "ks_stat3",
    "p_dir",
    "p_rate",
    "p_strength",
    "patterns",
    "pid",
    "y",
    # gridSearch : parallel_grid:
    "deviation",
    "p_dir",
    "p_n",
    "p_rate",
    "p_res",
    "p_strength",
    "p_t",
    "reference",
    "sim_counts",
    # simulateChangers:
    "confused",
    "confusion",
    "correspondence",
    "p_dir",
    "p_n",
    "p_rate",
    "p_rel",
    "p_res",
    "p_strength",
    "p_t",
    "sims",
    "slopes",
    # varyingSlopes:
    ".SD",
    "estimate",
    "pid"
  )
))
