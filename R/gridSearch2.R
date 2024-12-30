
# - `gridsearch` ------------------------------------------------------------- #
# - function: grid search ---------------------------------------------------- #

# FUNCTION CALL -------------------------------------------------------------- #

gridSearch2 <- function(data,
                              n_samples = 1000, # Number of samples from priors
                              caliper = .1,
                              ic_min = 0,
                              ic_max = 1,
                              pc_min = 0,
                              pc_max = 1,
                              bal_min = 0,
                              bal_max = 1,
                              rel_min = 0,
                              rel_max = 1,
                              fix = "none", # Fix any parameter? (values: "ic_sample", "pc_sample"...)
                              fix_at = 1, # At what value to fix it on?
                              n_sim,
                              t,
                        slopes = FALSE,
                        ks = FALSE,
                        verbose = TRUE) {

  ts <- paste0("V",1:t)
  j <- ts
  accepted <- list()
  counter <- 0

  if(slopes == TRUE) {
    ref <- varyingSlopes2(t)
    j <- "estimate" # column to join on
    data <- collapse::join(data, ref, on = ts,  how = "left", verbose = FALSE)
    data <- aggregate(N ~ estimate, data, sum)
  }

  ###---
  for (i in 1:n_samples) {

    ###--- Sample from prior
    ic_sample <- runif(n = 1, min = ic_min, max = ic_max)
    pc_sample <- runif(n = 1, min = pc_min, max = pc_max)
    bal_sample <- runif(n = 1, min = bal_min, max = bal_max)
    rel_sample <- runif(n = 1, min = rel_min, max = rel_max)
    assign(fix, fix_at)



    ###--- Simulate the data
    res <- buildDGP2(rate = pc_sample,
                     balance_dir = bal_sample,
                     strength = ic_sample,
                     reliable = rel_sample)


    ###--- Slopes
    if(slopes == TRUE) {
      res <- collapse::join(res, ref, on = ts,  how = "left", verbose = FALSE)
      res <- aggregate(N ~ estimate, res, sum)
    }
    ###--- Compute error
    res <- collapse::join(res, data, on = j,  how = "full", verbose = FALSE)
    res[is.na(res)] <- 0

    ###-- Ks test
    if(ks == TRUE && slopes == TRUE) {
      ks <- suppressWarnings(ks.test(x = rep(res$estimate, res$N * 100), y = rep(data$estimate, data$N * 100)))
      error <- ks$statistic
    }

    else{
      res$diff <- res$N - res$N_data
      error <- sum(abs(res$diff))
    }

    ###--- Accept/reject
    if (error < caliper) {
      accepted[[i]] <- data.table::data.table(ic_sample, pc_sample, bal_sample, rel_sample)
      counter <- counter + 1
    }
    if(i %% 1000 == 0 && verbose == TRUE) {
      print(paste0(round(i/n_samples*100, 1), "%"))
      print(paste0("AR: ", round(counter/i*100, 1), "%"))
    }



  }

  accepted <- data.table::rbindlist(accepted)

  return(accepted)

}


#'
#' @export


# ---------------------------------------------------------------------------- #
