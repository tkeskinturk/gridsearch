

varyingSlopes2 <- function(t) {


  ###--- Possible trajectories
  tx <- seq(0, 1, length.out = t)
  patterns <- do.call(expand.grid, rep(list(c(0,1)), t))


  estimate <- c()

  for(i in 1:nrow(patterns)) {

    fit <- suppressWarnings(glm(as.numeric(patterns[i,]) ~ tx, family = binomial()))
    prob1 <- plogis(fit$coefficients[1])
    prob2 <- plogis(fit$coefficients[1] + fit$coefficients[2])
    estimate[i] <- round(prob2 - prob1, 3)

  }

  patterns <-
    patterns |>
    cbind("estimate" = estimate)
  colnames(patterns) <- c(paste0("V",1:t), "estimate")

  return(patterns)
}
