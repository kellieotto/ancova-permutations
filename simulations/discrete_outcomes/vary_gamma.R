source("simulation_utils.R")
library(doParallel)

registerDoParallel(cores=30)

gamma_vec <- seq(0, 0.2, by=0.05)
gamma_res <- foreach(i = seq_along(gamma_vec)) %dopar% {
  set.seed(760682460) # Generated from random.org Timestamp: 2016-11-14 10:21:12 UTC
  tmp <- generate_simulated_data(gamma = gamma_vec[i], effect = "same effect", errors = "normal")
  pvalues <- replicate(10000, {
    tmp$Z <- permute_within_groups(tmp$Z, tmp$stratumID)
    tmp$Y1 <-gen_y(gamma = gamma_vec[i], v=tmp$v, error=tmp$delta, tmp$Z)
    generate_simulated_pvalues(tmp, reps=1e3)
  })
  return(pvalues)
}

gamma_res <- lapply(gamma_res, t)

tr_effect <- rep(0, 5)
for( i in seq_along(gamma_vec)){
  set.seed(760682460) # Generated from random.org Timestamp: 2016-11-14 10:21:12 UTC
  tmp <- generate_simulated_data(gamma = gamma_vec[i], effect = "same effect", errors = "normal")
  tr_effect[i] <- mean(gamma_vec[i] * exp(tmp$v))
}

save(tr_effect, gamma_res, file="vary_gamma_results.Rda")