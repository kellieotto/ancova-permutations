source("simulation_utils.R")
library(doParallel)

registerDoParallel(cores=30)

ntreat_vec <- list( c(8, 8, 8),
                 c(4, 4, 4),
                 c(4, 4, 8),
                 c(4, 8, 8),
                 c(4, 8, 12),
                 c(8, 8, 12),
                 c(8, 12, 12),
                 c(12, 12, 12))
  
  
gamma_res <- foreach(i = seq_along(ntreat_vec)) %dopar% {
  set.seed(760682460) # Generated from random.org Timestamp: 2016-11-14 10:21:12 UTC
  tmp <- generate_simulated_data(gamma = 0.2, effect = "same effect", errors = "normal", ntreat = ntreat_vec[[i]])
  pvalues <- replicate(10000, {
    tmp$Z <- permute_within_groups(tmp$Z, tmp$stratumID)
    tmp$Y1 <- gen_y(gamma = 0.2, v=tmp$v, error=tmp$delta, tmp$Z)
    generate_simulated_pvalues(tmp, reps=1e3)
  })
  return(pvalues)
}

gamma_res <- lapply(gamma_res, t)

tr_effect <- rep(0, length(ntreat_vec))
for( i in seq_along(ntreat_vec)){
  set.seed(760682460) # Generated from random.org Timestamp: 2016-11-14 10:21:12 UTC
  tmp <- generate_simulated_data(gamma = 0.2, effect = "same effect", errors = "normal", ntreat = ntreat_vec[[i]])
  tr_effect[i] <- mean(0.2 * exp(tmp$v))
}

save(tr_effect, gamma_res, file="vary_ntreat_results.Rda")
