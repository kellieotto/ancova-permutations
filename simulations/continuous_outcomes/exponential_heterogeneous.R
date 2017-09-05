source("simulation_utils.R")
library(doParallel)
registerDoParallel(cores=5)

set.seed(760682460) # Generated from random.org Timestamp: 2016-11-14 10:21:12 UTC
tmp <- generate_simulated_data(gamma = 0.2, effect = "heterogeneous", errors = "exponential")
pvalues <- foreach(i = 1:5) %dopar% {
  tmp$Z <- permute_within_groups(tmp$Z, tmp$stratumID)
  tmp$Y1 <-gen_y1(gamma = 0.15, v=tmp$v, error=tmp$delta, tmp$Z)
  generate_simulated_pvalues(tmp)
}
pvalues <- do.call(rbind, pvalues)
exponential_heterogeneous <- apply(pvalues, 2, compute_power)
save(exponential_heterogeneous, file="exponential_heterogeneous.Rda")