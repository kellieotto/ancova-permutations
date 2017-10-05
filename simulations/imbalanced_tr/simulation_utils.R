#devtools::install_github("statlab/permuter")
#require(permuter)
require(dplyr)
require(reshape2)
require(ggplot2)
require(xtable)
options(xtable.comment = FALSE, xtable.include.rownames = FALSE, xtable.align = "p{3.5cm}")
require(grid)
require(gridExtra)

report_theme <- theme_bw() + theme(
  axis.text.x = element_text(size = 14, angle = -45),
  axis.text.y = element_text(size = 14),
  axis.title = element_text(size = 20),
  title = element_text(size = 20),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  strip.text.x = element_text(size = 14)
)

permute_within_groups <- function(x, group) {
  if(is.vector(x)){
    for(g in unique(group)){
      gg <- (group == g)
      x[gg] <- sample(x[gg])
    }
  } else if(is.data.frame(x) | is.matrix(x)){
    for(g in unique(group)){
      gg <- which(group == g)
      x[gg, ] <- x[sample(gg), ]
    }
  } else {
    stop("x is an invalid data type")
  }
  return(x)
}



stratified_two_sample <- function(group, response, stratum, stat = "mean", reps = 1000) {
  if (!is.vector(group) | !is.vector(response) | !is.vector(stratum)) {
    stop("inputs must be vectors")
  }
  if (!is.numeric(response)) {
    stop("response must be numeric")
  }
  
  if (length(unique(group)) > 2) {
    stop("two samples only")
  }
  
  groups <- unique(group)
  strata <- unique(stratum)
  
  ordering <- order(group)
  response <- response[ordering]
  stratum <- stratum[ordering]
  group <- group[ordering]
  
  ntreat <- table(group)[1]
  N <- length(group)
  
  # If stat is callable, use it as the test function. Otherwise, look in the
  # dictionary
  stats <- list(mean = function(u) {
    mean(u[1:ntreat], na.rm = TRUE) - mean(u[(ntreat + 1):N], na.rm = TRUE)
  }, t = function(u) {
    t.test(u[1:ntreat], u[(ntreat + 1):N], var.equal = TRUE)$statistic
  }, mean_within_strata = function(u) {
    sum(abs(within_group_mean(group, u, stratum, groups, strata)))
  })
  if (is.function(stat)) {
    tst_fun <- stat
  } else {
    if (stat %in% names(stats)) {
      tst_fun <- stats[[stat]]
    } else {
      stop("stat must be in the dictionary of stats or a function")
    }
  }
  
  
  distr <- replicate(reps, {
    tst_fun(permute_within_groups(response, strata))
  })
  return(distr)
}


t2p <- function(tst, distr, alternative = c("greater", "less", "two-sided")) {
  
  # check that distr is a vector with appropriate size
  B <- sum(!is.na(distr))
  
  # check that t is just a single number
  
  p <- c()
  # pupper <- c('Upper' = mean(distr >= tst, na.rm = TRUE))
  pupper <- c(Upper = mean(distr >= tst, na.rm = TRUE))
  plower <- c(Lower = mean(distr <= tst, na.rm = TRUE))
  if ("greater" %in% alternative) {
    p <- c(p, pupper)
  }
  if ("less" %in% alternative) {
    p <- c(p, plower)
  }
  if ("two-sided" %in% alternative) {
    pboth <- c(`Two-sided` = 2 * min(c(pupper, plower)))
    p <- c(p, pboth)
  }
  return(p)
}



gen_y <- function(gamma, v, error, Z){
  y <- 0.5*((2*Z-1)*gamma*exp(v) + exp(v/2)) + error
  return(y)
}

gen_x <- function(gamma, v, error){
  x <- 0.5*(-1*gamma*exp(v) + exp(v/2)) + error
  return(x)
}

generate_simulated_data <- function(gamma, effect, errors, n = c(16, 16, 16), ntreat = c(8, 8, 8)){
  # Input:
  # gamma = multiplier for the magnitude of the treatment effect
  # effect = "same effect" or "heterogeneous"
  # errors = "normal" or "heavy"
  # n = number of individuals at each stratum
  # ntreat = number of treated individuals in each stratum
  # Returns: a dataframe containing columns named Y1 (response), Y0 (baseline), Z (treatment), gamma_vec (treatment effect per individual), stratumID (stratum), stratum_effect (beta coefficient per individual), and epsilon (errors)
  
  stratumID <- rep(1:3, times = n)
  N <- sum(n)
  
  # What is the treatment effect?
  if(effect == "same effect"){
    v <- runif(N, min=-4, max=4)
  } else if(effect == "heterogeneous"){
    v <- rep(0, N)
    v[stratumID == 1] <- runif(n[1], min=-4, max=-1)
    v[stratumID == 2] <- runif(n[2], min=-1, max=1)
    v[stratumID == 3] <- runif(n[3], min=1, max=4)
  } else {
    stop("invalid parameter effect")
  }
  
  # Generate errors
  if(errors == "normal"){
    epsilon <- rnorm(N)
    delta <- rnorm(N)
  } else if (errors == "t"){
    epsilon <- rt(N, df = 2)
    delta <- rt(N, df = 2)
  } else if (errors == "lognormal"){
    epsilon <- rlnorm(N)
    delta <- rlnorm(N)
  } else if (errors == "exponential"){
    epsilon <- rexp(N) - 1
    delta <- rexp(N) - 1
  } else {
    stop("invalid errors parameter")
  }
  
  # Generate covariates
  Z <- rep(0, N)
  for(i in 1:3){
    Z[stratumID == i] <- rep(1:0, times=c(ntreat[i], n[i] - ntreat[i]))
  }
  Y0 <- gen_x(gamma, v, epsilon)
  Y1 <- gen_y(gamma, v, delta, Z)
  return(data.frame(Y1, Y0, Z, v, stratumID, epsilon, delta))
}

generate_simulated_pvalues <- function(dataset, reps = 1e3){
  # Inputs:
  # dataset = a dataframe containing columns named Y1 (response), Y0 (baseline), Z (treatment), and stratumID (stratum)
  # Returns: a vector of p-values
  # first element is the p-value from the ANCOVA
  # second element is the p-value from the stratified two-sample permutation test
  # third element is the p-value from the linear model test, permuting treatment
  # fourth element is the p-value from the Freedman-Lane linear model test, permuting residuals
  
  # ANCOVA
  modelfit <- lm(Y1 ~ Y0 + Z + factor(stratumID), data = dataset)
  resanova <- summary(aov(modelfit))
  anova_pvalue <- resanova[[1]]["Z", "Pr(>F)"]
  
  # Stratified permutation test of Y1
  observed_diff_means <- mean(dataset$Y1[dataset$Z == 1]) - mean(dataset$Y1[dataset$Z == 0])
  diff_means_distr <- stratified_two_sample(group = dataset$Z, response = dataset$Y1, stratum = dataset$stratumID, reps = reps)
  perm_pvalue <- t2p(observed_diff_means, diff_means_distr, alternative = "two-sided")
  
  # Diffed permutation test of Y1-Y0
  dataset$diff <- dataset$Y1 - dataset$Y0
  observed_diff_means2 <- mean(dataset$diff[dataset$Z == 1]) - mean(dataset$diff[dataset$Z == 0])
  diff_means_distr2 <- stratified_two_sample(group = dataset$Z, response = dataset$diff, stratum = dataset$stratumID, reps = reps)
  perm_pvalue2 <- t2p(observed_diff_means2, diff_means_distr2, alternative = "two-sided")
  
  # Permutation of treatment in linear model
  observed_t1 <- summary(modelfit)[["coefficients"]]["Z", "t value"]
  
  # Freedman-Lane linear model residual permutation
  lm2_no_tr <- lm(Y1~Y0 + factor(stratumID), data = dataset)
  dataset$lm2_resid <- residuals(lm2_no_tr)
  lm2_yhat <- fitted(lm2_no_tr)
  
  lm1and2_t_distr <-  replicate(reps, {
    dataset[, c("Z_perm", "lm2_resid_perm")] <- permute_within_groups(dataset[, c("Z", "lm2_resid")], dataset$stratumID)
    lm1_perm <- lm(Y1 ~ Y0 + Z_perm + factor(stratumID), data = dataset)
    
    dataset$response_fl <-  lm2_yhat + dataset$lm2_resid_perm
    lm2_perm <- lm(response_fl ~ Y0 + Z + factor(stratumID), data = dataset)
    
    c(summary(lm1_perm)[["coefficients"]]["Z_perm", "t value"],
      summary(lm2_perm)[["coefficients"]]["Z", "t value"])
  })
  lm_pvalue <- t2p(observed_t1, lm1and2_t_distr[1, ], alternative = "two-sided")
  fl_pvalue <- t2p(observed_t1, lm1and2_t_distr[2, ], alternative = "two-sided")
  
  return(c("ANCOVA" = anova_pvalue, 
           "Stratified Permutation" = perm_pvalue, 
           "Differenced Permutation" = perm_pvalue2,
           "LM Permutation" = lm_pvalue, 
           "Freedman-Lane" = fl_pvalue))
}

compute_power <- function(pvalues){
  sapply((0:99)/100, function(p) mean(pvalues <= p, na.rm = TRUE))
}

plot_power_curves <- function(power_mat, title){
  melt(power_mat) %>% 
    mutate("pvalue" = Var1/100) %>%
    mutate("Method" = Var2) %>%
    ggplot(aes_string(x = "pvalue", y = "value", color = "Method")) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("P-value") +
    ylab("Power") +
    ggtitle(title) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 16),
      title = element_text(size = 16),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 14),
      strip.text.x = element_text(size = 12)
    ) 
}

plot_power_gamma <- function(powermat_list, gamma_vec, alpha, title){
  gamma_power_alpha <- t(sapply(powermat_list, function(x) x[floor(alpha*nrow(x)),]))
  colnames(gamma_power_alpha) <- c("ANCOVA", "Stratified Permutation", "Differenced Permutation", "LM Permutation", "Freedman-Lane")
  melt(gamma_power_alpha, value.name="Power") %>%
    mutate("Method" = Var2) %>%
    mutate("Effect" = rep(gamma_vec, 5)) %>% 
    ggplot(aes(x = Effect, y = Power)) +
    geom_line(aes(color = Method))+
    xlab("Effect size") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 16),
      title = element_text(size = 16),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 14),
      strip.text.x = element_text(size = 12)
    ) 
}