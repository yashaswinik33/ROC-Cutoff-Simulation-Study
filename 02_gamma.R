library(here)
library(cutpointr)


##### Function to estimate mu2 and sigma2 #####
estimate_gamma_parameters <- function(mu1, sigma1, d = 0.2, sigma2) {
  # Shape and scale for group 1
  shape1 <- (mu1 / sigma1)^2
  scale1 <- sigma1^2 / mu1
  
  # IQR function for gamma
  iqr_gamma <- function(shape, scale) {
    qgamma(0.75, shape, scale = scale) - qgamma(0.25, shape, scale = scale)
  }
  
  # Objective function: minimize difference between median2 and target
  objective <- function(mu2) {
    shape2 <- (mu2 / sigma2)^2
    scale2 <- sigma2^2 / mu2
    
    iqr1 <- iqr_gamma(shape1, scale1)
    iqr2 <- iqr_gamma(shape2, scale2)
    iqr_pooled <- (iqr1 + iqr2) / 2
    
    median1 <- qgamma(0.5, shape1, scale = scale1)
    target_median2 <- median1 + 1.35 * d * iqr_pooled
    median2 <- qgamma(0.5, shape2, scale = scale2)
    
    return((median2 - target_median2)^2)
  }
  
  # Initial guess for mu2
  init_mu2 <- mu1 + 1  # or something sensible
  
  # Optimize only mu2
  res <- optim(par = init_mu2, fn = objective, method = "L-BFGS-B", lower = 0.01, upper = 100)
  
  mu2 <- res$par
  shape2 <- (mu2 / sigma2)^2
  scale2 <- sigma2^2 / mu2
  
  return(list(
    mu2 = mu2,
    sigma2 = sigma2,
    shape2 = shape2,
    scale2 = scale2
  ))
}


##### Gamma ROC function #####
gamma_roc <- function(mu1, mu2, sigma1, sigma2,
                      n = 50, nsim = 10,
                      conf_level = 0.95) {
  
  # metrics & names
  metrics <- list(
    youden, accuracy, sum_sens_spec, sum_ppv_npv,
    prod_sens_spec, prod_ppv_npv, cohens_kappa,
    roc01, p_chisquared,
    odds_ratio, risk_ratio, F1_score
  )
  metric_names <- c(
    "youden", "accuracy", "sum_sens_spec", "sum_ppv_npv",
    "prod_sens_spec", "prod_ppv_npv", "cohens_kappa",
    "roc01", "p_chisquared",
    "odds_ratio", "risk_ratio", "F1_score"
  )
  K <- length(metrics)
  
  # storage
  sens_mat <- spec_mat <- acc_mat <- auc_mat <- matrix(NA_real_, nrow = nsim, ncol = K)
  
  # Monte Carlo loop
  s <- 1
  while (s <= nsim) {
    
    # Convert mu and sigma to shape and scale for Gamma
    shape1 <- (mu1 / sigma1)^2
    scale1 <- sigma1^2 / mu1
    shape2 <- (mu2 / sigma2)^2
    scale2 <- sigma2^2 / mu2
    
    
    x1 <- rgamma(n, shape = shape1, scale = scale1)
    x2 <- rgamma(n, shape = shape2, scale = scale2)
    
    D  <- data.frame(NLR = round(c(x1, x2), 2),
                     group = c(rep(1,n), rep(0,n)))
    
    success <- TRUE
    res_list <- mapply(function(mtf, nm) {
      tryCatch({
        cutpointr(D, x = NLR, class = group, metric = mtf)
      }, error = function(e) {
        success <<- FALSE
        return(NULL)
      })
    }, mtf = metrics, nm = metric_names,
    SIMPLIFY = FALSE)
    
    if (!success) next
    
    sens_mat[s, ] <- vapply(res_list, function(x) x$sensitivity, numeric(1))
    spec_mat[s, ] <- vapply(res_list, function(x) x$specificity, numeric(1))
    acc_mat [s, ] <- vapply(res_list, function(x) x$acc, numeric(1))
    auc_mat [s, ] <- vapply(res_list, function(x) x$AUC, numeric(1))
    
    s <- s + 1
  }
  
  summarise_metric <- function(mat, conf_level = 0.95) {
    m  <- colMeans(mat, na.rm = TRUE)
    sd <- apply(mat, 2, sd, na.rm = TRUE)
    
    pct <- (1 - conf_level) / 2
    ci_low  <- apply(mat, 2, quantile, probs = pct, na.rm = TRUE)
    ci_high <- apply(mat, 2, quantile, probs = 1 - pct, na.rm = TRUE)
    
    data.frame(
      mean = m,
      sd   = sd,
      ci_low = ci_low,
      ci_high = ci_high
    )
  }
  
  out <- data.frame(Metric = metric_names,
                    Sensitivity = summarise_metric(sens_mat),
                    Specificity = summarise_metric(spec_mat),
                    Accuracy    = summarise_metric(acc_mat),
                    AUC         = summarise_metric(auc_mat),
                    row.names = NULL,
                    check.names = FALSE)
  
  return(out)
}


##### Simulation 1:  k = 1, d = 0.2#####
nsim = 500
mu1 <- 2.5
sigma1 <- 0.6


k <- 1 # to be changed
d <- 0.2 # to be changed


sigma2 <- k*sigma1


e <- estimate_gamma_parameters(mu1, sigma1, d,sigma2)
mu2 <- e$mu2


ns <- c(30,seq(50,1000,50))

results <- list()

start_time <- Sys.time()
for (i in ns) {
  
  df <- gamma_roc(mu1, mu2, sigma1, sigma2,
                   n = i, nsim = nsim, conf_level = conf_level)
  df$n <- i
  results[[as.character(i)]] <- df
}

end_time <- Sys.time()
execution_time <- end_time - start_time
print(paste("Execution time:", execution_time))

results1 <- do.call(rbind, results)
write.csv(results1, "gamma_results/gamma_k1_d0.2.csv", row.names = FALSE)




##### Simulation 2:  k = 1, d = 0.5#####
nsim = 500
mu1 <- 2.5
sigma1 <- 0.6


k <- 1 # to be changed
d <- 0.5 # to be changed


sigma2 <- k*sigma1


e <- estimate_gamma_parameters(mu1, sigma1, d,sigma2)
mu2 <- e$mu2


ns <- c(30,seq(50,1000,50))

results <- list()

start_time <- Sys.time()
for (i in ns) {
  
  df <- gamma_roc(mu1, mu2, sigma1, sigma2,
                  n = i, nsim = nsim, conf_level = conf_level)
  df$n <- i
  results[[as.character(i)]] <- df
}

end_time <- Sys.time()
execution_time <- end_time - start_time
print(paste("Execution time:", execution_time))

results1 <- do.call(rbind, results)
write.csv(results1, "gamma_results/gamma_k1_d0.5.csv", row.names = FALSE)



##### Simulation 3:  k = 1, d = 0.8#####
nsim = 500
mu1 <- 2.5
sigma1 <- 0.6


k <- 1 # to be changed
d <- 0.8 # to be changed


sigma2 <- k*sigma1


e <- estimate_gamma_parameters(mu1, sigma1, d,sigma2)
mu2 <- e$mu2


ns <- c(30,seq(50,1000,50))

results <- list()

start_time <- Sys.time()
for (i in ns) {
  
  df <- gamma_roc(mu1, mu2, sigma1, sigma2,
                  n = i, nsim = nsim, conf_level = conf_level)
  df$n <- i
  results[[as.character(i)]] <- df
}

end_time <- Sys.time()
execution_time <- end_time - start_time
print(paste("Execution time:", execution_time))

results1 <- do.call(rbind, results)
write.csv(results1, "gamma_results/gamma_k1_d0.8.csv", row.names = FALSE)





##### Simulation 4:  k = 1.5, d = 0.2#####
nsim = 500
mu1 <- 2.5
sigma1 <- 0.6


k <- 1.5 # to be changed
d <- 0.2 # to be changed


sigma2 <- k*sigma1


e <- estimate_gamma_parameters(mu1, sigma1, d,sigma2)
mu2 <- e$mu2


ns <- c(30,seq(50,1000,50))

results <- list()

start_time <- Sys.time()
for (i in ns) {
  
  df <- gamma_roc(mu1, mu2, sigma1, sigma2,
                  n = i, nsim = nsim, conf_level = conf_level)
  df$n <- i
  results[[as.character(i)]] <- df
}

end_time <- Sys.time()
execution_time <- end_time - start_time
print(paste("Execution time:", execution_time))

results1 <- do.call(rbind, results)
write.csv(results1, "gamma_results/gamma_k1.5_d0.2.csv", row.names = FALSE)




##### Simulation 5:  k = 1.5, d = 0.5#####
nsim = 500
mu1 <- 2.5
sigma1 <- 0.6


k <- 1.5 # to be changed
d <- 0.5 # to be changed


sigma2 <- k*sigma1


e <- estimate_gamma_parameters(mu1, sigma1, d,sigma2)
mu2 <- e$mu2


ns <- c(30,seq(50,1000,50))

results <- list()

start_time <- Sys.time()
for (i in ns) {
  
  df <- gamma_roc(mu1, mu2, sigma1, sigma2,
                  n = i, nsim = nsim, conf_level = conf_level)
  df$n <- i
  results[[as.character(i)]] <- df
}

end_time <- Sys.time()
execution_time <- end_time - start_time
print(paste("Execution time:", execution_time))

results1 <- do.call(rbind, results)
write.csv(results1, "gamma_results/gamma_k1.5_d0.5.csv", row.names = FALSE)




##### Simulation 6:  k = 1.5, d = 0.8#####
nsim = 500
mu1 <- 2.5
sigma1 <- 0.6


k <- 1.5 # to be changed
d <- 0.8 # to be changed


sigma2 <- k*sigma1


e <- estimate_gamma_parameters(mu1, sigma1, d,sigma2)
mu2 <- e$mu2


ns <- c(30,seq(50,1000,50))

results <- list()

start_time <- Sys.time()
for (i in ns) {
  
  df <- gamma_roc(mu1, mu2, sigma1, sigma2,
                  n = i, nsim = nsim, conf_level = conf_level)
  df$n <- i
  results[[as.character(i)]] <- df
}

end_time <- Sys.time()
execution_time <- end_time - start_time
print(paste("Execution time:", execution_time))

results1 <- do.call(rbind, results)
write.csv(results1, "gamma_results/gamma_k1.5_d0.8.csv", row.names = FALSE)




##### Simulation 7:  k = 2, d = 0.2#####
nsim = 500
mu1 <- 2.5
sigma1 <- 0.6


k <- 2 # to be changed
d <- 0.2 # to be changed


sigma2 <- k*sigma1


e <- estimate_gamma_parameters(mu1, sigma1, d,sigma2)
mu2 <- e$mu2


ns <- c(30,seq(50,1000,50))

results <- list()

start_time <- Sys.time()
for (i in ns) {
  
  df <- gamma_roc(mu1, mu2, sigma1, sigma2,
                  n = i, nsim = nsim, conf_level = conf_level)
  df$n <- i
  results[[as.character(i)]] <- df
}

end_time <- Sys.time()
execution_time <- end_time - start_time
print(paste("Execution time:", execution_time))

results1 <- do.call(rbind, results)
write.csv(results1, "gamma_results/gamma_k2_d0.2.csv", row.names = FALSE)



##### Simulation 8:  k = 2, d = 0.5#####
nsim = 500
mu1 <- 2.5
sigma1 <- 0.6


k <- 2 # to be changed
d <- 0.5 # to be changed


sigma2 <- k*sigma1


e <- estimate_gamma_parameters(mu1, sigma1, d,sigma2)
mu2 <- e$mu2


ns <- c(30,seq(50,1000,50))

results <- list()

start_time <- Sys.time()
for (i in ns) {
  
  df <- gamma_roc(mu1, mu2, sigma1, sigma2,
                  n = i, nsim = nsim, conf_level = conf_level)
  df$n <- i
  results[[as.character(i)]] <- df
}

end_time <- Sys.time()
execution_time <- end_time - start_time
print(paste("Execution time:", execution_time))

results1 <- do.call(rbind, results)
write.csv(results1, "gamma_results/gamma_k2_d0.5.csv", row.names = FALSE)




##### Simulation 9:  k = 2, d = 0.8#####
nsim = 500
mu1 <- 2.5
sigma1 <- 0.6


k <- 2 # to be changed
d <- 0.8 # to be changed


sigma2 <- k*sigma1


e <- estimate_gamma_parameters(mu1, sigma1, d,sigma2)
mu2 <- e$mu2


ns <- c(30,seq(50,1000,50))

results <- list()

start_time <- Sys.time()
for (i in ns) {
  
  df <- gamma_roc(mu1, mu2, sigma1, sigma2,
                  n = i, nsim = nsim, conf_level = conf_level)
  df$n <- i
  results[[as.character(i)]] <- df
}

end_time <- Sys.time()
execution_time <- end_time - start_time
print(paste("Execution time:", execution_time))

results1 <- do.call(rbind, results)
write.csv(results1, "gamma_results/gamma_k_d0.8.csv", row.names = FALSE)







































### rough work ####
mu1 <- 2.5
sigma1 <- 0.6
sigma2 <- 0.9
k <- 1.5
sigma2 <- k*sigma1
d <- 0.2

res <- estimate_gamma_parameters(mu1, sigma1, d,sigma2)
print(res)


shape1 <- (mu1 / sigma1)^2
scale1 <- sigma1^2 / mu1
shape2 <- res$shape2
scale2 <- res$scale2

median1 <- qgamma(0.5, shape1, scale = scale1)
median2 <- qgamma(0.5, shape2, scale = scale2)
iqr1 <- qgamma(0.75, shape1, scale = scale1) - qgamma(0.25, shape1, scale = scale1)
iqr2 <- qgamma(0.75, shape2, scale = scale2) - qgamma(0.25, shape2, scale = scale2)
iqr_pooled <- (iqr1 + iqr2)/2

# Compute observed d
observed_d <- (median2 - median1) / (1.35 * iqr_pooled)
print(observed_d)  

