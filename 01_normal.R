library(here)
library(cutpointr)


##### Normal ROC function #####
normal_roc <- function(mu1, mu2, sigma1, sigma2,
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
  
  # store simulation results
  sens_mat <- spec_mat <- acc_mat <- auc_mat <- matrix(NA_real_, nrow = nsim, ncol = K)
  
  # Monte‑Carlo loop 
  s <- 1
  while (s <= nsim) {
  
    x1 <- rnorm(n, mu1, sigma1)
    x2 <- rnorm(n, mu2, sigma2)
    D  <- data.frame(NLR = round(c(x1, x2), 2),
                     group = c(rep(1,n), rep(0,n)))
    
    success <- TRUE
    
    res_list <- mapply(function(mtf, nm) {
      tryCatch({
        cutpointr(D, x = NLR, class = group, metric = mtf)
      }, error = function(e) {
        success <<- FALSE # Set flag to regenerate
        return(NULL)
      })
    }, mtf = metrics, nm = metric_names,
    SIMPLIFY = FALSE)
    
    if (!success) {
      next  # Skip to next iteration of while loop, regenerate data
    }
    
    sens_mat[s, ] <- vapply(res_list, function(x) x$sensitivity, numeric(1))
    spec_mat[s, ] <- vapply(res_list, function(x) x$specificity, numeric(1))
    acc_mat [s, ] <- vapply(res_list, function(x) x$acc, numeric(1))
    auc_mat [s, ] <- vapply(res_list, function(x) x$AUC, numeric(1))
    
    s <- s + 1  # Proceed only if successful
  }
  
  
  # mean(SD) and bootstrap CI
  summarise_metric <- function(mat, conf_level = 0.95) {
    m  <- colMeans(mat, na.rm = TRUE)
    sd <- apply(mat, 2, sd, na.rm = TRUE)
    
    # direct empirical percentile CI
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



##### Simulation 1:  k = 1, d = 0.2 #####
nsim = 500
mu1 <- 2.5
sigma1 <- 0.6

k <- 1 # to be changed
d <- 0.2 # to be changed


sigma2 <- k * sigma1
mu2 <- mu1 + d*sqrt((sigma1^2+sigma2^2)/2)

ns <- c(30,seq(50,1000,50))

results <- list()

start_time <- Sys.time()
for (i in ns) {
 
  df <- normal_roc(mu1, mu2, sigma1, sigma2,
                   n = i, nsim = nsim, conf_level = conf_level)
  df$n <- i
  results[[as.character(i)]] <- df
}

end_time <- Sys.time()
execution_time <- end_time - start_time
print(paste("Execution time:", execution_time))

results1 <- do.call(rbind, results)
write.csv(results1, "normal_results/normal_k1_d0.2.csv", row.names = FALSE)



##### Simulation 2:  k = 1, d = 0.5 #####
nsim = 500
mu1 <- 2.5
sigma1 <- 0.6

k <- 1 # to be changed
d <- 0.5 # to be changed


sigma2 <- k * sigma1
mu2 <- mu1 + d*sqrt((sigma1^2+sigma2^2)/2)

ns <- c(30,seq(50,1000,50))

results <- list()

start_time <- Sys.time()
for (i in ns) {
  
  df <- normal_roc(mu1, mu2, sigma1, sigma2,
                   n = i, nsim = nsim, conf_level = conf_level)
  df$n <- i
  results[[as.character(i)]] <- df
}

end_time <- Sys.time()
execution_time <- end_time - start_time
print(paste("Execution time:", execution_time))

results1 <- do.call(rbind, results)
write.csv(results1, "normal_results/normal_k1_d0.5.csv", row.names = FALSE)




##### Simulation 3:  k = 1, d = 0.8 #####
nsim = 500
m <- 2.5
sigma1 <- 0.6

k <- 1 # to be changed
d <- 0.8 # to be changed


sigma2 <- k * sigma1
mu2 <- mu1 + d*sqrt((sigma1^2+sigma2^2)/2)

ns <- c(30,seq(50,1000,50))

results <- list()

start_time <- Sys.time()
for (i in ns) {
  
  df <- normal_roc(mu1, mu2, sigma1, sigma2,
                   n = i, nsim = nsim, conf_level = conf_level)
  df$n <- i
  results[[as.character(i)]] <- df
}

end_time <- Sys.time()
execution_time <- end_time - start_time
print(paste("Execution time:", execution_time))

results1 <- do.call(rbind, results)
write.csv(results1, "normal_results/normal_k1_d0.8.csv", row.names = FALSE)




##### Simulation 4:  k = 1.5, d = 0.2 #####
nsim = 500
mu1 <- 2.5
sigma1 <- 0.6

k <- 1.5 # to be changed
d <- 0.2 # to be changed


sigma2 <- k * sigma1
mu2 <- mu1 + d*sqrt((sigma1^2+sigma2^2)/2)

ns <- c(30,seq(50,1000,50))

results <- list()

start_time <- Sys.time()
for (i in ns) {
  
  df <- normal_roc(mu1, mu2, sigma1, sigma2,
                   n = i, nsim = nsim, conf_level = conf_level)
  df$n <- i
  results[[as.character(i)]] <- df
}

end_time <- Sys.time()
execution_time <- end_time - start_time
print(paste("Execution time:", execution_time))

results1 <- do.call(rbind, results)
write.csv(results1, "normal_results/normal_k1.5_d0.2.csv", row.names = FALSE)



##### Simulation 5:  k = 1.5, d = 0.5 #####
nsim = 500
mu1 <- 2.5
sigma1 <- 0.6

k <- 1.5 # to be changed
d <- 0.5 # to be changed


sigma2 <- k * sigma1
mu2 <- mu1 + d*sqrt((sigma1^2+sigma2^2)/2)

ns <- c(30,seq(50,1000,50))

results <- list()

start_time <- Sys.time()
for (i in ns) {
  
  df <- normal_roc(mu1, mu2, sigma1, sigma2,
                   n = i, nsim = nsim, conf_level = conf_level)
  df$n <- i
  results[[as.character(i)]] <- df
}

end_time <- Sys.time()
execution_time <- end_time - start_time
print(paste("Execution time:", execution_time))

results1 <- do.call(rbind, results)
write.csv(results1, "normal_results/normal_k1.5_d0.5.csv", row.names = FALSE)



##### Simulation 6:  k = 1.5, d = 0.8 #####
nsim = 500
mu1 <- 2.5
sigma1 <- 0.6

k <- 1.5 # to be changed
d <- 0.8 # to be changed


sigma2 <- k * sigma1
mu2 <- mu1 + d*sqrt((sigma1^2+sigma2^2)/2)

ns <- c(30,seq(50,1000,50))

results <- list()

start_time <- Sys.time()
for (i in ns) {
  
  df <- normal_roc(mu1, mu2, sigma1, sigma2,
                   n = i, nsim = nsim, conf_level = conf_level)
  df$n <- i
  results[[as.character(i)]] <- df
}

end_time <- Sys.time()
execution_time <- end_time - start_time
print(paste("Execution time:", execution_time))

results1 <- do.call(rbind, results)
write.csv(results1, "normal_results/normal_k1.5_d0.8.csv", row.names = FALSE)




##### Simulation 7:  k = 2, d = 0.2 #####
nsim = 500
mu1 <- 2.5
sigma1 <- 0.6

k <- 2 # to be changed
d <- 0.2 # to be changed


sigma2 <- k * sigma1
mu2 <- mu1 + d*sqrt((sigma1^2+sigma2^2)/2)

ns <- c(30,seq(50,1000,50))

results <- list()

start_time <- Sys.time()
for (i in ns) {
  
  df <- normal_roc(mu1, mu2, sigma1, sigma2,
                   n = i, nsim = nsim, conf_level = conf_level)
  df$n <- i
  results[[as.character(i)]] <- df
}

end_time <- Sys.time()
execution_time <- end_time - start_time
print(paste("Execution time:", execution_time))

results1 <- do.call(rbind, results)
write.csv(results1, "normal_results/normal_k2_d0.2.csv", row.names = FALSE)



##### Simulation 8:  k = 2, d = 0.5 #####
nsim = 500
mu1 <- 2.5
sigma1 <- 0.6

k <- 2 # to be changed
d <- 0.5 # to be changed


sigma2 <- k * sigma1
mu2 <- mu1 + d*sqrt((sigma1^2+sigma2^2)/2)

ns <- c(30,seq(50,1000,50))

results <- list()

start_time <- Sys.time()
for (i in ns) {
  
  df <- normal_roc(mu1, mu2, sigma1, sigma2,
                   n = i, nsim = nsim, conf_level = conf_level)
  df$n <- i
  results[[as.character(i)]] <- df
}

end_time <- Sys.time()
execution_time <- end_time - start_time
print(paste("Execution time:", execution_time))

results1 <- do.call(rbind, results)
write.csv(results1, "normal_results/normal_k2_d0.5.csv", row.names = FALSE)




##### Simulation 9:  k = 2, d = 0.8 #####
nsim = 500
mu1 <- 2.5
sigma1 <- 0.6

k <- 2 # to be changed
d <- 0.8 # to be changed


sigma2 <- k * sigma1
mu2 <- mu1 + d*sqrt((sigma1^2+sigma2^2)/2)

ns <- c(30,seq(50,1000,50))

results <- list()

start_time <- Sys.time()
for (i in ns) {
  
  df <- normal_roc(mu1, mu2, sigma1, sigma2,
                   n = i, nsim = nsim, conf_level = conf_level)
  df$n <- i
  results[[as.character(i)]] <- df
}

end_time <- Sys.time()
execution_time <- end_time - start_time
print(paste("Execution time:", execution_time))

results1 <- do.call(rbind, results)
write.csv(results1, "normal_results/normal_k2_d0.8.csv", row.names = FALSE)








##### Rough #####
n <- 50
mu1 <- 3.5
sigma1 <- 0.65
mu2 <- 2.5
sigma2 <- 0.66

# Testing cutpoinr()
x1 <- rnorm(n, mu1, sigma1)
x2 <- rnorm(n, mu2, sigma2)
D  <- data.frame(NLR = round(c(x1, x2), 2),
                 group = rep(c(1, 0), each = n))
x <- cutpointr(D, x = NLR, class = group, metric = youden)

# testing normal_roc()
normal_roc(mu1, mu2, sigma1, sigma2,
           n = 50, nsim = 10,
           conf_level = 0.95)

