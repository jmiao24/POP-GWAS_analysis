rm(list = ls())
library(dplyr)
library(doParallel)
library(ggplot2)
library(data.table)

source("./Fun.R")
set.seed(42)
maf <- 0.25
var_factor <- 2 * maf * (1 - maf)
mar_effect_vec <- c(0, 0.02, 0.04, 0.06)
h2_X_Y_vec <- mar_effect_vec^2 * var_factor
sim.times <- 10^3 # number of replicates
n.total <- 1.2 * 10^5 # total sample size
h2_X_z <- 0.00015
# var_res <- 0.2
cov_yz_vec <- c(seq(h2_X_z, 1 - h2_X_z, 0.05), 1 - h2_X_z)
n_lab <- 10^4
n_unlab <- 10^5

for (cov_yz in cov_yz_vec){
  # Basic simulation parameters
  # cov_yz <- 0.2
  print(cov_yz) # This should be the imputation r^2

  #-------------- Known Model ------
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  result <- foreach(
    i = 1:sim.times, .combine = rbind, .packages = c("doParallel"),
    .errorhandling = "pass"
  ) %dopar% {
    h2_X_Y <- h2_X_Y_vec[2]
    g <- rbinom(n.total, 2, maf)
    y_tmp <- scale(g) * sqrt(h2_X_Y)
    y <- y_tmp + grt_err(n.total, 1) * sqrt(c(1 - var(y_tmp)))
    z <- scale(g) * sqrt(h2_X_z) + y * sqrt(cov_yz - h2_X_z) + grt_err(n.total, 1) * sqrt(1 - cov_yz)

    full_data <- data.frame(y, g, z)
    train_data <- full_data[1:10000, ]
    lab_data <- full_data[10001:20000, ]
    unlab_data <- full_data[20001:(20000 + n_unlab), ]

    # Fit the maching learning model
    train_fit <- lm(y ~ z, data = train_data)
    lab_data$y_hat <- predict(train_fit, newdata = lab_data)
    unlab_data$y_hat <- predict(train_fit, newdata = unlab_data)

    # Standardized the data
    train_data <- as.data.frame(scale(train_data))
    lab_data <- as.data.frame(scale(lab_data))
    unlab_data <- as.data.frame(scale(unlab_data))

    coeff_unlab_yhat_naive <- summary(lm(y_hat ~ g, data = unlab_data))$coefficient
    beta_unlab_yhat_naive <- coeff_unlab_yhat_naive[2, 1]
    z_unlab_yhat_naive <- coeff_unlab_yhat_naive[2, 3]
    p_unlab_yhat_naive <- coeff_unlab_yhat_naive[2, 4]

    coeff_lab_yhat_naive <- summary(lm(y_hat ~ g, data = lab_data))$coefficient
    beta_lab_yhat_naive <- coeff_lab_yhat_naive[2, 1]
    z_lab_yhat_naive <- coeff_lab_yhat_naive[2, 3]
    p_lab_yhat_naive <- coeff_lab_yhat_naive[2, 4]

    coeff_lab_y_naive <- summary(lm(y ~ g, data = lab_data))$coefficient
    beta_lab_y_naive <- coeff_lab_y_naive[2, 1]
    z_lab_y_naive <- coeff_lab_y_naive[2, 3]
    p_lab_y_naive <- coeff_lab_y_naive[2, 4]

    r <- cor(lab_data$y, lab_data$y_hat)

    Z_table <- data.frame(Z1 = z_lab_yhat_naive, Z2 = z_unlab_yhat_naive, Z3 = z_lab_y_naive)
      
    # pop
    pop_output <- cal_qt_opt(z_df=Z_table, n=n_lab, N=n_unlab, r=r, eaf = maf)
    beta_pop <- pop_output$BETA
    p_pop <- pop_output$P

    out <- c(beta_unlab_yhat_naive/sqrt(var_factor), p_unlab_yhat_naive,
            beta_lab_yhat_naive/sqrt(var_factor), p_lab_yhat_naive,
            beta_lab_y_naive/sqrt(var_factor), p_lab_y_naive,
            beta_pop, p_pop,
            r)
    names(out) <- c(
            "beta_unlab_yhat_naive", "p_unlab_yhat_naive",
            "beta_lab_yhat_naive", "p_lab_yhat_naive",
            "beta_lab_y_naive", "p_lab_y_naive",
            "beta_pop", "p_pop",
            "r")
    out
  }
  stopCluster(cl)

  fwrite(as.data.frame(result), paste0("./results/tmp/margin_cov_yz.power.", cov_yz, ".txt.gz"),
  sep = "\t", quote = F, row.names = F, col.names = T)
}
