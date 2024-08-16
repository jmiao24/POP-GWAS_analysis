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
h2_X_z <- 0.0015
var_res <- 0.2
n_lab <- 10^4
n_unlab <- 10^5

for (h2_X_Y in h2_X_Y_vec){
  # Basic simulation parameters
  # h2_X_Y <- 0.01
  print(h2_X_Y)
  mar_effect <- sqrt(h2_X_Y/var_factor)

  #-------------- Known Model ------
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  result <- foreach(
    i = 1:sim.times, .combine = rbind, .packages = c("doParallel"),
    .errorhandling = "pass"
  ) %dopar% {

    g <- rbinom(n.total, 2, maf)
    # 
    log_odds <- -1 + g * sqrt(mar_effect)
    prob_y <- 1 / (1 + exp(-log_odds))
    y <- rbinom(n = n.total, size = 1, prob = prob_y)
    # y_tmp <- scale(g) * sqrt(h2_X_Y)
    # y_tmp <- y_tmp + grt_err(n.total, 1) * sqrt(c(1 - var(y_tmp)))
    # y <- ifelse(y_tmp > quantile(y_tmp, 0.9), 1, 0)
    z <- scale(g) * sqrt(h2_X_z) + y * sqrt(1 - h2_X_z - var_res) + grt_err(n.total, 1) * sqrt(var_res)

    full_data <- data.frame(y, g, z)
    train_data <- full_data[1:10000, ]
    lab_data <- full_data[10001:20000, ]
    unlab_data <- full_data[20001:(20000 + n_unlab), ]

    # Fit the maching learning model
    train_fit <- glm(y ~ z, data = train_data, family = binomial)
    lab_data$y_hat <- predict(train_fit, newdata = lab_data, type="response")
    lab_data$y_hat <- ifelse(lab_data$y_hat > 0.5, 1, 0)
    unlab_data$y_hat <- predict(train_fit, newdata = unlab_data, type="response")
    unlab_data$y_hat <- ifelse(unlab_data$y_hat > 0.5, 1, 0)

    coeff_unlab_yhat_naive <- summary(glm(y_hat ~ g, data = unlab_data, family = binomial))$coefficient
    beta_unlab_yhat_naive <- coeff_unlab_yhat_naive[2, 1]
    z_unlab_yhat_naive <- coeff_unlab_yhat_naive[2, 3]
    p_unlab_yhat_naive <- coeff_unlab_yhat_naive[2, 4]

    coeff_lab_yhat_naive <- summary(glm(y_hat ~ g, data = lab_data, family = binomial))$coefficient
    beta_lab_yhat_naive <- coeff_lab_yhat_naive[2, 1]
    z_lab_yhat_naive <- coeff_lab_yhat_naive[2, 3]
    p_lab_yhat_naive <- coeff_lab_yhat_naive[2, 4]

    coeff_unlab_y_naive <- summary(glm(y ~ g, data = unlab_data, family = binomial))$coefficient
    beta_unlab_y_naive <- coeff_unlab_y_naive[2, 1]
    z_unlab_y_naive <- coeff_unlab_y_naive[2, 3]
    p_unlab_y_naive <- coeff_unlab_y_naive[2, 4]

    coeff_lab_y_naive <- summary(glm(y ~ g, data = lab_data, family = binomial))$coefficient
    beta_lab_y_naive <- coeff_lab_y_naive[2, 1]
    z_lab_y_naive <- coeff_lab_y_naive[2, 3]
    p_lab_y_naive <- coeff_lab_y_naive[2, 4]

    r <- cor(lab_data$y, lab_data$y_hat)

    Z_table <- data.frame(Z1 = z_lab_yhat_naive, Z2 = z_unlab_yhat_naive, Z3 = z_lab_y_naive)
      
    # pop
    pop_output <- cal_bt_opt(z_df=Z_table, n=n_lab, N=n_unlab, r=r, eaf = maf, case_prop = case_prop)
    beta_pop <- pop_output$BETA
    p_pop <- pop_output$P

    out <- c(beta_unlab_yhat_naive, p_unlab_yhat_naive,
            beta_unlab_y_naive, p_unlab_y_naive,
            beta_lab_yhat_naive, p_lab_yhat_naive,
            beta_lab_y_naive, p_lab_y_naive,
            beta_pop, p_pop,
            r)
    names(out) <- c("beta_unlab_yhat_naive", "p_unlab_yhat_naive",
            "beta_unlab_y_naive", "p_unlab_y_naive",
            "beta_lab_yhat_naive", "p_lab_yhat_naive",
            "beta_lab_y_naive", "p_lab_y_naive",
            "beta_pop", "p_pop",
            "r")
    out
  }
  stopCluster(cl)

  fwrite(as.data.frame(result), paste0("./results/tmp/margin_eff.bt.", mar_effect, ".txt.gz"),
  sep = "\t", quote = F, row.names = F, col.names = T)
}
