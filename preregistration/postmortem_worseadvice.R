library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(readr)
library(stringr)
library(glmnet)
library(broom)
library(metR)
library(lme4)
library(performance)


tasks <- 200

w_est <- function(df) {
  num <- sum((df$f1 - df$f0) * (df$fa - df$f0), na.rm = TRUE)
  denominator <- sum((df$fa - df$f0)^2, na.rm = TRUE)
  return(num / denominator)
}

load_responses <- function(file_path, treatment_name, min_duration = 180) {
  responses_df <- read_csv(file_path, show_col_types = FALSE) %>%
    mutate_all(as.numeric) %>%
    filter(duration >= min_duration) %>%
    pivot_longer(
      cols = starts_with(c("priorQ", "updatedQ", "Advice", "True", "SlNoQ")),
      names_to = c(".value", "Question"),
      names_pattern = "(\\D+)(\\d+)"
    ) %>%
    rename(
      Judge = Respondent,
      Duration = duration,
      Code = code,
      f0 = priorQ,
      f1 = updatedQ,
      fa = Advice,
      theta = True,
      Task = SlNoQ
    ) %>%
    mutate(
      duration = Duration,
      error_f0 = f0 - theta,
      error_fa = fa - theta,
      error_f1 = f1 - theta,
      f1_minus_f0 = f1 - f0,
      fa_minus_f0 = fa - f0,
      error_f0_abs = abs(f0 - theta),
      error_fa_abs = abs(fa - theta),
      error_f1_abs = abs(f1 - theta),
      treatment = treatment_name
    )
  return(responses_df)
}

responses_all <- load_responses("qualtricsdata_tr_worse.csv", "tr_worse")

summarize_tasks <- function(responses_df) {
  tasks_df <- responses_df %>%
    group_by(Task) %>%
    mutate(w = w_est(cur_data())) %>%
    summarise(
      treatment = first(treatment),
      J = n(),
      theta = first(theta),
      f0_bar = mean(f0, na.rm = TRUE),
      fa = first(fa),
      f1_bar = mean(f1, na.rm = TRUE),
      w = first(w),
      .groups = "drop"
    ) %>%
    mutate(
      error_f1_bar = f1_bar - theta,
      error_f0_bar = f0_bar - theta,
      error_fa = fa - theta,
      adjustment = f1_bar - f0_bar
    )
  return(tasks_df)
}

tasks_all <- summarize_tasks(responses_all)

responses_all <- responses_all %>%
  left_join(tasks_all %>% select(Task, treatment, w),
            by = c("Task", "treatment")) %>%
  mutate(res = f1 - (w * fa + (1 - w) * f0))




n_reps <- 5000

parameter_estimates <- data.frame(
  Task = integer(), 
  treatment = character(),
  rho_ha_est = numeric(), rho_ha_min = numeric(), rho_ha_max = numeric(),
  rho_hh_est = numeric(), rho_hh_min = numeric(), rho_hh_max = numeric(),
  phi_0_est = numeric(), phi_0_min = numeric(), phi_0_max = numeric(),
  alpha_est = numeric(), alpha_min = numeric(), alpha_max = numeric()
)


sample_correlation_ha <- function(df) {
  sampled_rows <- df %>%
    group_by(Task) %>%
    slice_sample(n = 1) %>%
    ungroup()
  rho_ha <- cor(sampled_rows$error_f0, sampled_rows$error_fa)
  return(rho_ha)
}
sample_correlation_hh <- function(df) {
  model <- lmer(error_f0 ~ 1 + 
                  (1 | Task), data = df)
  s <- icc(model, ci = TRUE)
  est <- c(s$ICC_adjusted[1], s$ICC_adjusted[2], s$ICC_adjusted[3])
  return(est)
}
sample_ratio <- function(df) {
  sampled_rows <- df %>%
    group_by(Task) %>%
    slice_sample(n = 1) %>%
    ungroup()
  var_f0 <- mean(sampled_rows$error_f0^2, na.rm = TRUE)
  var_fa <- mean(sampled_rows$error_fa^2, na.rm = TRUE)
  est <- sqrt(var_fa / var_f0)
  return(est)
}

### To run the original preregistered version without residuals, modify the function given below in Line 134/135.
### For details, check postmortem report.

compute_sse <- function(alpha, df) {
  w_star <- df$w / alpha
  y_hat <- w_star * df$fa + (1 - w_star) * df$f0 + df$res ### RESIDUALS INCLUDED
  #y_hat <- w_star * df$fa + (1 - w_star) * df$f0 ### RESIDUALS NOT INCLUDED
  sum((df$theta - y_hat)^2)
}
sample_aversion <- function(df) {
  sampled_rows <- df %>%
    group_by(Task) %>%
    slice_sample(n = 1) %>%
    ungroup() 
  optimize(compute_sse, interval = c(1e-6, 1.5), df = sampled_rows)$minimum
}

seed <- 1
set.seed(seed)

for (tr in unique(tasks_all$treatment)) {
  for (i in 1:tasks) {
    holdout_data <- responses_all %>% filter(treatment == tr, Task != i)
    
    correlation_ha <- replicate(n_reps, sample_correlation_ha(holdout_data))
    q <- quantile(correlation_ha, probs = c(0.025, 0.975))
    m <- mean(correlation_ha)
    rho_ha_vals <- c(m, q[1], q[2])
    print(paste("Task", i, "rho_ha estimation completed"))
    
    correlation_hh <- sample_correlation_hh(holdout_data)
    rho_hh_vals <- c(correlation_hh[1], correlation_hh[2], correlation_hh[3])
    print(paste("Task", i, "rho_hh estimation completed"))
    
    phi_0 <- replicate(n_reps, sample_ratio(holdout_data))
    q <- quantile(phi_0, probs = c(0.025, 0.975))
    m <- mean(phi_0)
    phi_0_vals <- c(m, q[1], q[2])
    print(paste("Task", i, "phi_0 estimation completed"))
    
    alpha <- replicate(n_reps, sample_aversion(holdout_data))
    q <- quantile(alpha, probs = c(0.025, 0.975))
    m <- mean(alpha)
    alpha_vals <- c(m, q[1], q[2])
    print(paste("Task", i, "alpha estimation completed"))
    
    
    parameter_estimates <- rbind(parameter_estimates, data.frame(
      Task = i, treatment = tr,
      rho_ha_est = rho_ha_vals[1], rho_ha_min = rho_ha_vals[2], rho_ha_max = rho_ha_vals[3],
      rho_hh_est = rho_hh_vals[1], rho_hh_min = rho_hh_vals[2], rho_hh_max = rho_hh_vals[3],
      phi_0_est = phi_0_vals[1], phi_0_min = phi_0_vals[2], phi_0_max = phi_0_vals[3],
      alpha_est = alpha_vals[1], alpha_min = alpha_vals[2], alpha_max = alpha_vals[3]
    ))
    
    print(paste("Task", i, "completed for treatment ", tr))
  }
}

## parameter_estimates.csv is provided in the package.
## The user can directly read the postmortem_worseadvice_estimates.csv file and proceed to the next step.

tasks_all <- tasks_all %>%
  left_join(parameter_estimates, by = c("Task", "treatment" = "treatment"))

sigmaxsq <- function(phi_0, rho_ha, rho_hh) {
  phi_0/rho_ha - 1
}
sigmaysq <- function(phi_0, rho_ha, rho_hh) {
  rho_hh/(phi_0 * rho_ha) - 1
}
sigma0sq <- function(phi_0, rho_ha, rho_hh) {
  (1 - rho_hh) / (phi_0 * rho_ha)
}
w_opt <- function(sigmaxsq, sigmaysq, sigma0sq) {
  (sigmaysq + sigma0sq) / (sigmaxsq + sigmaysq + sigma0sq)
}
w <- function(alpha, w_opt){
  alpha * w_opt
}
error_var <- function(sigmaxsq, sigmaysq, sigma0sq, gamma, w, J) {
  sigmaxsq * (w * (1-gamma))^2 +
    (sigmaysq + sigma0sq/J) * (1 - w * (1-gamma))^2
}

error_var_min <- function(sigmaxsq, sigmaysq, sigma0sq, J) {
  sigmaxsq * (sigmaysq + sigma0sq/J) / (sigmaxsq + sigmaysq + sigma0sq/J)
}


gran <- 101

robust <- function(sub_df) {
  sigmaxsq   <- sub_df$sigmaxsq
  sigmaysq   <- sub_df$sigmaysq
  sigma0sq   <- sub_df$sigma0sq
  w_vec      <- sub_df$w
  J_vec      <- sub_df$J
  error_min  <- sub_df$error_var_min
  
  max_regret <- function(gamma) {
    regret_gamma <- error_var(sigmaxsq, sigmaysq, sigma0sq, gamma, w_vec, J_vec) - error_min
    max(regret_gamma)
  }
  
  opt <- optimize(f = max_regret, interval = c(-1, 1))
  return(opt$minimum)
}

get_gamma_robust <- function(rho_ha_min, rho_ha_max, rho_hh_min, rho_hh_max, 
                             phi_0_min, phi_0_max, alpha_min, alpha_max, J) {
  
  rho_ha_vals <- seq(rho_ha_min, rho_ha_max, length.out = gran)
  rho_hh_vals <- seq(rho_hh_min, rho_hh_max, length.out = gran)
  phi_0_vals <- seq(phi_0_min, phi_0_max, length.out = gran)
  alpha_vals   <- seq(alpha_min, alpha_max, length.out = gran)
  
  df <- expand.grid(rho_ha = rho_ha_vals,
                    rho_hh = rho_hh_vals,
                    phi_0 = phi_0_vals,
                    alpha  = alpha_vals,
                    J     = J)
  
  df <- df %>%
    filter(phi_0 > rho_ha, rho_hh > phi_0 * rho_ha) %>%
    mutate(
      sigmaxsq = sigmaxsq(phi_0, rho_ha, rho_hh),
      sigmaysq = sigmaysq(phi_0, rho_ha, rho_hh),
      sigma0sq = sigma0sq(phi_0, rho_ha, rho_hh),
      w_opt    = w_opt(sigmaxsq, sigmaysq, sigma0sq),
      w        = w(alpha, w_opt),
      error_var_min = error_var_min(sigmaxsq, sigmaysq, sigma0sq, J)
    )
  
  robust(df)
}

tasks_all <- tasks_all %>% 
  rowwise() %>%
  mutate(gamma_robust = get_gamma_robust(
    rho_ha_min, rho_ha_max,
    rho_hh_min, rho_hh_max,
    phi_0_min,  phi_0_max,
    alpha_min,   alpha_max,
    J
  )) %>%
  ungroup()


tasks_all <- tasks_all %>% 
  mutate(
    est_gpe = f1_bar + gamma_robust * (f0_bar - f1_bar),
    error_est_gpe = est_gpe - theta
  )

mae_gpe <- mean(abs(tasks_all$error_est_gpe))
mae_f0bar <- mean(abs(tasks_all$error_f0_bar))
mae_fa <- mean(abs(tasks_all$error_fa))
mae_f1bar <- mean(abs(tasks_all$error_f1_bar))

t.test(abs(tasks_all$error_est_gpe), abs(tasks_all$error_f0_bar), paired = TRUE)
t.test(abs(tasks_all$error_est_gpe), abs(tasks_all$error_fa), paired = TRUE)
t.test(abs(tasks_all$error_est_gpe), abs(tasks_all$error_f1_bar), paired = TRUE)

