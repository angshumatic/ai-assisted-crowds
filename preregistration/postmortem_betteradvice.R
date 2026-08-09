library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)
library(glmnet)
library(broom)
library(metR)


seed <- 123
set.seed(seed)

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

responses_all <- load_responses("qualtricsdata_tr_better.csv", "tr_better")


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


tasks_all$gamma_est <- NA_real_
tasks_all$gamma_est_lasso <- NA_real_

responses_all <- responses_all %>%
  left_join(tasks_all %>% select(Task, treatment, w),
            by = c("Task", "treatment")) %>%
  mutate(res = f1 - (w * fa + (1 - w) * f0))



empirical_pivoting_lasso <- function(df) {
  x <- as.matrix(cbind(adjustment = df$adjustment, dummy = 0))
  y <- df$error_f1_bar
  
  cv_fit <- cv.glmnet(x, y, alpha = 1, nfolds = 10, intercept = FALSE)
  coef(cv_fit, s = "lambda.min")["adjustment", 1] 
}

for (tr in unique(tasks_all$treatment)) {
  for (i in 1:tasks) {
    holdout_data_tasks <- tasks_all %>% 
      filter(treatment == tr, Task != i)
    holdout_data_responses <- responses_all %>% 
      filter(treatment == tr, Task != i)
    
    tasks_all$gamma_est_lasso[
      tasks_all$treatment == tr & tasks_all$Task == i
    ] <- empirical_pivoting_lasso(holdout_data_tasks)
    
  }
}
rm(i, tr, holdout_data_tasks, holdout_data_responses)


tasks_all <- tasks_all %>%
  mutate(gpe_est = f1_bar + gamma_est * (f0_bar - f1_bar)) %>% 
  mutate(gpe_est_lasso = f1_bar + gamma_est_lasso * (f0_bar - f1_bar)) %>% 
  mutate(error_gpe_est_lasso = gpe_est_lasso - theta)

rmse_gpe <- sqrt(mean((tasks_all$error_gpe_est_lasso)^2))
rmse_f0bar <- sqrt(mean((tasks_all$error_f0_bar)^2))
rmse_fa <- sqrt(mean((tasks_all$error_fa)^2))
rmse_f1bar <- sqrt(mean((tasks_all$error_f1_bar)^2))

t.test((tasks_all$error_gpe_est_lasso)^2, (tasks_all$error_f0_bar)^2, paired = TRUE)
t.test((tasks_all$error_gpe_est_lasso)^2, (tasks_all$error_fa)^2, paired = TRUE)
t.test((tasks_all$error_gpe_est_lasso)^2, (tasks_all$error_f1_bar)^2, paired = TRUE)
