library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)
library(broom)
library(metR)
library(extraDistr)


seed <- 123
set.seed(seed)

tasks <- 200


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
      treatment = treatment_name
    )
  return(responses_df)
}

responses_all <- load_responses("qualtricsdata_baseline.csv", "baseline")


summarize_tasks <- function(responses_df) {
  tasks_df <- responses_df %>%
    group_by(Task) %>%
    summarise(
      treatment = first(treatment),
      J = n(),
      theta = first(theta),
      f0_bar = mean(f0, na.rm = TRUE),
      fa = first(fa),
      f1_bar = mean(f1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      error_f1_bar = f1_bar - theta,
      error_f0_bar = f0_bar - theta,
      error_fa = fa - theta
    )
  return(tasks_df)
}

tasks_all <- summarize_tasks(responses_all)

tasks_all$sigma_0_sq <- NA
tasks_all$sigma_x_sq <- NA
tasks_all$w_est <- NA
tasks_all$w_std_err <- NA
tasks_all$w_df <- NA

for (i in 1:tasks) {
  filtered_rows <- responses_all %>% 
    filter(Task == i)
  holdout_tasks <- tasks_all %>% 
    filter(Task != i)
  
  model_ols <- lm(f1_minus_f0 ~ fa_minus_f0 - 1, data = filtered_rows)
  model_ols_summary <- summary(model_ols)
  
  tasks_all$w_est[i] <- model_ols_summary$coefficients[1, "Estimate"]
  tasks_all$w_std_err[i] <- model_ols_summary$coefficients[1, "Std. Error"]
  tasks_all$w_df[i] <- model_ols_summary$df[2]
  
  tasks_all$sigma_x_sq[i] <- var(holdout_tasks$error_fa)
  tasks_all$sigma_0_sq[i] <- var(filtered_rows$f0)

}

tasks_all <- tasks_all %>% 
  mutate(
    L_w = ((J - 1)/J * (sigma_0_sq/sigma_x_sq)) / 
      (1 + (J - 1)/J * sigma_0_sq/sigma_x_sq),
    U_w = 1
  ) %>% 
  mutate(
    L_w_tstat = (L_w - w_est) / w_std_err,
    U_w_tstat = (U_w - w_est) / w_std_err
  ) %>% 
  mutate(
    area_w = pt(U_w_tstat, df = w_df) - pt(L_w_tstat, df = w_df)
  )

pdf_tdist <- function(x, location, scale, df) {
  dlst(x, df, mu = location, sigma = scale, log = FALSE)
}

function_gammastar <- function(w, sigma_0_sq, sigma_x_sq, J){
  1 - 
    (w/(1-w) - (J - 1)/J * (sigma_0_sq/sigma_x_sq)) / 
    w / (1 + w/(1-w) - (J - 1)/J * sigma_0_sq/sigma_x_sq)
}

integrand_gamma <- function(w, location_w, scale_w, df_w, sigma_0_sq, sigma_x_sq, J) {
  pdf_tdist(w, location_w, scale_w, df_w) * function_gammastar(w, sigma_0_sq, sigma_x_sq, J)
}

tasks_all <- tasks_all %>%
  rowwise() %>%
  mutate(
    gamma_est = integrate(
      integrand_gamma,
      lower = L_w, 
      upper = U_w, 
      location_w = w_est, 
      scale_w = w_std_err, 
      df_w = w_df,
      sigma_0_sq = sigma_0_sq,
      sigma_x_sq = sigma_x_sq,
      J = J
    )$value / area_w
  ) %>%
  ungroup()


tasks_all <- tasks_all %>% 
  mutate(
    est_gpe = f1_bar + gamma_est * (f0_bar - f1_bar),
    error_est_gpe = est_gpe - theta
  )


mae_gpe <- mean(abs(tasks_all$error_est_gpe))
mae_f0bar <- mean(abs(tasks_all$error_f0_bar))
mae_fa <- mean(abs(tasks_all$error_fa))
mae_f1bar <- mean(abs(tasks_all$error_f1_bar))

t.test(abs(tasks_all$error_est_gpe), abs(tasks_all$error_f0_bar), paired = TRUE)
t.test(abs(tasks_all$error_est_gpe), abs(tasks_all$error_fa), paired = TRUE)
t.test(abs(tasks_all$error_est_gpe), abs(tasks_all$error_f1_bar), paired = TRUE)

