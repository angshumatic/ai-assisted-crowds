# Run this script from the repository root.
# It creates the main analysis objects and writes the paper figures as PDFs.

suppressPackageStartupMessages({
  library(broom)
  library(dplyr)
  library(ggplot2)
  library(glmnet)
  library(readr)
  library(scales)
  library(tidyr)
})

set.seed(123)

n_tasks <- 200
min_duration <- 180

treatment_files <- c(
  baseline = "data/qualtricsdata_baseline.csv",
  tr_worse = "data/qualtricsdata_tr_worse.csv",
  tr_better = "data/qualtricsdata_tr_better.csv"
)

figures_dir <- "figures"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

treatment_levels <- c("tr_worse", "baseline", "tr_better")
treatment_labels <- c(
  tr_worse = "Worse Advice",
  baseline = "Baseline",
  tr_better = "Better Advice"
)


################
### Functions
################

w_est <- function(df) {
  num <- sum((df$f1 - df$f0) * (df$fa - df$f0), na.rm = TRUE)
  den <- sum((df$fa - df$f0)^2, na.rm = TRUE)
  num / den
}

load_responses <- function(file, treatment) {
  read_csv(file, show_col_types = FALSE) %>%
    mutate(across(everything(), as.numeric)) %>%
    filter(duration >= min_duration) %>%
    pivot_longer(
      cols = matches("^(priorQ|updatedQ|Advice|True|SlNoQ)[0-9]+$"),
      names_to = c(".value", "Question"),
      names_pattern = "(priorQ|updatedQ|Advice|True|SlNoQ)([0-9]+)"
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
      error_f0_abs = abs(error_f0),
      error_fa_abs = abs(error_fa),
      error_f1_abs = abs(error_f1),
      treatment = treatment
    )
}

summarize_tasks <- function(df) {
  df %>%
    group_by(Task) %>%
    summarise(
      treatment = first(treatment),
      J = n(),
      theta = first(theta),
      f0_bar = mean(f0, na.rm = TRUE),
      fa = first(fa),
      f1_bar = mean(f1, na.rm = TRUE),
      w = w_est(pick(f0, f1, fa)),
      .groups = "drop"
    ) %>%
    mutate(
      error_f1_bar = f1_bar - theta,
      error_f0_bar = f0_bar - theta,
      error_fa = fa - theta,
      adjustment = f1_bar - f0_bar
    )
}

participant_summary <- function(df) {
  df %>%
    group_by(treatment, Judge) %>%
    summarise(
      duration = first(duration),
      rmse_f0 = sqrt(mean((f0 - theta)^2, na.rm = TRUE)),
      rmse_f1 = sqrt(mean((f1 - theta)^2, na.rm = TRUE)),
      w = w_est(pick(f0, f1, fa)),
      .groups = "drop"
    )
}

empirical_pivoting <- function(df) {
  unname(coef(lm(error_f1_bar ~ adjustment - 1, data = df))[["adjustment"]])
}

empirical_pivoting_lasso <- function(df) {
  x <- as.matrix(cbind(adjustment = df$adjustment, dummy = 0))
  fit <- cv.glmnet(x, df$error_f1_bar, alpha = 1, nfolds = 10, intercept = FALSE)
  unname(as.numeric(coef(fit, s = "lambda.min")["adjustment", 1]))
}

compute_sse_alpha <- function(alpha, df) {
  w_star <- df$w / alpha
  y_hat <- w_star * df$fa + (1 - w_star) * df$f0 + df$res
  sum((df$theta - y_hat)^2)
}

alpha_estimation <- function(df) {
  optimize(compute_sse_alpha, interval = c(1e-6, 1.5), df = df)$minimum
}

compute_sse_w_star <- function(w_star, df) {
  y_hat <- w_star * df$fa + (1 - w_star) * df$f0
  sum((df$theta - y_hat)^2)
}

w_star_estimation <- function(df) {
  optimize(compute_sse_w_star, interval = c(1e-6, 1.5), df = df)$minimum
}

compute_sse_beta <- function(beta, w_star, df) {
  y_hat <- w_star * beta * df$fa + (1 - w_star * beta) * df$f0_bar
  sum((df$theta - y_hat)^2)
}

beta_estimation <- function(w_star, df) {
  optimize(compute_sse_beta, interval = c(1e-6, 1.5), w_star = w_star, df = df)$minimum
}

estimate_parameters <- function(tasks_df, responses_df) {
  tasks_df$gamma_est <- NA_real_
  tasks_df$gamma_est_lasso <- NA_real_
  tasks_df$alpha_est <- NA_real_
  tasks_df$beta_est <- NA_real_

  for (tr in unique(tasks_df$treatment)) {
    for (task in seq_len(n_tasks)) {
      holdout_tasks <- filter(tasks_df, treatment == tr, Task != task)
      holdout_responses <- filter(responses_df, treatment == tr, Task != task)
      row <- tasks_df$treatment == tr & tasks_df$Task == task

      tasks_df$gamma_est[row] <- empirical_pivoting(holdout_tasks)
      tasks_df$gamma_est_lasso[row] <- empirical_pivoting_lasso(holdout_tasks)

      w_star <- w_star_estimation(holdout_responses)
      tasks_df$alpha_est[row] <- alpha_estimation(holdout_responses)
      tasks_df$beta_est[row] <- beta_estimation(w_star, holdout_tasks)
    }
  }

  tasks_df
}

forecast_tests <- function(tasks_df) {
  pairs <- list(
    c("error_gpe_est_lasso", "error_f0_bar"),
    c("error_gpe_est_lasso", "error_f1_bar"),
    c("error_gpe_est_lasso", "error_fa")
  )

  out <- list()

  for (tr in unique(tasks_df$treatment)) {
    tr_df <- filter(tasks_df, treatment == tr)

    for (pair in pairs) {
      test <- tidy(t.test(tr_df[[pair[1]]]^2, tr_df[[pair[2]]]^2, paired = TRUE))
      out[[length(out) + 1]] <- test %>%
        mutate(treatment = tr, var1 = pair[1], var2 = pair[2]) %>%
        select(treatment, var1, var2, estimate, statistic, p.value)
    }
  }

  bind_rows(out)
}

make_segments <- function() {
  data.frame(x = c(1, 0), y = c(0, 1), xend = c(1, 1), yend = c(1, 1))
}

make_contours <- function(gammas) {
  contours <- do.call(rbind, lapply(gammas, function(g) {
    if (g > 0) {
      data.frame(x = 0, y = 0, xend = 1 - g, yend = 1)
    } else {
      data.frame(x = 0, y = 0, xend = 1, yend = 1 / (1 - g))
    }
  }))

  contours$gamma <- gammas
  contours
}

make_contour_labels <- function(contours, label_type = c("gamma", "kappa")) {
  label_type <- match.arg(label_type)
  x_mult <- ifelse(label_type == "gamma", 0.85, 0.75)
  label <- if (label_type == "gamma") {
    paste0("tilde(gamma) == ", contours$gamma)
  } else {
    benefit <- 1 - contours$gamma^2 / (1 - contours$gamma)^2
    paste0("tilde(kappa)~'='~", sprintf("%.0f", 100 * benefit), "*'%'")
  }

  contours %>%
    mutate(
      angle = atan2(yend - y, xend - x) * 180 / pi,
      label_x = ifelse(gamma >= 0, x_mult * (1 - gamma) - 0.04, x_mult),
      label_y = ifelse(gamma >= 0, x_mult, x_mult / (1 - gamma) + 0.03),
      label = label,
      x = label_x,
      y = label_y
    )
}

geom_abline_clipped <- function(slope, intercept,
                                xlim = c(0, 1), ylim = c(0, 1),
                                color = "black", linetype = "solid",
                                linewidth = 1) {
  pts <- data.frame(
    x = c(xlim, (ylim - intercept) / slope),
    y = c(slope * xlim + intercept, ylim)
  )
  pts <- subset(pts, x >= xlim[1] & x <= xlim[2] & y >= ylim[1] & y <= ylim[2])

  line_df <- data.frame(
    x = pts$x[1],
    y = pts$y[1],
    xend = pts$x[2],
    yend = pts$y[2]
  )

  geom_segment(
    data = line_df,
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    color = color,
    linetype = linetype,
    linewidth = linewidth
  )
}

make_gamma_contourplot <- function(segments, contours) {
  ggplot() +
    geom_segment(data = segments, aes(x = x, y = y, xend = xend, yend = yend),
                 color = "black", linetype = "dashed", linewidth = 0.5) +
    geom_segment(data = contours, aes(x = x, y = y, xend = xend, yend = yend),
                 color = "black", linewidth = 0.4) +
    geom_abline_clipped(1, 0, linewidth = 0.8) +
    geom_text(data = make_contour_labels(contours, "gamma"),
              aes(x = x, y = y, label = label, angle = angle),
              parse = TRUE, hjust = 0.5, vjust = 0.5, size = 4) +
    annotate("label", x = 0.3, y = 0.65, label = "Pivot back toward\ninitial forecast.",
             hjust = 0.5, vjust = 0, size = 4, label.size = 0.3,
             label.r = grid::unit(0.15, "lines")) +
    annotate("label", x = 0.6, y = 0.25, label = "Pivot away from\ninitial forecast.",
             hjust = 0.5, vjust = 0, size = 4, label.size = 0.3,
             label.r = grid::unit(0.15, "lines")) +
    annotate("text", x = 0.58, y = 0.62, label = "Two biases cancel out.",
             angle = 45, hjust = 0.5, vjust = 0.5, size = 4) +
    scale_x_continuous(limits = c(0, 1.01), breaks = pretty_breaks(n = 5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.05), breaks = pretty_breaks(n = 5), expand = c(0, 0)) +
    labs(x = expression(beta), y = expression(alpha)) +
    annotate("text", x = 0.75, y = 0.03, label = "\u2190 Increasing Borg Effect",
             hjust = 0.5, vjust = 0.5, size = 4, fontface = "bold") +
    annotate("text", x = 0.02, y = 0.65, label = "\u2190 Increasing Algorithm Aversion",
             angle = 90, hjust = 0.5, vjust = 0.5, size = 4, fontface = "bold") +
    annotate("text", x = 0.97, y = 0.5, label = "No Borg Effect",
             angle = 90, hjust = 0.5, vjust = 0.5, size = 4, fontface = "italic") +
    annotate("text", x = 0.5, y = 1.03, label = "Appropriate Use of AI Advice",
             hjust = 0.5, vjust = 0.5, size = 4, fontface = "italic") +
    theme_bw(base_size = 12) +
    theme(
      aspect.ratio = 1,
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank(),
      axis.ticks = element_line(linewidth = 0.4),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(face = "bold"),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.line.x.bottom = element_line(),
      axis.line.y.left = element_line()
    )
}

make_benefit_contourplot <- function(segments, contours) {
  shade_df <- data.frame(x = c(0, 0, 0.5), y = c(0, 1, 1))

  ggplot() +
    geom_polygon(data = shade_df, aes(x = x, y = y), fill = "grey90", color = NA) +
    geom_segment(data = segments, aes(x = x, y = y, xend = xend, yend = yend),
                 color = "black", linetype = "dashed", linewidth = 0.5) +
    geom_segment(data = contours, aes(x = x, y = y, xend = xend, yend = yend),
                 color = "black", linewidth = 0.4) +
    geom_abline_clipped(1, 0, linewidth = 0.8) +
    geom_text(data = make_contour_labels(contours, "kappa"),
              aes(x = x, y = y, label = label, angle = angle),
              parse = TRUE, hjust = 0.5, vjust = 0.5, size = 4) +
    scale_x_continuous(limits = c(0, 1.01), breaks = pretty_breaks(n = 5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.05), breaks = pretty_breaks(n = 5), expand = c(0, 0)) +
    labs(x = expression(beta), y = expression(alpha)) +
    annotate("text", x = 0.75, y = 0.03, label = "\u2190 Increasing Borg Effect",
             hjust = 0.5, vjust = 0.5, size = 4, fontface = "bold") +
    annotate("text", x = 0.02, y = 0.65, label = "\u2190 Increasing Algorithm Aversion",
             angle = 90, hjust = 0.5, vjust = 0.5, size = 4, fontface = "bold") +
    annotate("text", x = 0.97, y = 0.5, label = "No Borg Effect",
             angle = 90, hjust = 0.5, vjust = 0.5, size = 4, fontface = "italic") +
    annotate("text", x = 0.5, y = 1.03, label = "Appropriate Use of AI Advice",
             hjust = 0.5, vjust = 0.5, size = 4, fontface = "italic") +
    theme_bw(base_size = 12) +
    theme(
      aspect.ratio = 1,
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank(),
      axis.ticks = element_line(linewidth = 0.4),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(face = "bold"),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.line.x.bottom = element_line(),
      axis.line.y.left = element_line()
    )
}

make_scatterplot <- function(tasks_df, segments) {
  diagonal <- data.frame(x = 0, y = 0, xend = 1, yend = 1)

  plot <- ggplot(tasks_df, aes(x = beta_est, y = alpha_est)) +
    geom_point(aes(shape = treatment), color = "black", size = 2.2, stroke = 0.4) +
    geom_segment(data = segments, aes(x = x, y = y, xend = xend, yend = yend),
                 color = "black", linetype = "dashed", linewidth = 0.5) +
    geom_segment(data = diagonal, aes(x = x, y = y, xend = xend, yend = yend),
                 inherit.aes = FALSE, color = "black", linetype = "dashed",
                 linewidth = 0.5) +
    scale_shape_manual(values = c(tr_better = 24, baseline = 21, tr_worse = 25),
                       labels = treatment_labels, name = NULL) +
    scale_x_continuous(limits = c(0, 1.01), breaks = pretty_breaks(n = 5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.05), breaks = pretty_breaks(n = 5), expand = c(0, 0)) +
    labs(x = expression(hat(beta)), y = expression(hat(alpha))) +
    annotate("label", x = 0.25, y = 0.8, label = "tilde(gamma) > 0",
             parse = TRUE, hjust = 0.5, vjust = 0, size = 4, fontface = "bold",
             label.size = 0.3, label.r = grid::unit(0.15, "lines")) +
    annotate("label", x = 0.75, y = 0.2, label = "tilde(gamma) < 0",
             parse = TRUE, hjust = 0.5, vjust = 0, size = 4, fontface = "bold",
             label.size = 0.3, label.r = grid::unit(0.15, "lines")) +
    annotate("text", x = 0.47, y = 0.53,
             label = "alpha == beta ~ (tilde(gamma) == 0)",
             parse = TRUE, angle = 45, hjust = 0.5, vjust = 0.5,
             size = 4, fontface = "italic") +
    annotate("text", x = 0.73, y = 0.03, label = "\u2190 Increasing Borg Effect",
             hjust = 0.5, vjust = 0.5, size = 4, fontface = "bold") +
    annotate("text", x = 0.02, y = 0.63, label = "\u2190 Increasing Algorithm Aversion",
             angle = 90, hjust = 0.5, vjust = 0.5, size = 4, fontface = "bold") +
    annotate("text", x = 0.97, y = 0.5, label = "No Borg Effect",
             angle = 90, hjust = 0.5, vjust = 0.5, size = 4, fontface = "italic") +
    annotate("text", x = 0.5, y = 1.03, label = "Appropriate Use of AI Advice",
             hjust = 0.5, vjust = 0.5, size = 4, fontface = "italic") +
    theme_bw(base_size = 12) +
    theme(
      aspect.ratio = 1,
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank(),
      axis.ticks = element_line(linewidth = 0.4),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(face = "bold"),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.line.x.bottom = element_line(),
      axis.line.y.left = element_line(),
      #legend.position = "bottom"
      legend.position = "inside",
      legend.position.inside = c(0.10, 0.3),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = "white", colour = "black"),
      legend.box.background = element_blank(),
      legend.text = element_text(size = 10)
    )
}

make_gamma_hist <- function(tasks_df, stats_df) {
  ggplot(tasks_df, aes(x = gamma_est)) +
    geom_histogram(bins = 180, fill = "white", color = "black", alpha = 0.6) +
    scale_x_continuous(limits = c(-1.1, 1.1), breaks = pretty_breaks(n = 5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 200), expand = c(0, 0)) +
    labs(x = expression(hat(gamma)), y = "Frequency") +
    geom_segment(x = 0, y = 0, xend = 0, yend = 200,
                 color = "black", linetype = "dashed", linewidth = 0.5) +
    annotate("text", x = stats_df$gamma_mean[stats_df$treatment == "tr_worse"],
             y = 150, label = "Worse\nAdvice", hjust = 0.5, size = 4, lineheight = 0.7) +
    annotate("text", x = stats_df$gamma_mean[stats_df$treatment == "baseline"] + 0.05,
             y = 150, label = "Baseline", hjust = 0.5, size = 4) +
    annotate("text", x = stats_df$gamma_mean[stats_df$treatment == "tr_better"],
             y = 150, label = "Better\nAdvice", hjust = 0.5, size = 4, lineheight = 0.7) +
    theme_bw(base_size = 12) +
    theme(
      aspect.ratio = 0.5,
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.line.x.bottom = element_line(),
      axis.line.y.left = element_line(),
      panel.grid = element_blank(),
      axis.ticks = element_line(linewidth = 0.4),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(face = "bold"),
      legend.position = "none"
    )
}

make_kappa_hist <- function(tasks_df, stats_df) {
  neg10_trans <- trans_new(
    "neg10",
    transform = function(x) ifelse(x < 0, x / 1000, x / 100),
    inverse = function(x) ifelse(x < 0, x * 1000, x * 100)
  )

  plot_df <- tasks_df %>%
    mutate(kappa_percent = 100 * kappa_est)

  ggplot(plot_df, aes(x = kappa_percent)) +
    geom_histogram(bins = 180, fill = "white", color = "black", alpha = 0.6) +
    scale_x_continuous(
      trans = neg10_trans,
      breaks = c(-1000, 0, 100),
      labels = function(x) paste0(x, "%"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(limits = c(0, 200), expand = c(0, 0)) +
    coord_cartesian(xlim = c(-1100, 130)) +
    labs(x = expression(hat(kappa) ~ "(log scale)"), y = "Frequency") +
    geom_segment(x = 100, y = 0, xend = 100, yend = 200,
                 color = "black", linetype = "dashed", linewidth = 0.5) +
    annotate("text", x = 100 * stats_df$kappa_mean[stats_df$treatment == "tr_worse"],
             y = 140, label = "Worse\nAdvice", hjust = 0.5, size = 4, lineheight = 0.7) +
    annotate("text", x = 100 * stats_df$kappa_mean[stats_df$treatment == "baseline"] + 5,
             y = 190, label = "Baseline", hjust = 0.5, size = 4) +
    annotate("text", x = 100 * stats_df$kappa_mean[stats_df$treatment == "tr_better"],
             y = 165, label = "Better\nAdvice", hjust = 0.5, size = 4, lineheight = 0.7) +
    theme_bw(base_size = 12) +
    theme(
      aspect.ratio = 0.5,
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.line.x.bottom = element_line(),
      axis.line.y.left = element_line(),
      panel.grid = element_blank(),
      axis.ticks = element_line(linewidth = 0.4),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(face = "bold"),
      legend.position = "none"
    )
}


################
### Analysis
################

responses_all <- bind_rows(lapply(names(treatment_files), function(tr) {
  load_responses(treatment_files[[tr]], tr)
}))

participant_summary <- participant_summary(responses_all)

tasks_all <- bind_rows(lapply(names(treatment_files), function(tr) {
  summarize_tasks(filter(responses_all, treatment == tr))
}))

responses_all <- responses_all %>%
  left_join(select(tasks_all, Task, treatment, w), by = c("Task", "treatment")) %>%
  mutate(res = f1 - (w * fa + (1 - w) * f0))

tasks_all <- estimate_parameters(tasks_all, responses_all) %>%
  mutate(
    kappa_est = 1 - gamma_est^2 / (1 - gamma_est)^2,
    gpe_est = f1_bar + gamma_est * (f0_bar - f1_bar),
    gpe_est_lasso = f1_bar + gamma_est_lasso * (f0_bar - f1_bar),
    error_gpe_est = gpe_est - theta,
    error_gpe_est_lasso = gpe_est_lasso - theta
  )

desc_stats <- tasks_all %>%
  group_by(treatment) %>%
  summarise(
    alpha_mean = mean(alpha_est, na.rm = TRUE),
    alpha_min = min(alpha_est, na.rm = TRUE),
    alpha_max = max(alpha_est, na.rm = TRUE),
    alpha_sd = sd(alpha_est, na.rm = TRUE),
    beta_mean = mean(beta_est, na.rm = TRUE),
    beta_min = min(beta_est, na.rm = TRUE),
    beta_max = max(beta_est, na.rm = TRUE),
    beta_sd = sd(beta_est, na.rm = TRUE),
    gamma_mean = mean(gamma_est, na.rm = TRUE),
    gamma_min = min(gamma_est, na.rm = TRUE),
    gamma_max = max(gamma_est, na.rm = TRUE),
    gamma_sd = sd(gamma_est, na.rm = TRUE),
    kappa_mean = mean(kappa_est, na.rm = TRUE),
    kappa_min = min(kappa_est, na.rm = TRUE),
    kappa_max = max(kappa_est, na.rm = TRUE),
    kappa_sd = sd(kappa_est, na.rm = TRUE),
    gamma_LASSO_mean = mean(gamma_est_lasso, na.rm = TRUE),
    w_mean = mean(w, na.rm = TRUE),
    .groups = "drop"
  )

forecast_rmse <- tasks_all %>%
  select(treatment, error_f1_bar, error_f0_bar, error_fa, error_gpe_est_lasso) %>%
  pivot_longer(-treatment, names_to = "metric", values_to = "error") %>%
  group_by(treatment, metric) %>%
  summarise(rmse = sqrt(mean(error^2, na.rm = TRUE)), .groups = "drop") %>%
  pivot_wider(names_from = treatment, values_from = rmse)

forecast_ttests <- forecast_tests(tasks_all)

aov_f0_treatment_task <- summary(
  aov(se_f0 ~ treatment + factor(Task),
      data = mutate(tasks_all, se_f0 = error_f0_bar^2))
)

friedman_f0_by_task <- friedman.test(
  se_f0 ~ treatment | Task,
  data = mutate(tasks_all, se_f0 = error_f0_bar^2)
)

aov_alpha_treatment_task <- summary(aov(alpha_est ~ treatment + factor(Task), data = tasks_all))
aov_beta_treatment_task <- summary(aov(beta_est ~ treatment + factor(Task), data = tasks_all))


################
### Figures
################

segments <- make_segments()

gamma_values <- c(0.8, 0.6, 0.4, 0.2, 0, -0.5, -1, -2, -4)
gamma_contours <- make_contours(gamma_values)
gamma_contourplot <- make_gamma_contourplot(segments, gamma_contours)

benefit_gamma_values <- c(0.8, 0.6, 0.5, 0.4, 0.2, 0, -0.5, -1, -2, -4)
benefit_contours <- make_contours(benefit_gamma_values)
benefit_contourplot <- make_benefit_contourplot(segments, benefit_contours)

tasks_all <- tasks_all %>%
  mutate(treatment = factor(treatment, levels = treatment_levels))

scatterplot <- make_scatterplot(tasks_all, segments)
gamma_hist <- make_gamma_hist(tasks_all, desc_stats)
kappa_hist <- make_kappa_hist(tasks_all, desc_stats)

ggsave(file.path(figures_dir, "gamma_contourplot.pdf"), gamma_contourplot, width = 4.5, height = 4.5, device = cairo_pdf)
ggsave(file.path(figures_dir, "benefit_contourplot.pdf"), benefit_contourplot, width = 4.5, height = 4.5, device = cairo_pdf)
ggsave(file.path(figures_dir, "scatterplot.pdf"), scatterplot, width = 4.5, height = 4.5, device = cairo_pdf)
ggsave(file.path(figures_dir, "gamma_hist.pdf"), gamma_hist, width = 4.5, height = 2.25, device = cairo_pdf)
ggsave(file.path(figures_dir, "kappa_hist.pdf"), kappa_hist, width = 4.5, height = 2.25, device = cairo_pdf)

main_paper_results <- list(
  responses_all = responses_all,
  participant_summary = participant_summary,
  tasks_all = tasks_all,
  desc_stats = desc_stats,
  forecast_rmse = forecast_rmse,
  forecast_ttests = forecast_ttests,
  aov_f0_treatment_task = aov_f0_treatment_task,
  friedman_f0_by_task = friedman_f0_by_task,
  aov_alpha_treatment_task = aov_alpha_treatment_task,
  aov_beta_treatment_task = aov_beta_treatment_task,
  gamma_contourplot = gamma_contourplot,
  benefit_contourplot = benefit_contourplot,
  scatterplot = scatterplot,
  gamma_hist = gamma_hist,
  kappa_hist = kappa_hist
)

invisible(main_paper_results)
