###############################################################################
# PROGRAM NAME        : ex2_08_logit_reg_matched.R
# DESCRIPTION         : To conduct logistic regression analysis with covariate matching
# COMPOUND NAME       : Not Applicable
# PROGRAM VERSION     : 02
# CREATION DATE       : 31-July-2025
# PROGRAMMER          : pkpd_seminor
# EXTERNAL FILES USED : erdat.csv (er analysis dataset)
###############################################################################

## clear workspace
rm(list=ls())

# specify the path for analysis
path = ""

# specify the name of analysis dataset
dat = "erdat.csv"

# specify the endpoint for regression
endpoint = ""  #"CRS_onset_1", "CRS_onset_2", "CRS_onset_3",  "response"

# specify the name of exposure metrics for regression
exposure = ""      #"Cmax_1", "Cmax_2", "Cmax_3", "AUC_1", "AUC_2", "AUC_3", "Cmin_1", "Cmin_2", "Cmin_3"

# specify the covariate for matching
covariate = ""     # "age", "bsa", "ethnicity", "initial_tumor", "gender"

### program code ####################

  source(paste0(path, "/Source.R"))
  
  if (!dir.exists(file.path(path, "output"))) dir.create(file.path(path, "output"))
  
  # Read data
  df <- read_csv(file.path(path, "data", dat)) %>% filter(route==1)
  
  # matching on the covariate
  # Create a binary group for high/low exposure (median split, numeric 1/0)
  median_exposure <- median(df[[exposure]], na.rm = TRUE)
  df <- df %>%
    mutate(exposure_group = ifelse(get(exposure) > median_exposure, 1, 0))
  # Matching: High (1) vs Low (0) exposure, matched on initial_tumor
  m.out <- matchit(
    exposure_group ~ covariate,
    data     = df,
    method   = "nearest",
    distance = df[[covariate]], 
    caliper  = 0.1,
    std.caliper = TRUE
  )
  df_matched <- match.data(m.out)
  
  
  # Prepare formula for logistic regression
  formula <- as.formula(paste0(endpoint, "~", exposure))
  
  # Conduct logistic regression (assume endpoint is binary: 0/1 or factor)
  fit <- glm(formula, data = df_matched, family = binomial())
  
  # Save model object
  saveRDS(fit, file = file.path(path, "output", paste0("reg_", exposure, "_", endpoint, "_matched-", covariate, ".rds")))
  
  # Show summary
  summary(fit)
  
  # print odds ratio and 95% CI
  exp_coef <- exp(coef(fit))
  confint_fit <- exp(confint(fit))
  result <- data.frame(
    OR = exp_coef,
    lower95 = confint_fit[,1],
    upper95 = confint_fit[,2]
  )
  print(result)
  
  # Output summary to txt file
  sink(file.path(path, "output", paste0("reg_", exposure, "_", endpoint, "_matched-", covariate, ".txt")))
  cat("Logistic Regression Summary (Matched Population)\n")
  print(summary(fit))
  cat("\nOdds Ratios and 95% CI\n")
  print(result)
  sink()
  
  # Calculate observed response rate by exposure quantiles
  n_quantiles <- 4
  df_ <- df_matched %>%
    mutate(
      exposure_quantile = cut(get(exposure),
                              breaks = quantile(get(exposure), probs = seq(0, 1, length.out = n_quantiles + 1), na.rm = TRUE),
                              labels = paste0("Q", 1:n_quantiles),
                              include.lowest = TRUE)
    )
  
  # Calculate 90% confidence intervals for both observed quantile rates and predicted curve
  # 1. Observed quantile rates: 90% CI using normal approximation
  obs_rate <- df_ %>%
    group_by(exposure_quantile) %>%
    summarise(
      n = n(),
      n_resp = sum(get(endpoint)),
      obs_rate = mean(get(endpoint)),
      se = sqrt(obs_rate * (1 - obs_rate) / n),
      mean_exposure = mean(get(exposure), na.rm = TRUE),
      lower90 = obs_rate - 1.645 * se,
      upper90 = obs_rate + 1.645 * se,
      .groups = "drop"
    )
  
  # 2. Predicted curve: 90% CI using standard error of prediction
  exposure_range <- range(df_[[exposure]], na.rm = TRUE)
  exposure_seq <- seq(exposure_range[1], exposure_range[2], length.out = 100)
  newdata <- df_[rep(1, 100), ]
  newdata[[exposure]] <- exposure_seq
  
  pred_link <- predict(fit, newdata = newdata, type = "link", se.fit = TRUE)
  crit <- qnorm(0.95) # 90% CI (two-sided)
  fit_link <- pred_link$fit
  se_link <- pred_link$se.fit
  lower_link <- fit_link - crit * se_link
  upper_link <- fit_link + crit * se_link
  
  pred_curve <- data.frame(
    exposure_value = exposure_seq,
    predicted_prob = plogis(fit_link),
    lower90 = plogis(lower_link),
    upper90 = plogis(upper_link)
  )
  
  # 3. Plot with 90% CI for both observed and predicted, with reduced font size for title
  p <- ggplot() +
    geom_pointrange(
      data = obs_rate,
      aes(x = mean_exposure, y = obs_rate, ymin = lower90, ymax = upper90),
      color = "red", size = 1.2
    ) +
    geom_ribbon(
      data = pred_curve,
      aes(x = exposure_value, ymin = lower90, ymax = upper90),
      fill = "blue", alpha = 0.2
    ) +
    geom_line(
      data = pred_curve,
      aes(x = exposure_value, y = predicted_prob),
      color = "blue", size = 1
    ) +
    labs(
      title = paste("Observed Rate by Quantile and Predicted Probability (", exposure, ")"),
      x = exposure,
      y = "Observed Rate / Predicted Probability"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 12))
  
  print(p)
  ggsave(
    filename = file.path(path, paste0("output/logit_fit_", exposure, "_", endpoint, "_matched-", covariate, ".png")),
    plot = p,
    width = 7, height = 5, dpi = 300
  )

###############################################################################
#
# End of program
#
###############################################################################