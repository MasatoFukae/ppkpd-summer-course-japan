###############################################################################
# PROGRAM NAME        : ex2_07_logit_reg_cov.R
# DESCRIPTION         : To conduct logistic regression analysis with covariate effect
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

# specify the covariate for regression
covariate = ""     # "age", "bsa", "ethnicity", "initial_tumor", "gender"

### program code ####################

  source(paste0(path, "/Source.R"))
  
  if (!dir.exists(paste0(path, "/output"))) dir.create(paste0(path, "/output"))
  
  # Read data
  df <- readr::read_csv(file.path(path, paste0("/data/", dat))) %>% 
    filter(route == 1) %>%  
    mutate(gender = ifelse(gender == 0, "Male", "Female")) %>% 
    mutate(ethnicity = ifelse(ethnicity == 1, "Caucasian",
                              ifelse(ethnicity == 2, "African American",
                                     ifelse(ethnicity == 3, "Asian",
                                            ifelse(ethnicity == 4, "Hispanic", "Other")))))
  
  # Prepare formula for logistic regression including covariate
  formula <- as.formula(paste0(endpoint, " ~ ", exposure, " + ", covariate))
  
  # Conduct logistic regression (assume endpoint is binary: 0/1 or factor)
  fit <- glm(formula, data = df, family = binomial())
  
  # Save model object
  saveRDS(fit, file = file.path(path, paste0("/output/reg_", exposure, "_", endpoint, "_with_", covariate, ".rds")))
  
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
  sink(file.path(path, paste0("/output/reg_", exposure, "_", endpoint, "_with_", covariate, ".txt")))
  cat("Logistic Regression Summary (with covariate)\n")
  print(summary(fit))
  cat("\nOdds Ratios and 95% CI\n")
  print(result)
  sink()
  
  # Stratify population for fit check
  if (is.numeric(df[[covariate]])) {
    # Continuous: split by median
    median_cov <- median(df[[covariate]], na.rm = TRUE)
    df <- df %>%
      mutate(cov_group = ifelse(get(covariate) <= median_cov, 
                                paste0(covariate, "<=median"), 
                                paste0(covariate, ">median")))
    cov_group_levels <- unique(df$cov_group)
  } else {
    # Categorical: use as is
    df <- df %>% mutate(cov_group = as.factor(get(covariate)))
    cov_group_levels <- levels(as.factor(df$cov_group))
  }
  
  n_quantiles <- 4
  obs_rate_all <- list()
  pred_curve_all <- list()
  
  for (group in cov_group_levels) {
    df_sub <- df %>% filter(cov_group == group)
    # Observed quantile rates
    df_sub <- df_sub %>%
      mutate(
        exposure_quantile = cut(get(exposure),
                                breaks = quantile(get(exposure), probs = seq(0, 1, length.out = n_quantiles + 1), na.rm = TRUE),
                                labels = paste0("Q", 1:n_quantiles),
                                include.lowest = TRUE)
      )
    obs_rate <- df_sub %>%
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
    obs_rate$cov_group <- group
    
    # Predicted curve: use the same model, set covariate to group value (median for continuous, level for categorical)
    exposure_range <- range(df_sub[[exposure]], na.rm = TRUE)
    exposure_seq <- seq(exposure_range[1], exposure_range[2], length.out = 100)
    newdata <- data.frame(exposure_seq)
    names(newdata) <- exposure
    if (is.numeric(df[[covariate]])) {
      # Use median value for this group
      covariate_value <- if (grepl("<=median", group)) median(df[[covariate]][df[[covariate]] <= median_cov], na.rm = TRUE)
      else median(df[[covariate]][df[[covariate]] > median_cov], na.rm = TRUE)
    } else {
      covariate_value <- group
    }
    newdata[[covariate]] <- covariate_value
    
    pred_link <- predict(fit, newdata = newdata, type = "link", se.fit = TRUE)
    crit <- qnorm(0.95)
    fit_link <- pred_link$fit
    se_link <- pred_link$se.fit
    lower_link <- fit_link - crit * se_link
    upper_link <- fit_link + crit * se_link
    pred_curve <- data.frame(
      exposure_value = exposure_seq,
      predicted_prob = plogis(fit_link),
      lower90 = plogis(lower_link),
      upper90 = plogis(upper_link),
      cov_group = group
    )
    
    obs_rate_all[[group]] <- obs_rate
    pred_curve_all[[group]] <- pred_curve
  }
  
  obs_rate_all <- bind_rows(obs_rate_all)
  pred_curve_all <- bind_rows(pred_curve_all)
  
  p <- ggplot() +
    geom_pointrange(
      data = obs_rate_all,
      aes(x = mean_exposure, y = obs_rate, ymin = lower90, ymax = upper90),
      color = "red", size = 1.2
    ) +
    geom_ribbon(
      data = pred_curve_all,
      aes(x = exposure_value, ymin = lower90, ymax = upper90),
      fill = "blue", alpha = 0.2
    ) +
    geom_line(
      data = pred_curve_all,
      aes(x = exposure_value, y = predicted_prob),
      color = "blue", size = 1
    ) +
    facet_wrap(~cov_group) +
    labs(
      title = paste("Observed Rate by Quantile and Predicted Probability (", exposure, "+", covariate, ")"),
      x = exposure,
      y = "Observed Rate / Predicted Probability"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 12))
  
  print(p)
  ggsave(
    filename = file.path(path, paste0("output/logit_fit_", exposure, "_", endpoint, "_with_", covariate, ".png")),
    plot = p,
    width = 10, height = 5, dpi = 300
  )

###############################################################################
#
# End of program
#
###############################################################################