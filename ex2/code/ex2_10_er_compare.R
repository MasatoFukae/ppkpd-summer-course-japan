###############################################################################
# PROGRAM NAME        : ex2_10_er_compare.R
# DESCRIPTION         : To compare ER curve for multiple models
# COMPOUND NAME       : Not Applicable
# PROGRAM VERSION     : 02
# CREATION DATE       : 31-July-2025
# PROGRAMMER          : pkpd_seminor
# EXTERNAL FILES USED : model rds file(s) to compare
###############################################################################

## clear workspace
rm(list=ls())

# specify the path for analysis
path = ""

# specify the model files to compare
rds_files = c(".rds", 
              ".rds")

# Specify exposure variable and range for simulation
exposure_var <- ""  # Change as needed
exposure_range <- c(0, 25)  # Set your desired range

### program code ####################

  source(paste0(path, "/Source.R"))
  
  if (!dir.exists(paste0(path, "/output"))) dir.create(paste0(path, "/output"))

  # Read multiple RDS files (logistic regression models) and simulate exposure-response curves for a pre-defined range
  
  exposure_seq <- seq(exposure_range[1], exposure_range[2], length.out = 100)
  
  sim_results <- list()
  
  for (file in rds_files) {
    fit <- readRDS(file.path(path, "output", file))
    model_vars <- names(coef(fit))
    # Only use main effect covariates (not dummy variables created by factors)
    mf_vars <- attr(terms(fit), "term.labels")
    other_vars <- setdiff(mf_vars, exposure_var)
    
    # Prepare newdata for prediction
    newdata <- setNames(data.frame(exposure_seq), exposure_var)
    if (length(other_vars) > 0) {
      for (v in other_vars) {
        # Use most frequent for factor, median for numeric
        if (is.factor(fit$data[[v]])) {
          newdata[[v]] <- names(sort(table(fit$data[[v]]), decreasing = TRUE))[1]
        } else if (is.numeric(fit$data[[v]])) {
          newdata[[v]] <- median(fit$data[[v]], na.rm = TRUE)
        } else {
          newdata[[v]] <- unique(fit$data[[v]])[1]
        }
      }
    }
    
    pred_prob <- predict(fit, newdata = newdata, type = "response")
    sim_results[[basename(file)]] <- data.frame(
      exposure = exposure_seq,
      predicted_prob = pred_prob,
      model = basename(file)
    )
  }
  
  sim_df <- dplyr::bind_rows(sim_results)
  
  # Plot all exposure-response curves
  p = ggplot(sim_df, aes(x = exposure, y = predicted_prob, color = model)) +
    geom_line(size = 1) +
    labs(
      title = paste("Simulated Exposure-Response Curves for", exposure_var),
      x = exposure_var,
      y = "Predicted Probability"
    ) +
    theme_minimal() + theme(legend.position = "top") + 
    guides(color = guide_legend(nrow = length(rds_files), byrow = TRUE))
  
  print(p)
  ggsave(
    filename = file.path(path, paste0("/output/logit_fit_", exposure_var, "_compare_models.png")),
    plot = p,
    width = 7, height = 5, dpi = 300
  )

###############################################################################
#
# End of program
#
###############################################################################