###############################################################################
# PROGRAM NAME        : ex2_01_er_graphic.R
# DESCRIPTION         : To run exploatory data analysis
# COMPOUND NAME       : Not Applicable
# PROGRAM VERSION     : 01
# CREATION DATE       : 28-July-2025
# PROGRAMMER          : pkpd_seminor
# EXTERNAL FILES USED : erdat.csv (er analysis dataset)
###############################################################################

## clear workspace
rm(list=ls())

# specify the path for analysis
path = ""

# specify the name of analysis dataset
dat = "erdat.csv"

# specify the endpoint for eda
endpoint = ""  #"CRS_onset_1", "CRS_onset_2", "CRS_onset_3",  "response"

# specify the name of exposure metrics for eda
exposure = ""      #"Cmax_1", "Cmax_2", "Cmax_3", "AUC_1", "AUC_2", "AUC_3", "Cmin_1", "Cmin_2", "Cmin_3"


### program code ####################

  source(paste0(path, "/Source.R"))

  if (!dir.exists(paste0(path, "/output"))) dir.create(paste0(path, "/output"))

  # Function to calculate response rates by quantiles
  calc_response_rates <- function(data, exposure_var, n_quantiles = 4) {
    # Add quantile groups
    data <- data %>%
      mutate(
        quantile_group = cut(get(exposure_var), 
                             breaks = quantile(get(exposure_var), 
                                               probs = seq(0, 1, length.out = n_quantiles + 1)),
                             labels = paste0("Q", 1:n_quantiles),
                             include.lowest = TRUE)
      )
    
    # Calculate response rates and confidence intervals
    response_rates <- data %>%
      group_by(quantile_group) %>%
      summarise(
        n = n(),
        n_resp = sum(resp_binary),
        resp_rate = mean(resp_binary),
        se = sqrt((resp_rate * (1 - resp_rate)) / n),
        lower = resp_rate - 1.96 * se,
        upper = resp_rate + 1.96 * se,
        mean_exposure = mean(get(exposure_var))
      )
    
    return(response_rates)
  }
  
  # Function to create plot with quantiles
  fit_logistic_er_with_quantiles <- function(data, exposure_var, cycle, title) {
   
     # Fit logistic regression
    model <- glm(resp_binary ~ get(exposure_var), data = data, family = binomial())

    
    # Calculate response rates by quantiles
    resp_rates <- calc_response_rates(data, exposure_var)
    
    # Create plot
    p <- ggplot() +
      # Add observed data points
      geom_point(data = data, aes_string(x = exposure_var, y = "resp_binary"), 
                 alpha = 0.5) +
      # Add quantile-based response rates
      geom_pointrange(data = resp_rates,
                      aes(x = mean_exposure, y = resp_rate,
                          ymin = lower, ymax = upper),
                      color = "red", size = 1) +
      theme_minimal() +
      labs(title = paste(title, "- Cycle", cycle),
           x = paste(exposure_var, "(μg/L)"),
           y = "Probability of Response")
    
    # Extract parameters
    coef_summary <- summary(model)$coefficients
    
    return(list(plot = p, model = model, parameters = coef_summary))
  }
  
  
  # Prepare cycle-specific data for exposure-response analysis
  er_data <- read_csv(paste0(path, "/data/", dat)) %>% filter(route == 1) %>%
    rename(resp_binary := !!sym(endpoint))

  # Analyze exposure-response relationships by cycle
  exposure_metrics <- strsplit(exposure, "_")[[1]][1]
  cycle = strsplit(exposure, "_")[[1]][2]
  er_models <- list()
  er_plots <- list()
  
  for (metric in exposure_metrics) {
      metric_name <- paste0(metric, "_", cycle)
      result <- fit_logistic_er_with_quantiles(
        data = er_data,
        exposure_var = metric_name,
        cycle = cycle,
        title = paste(endpoint,"vs", metric)
      )
      
      er_models[[metric_name]] <- result$model
      er_plots[[metric_name]] <- result$plot
      
      # Print parameter estimates
      cat("\nParameter estimates for", metric_name, ":\n")
      # print(result$parameters)
    }
  
  # Create combined plots for each exposure metric
  for (metric in exposure_metrics) {
    metric_plots <- er_plots[grep(metric, names(er_plots))]
    # Arrange plots in a grid
    combined_plot <- do.call(gridExtra::grid.arrange, 
                             c(metric_plots, 
                               top = paste(metric, "Exposure-Response Relationships by Cycle"),
                               ncol = 1))
    # Print the combined plot
    print(combined_plot)
    
    # Save as PNG
    ggsave(
      filename = paste0(path, "/output/eda_", endpoint, "_", exposure, ".png"),
      plot = combined_plot,
      width = 6, height = 6, dpi = 300
    )
  }
  
###############################################################################
#
# End of program
#
###############################################################################






  
  