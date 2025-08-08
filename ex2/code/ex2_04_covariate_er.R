###############################################################################
# PROGRAM NAME        : ex2_04_covariate_er.R
# DESCRIPTION         : To explore the impact of covariate on ER
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

# specify the covaraite(s)
covariates = c() #"ethnicity", "initial_tumor", "gender"

### program code ####################

  source(paste0(path, "/Source.R"))

  if (!dir.exists(paste0(path, "/output"))) dir.create(paste0(path, "/output"))
  
  # Prepare data: binary outcome for CRS onset at Cycle 1
  er_data <- read_csv(paste0(path, "/data/", dat)) %>% filter(route == 1) %>%
    rename(resp_binary := !!sym(endpoint)) %>% 
    mutate(gender = ifelse(gender == 0, "Male", "Female")) %>% 
    mutate(ethnicity = ifelse(ethnicity == 1, "Caucasian",
                              ifelse(ethnicity == 2, "African American",
                                     ifelse(ethnicity == 3, "Asian",
                                            ifelse(ethnicity == 4, "Hispanic", "Other")))))
  
  # Define covariate groups (example: age quartile, BSA quartile, initial tumor quartile)
  data <- er_data %>%
    mutate(
      age_group = cut(age, breaks = quantile(age, probs = seq(0, 1, 0.25), na.rm = TRUE), labels = paste0("Q", 1:4), include.lowest = TRUE),
      bsa_group = cut(bsa, breaks = quantile(bsa, probs = seq(0, 1, 0.25), na.rm = TRUE), labels = paste0("Q", 1:4), include.lowest = TRUE),
      tumor_group = cut(initial_tumor, breaks = quantile(initial_tumor, probs = seq(0, 1, 0.25), na.rm = TRUE), labels = paste0("Q", 1:4), include.lowest = TRUE)
    )
  
  covariates_ <- c("age_group", "gender", "ethnicity", "bsa_group", "tumor_group")
  covariates_ <- covariates_[grepl(paste(covariates, collapse="|"), covariates_)]
  exposure_metrics <- exposure
  
  for (exposure_var in exposure_metrics) {
    for (cov in covariates_) {
      # Assign quantile groups for exposure metric
      data_ <- data %>%
        mutate(
          exposure_quantile = cut(get(exposure_var),
                                  breaks = quantile(get(exposure_var), probs = seq(0, 1, 0.25), na.rm = TRUE),
                                  labels = paste0("Q", 1:4),
                                  include.lowest = TRUE)
        )
      
      # Calculate observed rates and CIs by exposure quantile and covariate group
      rates <- data_ %>%
        group_by(exposure_quantile, !!sym(cov)) %>%
        summarise(
          n = n(),
          n_ = sum(resp_binary),
          rate = mean(resp_binary ),
          se = sqrt((rate * (1 - rate)) / n),
          lower = rate - 1.96 * se,
          upper = rate + 1.96 * se,
          mean_exposure = mean(get(exposure_var), na.rm = TRUE),
          .groups = "drop"
        )
      
      # Plot quantile-based observed CRS rates (x = mean_exposure)
      p <- ggplot(rates, aes(x = mean_exposure, y = rate, color = !!sym(cov), fill = !!sym(cov), group = !!sym(cov))) +
        geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.2), size = 1.2) +
        geom_point(size = 3, position = position_dodge(width = 0.2)) +
        geom_line(linetype=2) + 
        theme_minimal() +
        labs(
          title = paste(endpoint, "Rate by", exposure_var, "Quartile, Grouped by", cov),
          x = paste(exposure_var),
          y = paste(endpoint, "Rate"),
          color = cov,
          fill = cov
        ) +
        scale_y_continuous(limits = c(0, 1))
      
      print(p)
      ggsave(paste0(path, "/output/", endpoint, "_rate_by_", exposure_var, "_by_", cov, ".png"), plot = p, width = 7, height = 5, dpi = 300)
    }
  }
  
###############################################################################
#
# End of program
#
###############################################################################





