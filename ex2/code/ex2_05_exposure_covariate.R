###############################################################################
# PROGRAM NAME        : ex2_05_exposure_covariate.R
# DESCRIPTION         : To compare the exposure metrics by covaraite
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

# specify the name of exposure metrics for comparison
exposure = ""      #"Cmax_1", "Cmax_2", "Cmax_3", "AUC_1", "AUC_2", "AUC_3", "Cmin_1", "Cmin_2", "Cmin_3"

# specify a covariate
covariates = ""  #"age", "bsa", "ethnicity", "initial_tumor", "gender"


### program code ####################

  source(paste0(path, "/Source.R"))
  
  if (!dir.exists(paste0(path, "/output"))) dir.create(paste0(path, "/output"))
  
  # Read data
  df <- readr::read_csv(file.path(path, "data", dat))
  
  # Check if covariate is continuous or categorical
  if (is.numeric(df[[covariates]])) {
    # Continuous: split by median
    median_cov <- median(df[[covariates]], na.rm = TRUE)
    df <- df %>%
      mutate(cov_group = ifelse(get(covariates) > median_cov, 
                                paste0(covariates, ">median"), 
                                paste0(covariates, "<=median")))
    p <- ggplot(df, aes(x = cov_group, y = get(exposure), fill = cov_group)) +
      geom_boxplot() +
      labs(
        title = paste("Exposure (", exposure, ") by", covariates, "Group"),
        x = covariates,
        y = exposure
      ) +
      theme_minimal()
  } else {
    # Categorical: use as is
    df <- df %>% mutate(cov_group = as.factor(get(covariates)))
    p <- ggplot(df, aes(x = cov_group, y = get(exposure), fill = cov_group)) +
      geom_boxplot() +
      labs(
        title = paste("Exposure (", exposure, ") by", covariates),
        x = covariates,
        y = exposure
      ) +
      theme_minimal()
  }
  
  print(p)
  ggsave(
    filename = file.path(path, paste0("/output/", exposure, "_by_", covariates, ".png")),
    plot = p,
    width = 6, height = 4, dpi = 300
  )
  
###############################################################################
#
# End of program
#
###############################################################################





