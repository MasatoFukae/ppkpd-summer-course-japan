###############################################################################
# PROGRAM NAME        : ex2_03_exposure_matrix.R
# DESCRIPTION         : To explore the correlation of exposure metrics
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

# specify the covariates
covariates = c()  #"age", "bsa", "ethnicity", "initial_tumor", "gender"

### program code ####################

  source(paste0(path, "/Source.R"))
  
  if (!dir.exists(paste0(path, "/output"))) dir.create(paste0(path, "/output"))
  
  # Prepare data
  er_data <- read_csv(paste0(path, "/data/", dat)) %>% filter(route == 1)
  cycle_vars <- er_data %>%
      mutate(gender = as.factor(gender),
             ethnicity = as.factor(ethnicity)) %>% 
      dplyr::select(covariates)
    
  cycle_vars = tibble(cycle_vars)
  covariates_ = paste(covariates, collapse = "-")
  
    p_corr <- ggpairs(
      cycle_vars,
      #mapping = aes(color = Gender),
      lower = list(continuous = wrap("smooth", alpha = 0.3, size = 0.5)),
      diag = list(continuous = wrap("densityDiag")),
      upper = list(continuous = wrap("cor", size = 4)),
      title = paste("Covariate Correlation:", covariates_)
    ) +
      theme_minimal()
    
    print(p_corr)
    
    
    ggsave(paste0(path,"/output/", "covariates-corr-", covariates_, ".png"), p_corr, width = 12, height = 8, dpi = 300)
  
###############################################################################
#
# End of program
#
###############################################################################





