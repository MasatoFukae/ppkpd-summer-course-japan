###############################################################################
# PROGRAM NAME        : ex2_09_logit_reg_cov.R
# DESCRIPTION         : To compare AIC for multiple models
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

### program code ####################

  source(paste0(path, "/Source.R"))
  
  if (!dir.exists(paste0(path, "/output"))) dir.create(paste0(path, "/output"))
  
  # Compare AIC values of multiple logistic regression models saved as RDS files and export as txt

  # Specify the folder containing the RDS files
  model_dir <- file.path(path, "output")
  
  # List all RDS files in the folder (adjust pattern if needed)
  rds_files <- list.files(model_dir, pattern = "\\.rds$", full.names = TRUE)
  
  # Read models and calculate AIC
  aic_table <- lapply(rds_files, function(f) {
    model <- readRDS(f)
    data.frame(
      file = basename(f),
      AIC = AIC(model)
    )
  }) %>% bind_rows()
  
  # Sort by AIC
  aic_table <- aic_table %>% arrange(AIC)
  
  print(aic_table)
  
  # Export as txt file
  write.table(aic_table, file = file.path(model_dir, "logistic_regression_AIC_comparison.txt"), 
              row.names = FALSE, sep = "\t", quote = FALSE)

###############################################################################
#
# End of program
#
###############################################################################