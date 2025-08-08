###############################################################################
# PROGRAM NAME        : ex2_02_exposure_matrix.R
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

# specify the name of exposure metrics for eda
exposures = c("Cmax_1", "AUC_1", "Cmin_1")      #"Cmax_1", "Cmax_2", "Cmax_3", "AUC_1", "AUC_2", "AUC_3", "Cmin_1", "Cmin_2", "Cmin_3"

### program code ####################

  source(paste0(path, "/Source.R"))

  if (!dir.exists(paste0(path, "/output"))) dir.create(paste0(path, "/output"))
  
  # Prepare data: binary outcome for CRS onset at Cycle 1
  er_data <- read_csv(paste0(path, "/data/", dat)) %>% filter(route == 1) 
  
    # Add Route column: 
    cycle_vars <- er_data %>%
      mutate(Route = ifelse(route == 0, "IV", "SC")) %>%
      dplyr::select(ID, Route, exposures)
    
    # Custom color mapping: black for SC, red for IV
    my_colors <- c("SC" = "darkgray", "IV" = "red")
    
    cycle_vars = tibble(cycle_vars)
    
    # Only scatter, density, and correlation (no boxplot)
    p <- GGally::ggpairs(
      cycle_vars %>% dplyr::select(exposures),  # Route を除外
      mapping = ggplot2::aes(color = cycle_vars$Route),  # 色分けには使う
      lower = list(continuous = wrap("points", alpha = 0.7, size = 2)),
      diag = list(continuous = wrap("densityDiag")),
      upper = list(continuous = wrap("cor", size = 5)),
      title = paste("PK Parameter Correlation")
    ) +
      theme_minimal() +
      scale_color_manual(values = my_colors)+ scale_fill_manual(values=my_colors)
    
    print(p)
    
    exposures_ = paste0(exposures, collapse="-")
    
    ggsave(paste0(path,"/output/", "exposure-corr-", exposures_, ".png"), p, width = 8, height = 8, dpi = 300)
  
###############################################################################
#
# End of program
#
###############################################################################





