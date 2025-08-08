###############################################################################
# PROGRAM NAME        : ex2_11_simulation.R
# DESCRIPTION         : To simulate dose-response for efficacy (RR) and safety (>=CRS grade 2)
# COMPOUND NAME       : Not Applicable
# PROGRAM VERSION     : 02
# CREATION DATE       : 31-July-2025
# PROGRAMMER          : pkpd_seminor
# EXTERNAL FILES USED : model rds files, erdat (er analysis dataset)
###############################################################################

## clear workspace
rm(list=ls())

# specify the path for analysis
path = ""

# specify dataset name
dat = "erdat.csv"

# specify the model files for response rate and CRS
model_rr  = ".rds"
model_crs = ".rds"

# Define designs to simulate
maintenance_doses <- seq(10, 60, by = 5)
n_subjects <- 100
n_sim <- 100

### program code ####################

  source(paste0(path, "/Source.R"))
  
  if (!dir.exists(paste0(path, "/output"))) dir.create(paste0(path, "/output"))

  erdat <- readr::read_csv(file.path(path, "data", dat))
  
  # Prepare storage for results
  rr_mat <- matrix(NA, nrow = n_sim, ncol = length(maintenance_doses))
  crs_arr <- array(NA, dim = c(n_sim, length(maintenance_doses), 3)) # [sim, dose, cycle]
  
  # Load model objects correctly
  model_rr  <- readRDS(file.path(path, "output", model_rr))
  model_crs <- readRDS(file.path(path, "output", model_crs))
  
  # Create model object
  mod <- mcode("ex1_plus_pkpd", code) %>%
    Req(CP, CP_OBS, CYK_IPRED, CYK_OBS, TUMOR_IPRED, TUMOR_OBS, CL, BSA, BASE_TUMOR ,Dose, Cycle)
  
  for (j in 1:n_sim) {
    # Simulate BSA for each subject
    set.seed(1000 + j)
    bsa_values <- rnorm(n_subjects, mean = 1.88, sd = 0.25)
    bsa_values <- pmax(pmin(bsa_values, 2.5), 1.3)
    
    sim_results <- list()
    for (i in seq_along(maintenance_doses)) {
      dose <- maintenance_doses[i]
      # Dosing events for 3 cycles (step-up + maintenance)
      dosing <- tibble(
        ID = rep(1:n_subjects, each = 3),
        time = rep(c(0, 24*3, 24*7), times = n_subjects),
        evid = 1,
        amt = rep(c(3, 6, dose), times = n_subjects),
        cmt = 1,
        ii = 0,
        addl = 0
      )
      # Add maintenance dosing for cycles 2 and 3 (QW for 2 more cycles)
      dosing_maint <- tibble(
        ID = rep(1:n_subjects, each = 8),
        time = rep(seq(24*7*2, 24*7*9, by = 24*7), times = n_subjects),
        evid = 1,
        amt = dose,
        cmt = 1,
        ii = 168,
        addl = 0
      )
      dosing_all <- bind_rows(dosing, dosing_maint) %>% arrange(ID, time)
      # Observation times (dense for 3 cycles)
      obs_times <- expand.grid(
        ID = 1:n_subjects,
        time = seq(0, 24*28*3, by = 12)
      ) %>%
        mutate(evid = 0, amt = 0, cmt = 2, ii = 0, addl = 0)
      sim_df <- bind_rows(dosing_all, obs_times) %>% arrange(ID, time)
      idata <- tibble(ID = 1:n_subjects, BSA = bsa_values)
      sim <- mod %>%
        data_set(sim_df) %>%
        idata_set(idata) %>%
        mrgsim_df(atol=1e-100, rtol=1e-4, tad = TRUE, output='silent') %>%
        as_tibble() %>%
        group_by(ID) %>%
        mutate(
          Cycle = 1 + floor(time/(28*24)),
          time_in_cycle = time - (Cycle - 1) * 28 * 24
        ) %>%
        filter(Cycle <= 3) %>%
        ungroup() %>%
        mutate(mDose = dose)
      sim_results[[i]] <- sim
    }
    sim_all <- bind_rows(sim_results)
    
    # Calculate AUC for each cycle and subject
    calc_auc <- function(df, cycle) {
      df %>%
        filter(Cycle == cycle) %>%
        group_by(ID, mDose) %>%
        summarise(AUC = pracma::trapz(time, CP_OBS), .groups = "drop") %>%
        mutate(Cycle = cycle)
    }
    
    # Calculate Cmax for each cycle and subject
    calc_cmax <- function(df, cycle) {
      df %>%
        filter(Cycle == cycle) %>%
        group_by(ID, mDose) %>%
        summarise(Cmax = max(CP_OBS), .groups = "drop") %>%
        mutate(Cycle = cycle)
    }
    
    # Calculate Cmin for each cycle and subject
    calc_cmin <- function(df, cycle) {
      df %>%
        filter(Cycle == cycle) %>%
        group_by(ID, mDose) %>%
        summarise(Cmin = CP_OBS[time==max(time)][1], .groups = "drop") %>%
        mutate(Cycle = cycle)
    }
    
    # Get model formula and variables for RR and CRS
    model_terms_rr  <- attr(model_rr$terms, "term.labels")
    exposure_var_rr <- model_terms_rr[grepl("AUC|Cmax|Cmin|AUC_", model_terms_rr)][1]
    covariate_vars_rr <- setdiff(model_terms_rr, exposure_var_rr)
    
    model_terms_crs  <- attr(model_crs$terms, "term.labels")
    exposure_var_crs <- model_terms_crs[grepl("AUC|Cmax|Cmin|AUC_", model_terms_crs)][1]
    covariate_vars_crs <- setdiff(model_terms_crs, exposure_var_crs)
    
    # Simulate only necessary exposure metrics
    auc_all <- NULL
    cmax_all <- NULL
    cmin_all  <- NULL
    
    if (!is.null(exposure_var_rr) && grepl("^AUC_", exposure_var_rr)) {
      auc_cycle <- as.numeric(gsub("AUC_", "", exposure_var_rr))
      auc_all <- calc_auc(sim_all, auc_cycle) %>%
        pivot_wider(names_from = Cycle, values_from = AUC, names_prefix = "AUC_")
    }
    if (!is.null(exposure_var_rr) && grepl("^Cmax_", exposure_var_rr)) {
      cmax_cycle <- as.numeric(gsub("Cmax_", "", exposure_var_rr))
      cmax_all <- calc_cmax(sim_all, cmax_cycle) %>%
        pivot_wider(names_from = Cycle, values_from = Cmax, names_prefix = "Cmax_")
    }
    if (!is.null(exposure_var_rr) && grepl("^Cmin_", exposure_var_rr)) {
      cmin_cycle <- as.numeric(gsub("Cmin_", "", exposure_var_rr))
      cmin_all <- calc_cmin(sim_all, cmin_cycle) %>%
        pivot_wider(names_from = Cycle, values_from = Cmin, names_prefix = "Cmin_")
    }
    
    
    if (!is.null(exposure_var_crs) && grepl("^Cmax_", exposure_var_crs)) {
      cmax_cycle_crs <- as.numeric(gsub("Cmax_", "", exposure_var_crs))
      cmax_all_crs <- calc_cmax(sim_all, cmax_cycle_crs) %>%
        pivot_wider(names_from = Cycle, values_from = Cmax, names_prefix = "Cmax_")
     if (is.null(cmax_all)) {
        cmax_all <- cmax_all_crs
      } else if (exposure_var_crs != exposure_var_rr){
        cmax_all <- left_join(cmax_all, cmax_all_crs, by = c("ID", "mDose"))
      }
    }
    
    if (!is.null(exposure_var_crs) && grepl("^AUC_", exposure_var_crs)) {
      auc_cycle_crs <- as.numeric(gsub("AUC_", "", exposure_var_crs))
      auc_all_crs <- calc_auc(sim_all, auc_cycle_crs) %>%
        pivot_wider(names_from = Cycle, values_from = AUC, names_prefix = "AUC_")
      if (is.null(auc_all)) {
        auc_all <- auc_all_crs
      } else if (exposure_var_crs != exposure_var_rr) {
        auc_all <- left_join(auc_all, auc_all_crs, by = c("ID", "mDose"))
      }
    }
    
    if (!is.null(exposure_var_crs) && grepl("^Cmin_", exposure_var_crs)) {
      cmin_cycle_crs <- as.numeric(gsub("Cmin_", "", exposure_var_crs))
      cmin_all_crs <- calc_cmin(sim_all, cmin_cycle_crs) %>%
        pivot_wider(names_from = Cycle, values_from = Cmin, names_prefix = "Cmin_")
      if (is.null(cmin_all)) {
        cmin_all <- cmin_all_crs
      } else if (exposure_var_crs != exposure_var_rr) {
        cmin_all <- left_join(cmin_all, cmin_all_crs, by = c("ID", "mDose"))
      }
    }
  
    
    # 
    # Prepare data for RR prediction
    if (!is.null(auc_all) && exposure_var_rr %in% names(auc_all)) {
      rr_data <- auc_all %>% dplyr::select(ID, mDose, !!exposure_var_rr)
    } else if (!is.null(cmax_all) && exposure_var_rr %in% names(cmax_all)) {
      rr_data <- cmax_all %>% dplyr::select(ID, mDose, !!exposure_var_rr)
    } else if  (!is.null(cmin_all) && exposure_var_rr %in% names(cmin_all)) {
      rr_data <- cmin_all %>% dplyr::select(ID, mDose, !!exposure_var_rr)
    } else { 
      stop("Exposure metric for RR not found in simulated data.")
    }
    
    # Prepare newdata for RR prediction
    if (length(covariate_vars_rr) > 0) {
      valid_covs <- covariate_vars_rr[covariate_vars_rr %in% names(erdat)]
      n_pred <- nrow(rr_data)
      cov_sample <- erdat[sample(nrow(erdat), n_pred, replace = TRUE), valid_covs, drop = FALSE]
      newdata_rr <- cbind(
        setNames(data.frame(rr_data[[exposure_var_rr]]), exposure_var_rr),
        cov_sample
      )
    } else {
      newdata_rr <- setNames(data.frame(rr_data[[exposure_var_rr]]), exposure_var_rr)
    }
    
    # Predict RR using model_rr
    rr_pred <- predict(model_rr, newdata = newdata_rr, type = "response")
    rr_bin <- rbinom(length(rr_pred), 1, rr_pred)
    
    # Compute mean RR for each maintenance dose using dplyr, then assign to rr_mat
    rr_df <- tibble(
      mDose = as.numeric(rr_data$mDose),
      rr_bin = rr_bin
    ) %>%
      group_by(mDose) %>%
      summarise(mean_rr = mean(rr_bin), .groups = "drop")
    dose_idx <- match(rr_df$mDose, maintenance_doses)
    rr_mat[j, dose_idx] <- rr_df$mean_rr
    
    # Prepare data for CRS prediction (no cycle iteration, just use the exposure metric required)
    if (!is.null(auc_all) && exposure_var_crs %in% names(auc_all)) {
      crs_data <- auc_all %>% dplyr::select(ID, mDose, !!exposure_var_crs)
      mDose_vec_crs <- auc_all$mDose
      exposure_vec <- auc_all[[exposure_var_crs]]
    } else if (!is.null(cmax_all) && exposure_var_crs %in% names(cmax_all)) {
      crs_data <- cmax_all %>% dplyr::select(ID, mDose, !!exposure_var_crs)
      mDose_vec_crs <- cmax_all$mDose
      exposure_vec <- cmax_all[[exposure_var_crs]]
    } else if (!is.null(cmin_all) && exposure_var_crs %in% names(cmin_all)) {
      crs_data <- cmin_all %>% dplyr::select(ID, mDose, !!exposure_var_crs)
      mDose_vec_crs <- cmin_all$mDose
      exposure_vec <- cmin_all[[exposure_var_crs]]
    }  else {
      stop("Exposure metric for CRS not found in simulated data.")
    }
    
    # Prepare newdata for CRS prediction
    if (length(covariate_vars_crs) > 0) {
      valid_covs <- covariate_vars_crs[covariate_vars_crs %in% names(erdat)]
      n_pred <- length(exposure_vec)
      cov_sample <- erdat[sample(nrow(erdat), n_pred, replace = TRUE), valid_covs, drop = FALSE]
      newdata_crs <- cbind(
        setNames(data.frame(exposure_vec), exposure_var_crs),
        cov_sample
      )
    } else {
      newdata_crs <- setNames(data.frame(exposure_vec), exposure_var_crs)
    }
    
    # Predict CRS using model_crs
    crs_pred <- predict(model_crs, newdata = newdata_crs, type = "response")
    crs_bin <- rbinom(length(crs_pred), 1, crs_pred)
    for (i in seq_along(maintenance_doses)) {
      idx <- which(mDose_vec_crs == maintenance_doses[i])
      if (length(idx) > 0) {
        crs_arr[j, i, 1] <- mean(crs_bin[idx])
      } else {
        crs_arr[j, i, 1] <- NA
      }
    }
  }
  
  # Summarize: median and 90% CI for each dose/cycle
  rr_summary <- tibble(
    mDose = maintenance_doses,
    RR_median = apply(rr_mat, 2, median, na.rm = TRUE),
    RR_lower = apply(rr_mat, 2, quantile, probs = 0.05, na.rm = TRUE),
    RR_upper = apply(rr_mat, 2, quantile, probs = 0.95, na.rm = TRUE)
  )
  
  crs_summary <- map_dfr(1:3, function(cycle) {
    tibble(
      mDose = maintenance_doses,
      CRS_median = apply(crs_arr[,,cycle], 2, median, na.rm = TRUE),
      CRS_lower = apply(crs_arr[,,cycle], 2, quantile, probs = 0.05, na.rm = TRUE),
      CRS_upper = apply(crs_arr[,,cycle], 2, quantile, probs = 0.95, na.rm = TRUE),
      Cycle = paste0("Cycle ", cycle)
    )
  })
  
  # Prepare data for combined plot
  crs_long <- crs_summary %>%
    mutate(Endpoint = Cycle) %>%
    rename(Median = CRS_median, Lower = CRS_lower, Upper = CRS_upper) %>%
    dplyr::select(mDose, Median, Lower, Upper, Endpoint)
  
  rr_long <- rr_summary %>%
    mutate(Endpoint = "RR") %>%
    rename(Median = RR_median, Lower = RR_lower, Upper = RR_upper) %>%
    dplyr::select(mDose, Median, Lower, Upper, Endpoint)
  
  plot_data <- bind_rows(
    crs_long,
    rr_long
  ) %>%
    mutate(
      Endpoint = factor(Endpoint, levels = c("RR", "Cycle 1", "Cycle 2", "Cycle 3"),
                        labels = c("Response Rate", "CRS Cycle 1", "CRS Cycle 2", "CRS Cycle 3"))
    )
  
  # Plot RR and CRS in one panel
  p_combined <- ggplot(plot_data, aes(x = mDose, y = Median, color = Endpoint, fill = Endpoint, group = Endpoint)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
    theme_minimal() +
    labs(
      x = "Maintenance Dose (mg)",
      y = "Probability",
      title = "Predicted Response Rate and CRS Probability vs Maintenance Dose (90% CI)",
      color = "Endpoint",
      fill = "Endpoint"
    ) +
    scale_color_manual(values = c("Response Rate" = "#377eb8", "CRS Cycle 1" = "#e41a1c")) +
    scale_fill_manual(values = c("Response Rate" = "#377eb8", "CRS Cycle 1" = "#e41a1c")) +
    scale_y_continuous(limits = c(0, 1), breaks=seq(0, 1, by=0.1)) + geom_hline(yintercept=0.448, linetype="dashed", size=1.5) + 
    scale_x_continuous(breaks=seq(10, 90, by=5))
  
  print(p_combined)
  
  ggsave(paste0(path, "/output/predicted_rr_crs_vs_maintenance_dose_90CI_combined.png"), 
         plot = p_combined, width = 7, height = 5, dpi = 300)


###############################################################################
#
# End of program
#
###############################################################################