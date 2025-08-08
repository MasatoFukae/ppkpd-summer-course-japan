###########################
# ↓↓↓↓ Do not change ↓↓↓↓ #
###########################
# delete all variables
rm(list=ls())

# check working directory
getwd()

# library
library(gridExtra)
library(tidyverse)
library(mrgsolve)
library(xpose)

# set path
path <- "./.."

# set seed
seed <- 0702
###########################
# ↑↑↑↑ Do not change ↑↑↑↑ #
###########################

######################################
# ↓↓↓↓ Please change this block ↓↓↓↓ #
######################################
# set run number of final PKPD model
run.no <- "001"

# ---- Simulation settings ----
target <- c() # target dose (mg)

target.interval <- 24 * 21 # dose interval (h), default is Q3W
target.addl <- 3 # additional dose of target dose, default is 2

nsim <- 100  # number of patients/regimen

doses.regimen <- expand.grid(target = target)
######################################
# ↑↑↑↑ Please change this block ↑↑↑↑ #
######################################

###########################
# ↓↓↓↓ Do not change ↓↓↓↓ #
###########################
# ---- import xpose data ----
xpdb_pd <- xpose_data(runno = run.no, dir = paste0(path, "/02_model"))  ## PKPD model
xpdb_pk <- xpose_data(runno = '000',dir = paste0(path, "/02_model"))  ## PPK model

obs <- read.csv(paste0(path, "/01_data/pkpd02.csv"), sep=",")

# ---- read model file (mrgsolve) ----
mod_est <- readRDS(paste0("./00_mrgsolve_model/mod_pkpd_run", run.no, ".rds"))

# ---- create template dataset for simulation ----
d.target <- expand.ev(ID = 1:nsim,
                     cmt = 1,
                     amt = doses.regimen$target,
                     addl = target.addl,
                     ii = target.interval) %>%
  mutate(rate = amt)

d0 <- d.target 

d0.regimen <- d0 %>%
  distinct(ID, amt) %>%
  dplyr::rename(Regimen = amt)

# ---- simulation ----
# Resampling baseline IL6
obs1 <- obs %>%
  distinct(ID, .keep_all = T)
bsl.value = obs1$IL6BSL

d1 <- d0 %>%
  group_by(ID) %>%
  mutate(BSL = sample(bsl.value, 1))

# Simulation
set.seed(seed)

sim1 <- mod_est %>%
  env_eval(seed = seed) %>%
  mrgsim(data=d1, end= target.interval * (target.addl+1), delta=1) %>%
  as.data.frame() %>%
  left_join(d0.regimen, by = "ID")

df <- sim1 %>%
  dplyr::select(ID, Regimen, time, DV_CP, DV_IL6) %>%
  pivot_longer(cols = -c("ID", "Regimen", "time")) %>%
  group_by(Regimen, time, name) %>%
  summarise(Median = median(value, na.rm = T),
            Lower = quantile(value, 0.05, na.rm = T),
            Upper = quantile(value, 0.95, na.rm = T),
            .groups = "drop") %>%
  mutate(
    day = time / 24,
    Regimen = str_replace(Regimen, ",", ", ")
  )

# ---- Save simulation result ----
write_csv(df, "simulation_flat-dosing.csv")

# ---- Output to PDF file ----
pdf("simulation_flat-dosing.pdf", paper = "a4r", width = 20)

for (i in 1:length(unique(df$Regimen))) {
  p <- ggplot(df %>% filter(Regimen==unique(df$Regimen)[i]), 
         aes(x=time, y=Median, group=Regimen)) +
    geom_line() +
    geom_ribbon(aes(ymin=Lower, ymax=Upper), alpha=0.5, fill="darkblue") +
    facet_wrap(~name+Regimen, scales = "free_y", ncol=1)
  print(p)
}


dev.off()
###########################
# ↑↑↑↑ Do not change ↑↑↑↑ #
###########################