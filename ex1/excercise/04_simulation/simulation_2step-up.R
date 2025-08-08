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
step1 <- c() # 1st step-up dose (mg)
step2 <- c() # 2nd step-up dose (mg)
target <- c() # target dose (mg)

step1.time <- 0 # dosing time of 1st step-up dose (h)
step2.time <- 24 * 21 # dosing time of 2nd step-up dose (h), default is Q3W
target.time <- 24 * 21 # dosing time of 1st target dose (h), default is Q3W
target.interval <- 24 * 21 # dose interval (h), default is Q3W
target.addl <- 1 # additional dose of target dose, default is 2

nsim <- 100  # number of patients/regimen

doses.regimen <- expand.grid(step1 = step1,
                            step2 = step2,
                            target = target)
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
mod_est = readRDS(paste0("./00_mrgsolve_model/mod_pkpd_run", run.no, ".rds"))

# ---- create template dataset for simulation ----
d.step1 <- expand.ev(ID = 1:nsim,
                    cmt = 1,
                    amt = doses.regimen$step1,
                    addl = 0,
                    ii = 0) %>%
  mutate(rate = amt) %>%
  mutate(time = step1.time)

d.step2 <- expand.ev(ID = 1:nsim,
                    cmt = 1,
                    amt = doses.regimen$step2,
                    addl = 0,
                    ii = 0) %>%
  mutate(rate = amt) %>%
  mutate(time = step2.time)

d.target <- expand.ev(ID = 1:nsim,
                     cmt = 1,
                     amt = doses.regimen$target,
                     addl = target.addl,
                     ii = target.interval) %>%
  mutate(rate = amt) %>%
  mutate(time = step2.time + target.time)

d0 <- bind_rows(d.step1, d.step2, d.target) %>%
  arrange(ID, time) %>%
  group_by(ID) %>%
  mutate(Regimen = paste(amt, collapse = ",")) %>%
  ungroup()

d0.regimen <- d0 %>%
  distinct(ID, Regimen)

# ---- simulation ----
# Resampling baseline IL6
obs1 = obs %>%
  distinct(ID, .keep_all = T)
bsl.value = obs1$IL6BSL

d1 = d0 %>%
  group_by(ID) %>%
  mutate(BSL = sample(bsl.value, 1))

# Simulation

set.seed(seed)

sim1 <- mod_est %>%
  env_eval(seed = seed) %>%
  mrgsim(data=d1, end=step2.time + target.time + target.interval * (target.addl+1), delta=1) %>%
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
write_csv(df, "simulation_2step-up.csv")

# ---- Output to PDF file ----
pdf("simulation_2step-up.pdf", paper = "a4r", width = 20)

for (i in 1:length(unique(df$Regimen))) {
  p = ggplot(df %>% filter(Regimen==unique(df$Regimen)[i]), 
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