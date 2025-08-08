###########################
# ↓↓↓↓ Do not change ↓↓↓↓ #
###########################
# delete all variables
rm(list=ls())

# check working directory
getwd()

# library
library(tidyverse)
library(mrgsolve)
library(xpose)

# set path for NONMEM run folder
path_run_nonmem <- "./../02_model" # use '/' or '//' instead of '/'

# set path for NONMEM data file (pkpd02.csv)
path_data_nonmem <- "./../01_data/pkpd02.csv"
###########################
# ↑↑↑↑ Do not change ↑↑↑↑ #
###########################

######################################
# ↓↓↓↓ Please change this block ↓↓↓↓ #
######################################
# set run number of final PKPD model
run.no <- "001"
######################################
# ↑↑↑↑ Please change this block ↑↑↑↑ #
######################################

###########################
# ↓↓↓↓ Do not change ↓↓↓↓ #
###########################
# ---- import xpose data ----
xpdb_pk <- xpose_data(runno = "000", dir = path_run_nonmem) # PK model
xpdb_pd <- xpose_data(runno = run.no, dir = path_run_nonmem) # PKPD model

obs <- read_csv(path_data_nonmem)

# ---- read template model file (mrgsolve) ----
mod <- mread("./00_mrgsolve_model/mod_template_pkpd.cpp")

# ---- update to parameter estimates ----
d_pd <- get_prm(xpdb_pd, transform = F) 
d_pk <- get_prm(xpdb_pk, transform = F)

mod_est <- update(mod, param = list(TVCL = d_pk$value[d_pk$label=="TVCL"],
                                    TVVC = d_pk$value[d_pk$label=="TVV1"],
                                    TVQ = d_pk$value[d_pk$label=="TVQ"],
                                    TVVP = d_pk$value[d_pk$label=="TVV2"],
                                    TVEMAX = d_pd$value[d_pd$label=="TVEMAX"],
                                    TVEC50 = d_pd$value[d_pd$label=="TVEC50"],
                                    TVHILL = d_pd$value[d_pd$label=="TVHILL"],
                                    TVIMAX = d_pd$value[d_pd$label=="TVIMAX"],
                                    TVIC50 = d_pd$value[d_pd$label=="TVIC50"],
                                    TVKDEG = d_pd$value[d_pd$label=="TVKDEG"],
                                    TVKPRIM = d_pd$value[d_pd$label=="TVKPRIM"],
                                    TVMTT = d_pd$value[d_pd$label=="TVMTT"]),
                  omega = list(dmat(d_pk$value[d_pk$label=="IIV_CL"],
                                    d_pk$value[d_pk$label=="IIV_V1"],
                                    d_pk$value[d_pk$label=="IIV_Q"],
                                    d_pk$value[d_pk$label=="IIV_V2"],
                                    d_pd$value[d_pd$label=="IIV_EMAX"],
                                    d_pd$value[d_pd$label=="IIV_EC50"],
                                    d_pd$value[d_pd$label=="IIV_HILL"],
                                    d_pd$value[d_pd$label=="IIV_IMAX"],
                                    d_pd$value[d_pd$label=="IIV_IC50"],
                                    d_pd$value[d_pd$label=="IIV_KDEG"],
                                    d_pd$value[d_pd$label=="IIV_KPRIM"],
                                    0,  # OMEGA IL6BASE is not estimated in B2 method.
                                    d_pd$value[d_pd$label=="IIV_MTT"])),
                  sigma = list(dmat(d_pk$value[d_pk$label=="PROP_ERR"]^2,  # change scale to variance
                                    d_pd$value[d_pd$label=="PROP_ERR_IL6"]^2))) # change scale to variance

# ---- save files ----
saveRDS(mod_est, paste0("./00_mrgsolve_model/mod_pkpd_run", run.no, ".rds"))
###########################
# ↑↑↑↑ Do not change ↑↑↑↑ #
###########################
