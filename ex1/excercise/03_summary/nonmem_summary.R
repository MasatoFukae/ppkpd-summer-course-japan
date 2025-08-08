
rm(list=ls())
library(tidyverse)
library(gridExtra)
library(xpose)

###############################################################################
# import xpose data
###############################################################################

path = ""

run.no <- ""  # specify run number
xpdb <- xpose_data(runno = run.no, dir = paste0(path, "/02_model"))

###############################################################################
# Parameter estimates
###############################################################################

tab = get_prm(xpdb, transform = F) %>%
  mutate(value = format(value, scientific=F, digits=2),
         se = format(se, scientific=F, digits=2),
         rse = format(rse, scientific=F, digits=2))

###############################################################################
# GOF plot
###############################################################################

p1 <- dv_vs_pred(xpdb) + labs(title = "DV vs PRED", subtitle = NULL, caption = NULL)
p2 <- dv_vs_ipred(xpdb) + labs(title = "DV vs IPRED", subtitle = NULL, caption = NULL)
p3 <- res_vs_pred(xpdb, res = "NPDE") + labs(title = "NPDE vs PRED", subtitle = NULL, caption = NULL)
p4 <- res_vs_idv(xpdb, res = "NPDE") + labs(title = "NPDE vs TIME", subtitle = NULL, caption = NULL)

###############################################################################
# Individual plot
###############################################################################

ind1 <- ind_plots(xpdb, facets = c("ID")) + labs(subtitle = NULL, caption = NULL)

###############################################################################
# Output to PDF file
###############################################################################

pdf(paste0(path, "/03_summary/nonmem_summary_run", run.no, ".pdf"), paper = "a4r", width = 20)

grid.table(tab)
grid.arrange(p1, p2, p3, p4, nrow=2)
print(ind1)

dev.off()
