library(devtools)
devtools::install_github("kaishumason/NicheDE")


Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
install.packages("devtools")
options(timeout=9999999)
devtools::install_github("kaishumason/NicheDE") # install
library(nicheDE)
library(Matrix)
library(doParallel)


NDE_obj_m <- readRDS('20241007NDEobj_m_postEffniche.rds')

NDE_obj_m = niche_DE(NDE_obj_m,num_cores = 8, outfile = "",C = 30, M = 10, gamma = 0.8,print = T, Int = T, batch = T,self_EN = F,G = 0.001)
saveRDS(NDE_obj, '20241007NDEobj_m_postNDE.rds')