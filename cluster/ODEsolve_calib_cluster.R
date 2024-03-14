##################################################################################
# An R script to solve ODE's of the HIV model using deSolve package.             #
#                                                                                #
# Author:  Jeong Min, Xiao Zang, Emanuel Krebs                                   #
# updated: March 19, 2018                                                        #
##################################################################################

#A pdf document that provides very detailed information about the package "deSolve"
#vignette("deSolve")

## Module set-up ##
rm(list=ls())
# 
# library(rstudioapi)
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dfoptim) #for calibration function nmkb
library(openxlsx)

source("ode_model_func_v11.R")   #ode function module
source("obj_func_v11.R")         #objective function module


## Read in data and all inputs ##

#### SET CITY ####
CITY = "NYC"
#### SET CITY ####

source ("Data_input_v6.R")

#### Random sampling to obtain multiple starting values ####
library(prevalence)
library(mc2d)
library(MCMCpack)

###### Latin hypercube sampling to get 1,000 simplexes for calibration######
#install.packages("FME")
library(FME)

nsample = 10            #currently set 10 to reduce computing time
top = 5
npar    = nrow (calpar)

set.seed(16574)
# generates Latin hypercube samples from 90% interval
lhs = Latinhyper(matrix(c(0.05,0.95), npar, 2, byrow=T), nsample)

randsp = matrix(0, nsample, npar)
for (i in 1:npar){
  if (calpar$no.par[i] != 3){
    if (calpar$dist[i] == "beta"){
      randsp [ ,i] = qbeta(lhs[ ,i], calpar[i, ]$par1, calpar[i, ]$par2)
    }
    else if (calpar$dist[i] == "pert") {
      randsp [ ,i] = qbetagen(lhs[ ,i], calpar[i, ]$par1, calpar[i, ]$par2, calpar[i, ]$lower, calpar[i, ]$upper)
    }
    else if (calpar$dist[i] == "1/pert") {
      randsp [ ,i] = 1/(qbetagen(lhs[ ,i], calpar[i, ]$par1, calpar[i, ]$par2, 1/calpar[i, ]$upper, 1/calpar[i, ]$lower))
    }
    else if (calpar$dist[i] == "gamma") {
      randsp [ ,i] = qgamma(lhs[ ,i], calpar[i, ]$par1, 1/calpar[i, ]$par2)
    }
    else if (calpar$dist[i] == "uni") {
      randsp [ ,i] = qunif(lhs[ ,i], calpar[i, ]$lower, calpar[i, ]$upper)
    }
    else if (calpar$dist[i] == "1/uni") {
      randsp [ ,i] = 1/(qunif(lhs[ ,i], 1/calpar[i, ]$upper, 1/calpar[i, ]$lower))
    }
    else if (calpar$dist[i] == "ln") {
      randsp [ ,i] = qlnorm(lhs[ ,i], calpar[i, ]$par1, calpar[i, ]$par2)
    }
    else if (calpar$dist[i] == "1-lognormal") {
      randsp [ ,i] = 1- (qlnorm(lhs[ ,i], calpar[i, ]$par1, calpar[i, ]$par2))
    }
    else if (calpar$dist[i] == "1/ln") {
      randsp [ ,i] = 1/(qlnorm(lhs[ ,i], calpar[i, ]$par1, calpar[i, ]$par2))
    }
    else if (calpar$dist[i] == "poisson") {
      randsp [ ,i] = qpois(lhs[ ,i], calpar[i, ]$par1) / calpar[i, ]$par2
    }
  }
  else if (calpar$dist[i] == "dirichlet_sp"){
    dirsm        = rdirichlet(nsample, as.numeric(c(calpar[i, ]$par1, calpar[i, ]$par2, calpar[i, ]$par3)))
    randsp [ ,i] = dirsm[,1]*0.3 + dirsm[ ,2]*0 + dirsm[ ,3]*1
  }
}


## Multiple calibrations ##
nmrs = as.list(0)

#### FOR PARALLEL ####
library(foreach)
library(doParallel)

ncores = Sys.getenv("SLURM_NTASKS")
registerDoParallel(cores=ncores)

# library(doFuture)
# doFuture::registerDoFuture()
# 
# library("future.batchtools")
# plan(batchjobs_slurm)

# cores=detectCores()
# plan(cluster, workers = cores[1]-3)

# cores = Sys.getenv("SLURM_NTASKS") 
# plan(cluster, workers = cores)

nmcl =  matrix(0, nrow = npar, ncol = nsample)
#.export = c("CITY", "names.e", "names.msm", "names.pwid", "randsp", "obj", "calpar", "calpar.info", "vparameters", "vlist", "calib.target", "valid.target", "names.gp", "names18", "state.name", "par_info", "gp18.gn", "x", "vt", "ode_model", "n", "end_yr_ind", "black", "msm", "rname", "group11", "midu", "all.msm", "moplot", "sigmaO.FM", "sigmaO.MF", "sigmaS", "tau", "mor_S", "mor_I2", "mor_I3", "mor_T1", "mor_T2", "mor_T3", "O_T", "T1_T2", "T1_T3", "T2_T1", "T2_T3", "T3_T1", "T3_T2", "T1_O", "T2_O", "T3_O", "rho.m", "rho", "mat", "init.tot", "nOAT", "s", "FoI", "ass.eO", "ass.eS", "eta.m", "het", "ode_list", "oat", "off.oat", "ode_list_offOAT", "ode_list_OAT", "init.group.prop", "m", "het.m.l", "f", "het.f.l", "msm.l", "white", "hisp", "f.high", "m.high", "het.l", "all.idu", "plot", "yr", "nyr", "FoI", "gp.gn.fun", "gp18.gn.fun", "group11", "model.initial", "moplot", "msm.gn.fun", "ode_list", "ode_list_OAT", "ode_list_offOAT", "ode_model", "prepentry", "pwid.gn.fun")
nmcl <- foreach(i=1:nsample, .combine=cbind, .export = c("msm.h","inits", "CITY", "names.e", "names.msm", "names.pwid", "randsp", "obj", "calpar", "calpar.info", "vparameters", "vlist", "calib.target", "valid.target", "names.gp", "names18", "state.name", "par_info", "gp18.gn", "x", "vt", "ode_model", "n", "end_yr_ind", "black", "msm", "rname", "group11", "midu", "all.msm", "init.tot", "FoI", "ass.eO", "ass.eS", "het", "ode_list", "oat", "off.oat", "ode_list_offOAT", "ode_list_OAT", "init.group.prop", "m", "het.m.l", "f", "het.f.l", "msm.l", "white", "hisp", "f.high", "m.high", "het.l", "all.idu", "yr", "nyr", "gp.gn.fun", "gp18.gn.fun", "group11", "model.initial", "msm.gn.fun", "ode_list", "ode_list_OAT", "ode_list_offOAT", "ode_model", "pwid.gn.fun", "euler")) %dopar% {
  nmrs <- dfoptim::nmkb(par =randsp[i, ], fn =obj,
                        lower = calpar$lower, upper = calpar$upper,
                        calpar.info = calpar.info, fixed = vparameters, fixed.list =vlist,
                        calib.target = calib.target, valid.target = valid.target)
  
  nmrs$par
}

# write.excel(nmcl)

gof = numeric(nsample)

for (i in 1:nsample){
  gof[i] = obj(calib.par =nmcl[,i], calpar.info = calpar.info, fixed = vparameters, fixed.list =vlist, 
               calib.target = calib.target, valid.target = valid.target)
}

gof2000ord <- order(gof)[1:top]

nmrs2000 <- nmcl[ , gof2000ord]
gof2000 <- gof[gof2000ord]     

outcome <- list(nmcl = nmcl,
                gof  = gof,
                nmrs2000 = nmrs2000,
                gof2000  = gof2000)

save(outcome, file = "outputs.RData")

