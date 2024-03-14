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

library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dfoptim) #for calibration function nmkb
library(openxlsx)

source("ode_model_func_v11.R")   #ode function module
source("obj_func_v11_Seattle.R")         #objective function module

# # Set function to copy result to excel
write.excel <- function(tab, ...) write.table( tab, "clipboard", sep="\t", row.names=F)


## Read in data and all inputs ##

#### SET CITY ####
CITY = "SEA"
#### SET CITY ####

source ("Data_input_v6.R")
# output: vparameters, vt (time step), x (all initials in vector), calpar (calibration data), target data

#source("print.inc.R"); incplot(calib.par =calpar$pe, calpar.info =calpar.info, fixed =vparameters, fixed.list =vlist, calib.target =calib.target, valid.target =valid.target, plot =T)

obj (calib.par =calpar$pe, calpar.info =calpar.info, fixed =vparameters, fixed.list =vlist,
     calib.target =calib.target, valid.target =valid.target, plot =T)

# plot.par = calpar$pe
# source("full.plot.R")

###### Model calibration ######
#Nelder-Mead should not be used for high-dimensional optimization
# maximum # of parameters for nmk is 30
nm <- nmkb(par = calpar$pe, fn = obj,
          lower = calpar$lower, upper = calpar$upper,
          calpar.info = calpar.info, fixed = vparameters, fixed.list =vlist,
          calib.target = calib.target, valid.target = valid.target, plot=F)
# nm <- nmkb(par = calpar$pe, fn = obj,
#            lower = calpar$lower, upper = calpar$upper,
#            calpar.info = calpar.info, fixed = vparameters, fixed.list =vlist,
#            calib.target = calib.target, valid.target = valid.target)
# 
obj (calib.par = nm$par, calpar.info = calpar.info, fixed = vparameters, fixed.list =vlist,
     calib.target = calib.target, valid.target = valid.target, plot=T)
# # 
write.excel(nm$par)
# 
# plot.par = nm$par
# 
# source("full.plot.R")
# #source("print.inc.R"); incplot(calib.par =nm$par, calpar.info =calpar.info, fixed =vparameters, fixed.list =vlist, calib.target =calib.target, valid.target =valid.target, plot =T)

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

# any(randsp - matrix(rep(calpar$lower, each=nsample), nrow = nsample, ncol=npar) < 0 )  # random sample check
# any(matrix(rep(calpar$upper, each=nsample), nrow = nsample, ncol=npar) - randsp < 0 )  # random sample check
# 
# print(randsp)
# write.excel(randsp)
# 
# objval = numeric(nsample)
# for (i in 1:nsample){
#   objval[i] = obj(calib.par =randsp[i, ], calpar.info = calpar.info, fixed = vparameters, fixed.list =vlist,
#                   calib.target = calib.target, valid.target = valid.target, plot=F)
# }
# write.excel(objval)

## Multiple calibrations ##
nmrs = as.list(0)

#### FOR PARALLEL ####
library(foreach)
library(doParallel)

#ncores = Sys.getenv("SLURM_NTASKS")
registerDoParallel(cores=5)


# library(doFuture)
# doFuture::registerDoFuture()
# 
# library("future.batchtools")
# plan(batchjobs_slurm)

# cores=detectCores()
# plan(cluster, workers = cores[1]-3)

# library(doParallel)
# cores = Sys.getenv("SLURM_NTASKS")
# plan(cluster, workers = cores)



nmcl =  matrix(0, nrow = npar, ncol = nsample)
#.export = c("CITY", "names.e", "names.msm", "names.pwid", "randsp", "obj", "calpar", "calpar.info", "vparameters", "vlist", "calib.target", "valid.target", "names.gp", "names18", "state.name", "par_info", "gp18.gn", "x", "vt", "ode_model", "n", "end_yr_ind", "black", "msm", "rname", "group11", "midu", "all.msm", "moplot", "sigmaO.FM", "sigmaO.MF", "sigmaS", "tau", "mor_S", "mor_I2", "mor_I3", "mor_T1", "mor_T2", "mor_T3", "O_T", "T1_T2", "T1_T3", "T2_T1", "T2_T3", "T3_T1", "T3_T2", "T1_O", "T2_O", "T3_O", "rho.m", "rho", "mat", "init.tot", "nOAT", "s", "FoI", "ass.eO", "ass.eS", "eta.m", "het", "ode_list", "oat", "off.oat", "ode_list_offOAT", "ode_list_OAT", "init.group.prop", "m", "het.m.l", "f", "het.f.l", "msm.l", "white", "hisp", "f.high", "m.high", "het.l", "all.idu", "plot", "yr", "nyr", "FoI", "gp.gn.fun", "gp18.gn.fun", "group11", "model.initial", "moplot", "msm.gn.fun", "ode_list", "ode_list_OAT", "ode_list_offOAT", "ode_model", "prepentry", "pwid.gn.fun")
nmcl <- foreach(i=1:nsample, .combine=cbind, .export = c("msm.h","inits", "CITY", "names.e", "names.msm", "names.pwid", "randsp", "obj", "calpar", "calpar.info", "vparameters", "vlist", "calib.target", "valid.target", "names.gp", "names18", "state.name", "par_info", "gp18.gn", "x", "vt", "ode_model", "n", "end_yr_ind", "black", "msm", "rname", "group11", "midu", "all.msm", "init.tot", "s", "FoI", "ass.eO", "ass.eS", "het", "ode_list", "oat", "off.oat", "ode_list_offOAT", "ode_list_OAT", "init.group.prop", "m", "het.m.l", "f", "het.f.l", "msm.l", "white", "hisp", "f.high", "m.high", "het.l", "all.idu", "plot", "yr", "nyr", "FoI", "gp.gn.fun", "gp18.gn.fun", "group11", "model.initial", "msm.gn.fun", "ode_list", "ode_list_OAT", "ode_list_offOAT", "ode_model", "pwid.gn.fun", "euler")) %dopar% {
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

outcome

# obj (calib.par = nm$par, calpar.info = calpar.info, fixed = vparameters, fixed.list =vlist,
#      calib.target = calib.target, valid.target = valid.target, plot=T)

#stopCluster(cl)

#### END PARALLEL ####


# for (i in 1:nsample){
#   nmrs[[i]] <- nmkb(par =randsp[i, ], fn =obj,
#                     lower = calpar$lower, upper = calpar$upper,
#                     calpar.info = calpar.info, fixed = vparameters, fixed.list =vlist,
#                     calib.target = calib.target, valid.target = valid.target, plot=F)
# }


# #### Multiple calibration results
# source("calOut.R")
# 
# parrs = matrix(0, npar, nsample)
# objrs = numeric(nsample)
# 
# calOut10 = matrix(0, nsample, 164)   # ((11 groups * 3 targets) + 6 calib targets + 2 valid targets) * 4 years
# for (i in 1:nsample){
#   parrs[ ,i] = nmcl[,i]
#   calOut10[i, ] = calOut(calib.par =nmcl[,i], calpar.info = calpar.info, fixed = vparameters, fixed.list =vlist)
# }
# 
# diag11.model = calOut10[ , 1 : 44]
# ndiag.all    = calOut10[ ,133:136]
# ndiag.b      = calOut10[ ,137:140]
# ndiag.m      = calOut10[ ,141:144]
# death.all    = calOut10[ ,145:148]
# death.b      = calOut10[ ,149:152]
# death.m      = calOut10[ ,153:156]
# 
# ## Target operationalization
# diag.obs  = calib.target$diag18.obs [-1, ]
# ndiag.obs = calib.target$ndiag18.obs[-1, ]
# death.obs = calib.target$death18.obs[-1, ]
# inc.all   = valid.target$obs.inc.all
# inc.msm   = valid.target$obs.inc.msm
# 
# # observed total diagnosed PLHIV
# diag11.obs   = group11(diag.obs)
# names11 = c(names18[1:3], "MSM/PWID", "PWID", names18[13:18])
# 
# # observed new diagnosis
# ndiag.all.obs= rowSums(ndiag.obs)
# ndiag.b.obs  = rowSums(ndiag.obs[ ,grep("black", names18)]) # black
# ndiag.m.obs  = rowSums(ndiag.obs[ ,1:3]) #msm
# ndiag11.obs  = group11(ndiag.obs)
# 
# # observed mortality
# if (CITY != "MIA"){
#   death.all.obs= rowSums(death.obs)
#   death.b.obs  = rowSums(death.obs[ ,grep("black", names18)]) # black
#   death.m.obs  = rowSums(death.obs[ ,1:3]) #msm
# } else {
#   death.all.obs= death.obs[,1]
#   death.b.obs  = death.obs[,3] # black
# }
# 
# if (CITY != "MIA"){
#   obs   = c(diag11.obs,   ndiag.all.obs, ndiag.b.obs, ndiag.m.obs, death.all.obs, death.b.obs, death.m.obs)
# } else {
#   obs   = c(diag11.obs,   ndiag.all.obs, ndiag.b.obs, ndiag.m.obs, death.all.obs, death.b.obs)
# }
# 
# 
# ##GoF results
# for (i in 1: nsample){
#   if (CITY != "MIA"){
#     model = c(diag11.model[i, ], ndiag.all[i, ], ndiag.b[i, ], ndiag.m[i, ], death.all[i, ], death.b[i, ], death.m[i, ])
#   } else {
#     model = c(diag11.model[i, ], ndiag.all[i, ], ndiag.b[i, ], ndiag.m[i, ], death.all[i, ], death.b[i, ])
#   }
#   dev   = abs((model-obs)/obs)
#   objrs[i] = weighted.mean(dev, rep(vparameters$w, each=4), na.rm=T)
# }
# write.excel(parrs)
# print(objrs)
# write.excel(objrs)


# #### plots for the multiple calibrated parameter sets
# 
# # function to plot model outcomes compared with observed values
# moplot = function(model, obs, title, low =0, high =0){
#   ymin = min(c(model, obs), na.rm =T)*0.8
#   ymax = max(c(model, obs), na.rm =T)*1.2
#   if (any(low>0)){
#     ymin = min(c(model, low),  na.rm=T)*0.8
#     ymax = max(c(model, high), na.rm=T)*1.2
#   }
#   # observed output
#   plot(yr, obs, xlab ="Year", main =title,
#        ylab = 'Number of individuals', type ='l', ylim =c(ymin,ymax), lwd =2,
#        xaxt = "n") #xaxis label not shown
#   if (nyr <= 4) axis(side = 1, at = yr)
#   else axis(side =1)
#   # from model output
#   for (i in 1:nsample){
#     lines(yr, model[i, ], lty =2, lwd =1.5, col="blue")
#   }
#   # CI for incidence
#   if (any(low >0)){
#     lines(yr[1:4], low,  lty=3, lwd=2)
#     lines(yr[1:4], high, lty=3, lwd=2)
#   }
# }
# 
# # ####Result plots for all groups####
# # par(oma = c(3, 3, 3, 3))
# # par(mar =  c(2, 2, 2, 1))
# # par(mfrow=c(5, 7))
# # 
# # for (i in 1:11){
# #   moplot(calOut10[ ,(4*i-3):(4*i)],     diag11.obs[ ,i],  paste("Diag:", names11[i]))
# # }
# # for (i in 1:11){
# #   moplot(calOut10[ ,(4*i+41):(4*i+44)], ndiag11.obs[ ,i], paste("New diag:", names11[i]))
# # }
# # for (i in 1:11){
# #   moplot(calOut10[ ,(4*i+85):(4*i+88)], death11.obs[ ,i], paste("Death:", names11[i]))
# # }
# # #incidence
# # moplot(calOut10[ ,157:160], obs.inc.all$value, "Total HIV incidence",           obs.inc.all$low, obs.inc.all$high)
# # moplot(calOut10[ ,161:164], obs.inc.msm$value, "HIV incidence: MSM & MSM/PWID", obs.inc.msm$low, obs.inc.msm$high)
# # #### legend ####
# # par(fig=c(0, 1, 0, 0.3), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
# # plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
# # legend("bottom", c("Model","Obs"), lty=c(2,1), cex=1, lwd=2, col=c("blue","black"))
# 
# 
# ####Calibration/validation plots####
# par(mfrow=c(5, 4))
# par(oma =  c(3, 3, 3, 3))
# par(mar =  c(2, 2, 2, 1))
# for (i in 1:11){
#   moplot(calOut10[ ,(4*i-3):(4*i)],     diag11.obs[ ,i],  paste("Diag:", names11[i]))
# }
# # new diagnosis
# moplot(ndiag.all, ndiag.all.obs, "New diag: total")
# moplot(ndiag.b,   ndiag.b.obs,   "New diag: black")
# moplot(ndiag.m,   ndiag.m.obs,   "New diag: MSM")
# # death
# moplot(death.all, death.all.obs, "Death: total")
# moplot(death.b,   death.b.obs,   "Death: black")
# 
# if (CITY != "MIA"){
#   moplot(death.m,   death.m.obs,   "Death: MSM")
#   }
# #incidence
# moplot(calOut10[ ,157:160], obs.inc.all$value, "Total HIV incidence",                  obs.inc.all$low, obs.inc.all$high)
# moplot(calOut10[ ,161:164], obs.inc.msm$value, "HIV incidence: MSM & MSM/PWID", obs.inc.msm$low, obs.inc.msm$high)
# #### legend ####
# par(fig=c(0, 1, 0, 0.3), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
# plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
# legend("bottom", c("Model","Obs"), lty=c(2,1), cex=1, lwd=2, col=c("blue","black"))