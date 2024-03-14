
randsp<-randsp[104:106,]
nsample=3

nmrs = as.list(0)
library(foreach)
library(doFuture)
doFuture::registerDoFuture()
plan(cluster, workers = 3)


nmcl =  matrix(0, nrow = npar, ncol = nsample)
#.export = c("CITY", "names.e", "names.msm", "names.pwid", "randsp", "obj", "calpar", "calpar.info", "vparameters", "vlist", "calib.target", "valid.target", "names.gp", "names18", "state.name", "par_info", "gp18.gn", "x", "vt", "ode_model", "n", "end_yr_ind", "black", "msm", "rname", "group11", "midu", "all.msm", "moplot", "sigmaO.FM", "sigmaO.MF", "sigmaS", "tau", "mor_S", "mor_I2", "mor_I3", "mor_T1", "mor_T2", "mor_T3", "O_T", "T1_T2", "T1_T3", "T2_T1", "T2_T3", "T3_T1", "T3_T2", "T1_O", "T2_O", "T3_O", "rho.m", "rho", "mat", "init.tot", "nOAT", "s", "FoI", "ass.eO", "ass.eS", "eta.m", "het", "ode_list", "oat", "off.oat", "ode_list_offOAT", "ode_list_OAT", "init.group.prop", "m", "het.m.l", "f", "het.f.l", "msm.l", "white", "hisp", "f.high", "m.high", "het.l", "all.idu", "plot", "yr", "nyr", "FoI", "gp.gn.fun", "gp18.gn.fun", "group11", "model.initial", "moplot", "msm.gn.fun", "ode_list", "ode_list_OAT", "ode_list_offOAT", "ode_model", "prepentry", "pwid.gn.fun")
nmcl <- foreach(i=1:nsample, .combine=cbind, .errorhandling = 'remove', .export = c("msm.h","inits", "CITY", "names.e", "names.msm", "names.pwid", "randsp", "obj", "calpar", "calpar.info", "vparameters", "vlist", "calib.target", "valid.target", "names.gp", "names18", "state.name", "par_info", "gp18.gn", "x", "vt", "ode_model", "n", "end_yr_ind", "black", "msm", "rname", "group11", "midu", "all.msm", "init.tot", "s", "FoI", "ass.eO", "ass.eS", "het", "ode_list", "oat", "off.oat", "ode_list_offOAT", "ode_list_OAT", "init.group.prop", "m", "het.m.l", "f", "het.f.l", "msm.l", "white", "hisp", "f.high", "m.high", "het.l", "all.idu", "plot", "yr", "nyr", "FoI", "gp.gn.fun", "gp18.gn.fun", "group11", "model.initial", "msm.gn.fun", "ode_list", "ode_list_OAT", "ode_list_offOAT", "ode_model", "pwid.gn.fun", "euler")) %dopar% {
  nmrs <- dfoptim::nmkb(par =randsp[i, ], fn =obj,
                        lower = calpar$lower, upper = calpar$upper,
                        calpar.info = calpar.info, fixed = vparameters, fixed.list =vlist,
                        calib.target = calib.target, valid.target = valid.target)
  
  nmrs$par
}

write.excel(nmcl)