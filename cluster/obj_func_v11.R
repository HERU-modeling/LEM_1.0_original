################################################################################
## Objectove function                                                         ##
################################################################################
# objective function considering new diagnosis, PLHIV observed, number of death.
# return objective function to minimize

## Parameters ##
  # calib.par is a vector of parameters to be calibrated
  # leng: length of each parameter (i.e. noG has a length of 36)
  # names: parameter names for calib.par
  # (names of parameters for calibration should be assigned for every iteration)
  # groups: specific group names for a parameter
  # fixed is a list of fixed parameters (not going to be calibrated)
  # diag.obs: cumulative diagnosis for 18 groups
  # ndiag.obs: observed new diagnosis for 18 groupos
  # death.obs: death observed for 18 groups
  # obs.msm is the number of msm observed, before adjusting missing/other risk grouops
  # inc.all: estimated incidence for overall
  # inc.msm: estimated incidence for MSM

require(deSolve)

# obj <- function(calib.par, calpar.info, fixed, fixed.list, calib.target, valid.target, plot=F) { 
obj <- function(calib.par, calpar.info, fixed, fixed.list, calib.target, valid.target) { 

  # Attain inputs for objective function
  names.gp   = as.vector(fixed$names.gp)
  names18    = as.vector(fixed$names18)
  n.gp       = length(fixed$names.gp)
  state.name = fixed$state.name
  leng   = calpar.info$plength
  names  = calpar.info$names
  
  diag.obs  = calib.target$diag18.obs [-1, ]
  ndiag.obs = calib.target$ndiag18.obs[-1, ]
  death.obs = calib.target$death18.obs[-1, ]
  inc.all   = valid.target$obs.inc.all
  inc.msm   = valid.target$obs.inc.msm
  
  ## Assigning parameter values for Morris method ##
  vlist = fixed.list
  vparameters = fixed
  
  par2   = as.list(leng)
  
  for(l in 1:length(leng)){
    if (l==1) {
      par2[[l]]   = calib.par[1:leng[l]]
    }
    else {
      par2[[l]]   = calib.par[(sum(leng[1:(l-1)])+1):sum(leng[1:l])]
    }
  }
  
  names(par2)   = names
  
  # update parameter list with new free parameter values and save in par3
  
  for (l in 1:length(names)){
    if (length(unlist(vparameters[names[l]])) == 1) { # single parameter
      vparameters[names[l]][[1]] = par2[names[l]][[1]]
    }
    if (par_info[par_info$parameter == names[l], ]$stratification == 0) {
      vlist[[names[l]]]$pe = par2[names[l]][[1]]
    }
    if (par_info[par_info$parameter == names[l], ]$stratification == 3) { #parameters with 3 dimensions: gender, ethnicity, risk
      for (j in 1:leng[l]) {
        vlist[[names[l]]][vlist[[names[l]]]$gender == calpar[calpar$par == names[l], ][j, ]$gender & vlist[[names[l]]]$ethnicity == calpar[calpar$par == names[l], ][j,]$ethnicity & vlist[[names[l]]]$risk == calpar[calpar$par == names[l], ][j,]$risk, ]$pe = par2[names[l]][[1]][j]
      }}
    if (par_info[par_info$parameter == names[l], ]$stratification == 4) { #parameters with 4 dimensions: gender, ethnicity, risk, sexual intensity
      for (j in 1:leng[l]) {
        vlist[[names[l]]][vlist[[names[l]]]$gender == calpar[calpar$par == names[l], ][j, ]$gender & vlist[[names[l]]]$ethnicity == calpar[calpar$par == names[l], ][j,]$ethnicity & vlist[[names[l]]]$risk == calpar[calpar$par == names[l], ][j,]$risk & vlist[[names[l]]]$sexual.intensity == calpar[calpar$par == names[l], ][j,]$sexual.intensity, ]$pe = par2[names[l]][[1]][j]
      }}
    if (par_info[par_info$parameter == names[l], ]$stratification == 1) { #parameters with 1 dimension: CD4
      for (j in 1:leng[l]) {
        vlist[[names[l]]][vlist[[names[l]]]$cd4 == calpar[calpar$par == names[l], ][j, ]$cd4, ]$pe = par2[names[l]][[1]][j]
      }}
  }
  
  source("parameterization.R", local =T)
  # if ("v.ssp" %in% names){
  #   vparameters$v.ssp = vparameters$v.ssp * (vparameters$init.tot[gp18.gn$all.pwid]/ sum(vparameters$init.tot[gp18.gn$all.pwid])) 
  # }
  
  
  par3 = vparameters

  #parameters in "par3" overwrite ones in "fixed" if names are the same.
  out_euler <- euler(x, vt, ode_model, par3)[ ,-1]
  outa = array(out_euler[-1, ], dim = c(n, length(state.name), n.gp)) # initial value deleted


  dimnames(outa)[[2]] = state.name
  #number of individuals in each states at time t, excluding incidence
  outn = outa[ ,1:19, ] 
  
  #total diagnoses
  diag         = apply(outn[end_yr_ind, 10:19, ], c(1,3), sum)
  diag.all     = rowSums(diag)
  diag.b       = rowSums(diag[ , black])
  #print(diff(apply(outn[,10:19,msm],1,sum)))
  diag.m       = rowSums(diag[ , msm])
  diag18.model = matrix(0, 4, 18)
  for (i in 1:18){
    ind = which(rname %in% names18[i])
    diag18.model[ ,i] = rowSums(diag[ ,ind])
  }
  # total diagnoses in 11 groups (combining MSM/PWID into one; PWID into one)
  diag11.model = group11(diag18.model)
  diag11.obs   = group11(diag.obs)
  names11 = c(names18[1:3], "MSM/PWID", "PWID", names18[13:18])

  # new diagnosis from model
  ndiag        = rbind(outa[12, "diag", ],  apply(outa[end_yr_ind, "diag", ], 2, diff))
  ndiag.all    = rowSums(ndiag) # total
  ndiag.b      = rowSums(ndiag[ ,black]) # black
  ndiag.m      = rowSums(ndiag[ ,msm]) #msm
  # observed new diagnosis
  ndiag.all.obs= rowSums(ndiag.obs)
  ndiag.b.obs  = rowSums(ndiag.obs[ ,grep("black", names18)]) # black
  ndiag.m.obs  = rowSums(ndiag.obs[ ,1:3]) #msm
  
  # death from model
  death        = rbind(outa[12,"death", ],  apply(outa[end_yr_ind, "death", ], 2, diff))
  death.all    = rowSums(death) # total
  death.b      = rowSums(death[ ,black]) # black
  death.m      = rowSums(death[ ,msm]) #msm
  # observed death
  
  if (CITY != "MIA"){
      death.all.obs= rowSums(death.obs)
      death.b.obs  = rowSums(death.obs[ ,grep("black", names18)]) # black
      death.m.obs  = rowSums(death.obs[ ,1:3]) #msm
  } else {
      death.all.obs= death.obs[,1]
      death.b.obs  = death.obs[,3] # black
  }

  
  # incidence from model
  cum.inc    = apply(outa[ ,c("inc_bo","inc_bs","inc_g"), ], c(1,3), sum)
  inc        = rbind(cum.inc[12, ],  apply(cum.inc[end_yr_ind, ], 2, diff))
  tot.inc    = rowSums(inc)
  msm.inc    = rowSums(inc[ ,msm])
  midu.inc   = rowSums(inc[ ,midu])
  all.msm.inc= rowSums(inc[ ,all.msm])
  # observed incidence
  inc.all.obs= inc.all$value
  inc.msm.obs= inc.msm$value
  
  # #Create the plots for comapring estimated outcomes with the observed values
  # if (plot==T){
  #   par(mfrow = c(5, 4))
  #   par(oma =  c(3, 3, 3, 3))
  #   par(mar =  c(2, 2, 2, 1))
  #   # cumulative diagnosis
  #   for (i in 1:11){
  #     moplot(diag11.model[ ,i], diag11.obs[ ,i], paste("Diag:", names11[i]))
  #   }
  # 
  #   # new diagnosis
  #   moplot(ndiag.all, ndiag.all.obs, "New diag: total")
  #   moplot(ndiag.b,   ndiag.b.obs,   "New diag: black")
  #   moplot(ndiag.m,   ndiag.m.obs,   "New diag: MSM")
  # 
  #   # death
  #   moplot(death.all, death.all.obs, "Death: total")
  #   moplot(death.b,   death.b.obs,   "Death: black")
  #   moplot(death.m,   death.m.obs,   "Death: MSM")
  # 
  #   # incidence
  #   moplot(tot.inc, inc.all$value,   "Incidence: total",
  #          inc.all$low, inc.all$high)
  #   moplot(all.msm.inc, inc.msm$value, "Incidence: MSM & MSM/PWID",
  #          inc.msm$low, inc.msm$high)
  # 
  #   #### legend ####
  #   par(fig=c(0, 1, 0, 0.3), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
  #   plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  #   legend("bottom", c("Model","Obs"), lty=c(2,1), cex=1, lwd=2, col=c("blue","black"))
  # }
  
  if (CITY != "MIA"){
    model = c(diag11.model, ndiag.all,     ndiag.b,     ndiag.m,     death.all,     death.b,     death.m)
    obs   = c(diag11.obs,   ndiag.all.obs, ndiag.b.obs, ndiag.m.obs, death.all.obs, death.b.obs, death.m.obs)
  } else {
    model = c(diag11.model, ndiag.all,     ndiag.b,     ndiag.m,     death.all,     death.b)
    obs   = c(diag11.obs,   ndiag.all.obs, ndiag.b.obs, ndiag.m.obs, death.all.obs, death.b.obs)
  }
  
  dev = abs((model-obs)/obs)
  return(weighted.mean(dev, rep(fixed$w, each=4), na.rm=T))
}


# Function to combine MSM/PWID into one group; PWID into one group:
group11 = function(group18){
  g11 = matrix(0,4,11)
  g11[ ,4] = rowSums(group18[ ,4:6])
  g11[ ,5] = rowSums(group18[ ,7:12])
  g11[ ,c(1:3, 6:11)] = group18[ ,c(1:3, 13:18)]
  return(g11)
}
# 
# # Function for comparison plots
# moplot = function(model, obs, title, low=0, high=0){
#   ymin = min(c(model, obs), na.rm=T)*0.8
#   ymax = max(c(model, obs), na.rm=T)*1.2
#   if (any(low>0)){
#     ymin = min(c(model, low),  na.rm=T)*0.8
#     ymax = max(c(model, high), na.rm=T)*1.2
#   }
#   plot(yr, model, xlab="Year", main=title,
#        ylab ='Number of individuals', type='l', lty=2, ylim=c(ymin, ymax), lwd=2, col="blue",
#        xaxt = "n") #xaxis label not shown
#   if (nyr<=4) axis(side = 1, at = yr)
#   else axis(side=1)
#   lines(yr[1:4], obs, lty=1, lwd=2)
#   if (any(low>0)){
#     lines(yr[1:4], low,  lty=3, lwd=2)
#     lines(yr[1:4], high, lty=3, lwd=2)
#   }
# }
