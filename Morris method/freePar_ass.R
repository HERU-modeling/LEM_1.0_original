## Assigning parameter values for Morris method ##
par2   = as.list(leng)

if (pj ==0){
  for(l in 1:length(leng)){
    if (l==1) {
      par2[[l]]   = xarray[1, ,i][1:leng[l]]
    }
    else {
      par2[[l]]   = xarray[1, ,i][(sum(leng[1:(l-1)])+1):sum(leng[1:l])]
    }
  }
} else {
  for(l in 1:length(leng)){
    if (l==1) {
      par2[[l]]   = xarray[pj+1, ,i][1:leng[l]]
    }
    else {
      par2[[l]]   = xarray[pj+1, ,i][(sum(leng[1:(l-1)])+1):sum(leng[1:l])]
    }
  }
}
names(par2)   = names

# update parameter list with new free parameter values and save in par3

for (l in 1:length(names)){
  if (length(unlist(vparameters[names[l]])) == 1) { # single parameter
    vparameters[names[l]][[1]] = par2[names[l]][[1]]
  }
  if (par_info[par_info$parameter == names[l], ]$stratification == 0) { #parameters with 0 stratification
    vlist[[names[l]]]$pe = par2[names[l]][[1]]
  }
  if (par_info[par_info$parameter == names[l], ]$stratification == 3) { #parameters with 3 stratifications: gender, ethnicity, risk
    for (j in 1:leng[l]) {
      vlist[[names[l]]][vlist[[names[l]]]$gender == mmparn[mmparn$par == names[l], ][j, ]$gender & vlist[[names[l]]]$ethnicity == mmparn[mmparn$par == names[l], ][j,]$ethnicity & vlist[[names[l]]]$risk == mmparn[mmparn$par == names[l], ][j,]$risk, ]$pe = par2[names[l]][[1]][j]
    }}
  if (par_info[par_info$parameter == names[l], ]$stratification == 4) { #parameters with 4 stratifications: gender, ethnicity, risk, sexual intensity
    for (j in 1:leng[l]) {
      vlist[[names[l]]][vlist[[names[l]]]$gender == mmparn[mmparn$par == names[l], ][j, ]$gender & vlist[[names[l]]]$ethnicity == mmparn[mmparn$par == names[l], ][j,]$ethnicity & vlist[[names[l]]]$risk == mmparn[mmparn$par == names[l], ][j,]$risk & vlist[[names[l]]]$sexual.intensity == mmparn[mmparn$par == names[l], ][j,]$sexual.intensity, ]$pe = par2[names[l]][[1]][j]
    }}
  if (par_info[par_info$parameter == names[l], ]$stratification == 1) { #parameters with 1 stratification: CD4
    for (j in 1:leng[l]) {
      vlist[[names[l]]][vlist[[names[l]]]$cd4 == mmparn[mmparn$par == names[l], ][j, ]$cd4, ]$pe = par2[names[l]][[1]][j]
    }}
}

source("parameterization.R", local = T)
if ("v.ssp" %in% names){
  vparameters$v.ssp = vparameters$v.ssp * (vparameters$init.tot[gp18.gn$all.pwid]/ sum(vparameters$init.tot[gp18.gn$all.pwid])) 
}

par3 = vparameters
