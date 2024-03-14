#### Multiple calibration results
source("calOut.R")

calOut10 = numeric(164)   # ((11 groups * 3 targets) + 6 calib targets + 2 valid targets) * 4 years
calOut10 = calOut(calib.par = plot.par, calpar.info = calpar.info, fixed = vparameters, fixed.list =vlist)


diag11.model = calOut10[1 : 44]
ndiag.all    = calOut10[133:136]
ndiag.b      = calOut10[137:140]
ndiag.m      = calOut10[141:144]
death.all    = calOut10[145:148]
death.b      = calOut10[149:152]
death.m      = calOut10[153:156]

## Target operationalization
diag.obs  = calib.target$diag18.obs [-1, ]
ndiag.obs = calib.target$ndiag18.obs[-1, ]
death.obs = calib.target$death18.obs[-1, ]
inc.all   = valid.target$obs.inc.all
inc.msm   = valid.target$obs.inc.msm

# observed total diagnosed PLHIV
diag11.obs   = group11(diag.obs)
names11 = c(names18[1:3], "MSM/PWID", "PWID", names18[13:18])

# observed new diagnosis
ndiag.all.obs= rowSums(ndiag.obs)
ndiag.b.obs  = rowSums(ndiag.obs[ ,grep("black", names18)]) # black
ndiag.m.obs  = rowSums(ndiag.obs[ ,1:3]) #msm
ndiag11.obs  = group11(ndiag.obs)

# observed mortality
death.all.obs= rowSums(death.obs)
death.b.obs  = rowSums(death.obs[ ,grep("black", names18)]) # black
death.m.obs  = rowSums(death.obs[ ,1:3]) #msm
death11.obs  = group11(death.obs)

obs   = c(diag11.obs,   ndiag.all.obs, ndiag.b.obs, ndiag.m.obs, death.all.obs, death.b.obs, death.m.obs)

model = c(diag11.model, ndiag.all, ndiag.b, ndiag.m, death.all, death.b, death.m)


#### plots for the multiple calibrated parameter sets

# function to plot model outcomes compared with observed values
moplot = function(model, obs, title, low =0, high =0){
  ymin = min(c(model, obs), na.rm =T)*0.8
  ymax = max(c(model, obs), na.rm =T)*1.1
  # observed output
  plot(yr, obs, xlab ="Year", main =title,
       ylab = 'Number of individuals', type ='l', ylim =c(ymin,ymax), lwd =2,
       xaxt = "n") #xaxis label not shown
  if (nyr <= 4) axis(side = 1, at = yr)
  else axis(side =1)
  # from model output
  lines(yr, model, lty =2, lwd =1.5, col="blue")  
  
  # CI for incidence
  if (any(low >0)){
    lines(yr[1:4], low,  lty=3, lwd=2)
    lines(yr[1:4], high, lty=3, lwd=2)
  }
}

####Result plots for all groups####
par(oma = c(3, 3, 3, 3))
par(mar =  c(2, 2, 2, 1))
par(mfrow = c(5, 7))

for (i in 1:11){
  moplot(calOut10[(4*i-3):(4*i)],     diag11.obs[ ,i],  paste("Diag:", names11[i]))
}
for (i in 1:11){
  moplot(calOut10[(4*i+41):(4*i+44)], ndiag11.obs[ ,i], paste("New diag:", names11[i]))
}
for (i in 1:11){
  moplot(calOut10[(4*i+85):(4*i+88)], death11.obs[ ,i], paste("Death:", names11[i]))
}
#incidence
moplot(calOut10[157:160], obs.inc.all$value, "Total HIV incidence",           obs.inc.all$low, obs.inc.all$high)
moplot(calOut10[161:164], obs.inc.msm$value, "HIV incidence: MSM & MSM/PWID", obs.inc.msm$low, obs.inc.msm$high)
#### legend ####
par(fig=c(0, 1, 0, 0.3), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
legend("bottom", c("Model","Obs"), lty=c(2,1), cex=1, lwd=2, col=c("blue","black"))
