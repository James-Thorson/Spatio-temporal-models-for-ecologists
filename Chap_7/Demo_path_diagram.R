

setwd( R'(C:\Users\James.Thorson\Desktop\Git\Spatio-temporal-models-for-ecologists\Chap_7)' )
source( "../Shared_functions/build_ram.R" )

library(mvtnorm)
library(stars)
library(diagram)

###############
# Basic description
###############

# Visualize graph
names <- c("Tree cover", "Temperature", "log(Density)")
box.type = rep("square",length(names))

# Construct M
M2 = M1 = array(0, dim=c(length(names),length(names)), dimnames=list(names,names))
M1['Temperature','Tree cover'] = -1
M1['log(Density)','Temperature'] = 1

#
set.seed(101)
n_obs = 100
C = rnorm(n_obs, mean=0, sd=1)
T = M1['Temperature','Tree cover'] * C + rnorm(n_obs, mean=0, sd=0.1)
logD = 0 + M1['log(Density)','Temperature'] * T
N = rpois( n_obs, exp(logD) )

# Fit as GLM
fit = glm( N ~ C + T, family="poisson" )

#
M2['log(Density)','Tree cover'] = round(fit$coef['C'],2)
M2['log(Density)','Temperature'] = round(fit$coef['T'],2)

# Plot
png( "graph.png", width=6, height=3, res=200, un="in")
  par(  mfrow=c(1,2), mar=c(0,0,2,0))
  plotmat( M1, pos=cbind(1-0.5,(rev(seq_along(names))-0.5)/length(names)), curve = 0, name = names, lwd = 3, arr.pos=0.7,
           box.lwd = 2, box.type = box.type, box.size=0.2, box.prop = 0.3, cex.txt=0.8, main="Simulated dynamics" )
  plotmat( M2, pos=cbind(c(0.25,0.75,0.5),c(0.75,0.75,0.25)), curve = 0, name = names, lwd = 3, arr.pos=0.5,
           box.lwd = 2, box.type = box.type, box.size=0.2, box.prop = 0.3, cex.txt=0.8, main="Regression model" )
dev.off()

#############
# Simulation experiment
#############

#
Results = array( NA, dim=c(100,6), dimnames=list(NULL,c("Intercept","C","T","Collinearity","Extrapolation","Counterfactual")) )
f = function(x,y) 1 - mean((x-y)^2) / mean(x^2)
for( r in 1:100 ){
  # Simulate data
  C = rnorm(n_obs, mean=0, sd=1)
  T = M1['Temperature','Tree cover'] * C + rnorm(n_obs, mean=0, sd=0.1)
  logD = 0 + M1['log(Density)','Temperature'] * T
  N = rpois( n_obs, exp(logD) )

  # In-sample fiit
  fit = glm( N ~ C + T, data=data.frame(C,T,N), family="poisson" )
  Results[r,c("Intercept","C","T")] = fit$coef
  Results[r,"Collinearity"] = cor(C,T)

  # Extrapolation
  Newdata = data.frame("C"=C, "T"=T)
  logmu_true = 0 + M1['log(Density)','Temperature'] * Newdata$T
  logmu_pred = predict( fit, newdata=Newdata, type="link" )
  Results[r,"Extrapolation"] = f(logmu_true, logmu_pred)

  # Counter-factual change
  Newdata = data.frame("C"=C, "T"=T+1)
  logmu_true = 0 + M1['log(Density)','Temperature'] * Newdata$T
  logmu_pred = predict( fit, newdata=Newdata, type="link" )
  Results[r,"Counterfactual"] = f(logmu_true, logmu_pred)
}
Results[,"Counterfactual"] = ifelse( Results[,"Counterfactual"]<=0, 1e-10, Results[,"Counterfactual"] )

# Results as histograms
png( file="GLM_performance.png", width=4, height=4, res=200, units="in" )
  par( mfcol=c(2,2), yaxs='i', mar=c(2,2,2,0) )
  hist( Results[,'C'], xlab="Tree cover", main="Tree cover effect", ylab="" )
  abline( v=M1['log(Density)','Tree cover'], lwd=3, col="blue" )
  hist( Results[,'T'], xlab="Temperature", main="Temperature effect", ylab="" )
  abline( v=M1['log(Density)','Temperature'], lwd=3, col="blue" )
  hist( Results[,'Extrapolation'], main="Interpolation", xlab="Correlation", ylab="", breaks=seq(0,1,by=0.05) )
  hist( Results[,'Counterfactual'], main="Counterfactual", xlab="Correlation", ylab="", breaks=seq(0,1,by=0.05) )
dev.off()

#################
# SEM
#################

library(sem)
library(TMB)

# compile
file.copy( from="../Chap_4/make_covariance.hpp", to="../Chap_7/make_covariance.hpp", overwrite=TRUE )
compile("SEMGLM.cpp")
dyn.load( dynlib("SEMGLM") )

# Re-simulate original data
set.seed(101)
n_obs = 100
C = rnorm(n_obs, mean=0, sd=1)
T = M1['Temperature','Tree cover'] * C + rnorm(n_obs, mean=0, sd=0.1)
logD = 0 + M1['log(Density)','Temperature'] * T
N = rpois( n_obs, exp(logD) )

########## START IN-TEXT SNIPPET
# Define model and convert to RAM
text = "
  C -> T, b1
  T -> logD, b2
  logD <-> logD, NA, 0.01
"
SEM_model = sem::specifyModel( text=text, exog.variances=TRUE,
                               endog.variances=TRUE, covs=c("C","T","logD") )
RAM = build_ram( SEM_model, c("C","T","logD") )
########## END IN-TEXT SNIPPET
capture.output( RAM, file="RAM.txt" )

#
Data = list( y_iz = cbind(C,T,N),
             RAM = as.matrix(RAM[,1:4]),
             familycode_z = c(0,0,1),
             RAMstart = ifelse( is.na(RAM[,5]), 0, as.numeric(RAM[,5])) )
Params = list( x_iz = Data$y_iz,
               beta_j = rep(1,max(Data$RAM[,4])) )
Map = list( "x_iz" = factor(cbind(NA,NA,1:nrow(Data$y_iz))) )

# Make TMB object
Obj = MakeADFun( data=Data, parameters=Params, map=Map, random="x_iz" )
Opt = nlminb( objective=Obj$fn, grad=Obj$gr, start=Obj$par )
Opt$SD = sdreport( Obj )

# Coerce to SEM object
Sprime = Obj$report()$V_zz
  rownames(Sprime) = colnames(Sprime) = c("C","T","logD")
mysem = sem::sem( SEM_model,
           S = Sprime,
           N = nrow(Data$y_iz) )

# Use package-SEM to plot path diagram
sem::pathDiagram( model = mysem,
                  style = "traditional",
                  edge.labels = "values",
                  file = tempdir(),
                  graphics.fmt = "png",
                  output.type = "graphics" )
file.rename( from=paste0(tempdir(),".png"), to="path.png" )

# Run simulation experiment
Results = array( NA, dim=c(100,6), dimnames=list(NULL,c("Intercept","C","T","Collinearity","Extrapolation","Counterfactual")) )
f = function(x,y) 1 - mean((x-y)^2) / mean(x^2)
for( r in 1:100 ){
  # Simulate data
  C = rnorm(n_obs, mean=0, sd=1)
  T = M1['Temperature','Tree cover'] * C + rnorm(n_obs, mean=0, sd=0.1)
  logD = 0 + M1['log(Density)','Temperature'] * T
  N = rpois( n_obs, exp(logD))

  # SEM fit
  Data$y_iz = cbind(C,T,N)
  Params$x_iz = Data$y_iz
  Obj = MakeADFun( data=Data, parameters=Params, map=Map, random="x_iz" )
  Obj$env$beSilent()
  Opt = nlminb( objective=Obj$fn, grad=Obj$gr, start=Obj$par )
  Opt$SD = sdreport( Obj )
  Results[r,c("Intercept","C","T")] = c( NA, Opt$par[1:2] )
  Results[r,"Collinearity"] = cor(C,T)

  # Extrapolation
  Newdata = data.frame("C"=C, "T"=T)
  logmu_true = 0 + M1['log(Density)','Temperature'] * Newdata$T
  logmu_pred = 0 + Opt$par[2] * Newdata$T
  Results[r,"Extrapolation"] = f(logmu_true, logmu_pred)

  # Counter-factual change
  Newdata = data.frame("C"=C, "T"=T+1)
  logmu_true = 0 + M1['log(Density)','Temperature'] * Newdata$T
  logmu_pred = 0 + Opt$par[2] * Newdata$T
  Results[r,"Counterfactual"] = f(logmu_true, logmu_pred)
}

# Results as histograms
png( file="SEM_performance.png", width=4, height=4, res=200, units="in" )
  par( mfcol=c(2,2), yaxs='i', mar=c(2,2,2,0) )
  hist( Results[,'C'], xlab="Tree cover", main="Tree cover effect", ylab="" )
  abline( v=M1['Temperature','Tree cover'], lwd=3, col="blue" )
  hist( Results[,'T'], xlab="Temperature", main="Temperature effect", ylab="" )
  abline( v=M1['log(Density)','Temperature'], lwd=3, col="blue" )
  hist( Results[,'Extrapolation'], main="Extrapolation", xlab="Correlation", ylab="", breaks=seq(0,1,by=0.05) )
  hist( Results[,'Counterfactual'], main="Counterfactual", xlab="Correlation", ylab="", breaks=seq(0,1,by=0.05) )
dev.off()

