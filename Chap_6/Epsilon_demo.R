
setwd("C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_6")

# COmpile
TMB::compile( "epsilon_estimator.cpp" )
dyn.load( TMB:::dynlib("epsilon_estimator") )

########## START IN-LINE SNIPPET
# Define true distribution and transformation
Dist = 0  # 0=Normal; 1=Gamma
Trans = 2  # 0=Identity; 1=sqrt; 2=exp
if( Dist==0 ){
  mu = 1
  sigma = 0.5
}else{
  mu = 2
  sigma = 1
}

# Build object
Data = list( "options_z"=c(Dist, Trans), "mu"=mu, "sigma"=sigma )
Params = list( "epsilon"=1.2, "delta"=vector() )
Obj = TMB::MakeADFun( data=Data, parameters=Params, random="epsilon" )
Obj$fn(Obj$par)  # optimize random effects

# Compute SEs and apply epsilon estimator
SD = TMB::sdreport( Obj, bias.correct=TRUE )

# Extract various estimators
BiasCorr = summary(SD,"report")['Z','Est. (bias.correct)']
Plugin = summary(SD,"report")['Z','Estimate']

# Manually calculate epsilon estimator
parhat = Obj$env$parList()
parhat[["delta"]] = 0
Obj = TMB::MakeADFun( data=Data, parameters=parhat, random="epsilon" )
(BiasCorr_manual = Obj$gr(Obj$par)[1])
(BiasCorr_sampled = Obj$simulate()$Zmean_sampled)
########## END IN-LINE SNIPPET

# Known value in some cases
if(Dist==0 & Trans==0) True = mu
if(Dist==0 & Trans==1) True = NA
if(Dist==0 & Trans==2) True = exp(mu + sigma^2/2)
if(Dist==1 & Trans==0) True = mu*sigma
if(Dist==1 & Trans==1) True = NA
if(Dist==1 & Trans==2) True = NA

## Make lognormal plot
if( Dist==0 & Trans==2 ){
  xplot = seq(0, exp(6), by=0.01)
  yplot = dlnorm( xplot, meanlog=mu, sdlog=sigma )
  zplot = dnorm( log(xplot), mean=mu, sd=sigma )

  png( file="Retransformation_bias.png", width=6, height=4, res=200, units="in" )
    par( mfrow=c(2,1), yaxs="i", mar=c(3,3,2,1), mgp=c(1.5,0.25,0), tck=-0.02 )
    plot(x=log(xplot), y=zplot, xlim=c(-1,3), type="l", lwd=2, ylab="Density", xlab="X", main = "Normal distribution" )
    abline( v=c(log(Plugin)), col=viridisLite::viridis(3)[1], lwd=2 )
    plot(x=xplot, y=yplot, xlim=exp(c(-0.5,2)), type="l", lwd=2, ylab="Density", xlab=expression( Y == e^X ), log="", main = "Transformed (lognormal) distribution" )
    abline( v=c(True,Plugin,BiasCorr), col=viridisLite::viridis(3), lwd=2 )
    legend( "topright", bty="n", fill=viridisLite::viridis(3), legend=c("True value","Plug-in estimator","Epsilon estimator") )
  dev.off()
}

###############
# Check Eq. 6.8
#
# Y = exp(X) and X ~ Normal ... WORKS (3rd term = 0)
# Y = sqrt(X) and X ~ Normal ... WORKS (3rd term = 0)
# Y = X and X ~ Gamma() ... WORKS (2nd term = 0)
#
###############
library(calculus)
eps = Obj$env$parList()$epsilon
Params$epsilon = eps
Obj = TMB::MakeADFun( data=Data, parameters=Params, random=NULL )
f2x = derivative(f = function(x) Obj$fn(x), var=c(x=eps), order=2)
f3x = derivative(f = function(x) Obj$fn(x), var=c(x=eps), order=3)
phi1x = derivative(f = function(x) Obj$report(x)$Z, var=c(x=eps), order=1)
phi2x = derivative(f = function(x) Obj$report(x)$Z, var=c(x=eps), order=2)
# Eq. 6.8
BiasCorr_analytic = Plugin + 0.5*phi2x/f2x - 0.5*phi1x*f3x/f2x^2

# Show all
c( 'True'=True, 'Sampled'=BiasCorr_sampled, 'Plugin'=Plugin, 'Auto'=BiasCorr, 'Manual'=BiasCorr_manual, 'Analytic'=BiasCorr_analytic )
