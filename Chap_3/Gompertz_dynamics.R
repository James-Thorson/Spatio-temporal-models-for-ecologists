
setwd("C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_3")

# Illustrate Gompertz dynamics
K = 3000
beta = 0.4

# Convert parameters
rho = 1 - beta
alpha = log(K) * beta
# K = exp( alpha / beta )

x = seq( log(K/10), log(K*10), length=1000 )
y = exp( alpha - beta*x )

png( "Gompertz_production.png", width=3, height=3, res=200, unit="in" )
  par( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )
  matplot( x=exp(x), y=cbind(y,1), ylim=c(0.1,10), log="xy", lwd=2, type="l", col="black",
         lty=c("solid","dashed"), xlab=expression(n[t]), ylab=expression(lambda[t]) )
dev.off()

###############
# Download pollock EBS data
###############

#library(FishData)
#Data = download_catch_rates( localdir=getwd(), species_set="Hippoglossoides elassodon" )
#Index_t = tapply( Data$Wt/Data$AreaSwept_ha, INDEX=Data$Year, FUN=mean ) * 100 * 492989.9
#write.csv( Index_t, file="Biomass_index.csv" )

########## START IN-TEXT SNIPPET
# Get abundance in KG
Index_t = read.csv( file="Biomass_index.csv" )
Index_t = array(Index_t[,2], dimnames=list(Index_t[,1]))
Brange = range(Index_t) * c(0.8,1.2)
Prange = c(0.5,2)

# Fit Gompertz model with process errors
y = log(Index_t[-1])
x = log(Index_t[-length(Index_t)])
Lm = lm( y ~ 1 + x )
########## END IN-TEXT SNIPPET

# Plot
xpred = seq(log(Brange[1]),log(Brange[2]),length=1000)
ypred = predict.lm( Lm, newdata=data.frame("x"=xpred), se.fit=TRUE )
png( "gompertz_data.png", width=7, height=4, unit='in', res=200 )
  par( mfrow=c(1,2), mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02, xaxs="i", yaxs="i" )
  plot( x=names(Index_t), y=Index_t, log="y", type="p", main="Bimass timeseries", ylim=Brange, xlab=expression(t), ylab=expression(b[t]) )
  plot( y=exp(y-x), x=exp(x), log="xy", main="Production curve", ylim=Prange, xlim=Brange, xlab=expression(b[t]), ylab=expression(lambda[t]) )
  lines( x=exp(xpred), y=exp(ypred$fit-xpred), main="Production function" )
  polygon( x=exp(c(xpred,rev(xpred))), y=exp(c(ypred$fit-xpred-2*ypred$se,rev(ypred$fit-xpred+2*ypred$se))), col=rgb(0,0,1,0.2) ) #
  abline( h=1, lty="dotted" )
dev.off()

###############
# Correlation and semi-variance
###############

#
sigmaP = 0.4
sigmaM = 0.4
rho1 = 0.5
rho2 = 0.75

#
calc_semivariance = function(x,lag=0){
  diff = x[-seq_len(lag)] - rev(rev(x)[-seq_len(lag)])
  return(0.5*mean(diff^2))
}
sim = function( sigmaP, sigmaM, rho, n_t=1e6, K=3000 ){
  b_t = rep(0,n_t)
  b_t[1] = 10
  beta = 1 - rho
  alpha = log(K) * beta
  epsP_t = rnorm(n=n_t, mean=0, sd=sigmaP)
  for(t in 1:(n_t-1)) b_t[t+1] = b_t[t] * exp( alpha - beta * log(b_t[t]) + epsP_t[t] )

  #
  epsM_t = rnorm(n=n_t, mean=0, sd=sigmaM)
  bobs_t = b_t * exp(epsM_t)
  gamma_l = sapply(1:10, calc_semivariance, x=log(bobs_t) )

  #
  sill = sigmaP^2 / (1-rho^2) + sigmaM^2
  return( list("semivariance"=gamma_l, "correlation"=1-gamma_l/sill, "sill"=sill) )
}

png( file="Gompertz_semivariance.png", width=6, height=6, unit="in", res=200 )
  par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(2,0.5,0), oma=c(2,2,2,0), xaxs="i", yaxs="i")
  Y = cbind( sim(sigmaP=sigmaP,sigmaM=0,rho=rho1)$semivar, sim(sigmaP=sigmaP,sigmaM=0,rho=rho2)$semivar )
  matplot( x=1:10, y=Y, type="p", ylim=c(0,0.6), xlim=c(0,10.5), xlab="", ylab="", col=c('black',"red"), pch=20, cex=1.5 )
  mtext(side=3, text="Without measurement error")
  mtext(side=2, text="Semivariance", line=2)
  Y = cbind( sim(sigmaP=sigmaP,sigmaM=sigmaM,rho=rho1)$semivar, sim(sigmaP=sigmaP,sigmaM=sigmaM,rho=rho2)$semivar )
  matplot( x=1:10, y=Y, type="p", ylim=c(0,0.6), xlim=c(0,10.5), xlab="", ylab="", col=c('black',"red"), pch=20, cex=1.5 )
  mtext(side=3, text="With measurement error")
  Y = cbind( sim(sigmaP=sigmaP,sigmaM=0,rho=rho1)$corr, sim(sigmaP=sigmaP,sigmaM=0,rho=rho2)$corr )
  matplot( x=1:10, y=Y, type="p", ylim=c(0,1), xlim=c(0,10.5), xlab="", ylab="", col=c('black',"red"), pch=20, cex=1.5 )
  mtext(side=2, text="Correlation", line=2)
  Y = cbind( sim(sigmaP=sigmaP,sigmaM=sigmaM,rho=rho1)$corr, sim(sigmaP=sigmaP,sigmaM=sigmaM,rho=rho2)$corr )
  matplot( x=1:10, y=Y, type="p", ylim=c(0,1), xlim=c(0,10.5), xlab="", ylab="", col=c('black',"red"), pch=20, cex=1.5 )
  mtext( side=1, outer=TRUE, text=expression(Delta_t) )
dev.off()

# Calculate semi-variance for LM
gamma_l = sapply(1:10, calc_semivariance, x=Lm$residuals )


#############
# Visualize graph for CAR
#############

# Visualize graph
library(diagram)
n_groups = 4
names <- c("log(b0)", "alpha", "rho", paste0("x",1:n_groups), paste0("y",1:n_groups) )
M <- array(0, dim=c(length(names),length(names)), dimnames=list(names,names))
  M['x1','log(b0)'] = ""
  M[match(paste0("x",2:n_groups),names),'alpha'] = ""
  M[match(paste0("x",2:n_groups),names),'rho'] = ""
  M[cbind(grep("y",names),grep("x",names))] = ""
  M[cbind(match(paste0("x",2:n_groups),names),match(paste0("x",2:n_groups-1),names))] = ""
box.type = c( rep("circle",3), rep("diamond",n_groups), rep("square",n_groups) )
png( "graph.png", width=7, height=3, res=200, un="in")
  par( mar=c(0,0,0,0) )
  plotmat( M, pos = c(3,n_groups,n_groups), curve = 0, name = names, lwd = 3, arr.pos = 0.6,
           box.lwd = 2, box.type = box.type, box.prop = 0.5, cex.txt=0.8 )
dev.off()

####################
# Fit CAR in TMB
####################

########## START IN-TEXT SNIPPET
# Compile
library(TMB)
compile( "gompertz.cpp" )
dyn.load( dynlib("gompertz") )

# Build inputs
Data = list( "log_b_t"=log(Index_t), "log_bnew_z"=xpred, "simulate_t"=rep(0,length(Index_t)) )
Parameters = list( "log_d0"=0, "log_sigmaP"=1, "log_sigmaM"=1, "alpha"=0, "rho"=0, "log_d_t"=rep(0,length(Index_t)) )
Random = "log_d_t"

# Build and fit object
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random )
opt_CAR = nlminb( start=Obj$par, obj=Obj$fn, grad=Obj$gr )
opt_CAR$SD = sdreport( Obj )
########## END IN-TEXT SNIPPET

#
png( "gompertz_sparsity.png", width=4, height=4, unit='in', res=200 )
  par( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )
  Matrix::image( Obj$env$spHess(par=Obj$env$last.par.best,random=TRUE) )
dev.off()

# Plot fit
png( "gompertz_fit.png", width=7, height=4, unit='in', res=200 )
  par( mfrow=c(1,2), mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02, xaxs="i", yaxs="i" )
  # Time-series
  out = data.frame( names(Index_t), as.list(opt_CAR$SD,"Est")$log_d_t, as.list(opt_CAR$SD,"Std")$log_d_t )
  plot( x=out[,1], y=Index_t, log="y", type="p", main="Bimass timeseries", ylim=Brange, xlab=expression(t), ylab=expression(b[t]) )
  lines( x=out[,1], y=exp(as.list(opt_CAR$SD,"Est")$log_d_t), lwd=2 )
  polygon( x=c(out[,1],rev(out[,1])),
           y=exp(c(out[,2]-2*out[,3],rev(out[,2]+2*out[,3]))), col=rgb(1,0,0,alpha=0.2) )
  # Productio
  ypred = list("fit"=as.list(opt_CAR$SD,"Est",report=TRUE)$log_out_z, "se.fit"=as.list(opt_CAR$SD,"Std",report=TRUE)$log_out_z )
  plot( y=exp(out[-1,2]-out[-nrow(out),2]), x=exp(out[-nrow(out),2]), log="xy", main="Production curve", xlim=Brange, ylim=Prange, xlab=expression(b[t]), ylab=expression(lambda[t]) ) #, xlim=exp(range(xpred)), ylim=range(Ybounds) )
  lines( x=exp(xpred), y=exp(ypred$fit-xpred), main="Production function" )
  polygon( x=exp(c(xpred,rev(xpred))), y=exp(c(ypred$fit-xpred-2*ypred$se, rev(ypred$fit-xpred+2*ypred$se))), col=rgb(0,0,1,0.2) ) #
  abline( h=1, lty="dotted" )
dev.off()

####################
# Fit SAR in TMB
####################

# Check math
R = rho ^ as.matrix(dist(1:20))
Q = solve(R)

# Build inputs
compile( "gompertz_SAR.cpp" )
dyn.load( dynlib("gompertz_SAR") )
Data = list( "log_b_t"=log(Index_t), "log_bnew_z"=xpred )
Parameters = list( "log_delta"=0, "log_sigmaP"=1, "log_sigmaM"=1, "alpha"=0, "rho"=0.76, "eps_t"=rep(0,length(Index_t)) )
Random = c("eps_t")
#Map = list("rho"=factor(NA))

# Build and fit object
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random, DLL="gompertz_SAR" )  #
opt_SAR = nlminb( start=Obj$par, obj=Obj$fn, grad=Obj$gr )
opt_SAR$SD = sdreport(Obj)

# Compare two versions
out = data.frame( "CAR params"=names(opt_CAR$par), signif(summary(opt_CAR$SD,"fixed"),3),
  "SAR params"=names(opt_SAR$par), signif(summary(opt_SAR$SD,"fixed"),3) )
write.csv( out, file="TMB_comparison.csv", row.names=FALSE)

# Plot fit
png( "gompertz_fit_SAR.png", width=7, height=4, unit='in', res=200 )
  par( mfrow=c(1,2), mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02, xaxs="i", yaxs="i" )
  # Time-series
  out = data.frame( names(Index_t), as.list(opt_SAR$SD,"Est",report=TRUE)$log_d_t, as.list(opt_SAR$SD,"Std",report=TRUE)$log_d_t )
  plot( x=out[,1], y=Index_t, log="y", type="p", main="Bimass timeseries", ylim=Brange, xlab=expression(t), ylab=expression(b[t]) )
  lines( x=out[,1], y=exp(as.list(opt_SAR$SD,"Est",report=TRUE)$log_d_t), lwd=2 )
  polygon( x=c(out[,1],rev(out[,1])),
           y=exp(c(out[,2]-2*out[,3],rev(out[,2]+2*out[,3]))), col=rgb(1,0,0,alpha=0.2) )
  # Productio
  ypred = list("fit"=as.list(opt_SAR$SD,"Est",report=TRUE)$log_out_z, "se.fit"=as.list(opt_SAR$SD,"Std",report=TRUE)$log_out_z )
  plot( y=exp(out[-1,2]-out[-nrow(out),2]), x=exp(out[-nrow(out),2]), log="xy", main="Production curve", xlim=Brange, ylim=Prange, xlab=expression(b[t]), ylab=expression(lambda[t]) ) #, xlim=exp(range(xpred)), ylim=range(Ybounds) )
  lines( x=exp(xpred), y=exp(ypred$fit-xpred), main="Production function" )
  polygon( x=exp(c(xpred,rev(xpred))), y=exp(c(ypred$fit-xpred-2*ypred$se, rev(ypred$fit-xpred+2*ypred$se))), col=rgb(0,0,1,0.2) ) #
  abline( h=1, lty="dotted" )
dev.off()

####################
# Forecast CAR in TMB
####################

# Build inputs
Data = list( "log_b_t"=c(log(Index_t),rep(NA,10)), "log_bnew_z"=xpred, "simulate_t"=c(rep(0,length(Index_t)),rep(1,10)) )
Parameters = list( "log_d0"=0, "log_sigmaP"=1, "log_sigmaM"=1, "alpha"=0, "rho"=0, "log_d_t"=rep(0,length(Data$log_b_t)) )
Random = c("log_d_t")

# Build and fit object
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random, DLL="gompertz" )  #
opt_CAR = nlminb( start=Obj$par, obj=Obj$fn, grad=Obj$gr )
opt_CAR$SD = sdreport(Obj, par.fixed=opt_CAR$par)

years = c( as.numeric(names(Index_t)), max(as.numeric(names(Index_t))) + 1:10 )
######### START IN-TEXT SNIPPET
# Full uncertainty
out_full = data.frame( "year" = years,
                       "pred" = as.list(opt_CAR$SD,"Est")$log_d_t,
                       "SE" = as.list(opt_CAR$SD,"Std")$log_d_t )

# Samples from forecast variability
log_d_tz = sapply( 1:1000,
                   FUN=\(i) Obj$env$simulate(Obj$env$last.par.best)$log_d_t )
out_1 = data.frame( "year" = years,
                    "pred" = apply(log_d_tz,MARGIN=1,FUN=mean),
                    "SE" = apply(log_d_tz,MARGIN=1,FUN=sd) )

# Samples from forecast uncertainty + state uncertainty
u_zr = Obj$env$last.par.best %o% rep(1, 1000)
MC = Obj$env$MC( keep=TRUE, n=1000, antithetic=FALSE )
u_zr[Obj$env$random,] = attr(MC, "samples")
log_d_tz = sapply( 1:1000,
                   FUN=\(i) Obj$env$simulate(par=u_zr[,i])$log_d_t )
out_2 = data.frame( "year" = years,
                    "pred" = apply(log_d_tz,MARGIN=1,FUN=mean),
                    "SE" = apply(log_d_tz,MARGIN=1,FUN=sd) )
########### END IN-TEXT SNIPPET

# Plot fit
png( "gompertz_forecast.png", width=7, height=3, unit='in', res=200 )
  par( mfrow=c(1,3), mar=c(3,2,4,1), mgp=c(2,0.5,0), tck=-0.02, xaxs="i", yaxs="i" )
  for(z in 1:3 ){
    if(z==1) out = out_1
    if(z==2) out = out_2
    if(z==3) out = out_full
    main = c( "Future variability", "Future variability + \nstate uncertainty", "Future variability + \nstate uncertainty + \nparameter uncertainty" )[z]
    plot( x=out[,1], y=exp(Data$log_b_t), log="y", type="p", main=main, ylim=Brange, xlab=expression(t), ylab="" )
    lines( x=out[,1], y=exp(as.list(opt_CAR$SD,"Est")$log_d_t), lwd=2 )
    polygon( x=c(out[,1],rev(out[,1])),
             y=exp(c(out[,2]-2*out[,3],rev(out[,2]+2*out[,3]))), col=rgb(1,0,0,alpha=0.2) )
  }
dev.off()


