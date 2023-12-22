
setwd("C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_4")
source( "../Shared_functions/build_ram.R" )

##############
# Simulate central-place foraging trips
##############
library(TMB)
library(numDeriv)
library(mvtnorm)
library(pracma)

########## START IN-TEXT SNIPPET
# Parameters
n_t = 1000
beta_t = seq(0.5, -0.5, length=n_t)
gamma_t = rep(0, n_t)
get_preference = function( beta, gamma, s ){
  beta*sqrt(s[1]^2+s[2]^2) + gamma*(s[2])
}

# Simulation function
simulate_track = function( n_t, beta_t, gamma_t, get_preference ){
  s_t = array( NA, dim=c(n_t,2), dimnames=list(NULL,c("x","y")) )
  # Euler approximation
  for( t in seq_len(n_t) ){
    if(t==1) s_t[t,] = c(0,0)
    if(t>=2){
      gradient = grad( f=get_preference, x0=s_t[t-1,,drop=FALSE],
                       beta=beta_t[t-1], gamma=gamma_t[t-1])
      s_t[t,] = s_t[t-1,] + gradient + rmvnorm(n=1,sigma=diag(2))
    }
  }
  return(s_t)
}
########## END IN-TEXT SNIPPET

plot_vectors = function( beta, gamma, xlim=c(-200,200), ylim=c(-200,200), ...){
  loc_gz = expand.grid( "x"=seq(xlim[1],xlim[2],length=100),"y"=seq(ylim[1],ylim[2],length=100))
  h_g = apply(as.matrix(loc_gz), MARGIN=1, FUN=get_preference, beta=beta, gamma=gamma)
  out = sf::st_as_sf(cbind("h_g"=h_g-mean(h_g),loc_gz), coords = c("x","y"))
  out = stars::st_rasterize(out)
  image(out["h_g"], axes=TRUE, col=viridisLite::magma(10), ..., breaks=seq(-100,100,length=11) )
  loc_gz = expand.grid("x"=seq(xlim[1],xlim[2],length=10),"y"=seq(ylim[1],ylim[2],length=10))
  grad_gz = t( apply(as.matrix(loc_gz), MARGIN=1, FUN=grad, f=get_preference, beta=beta, gamma=gamma) )
  quiver( x=loc_gz[,1], y=loc_gz[,2], u=grad_gz[,1], v=grad_gz[,2], scale=25, lwd=2 )
}

set.seed(101)
png( file="Preference.png", units="in", res=200, width=6, height=6 )
  par( mfrow=c(2,2), mar=c(3,3,2,0), mgp=c(2,0.5,0) )
  plot_vectors( beta=beta_t[1], gamma=gamma_t[1], main="Initial (t=1)" )
  plot_vectors( beta=beta_t[n_t/2], gamma=gamma_t[n_t/2], main="Intermediate (t=500)" )
  plot_vectors( beta=beta_t[n_t], gamma=gamma_t[n_t], main="Final (t=1000)" )
  plot(1, type="n", xlim=c(-200,200), ylim=c(-200,200), xlab="x", ylab="y", main="Simulated tracks" )
  for( track in 1:10 ){
    s = simulate_track(n_t=n_t, beta_t=beta_t, gamma_t=gamma_t, get_preference=get_preference)
    s = s[round(seq(1,nrow(s),length=20)), ]
    quiver( x=s[-nrow(s),1], y=s[-nrow(s),2], u=diff(s[,1]), v=diff(s[,2]), scale=1, lwd=1,
            col=viridisLite::viridis(10)[track], length=0.05 )
  }
dev.off()

#############
# Fit using crawlTMB
#############

# Compile
compile( "crawlTMB.cpp" )
dyn.load( dynlib("crawlTMB") )

#
set.seed(4)
x_iz = simulate_track(n_t=n_t, beta_t=beta_t, gamma_t=gamma_t, get_preference=get_preference)
y_iz = x_iz + mvtnorm::rmvnorm(n=nrow(x_iz), sigma=diag(2) )
DeltaT_i = c(NA,rep(1,nrow(y_iz)-1))
t_i = c(0, cumsum(DeltaT_i[-1]))

########## START IN-TEXT SNIPPET
# Define drift as a function of time
formula = ~ 0
X_ij = as.matrix( model.matrix(formula, data=data.frame("t_i"=t_i)) )

# Build inputs
Data = list( "y_iz"=y_iz,
             "DeltaT_i"=DeltaT_i,
             "error2_i"=rep(1,nrow(y_iz)),
             "X_ij"=X_ij,
             "n_factors"=0,
             "RAM"=matrix(nrow=0,ncol=4) )
Parameters = list( "sigma2_z"=log(1),
                   "x_iz"=Data$y_iz,
                   "beta_jz"=matrix(0,nrow=ncol(Data$X_ij),ncol=2) )
Random = "x_iz"

# drop some data
which_include = seq(1, nrow(Data$y_iz), length=50 )
Data$y_iz[-which_include,] = NA

# Build and fit object
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random )
Opt = nlminb( start=Obj$par, obj=Obj$fn, grad=Obj$gr )
Opt$SD = sdreport( Obj )
########## END IN-TEXT SNIPPET

# Plot sparsity
png( "sparsity.png", width=4, height=4, unit='in', res=200 )
  par( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )
  Matrix::image( Obj$env$spHess(par=Obj$env$last.par.best,random=TRUE) )
dev.off()

# Plot smoother
plot_track = function( Opt, Data, x_iz, cols=1:2, which_include, ... ){
  tplot = seq(1, nrow(Data$y_iz), by=1 )
  G = as.list(Opt$SD,what="Estimate",report=TRUE)$Gsum_iz[,cols]
  X = as.list(Opt$SD,what="Estimate")$x_iz[,cols]
  Xsd = as.list(Opt$SD,what="Std. Error")$x_iz[,cols]
  plot( type="n", x_iz[,cols], pch=20, col="red", cex=1.5, ... )
  plotrix::draw.ellipse( x=X[tplot,1], y=X[tplot,2], a=Xsd[tplot,1], b=Xsd[tplot,2], col=rgb(0,0,1,alpha=0.1), border=NA )
  lines( G, pch=20, col="red", cex=1.5, lwd=2 )
  points( x_iz[,cols], pch=20, col="black", cex=1 )
  quiver( x=X[which_include[-length(which_include)],1], y=X[which_include[-length(which_include)],2],
    u=diff(X[which_include,1]), v=diff(X[which_include,2]), col="green", length=0.1, scale=0.8 )
}
png( "fitted_diffusion.png", width=4, height=4, res=200, units="in" )
  plot_track( Opt=Opt, Data=Data, x_iz=x_iz, xlab="Eastings", ylab="Northings", which_include=which_include )
dev.off()

########## START IN-TEXT SNIPPET
# Define drift as a function of time
formula = ~ splines::bs(t_i,5)
Data$X_ij = as.matrix( model.matrix(formula, data=data.frame("t_i"=t_i)) )
Parameters$beta_jz = matrix(0,nrow=ncol(Data$X_ij),ncol=2)

# Refit model
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random )
Opt = nlminb( start=Obj$par, obj=Obj$fn, grad=Obj$gr )
Opt$SD = sdreport( Obj )
########## START IN-TEXT SNIPPET
png( "fitted_drift.png", width=4, height=4, res=200, units="in" )
  plot_track( Opt=Opt, Data=Data, x_iz=x_iz, xlab="Eastings", ylab="Northings", which_include=which_include )
dev.off()

########################
# Northern fur seal demo
########################

library(sf)

DF = read.csv("FSdata_2016.csv")
DF = st_as_sf(DF, coords=c("longitude","latitude"), crs="+proj=longlat" )
DF = st_transform( DF, crs="+proj=utm +datum=WGS84 +units=km +zone=2" )

DF = subset(DF, tripno==1 & dbid==818)
#DF = subset(DF, dbid==818)
DF$duration = c(NA, diff(as.POSIXlt(DF$gmt)))
DF$duration[-1] = ifelse( DF$tripno[-1]==DF$tripno[-nrow(DF)] & DF$dbid[-1]==DF$dbid[-nrow(DF)], DF$duration[-1], NA  )
DF$t_i = c(0, cumsum(DF$duration[-1]))

# Build drift covariances
formula = ~ splines::bs(t_i,3)
X_ij = as.matrix( model.matrix(formula, data=data.frame("t_i"=t_i)) )

# Build inputs
Data = list( "y_iz"=st_coordinates(DF), "DeltaT_i"=DF$duration, "error2_i"=(DF$error_semi_major/1000)^2,
             "X_ij"=X_ij, "n_factors"=0, "RAM"=matrix(nrow=0,ncol=4) )
Parameters = list( "sigma2_z"=log(1), "x_iz"=Data$y_iz, "beta_jz"=matrix(0,nrow=ncol(Data$X_ij),ncol=2) )
Random = c("x_iz")

# drop some data
which_include = seq(1, nrow(Data$y_iz), length=20 )
Data$y_iz[-which_include,] = NA

# Build and fit object
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random )  #
Opt = nlminb( start=Obj$par, obj=Obj$fn, grad=Obj$gr )
Opt$SD = sdreport( Obj )
png( "fitted_NFS.png", width=4, height=4, res=200, units="in" )
  plot_track( Opt=Opt, Data=Data, x_iz=st_coordinates(DF), xlab="Eastings (km)", ylab="Northings (km)", which_include=which_include )
dev.off()

######################
# Simulate multiple tracks
######################

########## START IN-TEXT SNIPPET
# Parameters
gamma_t = beta_t = seq(0.5,-0.5,length=n_t)
beta_t = ifelse( beta_t>0, 0, beta_t )
gamma_t = ifelse( gamma_t<0, 0, gamma_t )

# Simulate three tracks
n_tracks = 4
set.seed(101)
x_iz = NULL
for( track in 1:n_tracks ){
  s = simulate_track( n_t=n_t, beta_t=beta_t,
                      gamma_t=gamma_t, get_preference=get_preference)
  colnames(s) = paste0( c("x","y"), track )
  x_iz = cbind(x_iz, s)
}
y_iz = x_iz + array(rnorm(prod(dim(x_iz))), dim=dim(x_iz))

# Build inputs
Data = list( "y_iz"=y_iz,
             "DeltaT_i"=c(NA,rep(1,nrow(y_iz)-1)),
             "error2_i"=rep(1,nrow(y_iz)),
             "X_ij"=array(dim=c(nrow(y_iz),0)),
             "n_factors"=0,
             "RAM"=matrix(nrow=0,ncol=4) )
Params = list( "sigma2_z"=log(1),
               "x_iz"=Data$y_iz,
               "beta_jz"=array(0,dim=c(ncol(Data$X_ij),ncol(Data$y_iz))) )
Random = c("x_iz")

# drop some data
which_include = seq(1, nrow(Data$y_iz), length=20*n_tracks )
Data$y_iz[-which_include,] = NA

# Build and fit object
Obj = MakeADFun(data=Data, parameters=Params, random=Random )
Opt = nlminb( start=Obj$par, obj=Obj$fn, grad=Obj$gr )
########## END IN-TEXT SNIPPET
Opt$SD = sdreport( Obj )
Opt$AIC = 2*Opt$objective + 2*length(Opt$par)
V_zz = Obj$report()$V_zz
dimnames(V_zz) = list( colnames(y_iz), colnames(y_iz) )

# Visualize habitat
set.seed(101)
png( file="Tracks-Correlated.png", units="in", res=200, width=6, height=6 )
  par( mfrow=c(2,2), mar=c(3,3,2,0), mgp=c(2,0.5,0) )
  plot_vectors( beta=beta_t[1], gamma=gamma_t[1], main="Initial (t=1)" )
  plot_vectors( beta=beta_t[n_t/2], gamma=gamma_t[n_t/2], main="Intermediate (t=500)" )
  plot_vectors( beta=beta_t[n_t], gamma=gamma_t[n_t], main="Final (t=1000)" )
  plot(1, type="n", xlim=c(-200,200), ylim=c(-200,200), xlab="x", ylab="y", main="Simulated tracks" )
  for( track in 1:10 ){
    s = simulate_track(n_t=n_t, beta_t=beta_t, gamma_t=gamma_t, get_preference=get_preference)
    s = s[round(seq(1,nrow(s),length=20)), ]
    quiver( x=s[-nrow(s),1], y=s[-nrow(s),2], u=diff(s[,1]), v=diff(s[,2]), scale=1, lwd=1,
            col=viridisLite::viridis(10)[track], length=0.05 )
  }
dev.off()

# Visualize tracks
png( "fitted_joint.png", width=8, height=8, res=200, units="in" )
  par( mfrow=c(2,2), mar=c(2,2,2,0) )
  for(track in 1:n_tracks) plot_track( Opt=Opt, Data=Data, x_iz=x_iz, cols=1:2+2*(track-1),
    main=paste("Animal",track), xlab="Eastings", ylab="Northings", which_include=which_include ) # , xlim=c(-200,200), ylim=c(-20,200)
dev.off()

# Switch to factor model
Data1 = Data
Data1$n_factors = 2
Params$sigma2_z = rep(0.1,2*ncol(y_iz))
Obj = MakeADFun(data=Data1, parameters=Params, random=Random )  #
Opt1 = nlminb( start=Obj$par, obj=Obj$fn, grad=Obj$gr )
Opt1$SD = sdreport( Obj )
Opt1$AIC = 2*Opt1$objective + 2*length(Opt1$par)
V1_zz = Obj$report()$V_zz
dimnames(V1_zz) = list( colnames(y_iz), colnames(y_iz) )

########## START IN-TEXT SNIPPET
# Specify SEM
text = "
  y1 -> y2, y
  y1 -> y3, y
  y1 -> y4, y
"
SEM_model = sem::specifyModel( text = text,
                               exog.variances = TRUE,
                               endog.variances = TRUE,
                               covs = colnames(y_iz) )
RAM = build_ram( SEM_model, colnames(y_iz) )
########## END IN-TEXT SNIPPET

# Build with RAM
Data2 = Data
Data2$n_factors = 0
Data2$RAM = as.matrix(RAM[,1:4])
Data2$RAM[Data2$RAM[,1]==2,4] = 2
Params$sigma2_z = rep(0.1, max(Data2$RAM[,4]) )
Obj = MakeADFun(data=Data2, parameters=Params, random=Random )
Opt2 = nlminb( start=Obj$par, obj=Obj$fn, grad=Obj$gr )
Opt2$SD = sdreport( Obj )
Opt2$AIC = 2*Opt2$objective + 2*length(Opt2$par)
V2_zz = Obj$report()$V_zz
dimnames(V2_zz) = list( colnames(y_iz), colnames(y_iz) )

png( "fitted_joint-tracks.png", width=8, height=4, res=200, units="in" )
  par( mfrow=c(1,3), oma=c(0,0,0,0), mar=c(2,2,2,0) )
  plot_track( Opt=Opt, Data=Data, x_iz=x_iz, cols=1:2, main="Independent-and-equal", which_include=which_include )
    legend("bottomright", legend=paste0("AIC = ",signif(Opt$AIC,4)), bty="n" )
  plot_track( Opt=Opt1, Data=Data1, x_iz=x_iz, cols=1:2, main="Factor model", which_include=which_include )
    legend("bottomright", legend=paste0("AIC = ",signif(Opt1$AIC,4)), bty="n" )
  plot_track( Opt=Opt2, Data=Data2, x_iz=x_iz, cols=1:2, main="Structural equations", which_include=which_include )
    legend("bottomright", legend=paste0("AIC = ",signif(Opt2$AIC,4)), bty="n" )
dev.off()

# Plot covariance
png( "fitted_joint-covariance.png", width=8, height=4, res=200, units="in" )
  par( mfrow=c(1,3), mar=c(0,0,0,0), oma=c(5,5,5,5) )
  corrplot::corrplot( (V_zz), is.corr=FALSE, tl.pos="lt", type="lower", col=viridisLite::viridis(20), cl.cex=0.8, cl.length=6 )
  mtext( side=3, text="Independent-and-equal", line=0 )
  corrplot::corrplot( (V1_zz), is.corr=FALSE, tl.pos="lt", type="lower", col=viridisLite::viridis(20), cl.cex=0.8, cl.length=6 )
  mtext( side=3, text="Factor model" )
  corrplot::corrplot( (V2_zz), is.corr=FALSE, tl.pos="lt", type="lower", col=viridisLite::viridis(20), cl.cex=0.8, cl.length=6 )
  mtext( side=3, text="Structural equations" )
dev.off()
