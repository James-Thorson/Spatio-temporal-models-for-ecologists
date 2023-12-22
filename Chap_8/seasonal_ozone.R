
# Data processed from: C:\Users\James.Thorson\Desktop\Work files\AFSC\2022-06 -- Breeding Bird survey
setwd("C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_8")

library(sf)
library(fmesher)
library(TMB)
library(splines2)
library(rnaturalearth)
source("../Shared_functions/rmvnorm_prec.R")
source("../Shared_functions/sample_var.R")
source("../Shared_functions/add_legend.R")

# Ozone data -- FULL
data_dir = R'(C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_6)'
ozone = st_read( file.path(data_dir,"2019_ozone.csv"), options=c("X_POSSIBLE_NAMES=X","Y_POSSIBLE_NAMES=Y"), crs=st_crs("+proj=longlat +datum=WGS84") )
  ozone$Daily.Max.8.hour.Ozone.Concentration = as.numeric(ozone$Daily.Max.8.hour.Ozone.Concentration)
  ozone$Date = as.POSIXlt( ozone$Date, tryFormats="%m/%d/%Y" )
  ozone$Julian = as.numeric(julian(ozone$Date, origin = as.POSIXct("01-01-2019", tryFormats="%m-%d-%Y"))) / 365
ozone = ozone[ which(ozone$Daily.Max.8.hour.Ozone.Concentration>0), ]

# Define states
state_set = c("Virginia","North Carolina","South Carolina","Georgia","Florida","Maryland","Delaware","District of Columbia","New Jersey","Pennsylvania","New York")
states_sf = ne_states( c("United States of America","Canada"), return="sf")
states_sf = states_sf[pmatch(state_set, states_sf$name),]
domain_sf = st_union( states_sf )

# Create extrapolation grid
cellsize = 0.25
grid = st_make_grid( states_sf, cellsize=cellsize )
grid = st_intersection( grid, domain_sf )

# create mesh
mesh = fm_mesh_2d( st_coordinates(st_centroid(grid)), refine=TRUE, cutoff=0.2)
# Create matrices in INLA
spde <- fm_fem(mesh, order=2)

#
season_df = 6
season_plots = 12

#
A_is = fm_evaluator( mesh, loc=st_coordinates(ozone) )$proj$A
A_gs = fm_evaluator( mesh, loc=st_coordinates(st_centroid(grid)) )$proj$A
Afull_gs = A_gs
for( r in 2:season_plots ) Afull_gs = rbind( Afull_gs, A_gs )

# formula
formula = ~ 0 + mSpline(Julian, knots=seq(0,1,length=season_df+1)[2:season_df], degree=3, intercept=TRUE, periodic=TRUE, Boundary.knots=c(0,1) )  #
options(na.action='na.pass')
X_ij = model.matrix( formula, data=data.frame("Julian"=c(ozone$Julian)) )
#X_gj = model.matrix( formula, data=data.frame("Julian"=rep(0.5,length(grid))) )
X_gj = model.matrix( formula, data=data.frame("Julian"=rep(seq(0,1,length=season_plots+1)[1:season_plots],each=length(grid))) )

#
png( file="Cyclic_spline_basis_functions.png", width=5, height=5, res=200, units="in")
  par( mar=c(2,2,1,1), mgp=c(2,0.5,0), tck=-0.02 )
  matplot( y=X_ij[order(ozone$Julian),], x=ozone$Julian[order(ozone$Julian)], type="l", lwd=2, lty="solid", col=viridisLite::viridis(6) )
dev.off()

#
plot_response = function( formula, name, f=identity, seed=101 ){
  set.seed(seed)
  xpred = seq( 0, 1, length=1000)
  X_ij = model.matrix( formula, data=data.frame(Julian=xpred) )
  beta_zj = mvtnorm::rmvnorm( n=10, mean=rep(0,ncol(X_ij)), sigma=diag(ncol(X_ij)) )
  ypred_iz = X_ij %*% t(f(beta_zj))
  png( name, width=6, height=6, res=200, units='in')
    par( mfrow=c(2,1), mar=c(3,3,3,1), mgp=c(2,0.5,0) )
    matplot( x=xpred, y=X_ij, type="l", lwd=2, lty="solid", main="f(x)", xlab="x" )
    matplot( x=xpred, y=ypred_iz, type="l", lwd=2, col=viridisLite::viridis(nrow(beta_zj)), lty="solid", main="Potential response functions", xlab="Elevation" )
  dev.off()
}
plot_response( formula=formula, name="Basis-cyclic.png" )

# COmpile
compile( "SPDE_SVC.cpp" )
dyn.load( dynlib("SPDE_SVC") )

# Build object
Data = list( "y_i" = ozone$Daily.Max.8.hour.Ozone.Concentration,
             "A_is" = A_is,
             "A_gs" = Afull_gs,
             "M0" = spde$c0,
             "M1" = spde$g1,
             "M2" = spde$g2,
             "X_ij" = X_ij,
             "X_gj" = X_gj )
Params = list( "ln_tau" = 0,
               "ln_kappa" = 0,
               "ln_sigma" = 0,
               "beta_j" = rep(0,ncol(X_ij)),
               "xi_sj" = matrix(0,ncol=ncol(X_ij),nrow=mesh$n) )
Obj = MakeADFun( data=Data, parameters=Params, random="xi_sj" )

# Optimize
Opt = nlminb( start=Obj$par, obj=Obj$fn, gr=Obj$gr )

## Plot
report = Obj$report()
grid_sf = st_sf(grid, Ozone=matrix(report$yhat_g,ncol=season_plots) )
names(grid_sf)[1:12] = paste0( month.abb," 1" )

# Plot seasonal dynamics
#png( "Seasonal_ozone.png", width=6, height=9, res=200, units="in")
#  plot(grid_sf, max.plot=12, border=NA )
#dev.off()

#
png( "Seasonal_ozone.png", width=7, height=9, res=200, units="in")
  par( mfrow=c(3,4), mar=c(1,1,1,1) )
  zlim = range( data.frame(grid_sf)[,1:12] )
  for( i in 1:season_plots ){
    plot( grid_sf[i], key.pos=NULL, reset=FALSE, border=NA, breaks=seq(zlim[1],zlim[2],length=100) )
    plot( st_geometry(states_sf), add=TRUE )
  }
  add_legend( legend=round(zlim,3) )
dev.off()

#
png( "Seasonal_ozone_V2.png", width=7, height=9, res=200, units="in")
  par( mfrow=c(3,4), mar=c(1,1,1,1) )
  for( i in 1:season_plots ){
    zlim = range( data.frame(grid_sf)[,i] )
    plot( grid_sf[i], key.pos=NULL, reset=FALSE, border=NA, breaks=seq(zlim[1],zlim[2],length=100) )
    plot( st_geometry(states_sf), add=TRUE )
    add_legend( legend=round(zlim,3) )
  }
dev.off()

#png( "Mapped_ozone.png", width=7, height=7, res=200, units="in")
#  par( mfrow=c(2,2), mar=c(1,1,1,1) )
#  plot( st_geometry(states_sf), main=expression(Mesh ~~ and ~~ data) )
#  plot( mesh, add=TRUE )
#  plot( ozone[,'Daily.Max.8.hour.Ozone.Concentration'], add=TRUE )
#  add_legend( legend=round(range(ozone$Daily.Max.8.hour.Ozone.Concentration),2) )
#  plot( plot_grid['Ozone'], key.pos=NULL, reset=FALSE, main=expression(Predicted ~~ ozone ~~ (ppm)) )
#  add_legend( legend=round(range(plot_grid$Ozone),2) )
#  plot( plot_grid['SE'], key.pos=NULL, reset=FALSE, main=expression(Standard ~~ Error ~~ (ppm)) )
#  add_legend( legend=round(range(plot_grid$SE),3) )
#  plot( plot_grid['Dens2020'], key.pos=NULL, reset=FALSE, main=expression(Density ~~ (Count/km^2)), logz=TRUE )
#  add_legend( legend=round(range(plot_grid$Dens2020),2) )
#dev.off()
