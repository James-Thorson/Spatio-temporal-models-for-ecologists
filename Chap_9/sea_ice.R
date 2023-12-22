
setwd( "C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_9/")

Ice = read.csv( "Ice.csv" )

library(sf)
library(fmesher)
library(TMB)
library(rnaturalearth)
source("C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Shared_functions/rmvnorm_prec.R")
source( "C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Shared_functions/rotate_pca.R")
source( "C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Shared_functions/add_legend.R")

# project data
sf_ice = st_as_sf( Ice, coords = c("Longitude","Latitude") )
st_crs(sf_ice) = "+proj=longlat +datum=WGS84"
sf_ice = st_transform( sf_ice, crs=st_crs("+proj=laea +lat_0=90 +lon_0=-30 +units=km") )

# 
sf_pole = st_point( c(0,90) )
sf_pole = st_sfc( sf_pole, crs="+proj=longlat +datum=WGS84" )
sf_pole = st_transform( sf_pole, crs=st_crs("+proj=laea +lat_0=90 +lon_0=-30 +units=km") )
sf_pole = st_buffer( sf_pole, dist=3000 )
sf_ice = st_intersection( sf_ice, sf_pole )

# Country shapefiles for plotting
sf_maps = ne_countries( return="sf", scale=10, continent=c("north america","europe","asia") )
sf_maps = st_transform( sf_maps, crs=st_crs("+proj=laea +lat_0=90 +lon_0=-30 +units=km") )
sf_maps = st_union( sf_maps )

# Shapefile for water
sf_water = st_difference( st_as_sfc(st_bbox(sf_maps)), sf_maps )

#
plot( sf_ice[,'Ice'], reset=FALSE )
plot( st_geometry(sf_maps), add=TRUE, border="brown", lwd=2 )

# create mesh
xy_i = st_coordinates(sf_ice)
mesh = fm_mesh_2d( xy_i, plot.delay=NULL, cutoff=200)
# Create matrices in INLA
spde <- fm_fem(mesh, order=2)
# create projection matrix from vertices to samples
A_is = fm_evaluator( mesh, loc=xy_i )$proj$A
# Create extrapolation grid
cellsize = 50
sf_grid = st_make_grid( sf_pole, cellsize=cellsize )
# Restrict to water
grid_i = st_intersects( sf_water, sf_grid )
sf_grid = sf_grid[ unique(unlist(grid_i)) ]
# Restrict to 3000 km from North Pole
grid_i = st_intersects( sf_pole, sf_grid )
sf_grid = sf_grid[ unique(unlist(grid_i)) ]
# create projection matrix from vertices to grid
A_gs = fm_evaluator( mesh, loc=st_coordinates(st_centroid(sf_grid)) )$proj$A
# Area for each grid cell
a_g = rep( cellsize^2, length(sf_grid) )
# sum(a_g) < pi*3000^2 due to missing land

# COmpile
compile( "EOF.cpp" )
dyn.load( dynlib("EOF") )

######## ######## SPDE-based
# Build object
rmatrix = function(nrow,ncol,sd=1) matrix(sd*rnorm(nrow*ncol), nrow=nrow)
n_factors = 2
n_years = max(sf_ice$Year)-min(sf_ice$Year)+1
Data = list( "y_i" = sf_ice$Ice,
             "t_i" = sf_ice$Year-min(Ice$Year),
             "A_is" = A_is,
             "A_gs" = A_gs,
             "a_g" = a_g,
             "M0" = spde$c0,
             "M1" = spde$g1,
             "M2" = spde$g2 )
Params = list( "beta0" = 0,
               "ln_tau" = 0,
               "ln_kappa" = log(1),
               "ln_sigma" = 0,
               "L_ft" = rmatrix(nrow=n_factors, ncol=n_years),
               "omega_s" = rnorm( nrow(spde$c0) ),
               "epsilon_sf" = rmatrix(nrow=nrow(spde$c0), ncol=n_factors) )

#
Map = list( "L_ft"=matrix(1:prod(dim(Params$L_ft)),nrow=n_factors) )
Map$L_ft[lower.tri(Map$L_ft)] = NA
Map$L_ft = factor(Map$L_ft)
Params$L_ft[lower.tri(Params$L_ft)] = 0

# Build
Obj = MakeADFun( data=Data, parameters=Params, map=Map, random=c("omega_s","epsilon_sf") )
Obj$env$beSilent()

# Optimize
Opt = nlminb( start=Obj$par, obj=Obj$fn, gr=Obj$gr, control=list(trace=1, eval.max=1e4, iter.max=1e4) )
Opt$SD = sdreport( Obj )
#Opt = TMBhelper::fit_tmb( obj=Obj, newtonsteps=0, getsd=TRUE, bias.correct=FALSE )
report = Obj$report()
parhat = Obj$env$parList()

#
L_tf = t(parhat$L_ft)
epsilon_gf = report$epsilon_gf
omega_g = report$omega_g

#
rotated_results = rotate_pca( L_tf=L_tf, x_sf=epsilon_gf, order="decreasing" )
L_tf = rotated_results$L_tf
epsilon_gf = rotated_results$x_sf

# Plot
png( file=paste0("EOF=",n_factors,".png"), width=8, height=8, res=200, units="in")
  sf_plot = st_sf( sf_grid, "omega"=omega_g, "epsilon"=epsilon_gf )
  par(mfrow=c(2,2), oma=c(2,2,0,0) )
  plot( sf_plot[,'omega'], reset=FALSE, key.pos=NULL, border=NA )
  plot( st_geometry(sf_maps), add=TRUE, border=NA, col="grey" ) # , reset=FALSE
  add_legend( legend=round(range(sf_plot[['omega']]),2), legend_x=c(0.95,1), legend_y=c(0.55,0.95) )
  for( f in 1:n_factors ){
    plot( sf_plot[,paste0("epsilon.",f)], reset=FALSE, key.pos=NULL, border=NA )
    plot( st_geometry(sf_maps), add=TRUE, border=NA, col="grey" ) # , reset=FALSE
    add_legend( legend=round(range(sf_plot[[paste0("epsilon.",f)]]),2), legend_x=c(0.95,1), legend_y=c(0.55,0.95) )
  }
  matplot( y=L_tf, x=min(Ice$Year):max(Ice$Year), type="l", col=viridisLite::viridis(n_factors), lwd=2, lty="solid" )
  legend( "top", ncol=n_factors, legend=1:n_factors, fill=viridisLite::viridis(n_factors) )
dev.off()

png( file=paste0("sea_ice_EOF=",n_factors,".png"), width=8, height=8, res=200, units="in")
  sf_plot = st_sf( sf_grid, "conc"=report$yhat_gt )
  par(mfrow=c(5,5), oma=c(2,2,0,0) )
  for( t in 1:ncol(report$yhat_gt) ){
    plot( sf_plot[,paste0("conc.",t)], reset=FALSE, key.pos=NULL, border=NA, main=(min(Ice$Year):max(Ice$Year))[t] )
    plot( st_geometry(sf_maps), add=TRUE, border=NA, col="grey" ) # , reset=FALSE
    if(t==1) add_legend( legend=round(range(st_drop_geometry(sf_plot)),2), legend_x=c(0.95,1), legend_y=c(0.55,0.95) )
  }
dev.off()

#############
# Index plot
#############

Data = list( "y_i" = sf_ice$Ice,
             "t_i" = sf_ice$Year-min(Ice$Year),
             "A_is" = A_is,
             "A_gs" = A_gs,
             "a_g" = a_g,
             "M0" = spde$c0,
             "M1" = spde$g1,
             "M2" = spde$g2 )
Params = list( "beta0" = 0,
               "ln_tau" = 0,
               "ln_kappa" = log(1),
               "ln_sigma" = 0,
               "L_ft" = rmatrix(nrow=n_years, ncol=n_years),
               "omega_s" = rnorm( nrow(spde$c0) ),
               "epsilon_sf" = rmatrix(nrow=nrow(spde$c0), ncol=n_years) )

#
Map = list( "L_ft" = factor(ifelse(diag(n_years)==1,1,NA)) )
Params$L_ft = diag(n_years)

# Build
Obj_full = MakeADFun( data=Data, parameters=Params, map=Map, random=c("omega_s","epsilon_sf") )
Obj_full$env$beSilent()

# Optimize
#Opt_full = TMBhelper::fit_tmb( obj=Obj_full, newtonsteps=0, getsd=FALSE, bias.correct=FALSE )
Opt_full = nlminb( start=Obj_full$par, obj=Obj_full$fn, gr=Obj_full$gr, control=list(trace=1, eval.max=1e4, iter.max=1e4) )
Opt_full$SD = sdreport( Obj_full )
report_full = Obj_full$report()
parhat_full = Obj_full$env$parList()

png( "sea_ice_full.png", width=8, height=8, res=200, units="in")
  sf_plot = st_sf( sf_grid, "conc"=report$yhat_gt )
  par(mfrow=c(5,5), oma=c(2,2,0,0) )
  for( t in 1:ncol(report$yhat_gt) ){
    plot( sf_plot[,paste0("conc.",t)], reset=FALSE, key.pos=NULL, border=NA, main=(min(Ice$Year):max(Ice$Year))[t] )
    plot( st_geometry(sf_maps), add=TRUE, border=NA, col="grey" ) # , reset=FALSE
    if(t==1) add_legend( legend=round(range(st_drop_geometry(sf_plot)),2), legend_x=c(0.95,1), legend_y=c(0.55,0.95) )
  }
dev.off()

png( "sea_ice_extent.png", width=4, height=4, res=200, units="in")
  par( mar=c(3,4,1,1), mgp=c(2,0.5,0), tck=-0.02 )
  matplot( x=min(Ice$Year):max(Ice$Year),
           y=cbind(report_full$a_t,report$a_t)/1000,
           type='l',
           lwd=2,
           lty="solid",
           col=c("black","red"),
           ylab="Sept. concentration-weighted sea ice\nextent (1000 km^2)",
           xlab="Year" )
  legend( "bottomleft", bty='n', fill=c("black","red"), legend=c("full", paste0("EOF: dim=",n_factors)) )
dev.off()
