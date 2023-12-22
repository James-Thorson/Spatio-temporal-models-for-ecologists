
setwd( R'(C:\Users\James.Thorson\Desktop\Git\Spatio-temporal-models-for-ecologists\Chap_7)' )
source( "../Shared_functions/marginaleffects.R" )
source( "../Shared_functions/sample_var.R" )
source( "../Shared_functions/rmvnorm_prec.R" )
source( "../Shared_functions/add_legend.R" )
Sys.setenv("_SP_EVOLUTION_STATUS_"=2)

library(sf)
library(TMB)
library(viridisLite)
library(rnaturalearth)
library(fmesher)
library(terra)

# Read data
Data = read.csv( "red_snapper_data.csv" )
sf_data = st_as_sf( Data, coords=c("Lon","Lat"), crs=st_crs("EPSG:4326") )
sf_data$Year = factor(sf_data$Year)
extent = st_read( "red_snapper_extent.shp", crs = st_crs("EPSG:4326") )
coast_sf = ne_coastline(scale=10, returnclass="sf")

# Download bathymetry and save copy for stability
# library(marmap)
# bathy = getNOAA.bathy( lon1 = -99, lon2 = -79,
#                         lat1 = 23, lat2 = 32,
#                         resolution = 4, keep=TRUE)
# bathy_raster = as.raster(bathy)
# writeRaster( bathy_raster, filename="bathy" )

#
bathy = rast( "bathy.grd" )
sf_bathy = st_as_sf( as.data.frame(bathy/1000,xy=TRUE), coords=c("x","y"), crs = st_crs("EPSG:4326") )

# Add bathymetry at sample points
nearest = RANN::nn2( st_coordinates(sf_bathy), st_coordinates(sf_data), k=1 )
sf_data$bathy = st_drop_geometry(sf_bathy)[nearest$nn.idx]

# create mesh
mesh = fm_mesh_2d( st_coordinates(sf_data), plot.delay=NULL, cutoff=0.5)
# Create matrices in INLA
spde <- fm_fem(mesh, order=2)
# create projection matrix from vertices to samples
A_is = fm_evaluator( mesh, loc=st_coordinates(sf_data) )$proj$A
# Create extrapolation grid
cellsize = 0.1
sf_fullgrid = st_make_grid( extent, cellsize=cellsize )
sf_grid = st_make_valid(st_intersection( sf_fullgrid, extent ))
# Add bathymetry
nearest = RANN::nn2( st_coordinates(sf_bathy), st_coordinates(st_centroid(sf_grid)), k=1 )
sf_grid = st_sf( sf_grid, "bathy"=st_drop_geometry(sf_bathy)[nearest$nn.idx] )
# create projection matrix from vertices to grid
A_gs = fm_evaluator( mesh, loc=st_coordinates(st_centroid(sf_grid)) )$proj$A
a_g = st_area( sf_grid )

#
plot( st_geometry(sf_grid) )
plot( sf_data[1], add=TRUE )

# Define covariates
Q_formula = ~ 0 + Data_type
X_formula = ~ 0 + poly( bathy, 2, raw=TRUE ) + Year

# Make Q-matrix
Q_ij = model.matrix( Q_formula, sf_data )
Q_ij = Q_ij[, c("Data_typeCount","Data_typeEncounter") ]

# Make X-matrix for data
X_ik = model.matrix( X_formula, sf_data )

# Make X-matrix for grid
sf_grid$Year = factor(2014,levels=sort(unique(sf_data$Year)))
frame0 = model.frame( formula=X_formula, data=sf_data )
terms0 = terms( frame0 )
xlevels = .getXlevels( terms0, frame0 )
terms1 = delete.response( terms0 )
frame1 = model.frame( terms1, sf_grid, xlev=xlevels )
X_gk = model.matrix( terms1, frame1 )

# COmpile
compile( "integrated_model.cpp" )
dyn.load( dynlib("integrated_model") )

# Build object
Data = list( "c_i" = sf_data$Response,
             "e_i" = as.numeric(factor(sf_data$Data_type,levels=c("Encounter","Count","Biomass_KG")))-1,
             "X_ik" = X_ik,
             "X_gk" = X_gk,
             "Q_ij" = Q_ij,
             "A_is" = A_is,
             "A_gs" = A_gs,
             "M0" = spde$c0,
             "M1" = spde$g1,
             "M2" = spde$g2 )
Params = list( "ln_tau" = 0,
               "ln_kappa" = 0,
               "ln_phi" = 0,
               "finv_power" = 0,
               "gamma_k" = rep(0,ncol(X_ik)),
               "eta_j" = rep(0,ncol(Q_ij)),
               "omega_s" = rnorm(nrow(spde$c0),sd=0) )
obj = MakeADFun( data=Data, parameters=Params, random="omega_s" )
# Optimize
opt = nlminb( objective=obj$fn, grad=obj$gr, start=obj$par )
opt$SD = sdreport( obj, getJointPrecision=TRUE )
report = obj$report()

# Plot response
library(marginaleffects)
library(ggplot2)
library(gridExtra)
fit = list( obj=obj,
            opt=opt,
            formula=X_formula,
            parhat=obj$env$parList(),
            data = sf_grid )
class(fit) = "custom_tmb"
quant = function(x) seq(min(x),max(x),length=21)
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}
options("marginaleffects_model_classes" = "custom_tmb")

# Get prediction -- bathymetry
newdata = datagrid( newdata=data.frame(sf_grid)[,c('bathy','Year')], bathy=quant, Year=getmode )
  pred = predictions( fit, newdata=newdata, param="gamma_k", center=TRUE )
g1 = ggplot( pred, aes(bathy, estimate)) +
  geom_line( aes(y=estimate), color="blue", size=1 ) +
  geom_ribbon( aes( x=bathy, ymin=conf.low, ymax=conf.high), fill=rgb(0,0,1,0.2) ) +
  labs(y="Predicted response")
# Get prediction -- Year
newdata = datagrid( newdata=data.frame(sf_grid)[,c('bathy','Year')], bathy=mean, Year=levels(sf_data$Year) )
  pred = predictions( fit, newdata=newdata, param="gamma_k", center=TRUE )
g2 = ggplot(pred, aes(x=Year, y=estimate)) +
  geom_point( ) +
  geom_errorbar(aes(ymin=conf.low , ymax=conf.high ), width=.2)
ggsave( plot=grid.arrange(g1,g2, nrow=1), paste0("covariate_response.png"), width=6, height=3 )

#
SE_logmu_g = sample_var( obj=obj, var_name="logmu_g", mu=obj$env$last.par.best, prec=opt$SD$jointPrecision )
plotdata = st_sf( st_geometry(sf_grid), "log_density"=report$logmu_g, "SE"=SE_logmu_g )
  plotdata$SE = ifelse( plotdata$log_density<(max(plotdata$log_density)-log(1000)), NA, plotdata$SE )
  plotdata$log_density = ifelse( plotdata$log_density<(max(plotdata$log_density)-log(1000)), NA, plotdata$log_density )
png( file="Predicted densities.png", width=6, height=7.5, res=200, units="in" )
  par( mfrow=c(3,1), mar=c(0,0,2,0), tck=-0.02 )
  #
  plot( st_geometry(extent), add=FALSE, col="grey", lwd=0.01, main="Sample locations" )
  plot( st_geometry(coast_sf), add=TRUE )
  plot( sf_data['Data_type'], border=NA, add=TRUE, pch=20, pal=sf.colors(n=3,alpha=0.3) )
  legend("topright", bty="n", legend=unique(sf_data$Data_type), fill=sf.colors(n=3) )
  #
  plot( st_geometry(extent), add=FALSE, col="grey", lwd=0.01, main="log(density)" )
  plot( st_geometry(coast_sf), add=TRUE )
  plot( plotdata["log_density"], border=NA, add=TRUE )
  add_legend( legend=round(range(plotdata$log_density,na.rm=TRUE),2), legend_x=c(0.95,1) )
  #
  plot( st_geometry(extent), add=FALSE, col="grey", lwd=0.01, main="SE[log(density)]" )
  plot( st_geometry(coast_sf), add=TRUE )
  plot( plotdata["SE"], border=NA, add=TRUE )
  add_legend( legend=round(range(plotdata$SE,na.rm=TRUE),2), legend_x=c(0.95,1) )
dev.off()

