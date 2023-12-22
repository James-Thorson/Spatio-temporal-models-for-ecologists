
data_dir = "C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_8/"
setwd( "C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_9/" )

# Libraries and functions
library(sf)
library(fmesher)
library(TMB)
library(rnaturalearth)
source("C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Shared_functions/rmvnorm_prec.R")
source("C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Shared_functions/sample_var.R")
source("C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Shared_functions/add_legend.R")

#
pollock = readRDS( file.path(data_dir,"pollock.rds") )
pollock = st_as_sf(pollock, coords = c("Long","Lat"), crs="+proj=longlat +datum=WGS84" )
pollock = st_transform( pollock, crs=st_crs("+proj=utm +zone=2 +datum=WGS84 +units=km") )
coldpool = readRDS( file.path(data_dir,"coldpool.rds") )
survey_domain = readRDS( file.path(data_dir,"survey_domain.rds") )
survey_domain = st_sfc( survey_domain, crs="+proj=longlat +datum=WGS84" )
survey_domain = st_transform( survey_domain, crs=st_crs("+proj=utm +zone=2 +datum=WGS84 +units=km") )

# Country shapefiles for plotting
sf_maps = ne_countries( return="sf", scale=10, country=c("russia","united states of america") )
sf_maps = st_transform( sf_maps, crs=st_crs("+proj=utm +zone=2 +datum=WGS84 +units=km") )
sf_maps = st_intersection( sf_maps, st_as_sfc(st_bbox(survey_domain)) )
sf_maps = st_union( sf_maps )

# Make triangulated mesh
mesh = fm_mesh_2d( st_coordinates(pollock), cutoff=100, refine=TRUE )
# Create matrices in INLA
spde = fm_fem(mesh, order=2)
# create projection matrix from vertices to samples
A_is = fm_evaluator( mesh, loc=st_coordinates(pollock) )$proj$A
# Create extrapolation grid
cellsize = 25
grid = st_make_grid( survey_domain, cellsize=cellsize )
grid = st_intersection( grid, survey_domain )
# create projection matrix from vertices to grid
A_gs = fm_evaluator( mesh, loc=st_coordinates(st_centroid(grid)) )$proj$A

# Compile
compile( "pollock_SVC.cpp" )
dyn.load( dynlib("pollock_SVC") )

# Make inputs
year_set = min(pollock$Year):max(pollock$Year)
Data = list( "n_t" = length(year_set),
             "a_g" = as.numeric(st_area(grid)),
             "z_g" = st_coordinates(st_centroid(grid))[,2],
             "b_i" = pollock$Wt,
             "a_i" = pollock$AreaSwept_ha / 100, # Convert hectares to km^2
             "t_i" = pollock$Year - min(pollock$Year),
             "x_t" = scale(coldpool[match(year_set,coldpool$YEAR), 'AREA_LTE2_KM2']), 
             "A_is" = A_is,
             "A_gs" = A_gs,
             "M0" = spde$c0,
             "M1" = spde$g1,
             "M2" = spde$g2 )
Params = list( "beta_t" = rep(0,Data$n_t),
               "ln_tauO" = log(1),
               "ln_tauE" = log(1),
               "ln_tauX" = log(1),
               "ln_kappa" = 1,
               "ln_phi" = log(1),
               "logit_rhoE" = 0,
               "finv_power" = 0,
               "omega_s" = rep(0, mesh$n),
               "xi_s" = rep(0, mesh$n),
               "epsilon_st" = matrix(0, nrow=mesh$n, ncol=Data$n_t) )
Random = c("omega_s", "xi_s", "epsilon_st")

# Build and run
Obj = MakeADFun( data=Data, parameters=Params, random=Random )
Obj$env$beSilent()

# Run
Opt = nlminb( start=Obj$par, obj=Obj$fn, gr=Obj$gr, control=list(trace=1, eval.max=1e4, iter.max=1e4) )
Opt$SD = sdreport( Obj, bias.correct=FALSE, getJointPrecision=TRUE )
Report = Obj$report()

# making plotting grids
grid_xi = st_sf(grid, xi_g=Report$xi_g, crs=st_crs(grid) )
SE_g = sample_var( var_name="xi_g", n_samples=500, obj=Obj, mu=Obj$env$last.par.best, prec=Opt$SD$jointPrecision )
grid_SE = st_sf(grid, SE_xi_g=SE_g, crs=st_crs(grid) )

# make plot
png( "SVC_response.png", width=7, height=3, res=200, units="in")
  par( mfrow=c(1,2) )
  plot( grid_xi, key.pos=NULL, reset=FALSE, main="SVC response", border=NA, breaks=seq(min(Report$xi_g),max(Report$xi_g),length=1000) ) #
  plot( sf_maps, col="grey", add=TRUE )
  add_legend( legend=round(range(Report$xi_g),3), legend_y = c(0.5,0.95), legend_x = c(0.85, 0.9) )
  plot( grid_SE, key.pos=NULL, reset=FALSE, main="SE[SVC response]", breaks=seq(min(SE_g),max(SE_g),length=1000), border=NA )
  plot( sf_maps, col="grey", add=TRUE )
  add_legend( legend=round(range(SE_g),3), legend_y = c(0.5,0.95), legend_x = c(0.85, 0.9) )
dev.off()

#############
# Repeat with random time-series
#############

set.seed(123)
Data_null = Data
Data_null$x_t = sample(Data$x_t)
Obj_null = MakeADFun( data=Data_null, parameters=Params, random=Random )
Obj_null$env$beSilent()
Opt_null = nlminb( start=Obj_null$par, obj=Obj_null$fn, gr=Obj_null$gr, control=list(trace=1, eval.max=1e4, iter.max=1e4) )
Opt_null$SD = sdreport( Obj_null, bias.correct=FALSE, getJointPrecision=TRUE )
Report_null = Obj_null$report()

# making plotting grids
grid_xi = st_sf(grid, xi_g=Report_null$xi_g, crs=st_crs(grid) )
SE_g = sample_var( var_name="xi_g", n_samples=500, obj=Obj, mu=Obj_null$env$last.par.best, prec=Opt_null$SD$jointPrecision )
grid_SE = st_sf(grid, SE_xi_g=SE_g, crs=st_crs(grid) )

# make plot
png( "SVC_response_randomized.png", width=7, height=3, res=200, units="in")
  par( mfrow=c(1,2) )
  vec = c(-0.1,0.1,Report_null$xi_g)
  plot( grid_xi, key.pos=NULL, reset=FALSE, main="SVC response", border=NA, breaks=seq(min(vec),max(vec),length=1000) ) # 
  plot( sf_maps, col="grey", add=TRUE )
  add_legend( legend=round(range(vec),3), legend_y = c(0.5,0.95), legend_x = c(0.85, 0.9) )
  plot( grid_SE, key.pos=NULL, reset=FALSE, main="SE[SVC response]", breaks=seq(min(SE_g),max(SE_g),length=1000), border=NA )
  plot( sf_maps, col="grey", add=TRUE )
  add_legend( legend=round(range(SE_g),3), legend_y = c(0.5,0.95), legend_x = c(0.85, 0.9) )
dev.off()

