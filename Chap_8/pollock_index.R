
data_dir = "C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_8/"
setwd( "C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_8/" )

# Libraries and functions
library( TMB )
library( fmesher )
library( sf )
library( rasterVis )
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
survey_domain = st_transform( survey_domain, crs=st_crs(pollock) )

# Country shapefiles for plotting
sf_maps = rnaturalearth::ne_countries( return="sf", scale=10, country="united states of america" )
sf_maps = st_transform( sf_maps, crs=st_crs(pollock) )
#sf_maps = st_intersection( sf_maps, st_as_sfc(st_bbox(survey_domain)) )
#sf_maps = st_union( sf_maps )

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
Version = "pollock_index"
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )

# Make inputs
year_set = min(pollock$Year):max(pollock$Year)
Data = list( "n_t" = length(year_set),
             "a_g" = as.numeric(st_area(grid)),
             "z_g" = st_coordinates(st_centroid(grid))[,2],
             "b_i" = pollock$Wt,
             "a_i" = pollock$AreaSwept_ha / 100, # Convert hectares to km^2
             "t_i" = pollock$Year - min(pollock$Year),
             "A_is" = A_is,
             "A_gs" = A_gs,
             "M0" = spde$c0,
             "M1" = spde$g1,
             "M2" = spde$g2 )
Params = list( "beta_t" = rep(0,Data$n_t),
               "ln_tauO" = log(1),
               "ln_tauE" = log(1),
               "ln_kappa" = 1,
               "ln_phi" = log(1),
               "logit_rhoE" = 0,
               "finv_power" = 0,
               "omega_s" = rep(0, mesh$n),
               "epsilon_st" = matrix(0, nrow=mesh$n, ncol=Data$n_t) )
Random = c("omega_s", "epsilon_st")

# Build and run
Obj = MakeADFun( data=Data, parameters=Params, random=Random )
Obj$env$beSilent()

# Run
Opt = nlminb( start=Obj$par, obj=Obj$fn, gr=Obj$gr, control=list(trace=1, eval.max=1e4, iter.max=1e4) )
Report = Obj$report()
Opt$SD = sdreport( Obj, bias.correct=FALSE, getJointPrecision=TRUE )

# making plotting grids
SE_gt = sample_var( var_name="ln_d_gt", n_samples=500, obj=Obj, mu=Obj$env$last.par.best, prec=Opt$SD$jointPrecision )
ln_d_gt = ifelse( Report$ln_d_gt<max(Report$ln_d_gt-log(1e3)), NA, Report$ln_d_gt )
grid_density = st_sf(grid, ln_d_gt=ln_d_gt, crs=st_crs(grid) )
SE_gt = ifelse( Report$ln_d_gt<max(Report$ln_d_gt-log(1e3)), NA, SE_gt )
grid_SE = st_sf(grid, SE_ln_d_gt=SE_gt, crs=st_crs(grid) )

# make plot
years_to_plot = c(1982,1989,1996,2003,2010, 2015,2016,2017,2018,2019)
png( "pollock_densities_and_SEs.png", width=7, height=8, res=200, units="in")
  par( mfrow=c(length(years_to_plot)/2,4) )
  for(t in seq_along(years_to_plot)){
    tI = which( years_to_plot[t] == 1982:2019 )
    plot( grid_density[tI], key.pos=NULL, reset=FALSE, main=paste0("log-density: ", years_to_plot[t]),
          breaks=seq(min(ln_d_gt,na.rm=TRUE),max(ln_d_gt,na.rm=TRUE),length=1000), border=NA, pal=viridisLite::viridis )
    plot( st_geometry(pollock[which(pollock$Year==years_to_plot[t]),]), add=TRUE, pch=20, border=NA, cex=0.1 )
    plot( st_geometry(sf_maps), add=TRUE, col="brown" )
    if(t==1) add_legend( legend=round(range(ln_d_gt,na.rm=TRUE),1), legend_y=c(0.5,0.95), text_col="white", cex=0.8, col=viridisLite::viridis(10) )
    plot( grid_SE[tI], key.pos=NULL, reset=FALSE, main=paste0("SE[log-density]: ", years_to_plot[t]),
          breaks=seq(min(SE_gt,na.rm=TRUE),max(SE_gt,na.rm=TRUE),length=1000), border=NA )
    plot( st_geometry(sf_maps), add=TRUE, col="brown" )
    if(t==1) add_legend( legend=round(range(SE_gt,na.rm=TRUE),2), legend_y=c(0.5,0.95), text_col="white", cex=0.8 )
  }
dev.off()

# Biomass
calculate_total = function( ln_d_gt, a_g ){
  calc_t = rep(NA, ncol(ln_d_gt))
  for( t in seq_along(calc_t)){
    calc_t[t] = sum(a_g * exp(ln_d_gt[,t]))
  }
  return(calc_t)
}
B_tz = sample_var( var_name = "ln_d_gt",
                     n_samples = 500,
                     obj = Obj,
                     mu = Obj$env$last.par.best,
                     prec = Opt$SD$jointPrecision,
                     fun1 = calculate_total,
                     fun2 = function(x){c(mean(x),sd(x))},
                     a_g = Data$a_g )

# Center-of-gravity
calculate_COG = function( ln_d_gt, a_g, z_g ){
  calc_t = rep(NA, ncol(ln_d_gt))
  for( t in seq_along(calc_t)){
    calc_t[t] = weighted.mean( z_g, w=a_g*exp(ln_d_gt[,t]) )
  }
  return(calc_t)
}
COG_tz = sample_var( var_name = "ln_d_gt",
                     n_samples = 500,
                     obj = Obj,
                     mu = Obj$env$last.par.best,
                     prec = Opt$SD$jointPrecision,
                     fun1 = calculate_COG,
                     fun2 = function(x){c(mean(x),sd(x))},
                     z_g = Data$z_g,
                     a_g = Data$a_g )

# Equivalent area
calculate_EA = function( ln_d_gt, a_g ){
  calc_t = rep(NA, ncol(ln_d_gt))
  for( t in seq_along(calc_t)){
    calc_t[t] = sum(a_g*exp(ln_d_gt[,t])) / weighted.mean( exp(ln_d_gt[,t]), w=a_g*exp(ln_d_gt[,t]) )
  }
  return(calc_t)
}
EA_tz = sample_var( var_name = "ln_d_gt",
                     n_samples = 500,
                     obj = Obj,
                     mu = Obj$env$last.par.best,
                     prec = Opt$SD$jointPrecision,
                     fun1 = calculate_EA,
                     fun2 = function(x){c(mean(x),sd(x))},
                     a_g = Data$a_g )
EA_tz = EA_tz / 1000

# Northward range edge
calculate_NRE = function( ln_d_gt, a_g, z_g, edge_quantile=0.95 ){
  order_g = order(z_g, decreasing=FALSE)
  calc_t = rep(NA, ncol(ln_d_gt))
  for( t in seq_along(calc_t)){
    prop_z = cumsum(exp(ln_d_gt[order_g,t])) / sum(exp(ln_d_gt[order_g,t]))
    calc_t[t] = z_g[order_g][ which.min((prop_z-edge_quantile)^2) ]
  }
  return(calc_t)
}
NRE_tz = sample_var( var_name = "ln_d_gt",
                     n_samples = 500,
                     obj = Obj,
                     mu = Obj$env$last.par.best,
                     prec = Opt$SD$jointPrecision,
                     fun1 = calculate_NRE,
                     fun2 = function(x){c(mean(x),sd(x))},
                     z_g = Data$z_g,
                     a_g = Data$a_g )


# Plot fit
png( "statistics.png", width=7, height=4, unit='in', res=200 )
  par( mfrow=c(2,2), mar=c(2,3,1,1), mgp=c(2,0.5,0), tck=-0.02, xaxs="i", yaxs="i", oma=c(1,0,0,0) )
  # Time-series
  for( z in 1:4 ){
    if( z==1 ) matrix_tz = B_tz
    if( z==2 ) matrix_tz = COG_tz
    if( z==3 ) matrix_tz = EA_tz
    if( z==4 ) matrix_tz = NRE_tz
    Y = c(matrix_tz[1,]-2*matrix_tz[2,], rev(matrix_tz[1,]+2*matrix_tz[2,]))
    plot( x=year_set, y=matrix_tz[1,], type="l", ylim=range(Y),
          main= c("Population size", "Northward center of gravity", "Equivalent area", "Northward range edge")[z],
          xlab="", ylab=c("Biomass (kg)", "position (km north)", "area (1000 km^2)", "position (km north)" )[z] )
    polygon( x=c(year_set,rev(year_set)),
             y=Y, col=rgb(0,0,0,alpha=0.2), border=NA )
  }
  mtext( side=1, text="Year", outer=TRUE )
dev.off()


####################
# Calculate local trends
####################

#
library(terra)
years = c(2010, 2019)

#
Y_gt = Report$ln_d_gt[,match(years,year_set)]
n_cells = nrow(Y_gt)

# Convert to Raster_stack
r_list = list()
for( tI in seq_len(ncol(Y_gt)) ){
  Points = st_sf(grid, "ln_d"=Y_gt[,tI], crs=st_crs(grid) )

  # Define new grid
  Terra_layer = rast( Points,
                      nrows = floor(sqrt(n_cells)),
                      ncols = floor(sqrt(n_cells)) )
  # Rasterize
  r_list[[tI]] = rasterize( x = Points,
                            y = Terra_layer,
                            field = "ln_d",
                            fun = mean )
}
# Combine into raster-brick
names(r_list) = years
Brick = rast(r_list)

# Calculate negative-gradient
# Aspect points downhill (negative of gradient) and is measured in radians clockwise from North
Raster_mean = mean( Brick )
negative_grad = terrain( Raster_mean,
                        v = c('slope', 'aspect'),
                        neighbors = 8,
                        unit = "radian" )

# Calculate trend
calc_trend = function(vec, times=seq_along(vec)){
  if( sum(!is.na(vec))<2 ){
    return(NA)
  }else{
    tvec = times - mean(times)
    Lm = lm( vec ~ 1 + tvec )
    return(coef(Lm)['tvec'])
  }
}
trend_vec = apply( values(Brick),
                   MARGIN = 1,
                   FUN = calc_trend,
                   times = years )
trend = Terra_layer
values(trend) = trend_vec

# Calculate climate velocity
# Shift = Positive trend * negative-gradient
vel = rast( list( "vocc_mag" = trend[[1]] / negative_grad[["slope"]],
                           "vocc_angle" = negative_grad[["aspect"]]) )

# flip vector for negative magnitude resulting from negative local trend
values(vel[["vocc_angle"]]) = ifelse( values(vel[["vocc_mag"]])<0,
                                     (values(vel[["vocc_angle"]])+pi) %% (2*pi),
                                     values(vel[["vocc_angle"]]) )
vel[["vocc_mag"]] = log(abs(vel[["vocc_mag"]]))

# Plot
png( file = "local_trends.png", width = 6, height = 6, res = 200, units = "in" )
  vectorplot( raster::stack(vel),
              isField = TRUE,
              scaleSlope = FALSE,
              col = "green",
              aspX = 10,
              axpY = 10,
              #narrows = 100,
              key.arrow = list(label='log(km/year)'),
              unit = "radian",
              add=TRUE )
dev.off()



