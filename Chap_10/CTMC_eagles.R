
# Data processed from: C:\Users\James.Thorson\Desktop\Work files\AFSC\2022-06 -- Breeding Bird survey
setwd("C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_10")
source( "../Shared_functions/add_legend.R" )

#
library(TMB)
library(sf)
library(Matrix)
library(raster)
library(elevatr)
library(rnaturalearth)
library(ggplot2)
library(gridExtra)

# Get spatial domain
sf_states = ne_states( c("United States of America","Canada"), return="sf")
sf_states = sf_states[pmatch(c("Brit","Alas","Yukon","Washington","Ore","Calif"), sf_states$name_en),]
p = st_polygon(list(matrix(c(-180,0,-180,80,-50,80,-50,0,-180,0),byrow=TRUE,ncol=2)))
sf_states <- st_intersection( st_sfc(p, crs="+proj=longlat +datum=WGS84"), sf_states )
sf_states = st_union(sf_states)
sf_states = st_transform( sf_states, crs=st_crs("+init=epsg:3338 +units=km") )

# Get coastline
sf_coast = ne_coastline( scale=50, return="sf")
sf_coast = st_transform( sf_coast, crs=st_crs("+init=epsg:3338 +units=km") )

# Keep mainland Alaska
sf_states = st_cast(sf_states,"POLYGON")
which_max = which.max(st_area(sf_states))
sf_states = sf_states[which_max]

# simplify geometry somewhat
sf_states = st_simplify( sf_states, dTolerance = 10 )

# Load copNDVI (saved from rasterdiv)
copNDVI = raster( "NDVI.tif" )

#
samples = st_read( paste0("../Chap_5/samples_3520.csv"), options=c("X_POSSIBLE_NAMES=X","Y_POSSIBLE_NAMES=Y") )
st_crs(samples) = "+proj=longlat +datum=WGS84"
samples = st_transform( samples, crs=st_crs("+init=epsg:3338 +units=km") )
samples$Year = as.numeric(samples$Year)
samples$SpeciesTotal = as.numeric(samples$SpeciesTotal)

#
cellsize = 200
sf_fullgrid = st_make_grid( sf_states, cellsize=cellsize)
sf_grid = st_make_valid(st_intersection( sf_fullgrid, sf_states ))
sf_grid = sf_grid[ st_area(sf_grid)>(0.01*max(st_area(sf_grid))) ]

# make data frame of covariates
df_grid = st_centroid(sf_grid)
df_grid = get_elev_point( df_grid, src = "aws" )
df_grid[,c('elevation','elev_units')] = data.frame(df_grid$elevation / 1000, "kilometers")
df_grid$elevation = ifelse( is.na(df_grid$elevation), 0, df_grid$elevation )
df_grid$NDVI = extract( x=copNDVI, y=as(df_grid,"Spatial") ) / 255
dist_to_coast = st_distance( sf_grid, st_cast(sf_coast,"LINESTRING") ) / 1000
df_grid$dist_to_coast = apply( dist_to_coast, MARGIN=1, FUN=min )
df_grid = data.frame(df_grid)

# Bin into grids
samples = st_intersection( samples, sf_grid )
grid_i = as.integer(st_intersects(samples, sf_grid))

# Get adjacency
st_rook = function(a, b = a, ...) st_relate(a, b, pattern = "F***1****", ... )
grid_A = st_rook( sf_grid, sparse=TRUE )
A_ss = as(grid_A, "sparseMatrix")
A_ss = as(A_ss, "TsparseMatrix")
I_ss = Matrix::sparseMatrix( i=1:length(sf_grid), j=1:length(sf_grid), x=rep(1,length(sf_grid)) )
I_ss = as(I_ss, "TsparseMatrix")
At_zz = cbind( attr(A_ss,"i"), attr(A_ss,"j") )
colsumA_s = colSums(A_ss)

#####################
# Demonstrate CAR with expm -- VERSION-2
# 80 km/year natal dispersal (https://doi.org/10.3356/JRR-13-00005.1)
# age-5 maturity
# Mean-squared distance MSD = (80)^2 per 5-years = 6400 km^2 / 5 years
# MSD = 2nDt ->    D = MSD / 2nt = 6400 / (2*2*5) = 320 km^2 / year
#####################

compile( "CAR_expm.cpp" ) # framework='TMBad'
dyn.load( dynlib("CAR_expm") )

# Covariates ... I(.) doesn't work right with scale
preference_formula = ~ 0 + poly(elevation,2,raw=TRUE) + poly(NDVI,2,raw=TRUE) + poly(dist_to_coast,2,raw=TRUE)
X_sz = model.matrix( preference_formula, data=df_grid )

# Assemble inputs
Data = list( #"movement_penalty" = 1000, 
             "CTMC_version" = 1,
             "DeltaD" = cellsize,
             "n_t" = diff(range(samples$Year))+1, 
             "n_s" = nrow(X_sz), 
             "c_i" = samples$SpeciesTotal, 
             "s_i" = grid_i-1, 
             "t_i" = samples$Year-min(samples$Year), 
             "X_sz" = X_sz, 
             "colsumA_s" = colsumA_s, 
             "At_zz" = At_zz,
             "rho_bounds" = 1/range(eigen(as.matrix(grid_A))$values),
             "I_ss" = I_ss,
             "A_ss" = A_ss )
Params = list( # "beta0" = 0, 
               "rho_prime" = 0, 
               "ln_sigmaO" = 0,
               "ln_sigmaB"  =  0,  
               "ln_D" = log( 80^2 / (2*2*5) ),
               "gamma_z" = rep(0,ncol(Data$X_sz)), 
               "beta_t" = rep(0,Data$n_t), 
               "ln_D_st" = matrix(0,nrow=Data$n_s,ncol=Data$n_t) ) # "gamma_z" = rep(0,ncol(X_sz)), 
Map = NULL
  Map$ln_D = factor(NA)

# Build object
obj = MakeADFun( data=Data, parameters=Params, random=c("beta_t","ln_D_st"), map=Map, DLL="CAR_expm" )
opt = nlminb( start=obj$par, obj=obj$fn, gr=obj$gr )
opt$SD = sdreport( obj )
report = obj$report()

#
ln_D_st = obj$env$parList()$ln_D_st
colnames(ln_D_st) = sort(unique(samples$Year))
plotgrid = st_sf( sf_grid, 
                  df_grid[,c("elevation","NDVI","dist_to_coast")],
                  preference = report$h_s,
                  log_D = ln_D_st, 
                  crs = st_crs(sf_grid) )
png( paste0("Eagle_results.png"), width=6, height=8, res=200, units="in" )
  par( mfrow=c(5,2), mar=c(0,0,0,0), mgp=c(2,0.5,0), tck=-0.02 )
  for( b in 1:4 ){
    plot( plotgrid[,b], key.pos=NULL, reset=FALSE, border=NA, main=c("elevation (m)","NDVI","distance to coast (km)","preference")[b] )
    add_legend( round(range(plotgrid[[b]]),2), legend_y=c(0.5,0.9) )
  }
  for( t in 4+seq(1,Data$n_t,len=6) ){
    plot( plotgrid[,t], key.pos=NULL, reset=FALSE, border=NA )
    add_legend( round(range(plotgrid[[t]]),2), legend_y=c(0.5,0.9) )
  }
dev.off()

# 
png( "Eagle_sparsity.png", width=4, height=4, unit='in', res=200 )
  par( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )
  Matrix::image( obj$env$spHess(par=obj$env$last.par.best,random=TRUE) )
dev.off()

# Plot covariates
png( paste0("Eagle_covariate_response.png"), width=6, height=2, res=200, units="in" )
  par( mfrow=c(1,3), mar=c(3,2,1,1), mgp=c(2,0.5,0), tck=-0.02 )
  for( cI in 1:3 ){
    Xval = seq( min(df_grid[[c("elevation","NDVI","dist_to_coast")[cI]]]),max(df_grid[[c("elevation","NDVI","dist_to_coast")[cI]]]),length=1000)
    X_sz = model.matrix( ~ 0 + poly(Xval,2,raw=TRUE), data=data.frame("Xval"=Xval) )
    Yval = X_sz %*% obj$env$parList()$gamma_z[2*(cI-1)+1:2]
    plot( x=Xval, y=Yval, type="l", lwd=2, ylab="Preference", xlab=c("elevation","NDVI","distance to coast")[cI] )
  }
dev.off()

# Use marginaleffects
library(marginaleffects)
library(ggplot2)
source( "../Shared_functions/marginaleffects.R" )
# Bundle materials
fit = list( obj = obj,
            opt = opt,
            data = df_grid,
            formula = preference_formula,
            parhat = obj$env$parList() )
class(fit) = "custom_tmb"
quant = function(x) seq(min(x),max(x),length=21)
# get_coef( fit, param="gamma_z" )
# get_vcov( fit, param="gamma_z" )
# set_coef( fit, param="gamma_z", newpar=rep(0,6) )
# get_predict( fit, param="gamma_z", newpar=rep(0,length(get_coef(fit,param="gamma_z"))), newdata=new_elev )

# Get prediction for partial dependence plots
new_elev = datagrid( newdata=data.frame(df_grid)[,c('elevation','NDVI','dist_to_coast')], elevation=quant, NDVI=mean, dist_to_coast=mean) #, NDVI=quant, dist_to_coast=quant )
  pred_elev = predictions( fit, newdata=new_elev, center=TRUE, param="gamma_z" )
new_NDVI = datagrid( newdata=data.frame(df_grid), elevation=mean, NDVI=quant, dist_to_coast=mean) #, NDVI=quant, dist_to_coast=quant )
  pred_NDVI = predictions( fit, newdata=new_NDVI, center=TRUE, param="gamma_z" )
new_dist = datagrid( newdata=data.frame(df_grid), elevation=mean, NDVI=mean, dist_to_coast=quant) #, NDVI=quant, dist_to_coast=quant )
  pred_dist = predictions( fit, newdata=new_dist, center=TRUE, param="gamma_z" )

# Make plot of marginal effects
p1 <- ggplot( as.data.frame(pred_elev) ) +
  geom_line( aes(y=estimate, x=elevation), color="blue", size=1 ) +
  geom_ribbon( aes( x=elevation, ymin=conf.low, ymax=conf.high), fill=rgb(0,0,1,0.2) )
p2 <- ggplot( as.data.frame(pred_NDVI) ) +
  geom_line( aes(y=estimate, x=NDVI), color="blue", size=1 ) +
  geom_ribbon( aes( x=NDVI, ymin=conf.low, ymax=conf.high), fill=rgb(0,0,1,0.2) )
p3 <- ggplot( as.data.frame(pred_dist) ) +
  geom_line( aes(y=estimate, x=dist_to_coast ), color="blue", size=1 ) +
  geom_ribbon( aes( x=dist_to_coast, ymin=conf.low, ymax=conf.high), fill=rgb(0,0,1,0.2) )
p <- grid.arrange(p1, p2, p3, nrow = 1)
ggsave( paste0("Eagle_covariate_response.png"), p, width=6, height=2 )

