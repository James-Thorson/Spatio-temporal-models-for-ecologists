
# Data processed from: C:\Users\James.Thorson\Desktop\Work files\AFSC\2022-06 -- Breeding Bird survey
setwd("C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_5")
source("../Shared_functions/add_legend.R")

library(sf)

#
library(rnaturalearth)
sf_states = ne_states( c("United States of America","Canada"), return="sf")
sf_states = sf_states[pmatch(c("Brit","Alas","Yukon"), sf_states$name_en),]
borders = st_polygon(list(matrix(c(-180,0,-180,80,-50,80,-50,0,-180,0),byrow=TRUE,ncol=2)))
sf_states <- st_intersection( st_sfc(borders, crs="+proj=longlat +datum=WGS84"), sf_states )
sf_states = st_union(sf_states)

#
samples = st_read( "samples_3520.csv", options=c("X_POSSIBLE_NAMES=X","Y_POSSIBLE_NAMES=Y") )
st_crs(samples) = "+proj=longlat +datum=WGS84"
samples$Year = as.numeric(samples$Year)
samples$SpeciesTotal = as.numeric(samples$SpeciesTotal)
samples = st_intersection( samples, sf_states )
samples = subset( samples, Year %in% 2019 )   #  &  Year==2019

###########
# Fit a GAM
###########

########## START IN-LINE CODE
library(mgcv)
Gam = gam( SpeciesTotal ~ 1 + s(X,Y,bs="gp",m=1), data=samples, family=poisson )
########## END IN-LINE CODE
plot(Gam)

cellsize = 1
sf_fullgrid = st_make_grid( sf_states, cellsize=cellsize)
sf_grid = st_make_valid(st_intersection( sf_fullgrid, sf_states ))
plotgrid = data.frame(st_coordinates(st_centroid(sf_grid)))
plotdata = st_sf( sf_grid, log_density=log(predict(Gam,type="response",newdata=plotgrid)) )
plotdata$log_density_trimmed = ifelse( plotdata$log_density<(max(plotdata$log_density)-log(1000)), NA, plotdata$log_density )

#plotdata = cbind( plotgrid, prediction=predict(Gam,type="response",newdata=plotgrid) )
#library(ggplot2)
#ggplot(plotdata) + stat_summary_2d(aes(z=log(prediction), y=Y, x=X), bins=20) + ggtitle("log(density)") + xlim(c(-180,-120))
#ggsave( filename="Mapped_GAM.png", width=4, height=3 )

png( "Mapped_GAM.png", width=6, height=3, res=200, units="in")
  par( mfrow=c(1,2) )
  # Full
  plot( plotdata[,'log_density'], key.pos=NULL, reset=FALSE, border=NA, main="Full results" )
  plot( sf_states, add=TRUE )
  add_legend( round(range(plotdata[["log_density"]],na.rm=TRUE),2), legend_y=c(0.55,0.95), legend_x=c(0.95,1), cex=0.7 )
  # Trimmed
  plot( plotdata[,'log_density_trimmed'], key.pos=NULL, reset=FALSE, border=NA, main="Trimmed results" )
  plot( sf_states, add=TRUE )
  add_legend( round(range(plotdata[["log_density_trimmed"]],na.rm=TRUE),2), legend_y=c(0.55,0.95), legend_x=c(0.95,1), cex=0.7 )
dev.off()


###############
# Equally spaced grid
###############

library(TMB)

# Bin into grids
grid_i = as.integer(st_intersects( samples, sf_fullgrid ))
n_x = 1 + round(diff(range(plotgrid$X)) / cellsize)
n_y = 1 + round(diff(range(plotgrid$Y)) / cellsize)
xy_i = as.matrix(expand.grid(1:n_x,1:n_y))[grid_i,]

# Compile
Params = list( "beta0"=0, "ln_sigma2"=0, "logit_rho"=0, "epsilon_xy"=array(0,dim=c(n_x,n_y)) )
# dyn.unload( dynlib("autoregressive_grid") )
compile( "autoregressive_grid.cpp" )
dyn.load( dynlib("autoregressive_grid") )

######## Version 1 -- Conditional in dimension-Y, joint in dimension-X
# Build object
Data = list("Options_vec"=1, "c_i"=samples$SpeciesTotal, "xy_i"=xy_i-1 )
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_xy", DLL="autoregressive_grid" )
# Optimize
Opt1 = nlminb( start=Obj$par, obj=Obj$fn, grad=Obj$gr )
Opt1$SD = sdreport( Obj )
report1 = Obj$report()

png( "Hessian.png", width=4, height=4, res=200, units="in")
  Matrix::image( Obj$env$spHess(random=TRUE), main="Sparsity pattern" )
dev.off()

######## Version 2 -- Kroenekcer product of precision in both dimensions
# Build object
Data$Options_vec = 2
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_xy", DLL="autoregressive_grid" )
# Optimize
Opt2 = nlminb( start=Obj$par, obj=Obj$fn, grad=Obj$gr )
Opt2$SD = sdreport( Obj )
report2 = Obj$report()

# Plot basis functions
E = eigen(report2$Q_zz)
V = E$vectors * rep(1,nrow(E$vectors))%o%(E$values^-0.5)
V = V * rep(1,nrow(E$vectors))%o%sign(colSums(V))
plotgrid = st_sf(sf_fullgrid, v=V[,(n_x*n_y):1], crs=st_crs(sf_fullgrid) )
plotgrid = st_intersection( plotgrid, sf_states )
png( "Mapped_AR_basis_functions.png", width=8, height=8, res=200, units="in")
  plot( plotgrid, breaks=seq(-1*max(abs(V)),1*max(abs(V)),length=11), max.plot=16, border=NA ) # , key.pos=NULL, reset=FALSE, main="Conditional in axis-Y" )
dev.off()

######## Version 3 -- Built-in function for AR process
# Build object
Data$Options_vec = 3
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_xy", DLL="autoregressive_grid" )
# Optimize
Opt3 = nlminb( start=Obj$par, obj=Obj$fn, grad=Obj$gr )
Opt3$SD = sdreport( Obj )
report3 = Obj$report()

# Compare results
Mat = cbind( "Conditional"=Opt1$par, "Simultaneous"=Opt2$par, "Density_package"=Opt3$par )
Mat = rbind(Mat, "Time_seconds"=c(Opt1$time_for_run, Opt2$time_for_run, Opt3$time_for_run) )
write.csv( round(Mat,3), file=paste0("Comparison.csv"))

# Plot
plotgrid = st_sf( sf_fullgrid, D1=log(as.vector(report1$D_xy)), D2=log(as.vector(report2$D_xy)),
          D3=log(as.vector(report3$D_xy)), crs=st_crs(sf_grid) )
plotgrid = st_intersection( plotgrid, sf_states )

png( "Mapped_densities.png", width=6, height=2, res=200, units="in")
  par( mfrow=c(1,3), mar=c(1,1,1,1) )
  Range = range( cbind(plotgrid[['D1']], plotgrid[['D2']], plotgrid[['D3']]) )
  plot( plotgrid['D1'], key.pos=NULL, reset=FALSE, main="Conditional in axis-Y", border=NA, breaks=seq(Range[1],Range[2],length=11) )
  #add_legend( round(range(plotgrid[['D1']]),2), legend_y=c(0.55,0.95), legend_x=c(0.95,1), cex=0.7 )
  plot( plotgrid['D2'], key.pos=NULL, reset=FALSE, main="Simultaneous", border=NA, breaks=seq(Range[1],Range[2],length=11) )
  #add_legend( round(range(plotgrid[['D2']]),2), legend_y=c(0.55,0.95), legend_x=c(0.95,1), cex=0.7 )
  plot( plotgrid['D3'], key.pos=NULL, reset=FALSE, main="Using density package", border=NA, breaks=seq(Range[1],Range[2],length=11) )
  add_legend( round(Range,2), legend_y=c(0.55,0.95), legend_x=c(0.95,1), cex=0.7 )
dev.off()

#####################
# Demonstrate CAR
#####################

# make grid
cellsize = 1
sf_fullgrid = st_make_grid( sf_states, cellsize=cellsize)
sf_grid = st_make_valid(st_intersection( sf_fullgrid, sf_states ))

# Get adjacency
grid_i = as.integer(st_intersects(samples, sf_grid))
st_rook = function(a, b = a, ...) st_relate(a, b, pattern = "F***1****", ... )
grid_A = st_rook( sf_grid, sparse=TRUE )
A_ss = as(grid_A, "sparseMatrix")
A_ss = as(A_ss, "TsparseMatrix")
I_ss = Matrix::sparseMatrix( i=1:length(sf_grid), j=1:length(sf_grid), x=rep(1,length(sf_grid)) )
I_ss = as(I_ss, "TsparseMatrix")

# Calculate min and max eigenvalues, range( 1 / eigen(as.matrix(grid_A))$values )
# Using igraph for sparse-matrix calculations
graph = igraph::graph_from_adjacency_matrix( A_ss, weighted=TRUE )
f2 = function(x, extra=NULL) { as.vector(x %*% A_ss) }
rho_min = Re(igraph::arpack(f2, sym=FALSE, options=list(n=nrow(A_ss), nev=3, ncv=8, which="SR"))$values[1]) 
rho_max = Re(igraph::arpack(f2, sym=FALSE, options=list(n=nrow(A_ss), nev=3, ncv=8, which="LR"))$values[1]) 

#
compile( "CAR.cpp" ) # framework='TMBad'
dyn.load( dynlib("CAR") )

# Compile
Params = list( "beta0" = 0, 
               "rho_prime" = 0, 
               "ln_sigma" = 0, 
               "omega_s" = rep(0,length(sf_grid)) )
Data = list( "c_i" = samples$SpeciesTotal, 
             "rho_bounds" = c( rho_min, rho_max ), 
             "s_i" = grid_i-1, 
             "I_ss" = I_ss, 
             "A_ss" = A_ss )

# Build and optimize
obj = MakeADFun( data=Data, parameters=Params, random="omega_s", DLL="CAR" )
opt = nlminb( start=obj$par, obj=obj$fn, grad=obj$gr )
opt$SD = sdreport( obj )
report = obj$report()

# Plot basis functions
E = eigen(report$Q_ss)
V = E$vectors * rep(1,nrow(E$vectors))%o%(E$values^-0.5)
V = V * rep(1,nrow(E$vectors))%o%sign(colSums(V))
V = V[,length(sf_grid):1][,1:16]
plotgrid = st_sf(sf_grid, v=V, crs=st_crs(sf_grid) )
png( "Mapped_CAR_basis_functions.png", width=8, height=8, res=200, units="in")
  plot( plotgrid, breaks=seq(-1*max(abs(V)),1*max(abs(V)),length=11), max.plot=16, border=NA ) # , key.pos=NULL, reset=FALSE, main="Conditional in axis-Y" )
dev.off()

plotgrid = st_sf(sf_grid, ln_N=report$logN_s, crs=st_crs(sf_grid) )
png( "Mapped_densities--CAR.png", width=3, height=3, res=200, units="in")
  par( mfrow=c(1,1), mar=c(1,1,1,1) )
  plot( plotgrid, border=NA, key.pos=NULL, main="CAR method" )
  add_legend( round(range(plotgrid[['ln_N']]),2), legend_y=c(0.55,0.95), legend_x=c(0.95,1), cex=0.7 )
dev.off()
