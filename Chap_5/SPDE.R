
# Data processed from: C:\Users\James.Thorson\Desktop\Work files\AFSC\2022-06 -- Breeding Bird survey
setwd("C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_5")
source("../Shared_functions/add_legend.R")

library(sf)
library(TMB)
library(rnaturalearth)
library(fmesher)

#
out = st_read( "samples_3520.csv", options=c("X_POSSIBLE_NAMES=X","Y_POSSIBLE_NAMES=Y") )
st_crs(out) = "+proj=longlat +datum=WGS84"
out$Year = as.numeric(out$Year)
out$SpeciesTotal = as.numeric(out$SpeciesTotal)
out = subset( out, State %in% c("BritCol","Alaska","Yukon") )   #  &  Year==2018
out = subset( out, Year %in% 2019 )   #  &  Year==2018

#
sf_states = ne_states( c("United States of America","Canada"), return="sf")
sf_states = sf_states[pmatch(c("Brit","Alas","Yukon"), sf_states$name_en),]
borders = st_polygon(list(matrix(c(-180,0,-180,80,-50,80,-50,0,-180,0),byrow=TRUE,ncol=2)))
sf_states <- st_intersection( st_sfc(borders, crs="+proj=longlat +datum=WGS84"), sf_states )
sf_states = st_union(sf_states)

############################
# Show SPDE approximation
############################

# Simulate locations
loc_xy = expand.grid( "x"=1:10, "y"=1:10 )
#loc_xy = cbind( "x"=runif(10), "y"=runif(10))

# create mesh
mesh = fm_mesh_2d( loc_xy, plot.delay=NULL, refine=TRUE )

# Create matrices in fmesher / INLA
spde <- fm_fem(mesh, order=2)

png(file="SPDE_explanation_pt1.png", width=8, height=4, res=200, units="in")
  par( mfrow=c(1,2), mar=c(0,0,2,0), mgp=c(2,0.5,0), tck=-0.02)
  # Plot samples
  plot( loc_xy, xlim=range(mesh$loc[,1]), ylim=range(mesh$loc[,2]), main="Sample locations")
  # Plot mesh
  plot(mesh, main="Mesh composed of triangles")
  text( x=mesh$loc[,1], y=mesh$loc[,2], labels=1:mesh$n, col=ifelse(1:mesh$n%in%mesh$idx$loc,"blue","black"))
  title("Mesh composed of triangles")
dev.off()

png(file="SPDE_explanation_pt2.png", width=9, height=3, res=200, units="in")
  par( mfrow=c(1,3), mar=c(2,2,2,0), mgp=c(2,0.5,0), tck=-0.02, xaxs="i", yaxs='i')
  # Visualize SPDE approx.
  Col = colorRampPalette(colors=c("blue","white","red"))
  Points = function(X,Y,Z){
    DF = cbind( expand.grid('X'=as.vector(X), 'Y'=as.vector(Y)), 'Z'=as.vector(Z) )
    DF = DF[which(DF[,'Z']!=0),]
    points( x=DF[,'X'], y=DF[,'Y'], pch=20)
  }
  image( z=as.matrix(spde$c0), x=1:mesh$n, y=1:mesh$n, main="M0", zlim=c(-1,1)*max(abs(spde$c0)), col=Col(11), xlab="", ylab=""); box()
  Points( X=1:mesh$n, Y=1:mesh$n, Z=as.matrix(spde$c0) )
  image(z=as.matrix(spde$g1), x=1:mesh$n, y=1:mesh$n, main="M1", zlim=c(-1,1)*max(abs(spde$g1)), col=Col(11), xlab="", ylab=""); box()
  Points( X=1:mesh$n, Y=1:mesh$n, Z=as.matrix(spde$g1) )
  image(z=as.matrix(spde$g2), x=1:mesh$n, y=1:mesh$n, main="M2", zlim=c(-1,1)*max(abs(spde$g2)), col=Col(11), xlab="", ylab=""); box()
  Points( X=1:mesh$n, Y=1:mesh$n, Z=as.matrix(spde$g2) )
dev.off()

png(file="SPDE_explanation_pt2.png", width=9, height=3, res=200, units="in")
  par( mfrow=c(1,3), mar=c(2,2,2,0), mgp=c(2,0.5,0), tck=-0.02, xaxs="i", yaxs='i')
  # Visualize SPDE approx.
  Col = colorRampPalette(colors=c("blue","white","red"))
  Points = function(X,Y,Z){
    DF = cbind( expand.grid('X'=as.vector(X), 'Y'=as.vector(Y)), 'Z'=as.vector(Z) )
    DF = DF[which(DF[,'Z']!=0),]
    points( x=DF[,'X'], y=DF[,'Y'], pch=20)
  }
  image(z=as.matrix(spde$c0), x=1:mesh$n, y=1:mesh$n, main="M0", zlim=c(-1,1)*max(abs(spde$c0)), col=Col(11), xlab="", ylab=""); box()
  Points( X=1:mesh$n, Y=1:mesh$n, Z=as.matrix(spde$c0))
  image(z=as.matrix(spde$g1), x=1:mesh$n, y=1:mesh$n, main="M1", zlim=c(-1,1)*max(abs(spde$g1)), col=Col(11), xlab="", ylab=""); box()
  Points( X=1:mesh$n, Y=1:mesh$n, Z=as.matrix(spde$g1))
  image(z=as.matrix(spde$g2), x=1:mesh$n, y=1:mesh$n, main="M2", zlim=c(-1,1)*max(abs(spde$g2)), col=Col(11), xlab="", ylab=""); box()
  Points( X=1:mesh$n, Y=1:mesh$n, Z=as.matrix(spde$g2))
dev.off()

#
ln_kappa = -3
Q_spde = exp(4*ln_kappa)*spde$c0 + 2*exp(2*ln_kappa)*spde$g1 + spde$g2
E = eigen(Q_spde)
V = E$vectors * rep(1,nrow(E$vectors))%o%(E$values^-0.5)
V = V * rep(1,nrow(E$vectors))%o%sign(colSums(V))

###################
# Show Matern correlation function
###################

Matern_Correlation = function( distance_vec, Nu, Scale ){
  Correlation = 2^(1-Nu) * gamma(Nu)^(-1) * (sqrt(2*Nu)*distance_vec/Scale)^Nu * besselK(sqrt(2*Nu) * distance_vec/Scale, nu=Nu)
  return( Correlation )
}

# Distance with 10% correlation is approximately Scale/2
Distance = seq(0,3, length=1e4)
Nu = c(0.5, 1, 2, 10, 100)
Corr_Matrix = sapply( Nu, FUN=Matern_Correlation, distance_vec=Distance, Scale=1)

# Plot
png(file="Matern_explanation.png", width=4, height=4, res=200, units="in")
  par( mar=c(2,2,1,0), mgp=c(2,0.5,0), tck=-0.02, xaxs="i", yaxs='i')
  matplot( y=Corr_Matrix, x=Distance, type="l", lwd=3, ylim=c(0,1.2), lty="solid", ylab="Correlation", xlab="Distance" )
  abline( h=0.1, lty="dotted" )
  legend( "topright", legend=Nu, fill=1:6, bty="n", title=expression(nu))
  title( "Matern correlation function (Scale=1)" )
dev.off()

###################
# Unequal distance 2D autoregressive
###################

########## START IN-LINE CODE
# create mesh
xy_i = st_coordinates(out)
mesh = fm_mesh_2d( xy_i, refine=TRUE, cutoff=0.5)

# Create matrices in fmesher / INLA
spde <- fm_fem(mesh, order=2)

# create projection matrix from vertices to samples
A_is = fm_evaluator( mesh, loc=xy_i )$proj$A

# Create extrapolation grid
cellsize = 1
grid = st_make_grid( sf_states, cellsize=cellsize )

# create projection matrix from vertices to grid
A_gs = fm_evaluator( mesh, loc=st_coordinates(st_centroid(grid)) )$proj$A
########## END IN-LINE CODE

# plot(xy_i); points(st_coordinates(grid)[,1:2],col="red")

png(file="SPDE_mesh.png", width=8, height=4, res=200, units="in")
  par( mfrow=c(1,2), mar=c(0,0,2,0), mgp=c(2,0.5,0), tck=-0.02)
  # Plot samples
  #plot( xy_i, xlim=range(mesh$loc[,1]), ylim=range(mesh$loc[,2]), main="Sample locations")
  plot(sf_states, xlim=c(-180,-110), main="Sample locations")
  plot(out,add=TRUE, col="blue")
  # Plot mesh
  plot(sf_states, xlim=c(-180,-110))
  plot(mesh, main="Finite element analysis mesh",add=TRUE)
  points( x=mesh$loc[,1], y=mesh$loc[,2], col=ifelse(1:mesh$n%in%mesh$idx$loc,"blue",NA) ) # , labels=1:mesh$n
  title("Mesh composed of triangles")
dev.off()

library(gridExtra)
png(file="SPDE_matrices.png", width=7, height=3, res=200, units="in")
  m0 = image(spde$c0, colorkey=FALSE, sub=paste0("Nonzero: ",prod(dim(spde$c0))-sum(spde$c0==0)), main=expression(M[0]) ) #, useRaster=TRUE)
  m1 = image(spde$g1, colorkey=FALSE, sub=paste0("Nonzero: ",prod(dim(spde$g1))-sum(spde$g1==0)), main=expression(M[1])) #, useRaster=TRUE)
  m2 = image(spde$g2, colorkey=FALSE, sub=paste0("Nonzero: ",prod(dim(spde$g2))-sum(spde$g2==0)), main=expression(M[2])) #, useRaster=TRUE)
  grid.arrange(m0, m1, m2, ncol=3)
dev.off()

# COmpile
compile( "SPDE.cpp" )
dyn.load( dynlib("SPDE") )

######## ######## SPDE-based
# Build object
Data = list( "c_i"=out$SpeciesTotal, "A_is"=A_is, "A_gs"=A_gs,
             "M0"=spde$c0, "M1"=spde$g1, "M2"=spde$g2 )
Params = list( "beta0"=0, "ln_tau"=0, "ln_kappa"=0,
               "omega_s"=rnorm(nrow(spde$c0)) )
Obj = MakeADFun( data=Data, parameters=Params, random="omega_s" )

# Optimize
Opt = nlminb( start=Obj$par, obj=Obj$fn, grad=Obj$gr )
Opt$SD = sdreport( Obj, bias.correct=TRUE )
report = Obj$report()

# Plot basis functions
E = eigen(report$Q)
V = E$vectors * rep(1,nrow(E$vectors))%o%(E$values^-0.5)
V = V * rep(1,nrow(E$vectors))%o%sign(colSums(V))
V = as.matrix(A_gs %*% V)[,ncol(V):1]
grid_sf = st_sf(grid, v=V, crs=st_crs(grid) )
plot_grid = st_intersection( grid_sf, sf_states )
png( "Mapped_SPDE_basis_functions.png", width=8, height=8, res=200, units="in")
  plot( plot_grid, breaks=seq(-1*max(abs(V)),1*max(abs(V)),length=11), max.plot=16, border=NA ) # , key.pos=NULL, reset=FALSE, main="Conditional in axis-Y" )
dev.off()

# Plot
grid_sf = st_sf(grid, log_density=report$logN_g )
plot_grid = st_intersection( grid_sf, sf_states )

png( "Mapped_densities-SPDE.png", width=3, height=3, res=200, units="in")
  par( mfrow=c(1,1), mar=c(1,1,1,1) )
  plot( plot_grid['log_density'], key.pos=NULL, main="SPDE method", border=NA )
  add_legend( round(range(plot_grid[["log_density"]]),2), legend_y=c(0.55,0.95), legend_x=c(0.95,1), cex=0.7 )
dev.off()


