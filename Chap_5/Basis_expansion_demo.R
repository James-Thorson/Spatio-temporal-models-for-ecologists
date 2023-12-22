
setwd("C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_5")
source("../Shared_functions/add_legend.R")

#
library(rnaturalearth)
library(elevatr)
library(sf)
library(viridisLite)
library(stars)
library(splines2)
library(mvtnorm)

# Download data using sf format ... seems slow
Map = ne_countries( country="China", returnclass = "sf" )
Map = st_geometry(Map)  # Reduce to just map
Grid = st_make_grid( Map, cellsize=1 )
Grid = st_intersection( Map, Grid )
Centroid = st_centroid(Grid)
Centroid = get_elev_point( Centroid, src = c("aws","epqs")[1] )
Grid = st_sf( Grid, "elevation" = Centroid$elevation )

#
plot_response = function( formula, name, f=identity, seed=101, sd=1 ){
  set.seed(seed)
  xpred = seq( min(Grid$elevation,na.rm=TRUE), max(Grid$elevation,na.rm=TRUE),length=1000) / 1000
  X_ij = model.matrix( formula, data=data.frame(elevation=xpred) )
  beta_zj = rmvnorm( n=10, mean=rep(0,ncol(X_ij)), sigma=sd^2*diag(ncol(X_ij)) )
  ypred_iz = X_ij %*% t(f(beta_zj))
  png( name, width=6, height=6, res=200, units='in')
    par( mfrow=c(2,1), mar=c(3,3,3,1), mgp=c(2,0.5,0) )
    matplot( x=xpred, y=X_ij, type="l", lwd=2, lty="solid", main="Basis functions", xlab="Elevation" )
    matplot( x=xpred, y=ypred_iz, type="l", lwd=2, col=viridis(nrow(beta_zj)), lty="solid", main="Potential response functions", xlab="Elevation" )
  dev.off()
}
plot_response( formula=~ 1 + elevation, name="Basis-linear.png" )
plot_response( formula=~ 1 + elevation + I(0.2*elevation^2), name="Basis-polynomial.png" )
plot_response( formula=~ 1 + bSpline(elevation,10), name="Basis-bs.png" )
plot_response( formula=~ 1 + iSpline(elevation,10), name="Basis-is.png" )
plot_response( formula=~ 1 + cSpline(elevation,10), f=exp, name="Basis-cs.png" )

# Map basis functions
png( file="Basis-Map.png", width=7, height=3, res=200, units="in" )
  X_ij = model.matrix( ~ 1 + bSpline(elevation,4,intercept=FALSE), data=as.data.frame(Grid) )
  par(mfrow=c(2,3), oma=c(0,0,0,0), mar=c(0,0,2,0), mgp=c(2,0.5,0))
  plot( Grid[,'elevation'], main="Elevation", reset=FALSE, key.pos=NULL, pch=20, cex=0.2, border=NA )
  add_legend( round(range(Grid[['elevation']]),2), legend_y=c(0.05,0.45), legend_x=c(1,1.05) )
  for( i in 1:ncol(X_ij) ){
    Grid$value = X_ij[,i]
    plot(Grid[,'value'], main=ifelse(i==1,"Intercept",paste0("Basis_",i-1)), reset=FALSE, key.pos=NULL, pch=20, cex=0.2, border=NA )
    add_legend( round(range(Grid[['value']]),2), legend_y=c(0.05,0.45), legend_x=c(1,1.05) )
  }
dev.off()

# Illustrate 2D basis
S_ix = bSpline( st_coordinates(st_centroid(Grid))[,1], 4, intercept=TRUE )
S_iy = bSpline( st_coordinates(st_centroid(Grid))[,2], 4, intercept=TRUE )
S_iz = t( sapply(1:nrow(S_ix), FUN=function(i){outer(S_ix[i,],S_iy[i,])}) )
samples = st_as_sf( st_geometry(Grid), S_iz )
#samples = st_as_sf( data.frame(st_coordinates(st_centroid(Grid)), S_iz), coords=c("X","Y") )

png( file="Basis-2D.png", width=7, height=7, res=200, units="in" )
  par( mfrow=c(5,5), oma=c(0,2,2,0), mar=c(1,1,1,1), mgp=c(2,0.5,0) )
  for( rowI in 1:5 ){
  for( colI in 1:5 ){
    if( rowI==1 & colI==1 ){
      plot.new()
    }else if( rowI==1 ){
      xpred = seq( min(st_coordinates(st_centroid(Grid))[,1]), max(st_coordinates(st_centroid(Grid))[,1]),length=1000)
      ypred = predict( S_ix, new=xpred )
      plot( x=xpred, y=ypred[,colI-1] )
    }else if( colI==1 ){
      xpred = seq( min(st_coordinates(st_centroid(Grid))[,2]), max(st_coordinates(st_centroid(Grid))[,2]),length=1000)
      ypred = predict( S_iy, new=xpred )
      plot( y=xpred, x=ypred[,rowI-1] )
    }else{
      zindex = (colI-1) + 4*(rowI-2)
      samples = st_sf( st_geometry(Grid), "spline"=S_iz[,zindex] )
      plot(samples, reset=FALSE, key.pos=NULL, main="", pch=20, cex=0.2, border=NA)
      add_legend( round(range(st_drop_geometry(samples)),2), legend_y=c(0.05,0.45), legend_x=c(1,1.05) )
      #plot(samples, reset=FALSE, key.pos=NULL, main="", pch=20, cex=0.2, border=NA, breaks=seq(0,0.5,length=11) )
      #if(rowI==5 & colI==5) add_legend( seq(0,0.5,length=3), legend_y=c(0.05,0.45), legend_x=c(1,1.05) )
    }
  }}
  mtext( side=c(2,3), c("Latitude spline", "Longitude spline"), outer=TRUE, line=0.5 )
dev.off()

