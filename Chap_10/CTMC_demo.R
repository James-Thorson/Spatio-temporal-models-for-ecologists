
# Data processed from: C:\Users\James.Thorson\Desktop\Work files\AFSC\2022-06 -- Breeding Bird survey
setwd("C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_10")
source( "../Shared_functions/matexp.R" )
source( "../Shared_functions/add_legend.R" )

# load shapefile
library(sf)
library(Matrix)
sf_area = rnaturalearth::ne_countries( scale=10, country="Madagascar", return="sf")

# extract main island
sf_area = st_cast(st_geometry(sf_area),"POLYGON")[[1]]
sf_area = st_sfc( sf_area, crs="+proj=longlat +datum=WGS84" )
plot(sf_area)

# START HERE FOR BOOK SNIPPET
# make grid and exclude small boundary cells
cellsize = 0.8
sf_fullgrid = st_make_grid( sf_area, cellsize=cellsize )
sf_grid = st_make_valid(st_intersection( sf_fullgrid, sf_area ))
sf_grid = sf_grid[ st_area(sf_grid)>(0.2*max(st_area(sf_grid))) ]

# Download elevation and truncate negative values
grid = elevatr::get_elev_point( locations=st_point_on_surface(sf_grid), src = "aws" )
grid$elevation = ifelse( grid$elevation<1, 1, grid$elevation)

# Get adjacency
st_rook = function(m, ...) st_relate(m, m, pattern="F***1****", ... )
grid_A = st_rook( sf_grid, sparse=TRUE )
A = as(grid_A,"sparseMatrix")

# Diffusion rate
diffusion_coefficient = 1 / cellsize^2
D = diffusion_coefficient * A
diag(D) = -1 * colSums(D)

# Movement matrix at for different steps sizes
M = NULL
for( p in 1:3 ){
  M[[p]] = matexp( D, log2steps=c(2,3,Inf)[p] )
}
# END HERE FOR BOOK SNIPPET

# Visualize Euler approximation for diffusion-only case
Image = NULL
png( file=paste0("euler_cellsize",gsub(pattern=".",replacement="",x=as.character(cellsize), fixed=TRUE),".png"), width=6, height=6, units="in", res=200 )
  library(gridExtra)
  Image[[1]] = image( D, colorkey=TRUE, main="Diffusion rate", sub="", xlab="", ylab="" )
  for(p in 1:3 ){
    main = ifelse( p==3, "Matrix exponential", paste0("Euler: nstep=",2^c(2,3,Inf)[p]) )
    Image[[p+1]] = image( M[[p]], colorkey=TRUE, main=main, sub="", xlab="", ylab="", lwd=0, col.regions=viridisLite::viridis(20) )
  }
  grid.arrange(Image[[1]], Image[[2]], Image[[3]], Image[[4]], ncol=2)
dev.off()

# Taxis rate
grid$preference = 20 / cellsize * dnorm( log(grid$elevation), mean=log(100), sd=3 )
Z = A * outer( grid$preference, grid$preference, FUN="-")
diag(Z) = -1 * colSums(Z)

# Movement matrix
Mrate = D + Z
range(Mrate - diag(diag(Mrate)))
M = matexp(Mrate, log2steps=Inf )

# Project abundance
N_gt = matrix( 0, ncol=9, nrow=nrow(grid), dimnames=list(NULL,paste0("Time ",1:9)) )
N_gt[as.numeric(st_within(st_point(c(48.5,-20)),sf_grid)),1] = 1000
for( t in 2:9 ) N_gt[,t] = as.numeric(M %*% N_gt[,t-1])
stationary = abs(eigen(M)$vectors[,1])
stationary = stationary / mean(stationary)

# Plot
png( paste0("Demo_cellsize",gsub(pattern=".",replacement="",x=as.character(cellsize), fixed=TRUE),".png"), width=6, height=6, res=200, units="in" )
  plotgrid = st_sf( sf_grid, data.frame(cbind("log_elevation"=grid$elevation, "log_preference"=grid$preference, "log_stationary"=stationary, N_gt)) )
  par( mfrow=c(3,4), mar=c(0,0,0,0), mgp=c(2,0.5,0), tck=-0.02 )
  for( p in 1:12 ){
    if(p<3){f=log}else{f=identity}
    plot( plotgrid[,p], key.pos=NULL, reset=FALSE, logz=ifelse(p<3,TRUE,FALSE), border=NA )
    add_legend( round(range(f(plotgrid[[p]])),2) )
  }
dev.off()
