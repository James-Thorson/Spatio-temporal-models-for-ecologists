
library(terra)
library(sf)
library(TMB)
library(Matrix)
library(rnaturalearth)
library(ggplot2)
library(marginaleffects)

setwd("C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_10")
source( "../Shared_functions/add_legend.R" )

# load data
bathy_terra = readRDS("ai_bathy_3km.Rds")
likelihood_terra = readRDS("likelihood_3km.Rds")

# Get land layer
sf_states = ne_states( "united states of america", return="sf")
sf_states = sf_states[pmatch("Alas", sf_states$name_en),]
sf_states = st_transform( sf_states, crs=st_crs("+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") )
sf_states = st_geometry(sf_states)

# change resolution
fact = 1
bathy_terra = aggregate( bathy_terra, fact = fact )
likelihood_terra = as.array( aggregate(likelihood_terra, fact=fact) )

# Convert bathymetry
bathy_sf = st_as_sf(as.polygons(bathy_terra, trunc=FALSE, dissolve=FALSE))
bathy_sf$northings = st_coordinates(st_centroid(bathy_sf))[,'Y']
bathy_sf$eastings = st_coordinates(st_centroid(bathy_sf))[,'X']

# Format likelihood
L_gt = apply( likelihood_terra, MARGIN=3, FUN=function(mat){as.vector(t(mat))} )  
sf_L_gt = st_sf( st_geometry(bathy_sf), L_gt )

# Get adjacency matrix using Raster
A_gg = adjacent( bathy_terra, cells=1:prod(dim(bathy_terra)), pairs=TRUE )
A_gg = Matrix::sparseMatrix( i=A_gg[,1], j=A_gg[,2], x=rep(1,nrow(A_gg)) )
A_gg = as(A_gg, "TsparseMatrix")

# Drop geometry and likelihood with depth <1 m
which_exclude = which( bathy_sf$ai_bathy_fill <= 1 )
sf_L_gt = sf_L_gt[-which_exclude,]
bathy_sf = bathy_sf[-which_exclude,]
A_gg = A_gg[-which_exclude,-which_exclude]
bathy_sf$ai_bathy_fill = bathy_sf$ai_bathy_fill / 1000

# Drop geometry and likelihood on land
#which_exclude = which( as.integer(st_intersects(sf_L_gt, sf_states))==TRUE )
#sf_L_gt = sf_L_gt[-which_exclude,]
#bathy_sf = bathy_sf[-which_exclude,]
#A_gg = A_gg[-which_exclude,-which_exclude]

#
#graph = igraph::graph_from_adjacency_matrix( A_gg )
#which_include = which( igraph::components(graph)$membership == 1 )
#sf_L_gt = sf_L_gt[which_include,]
#bathy_sf = bathy_sf[which_include,]
#A_gg = A_gg[which_include,which_include]

#
#graph = igraph::graph_from_adjacency_matrix( A_gg )
#igraph::is_connected( graph, mode = "weak" )
#igraph::is_connected( graph, mode = "strong" )

# Assemble inputs
colsumA_g = colSums(A_gg)
I_gg = Matrix::sparseMatrix( i=1:nrow(bathy_sf), j=1:nrow(bathy_sf), x=rep(1,nrow(bathy_sf)) )
  I_gg = as(I_gg, "TsparseMatrix")
D_gg = Matrix::sparseMatrix( i=1:nrow(bathy_sf), j=1:nrow(bathy_sf), x=colsumA_g )
  D_gg = as(D_gg, "TsparseMatrix")
At_zz = cbind( attr(A_gg,"i"), attr(A_gg,'j') )

# Make covariates
preference_formula = ~ 0 + poly(ai_bathy_fill, 2)
#preference_formula = ~ 0 + poly(scale_ai_bathy_fill,2)
#preference_formula = ~ 0 + poly(scale(ai_bathy_fill),2)
X_gz = model.matrix( preference_formula, data=bathy_sf )
#X_gz = model.matrix( ~ 0 + poly(scale(ai_bathy_fill),2), data=bathy_sf )

#####################
# Diffusion-only model
#####################

# Compile
compile( "hmm_TMB.cpp", framework='TMBad' )
dyn.load( dynlib("hmm_TMB") )

# Build inputs
Params = list( "ln_D" = 1,
               "gamma_z" = 0.01*rnorm(ncol(X_gz)) )
Data = list( "CTMC_version" = 1,  # 0=diffusion-taxis;  1=logspace diffusion-taxis
             "expm_version" = 0,  # 0=uniformization;  1=series
             "Nmax" = ifelse(mean(res(bathy_terra))==3000, 500, 250),
             "DeltaD" = mean(res(bathy_terra)) / 1000,
             "colsumA_g" = colsumA_g,
             "L_gt" = as.matrix(st_drop_geometry(sf_L_gt)),
             "X_gz" = X_gz,
             "At_zz" = At_zz ) 

# Build and optimize object
obj = MakeADFun( data=Data, parameters=Params )
opt = nlminb( start=obj$par, obj=obj$fn, gr=obj$gr )
opt$SD = sdreport( obj )
report = obj$report()

# Calculate stationary distributino, i.e., v such that report$Mrate_gg %*% v = 0
graph = igraph::graph_from_adjacency_matrix( report$Mrate_gg, weighted=TRUE )
f2 <- function(x, extra=NULL) as.vector(x %*% report$Mrate_gg)
baev <- igraph::arpack(f2, sym=FALSE, options=list(n=nrow(report$Mrate_gg), nev=3, ncv=8, which="LR", maxiter=1e5)) # SR is fast but gives negative eigenvalues;
stationary_g = Re(baev$vectors[,1]) * sign(median(Re(baev$vectors[,1])))

# Plot covariates
Xval = seq( min(bathy_sf$ai_bathy_fill),max(bathy_sf$ai_bathy_fill),length=21)
Xpred_xz = model.matrix( preference_formula, data=data.frame(ai_bathy_fill=Xval) )
matplot( y=Xpred_xz, x=Xval, type="l", lwd=2 )
Yval = Xpred_xz %*% obj$env$parList()$gamma_z
png( paste0("covariate_response-",Data$DeltaD,".png"), width=3, height=3, res=200, units="in" )
  plot( x=Xval, y=Yval, type="l", lwd=2, ylab="Preference", xlab="Depth" )
dev.off()

# Use marginaleffects
source( "../Shared_functions/marginaleffects.R" )
# Bundle materials
fit = list( obj = obj,
            opt = opt,
            data = bathy_sf,
            formula = preference_formula,
            parhat = obj$env$parList() )
class(fit) = "custom_tmb"
quant = function(x) seq(min(x),max(x),length=21)
# get_coef( fit, param="gamma_z" )
# get_vcov( fit, param="gamma_z" )
# set_coef( fit, param="gamma_z", newpar=rep(0,2) )
# get_predict( fit, param="gamma_z", newpar=rep(0,length(get_coef(fit,param="gamma_z"))), newdata=newdata )

# Get prediction for elevation
newdata = datagrid( newdata=data.frame(bathy_sf), ai_bathy_fill=quant )
  pred = predictions( fit, newdata=newdata, param="gamma_z" )
ggplot( pred, aes(ai_bathy_fill, estimate)) +
  geom_line( aes(y=estimate), color="blue", size=1 ) +
  geom_ribbon( aes( x=ai_bathy_fill, ymin=conf.low, ymax=conf.high), fill=rgb(0,0,1,0.2) ) +
  labs(y="Predicted response")
ggsave( paste0("covariate_response_",Data$DeltaD,".png"), width=3, height=3 )

# Plot data likelihood
colnames(sf_L_gt)[-ncol(sf_L_gt)] = as.character(as.Date(seq( as.POSIXlt("2019-02-21"), as.POSIXlt("2019-05-24"), by="days")))
dims = c(3,3)
year_set = round(seq(1,ncol(sf_L_gt)-1,length=prod(dims)))
png( paste0("Data_likelihood_",Data$DeltaD,".png"), width=dims[2]*4, height=dims[1]*2, res=200, units="in" )
  par( mfrow=dims )
  for( panelI in year_set ){
    plot( sf_L_gt[,panelI], border=NA, logz=FALSE, key.pos=NULL, reset=FALSE )
    add_legend( round(range((sf_L_gt[[panelI]]),na.rm=TRUE),2), text_col="white", legend_x=c(0.9,1.0), legend_y=c(0.45,0.85) )
  }
dev.off()

# Unconditional projection
plotgrid = st_sf( st_geometry(bathy_sf), 
                  "bathymetry" = bathy_sf$ai_bathy_fill, 
                  "preference" = report$h_g,
                  "stationary density" = stationary_g )
png( paste0("Small_results_",Data$DeltaD,".png"), width=3, height=4, res=200, units="in" )
  par( mfrow=c(3,1) )
  for( panelI in 1:3 ){
    plot( plotgrid[,panelI], border=NA, logz=FALSE, key.pos=NULL, reset=FALSE )
    add_legend( round(range((plotgrid[[panelI]]),na.rm=TRUE),2), text_col="white", legend_x=c(0.9,1.0), legend_y=c(0.45,0.85) )
  }
dev.off()

# Smoothed prediction
prob_gt = report$Forward_gt * report$Backward_gt
prob_gt = prob_gt / outer( rep(1,nrow(prob_gt)), colSums(prob_gt) )
colnames(prob_gt) = as.character(as.Date(seq( as.POSIXlt("2019-02-21"), as.POSIXlt("2019-05-24"), by="days")))

plotgrid = st_sf( st_geometry(bathy_sf), 
                  (prob_gt) )
colnames(plotgrid)[-ncol(plotgrid)] = as.character(as.Date(seq( as.POSIXlt("2019-02-21"), as.POSIXlt("2019-05-24"), by="days")))
dims = c(3,3)
year_set = round(seq(1,ncol(prob_gt),length=prod(dims)))
png( paste0("Smoothed_results_",Data$DeltaD,"_Small.png"), width=dims[2]*4, height=dims[1]*2, res=200, units="in" )
  par( mfrow=dims )
  for( panelI in year_set ){
    plot( plotgrid[,panelI], border=NA, logz=FALSE, key.pos=NULL, reset=FALSE )
    add_legend( round(range((plotgrid[[panelI]]),na.rm=TRUE),2), text_col="white", legend_x=c(0.9,1.0), legend_y=c(0.45,0.85) )
  }
dev.off()

dims = c(6,6)
year_set = round(seq(1,ncol(prob_gt),length=prod(dims)))
png( paste0("Smoothed_results_",Data$DeltaD,"_Large.png"), width=dims[2]*4, height=dims[1]*2, res=200, units="in" )
  par( mfrow=dims )
  for( panelI in year_set ){
    plot( plotgrid[,panelI], border=NA, logz=FALSE, key.pos=NULL, reset=FALSE )
    add_legend( round(range((plotgrid[[panelI]]),na.rm=TRUE),2), text_col="white", legend_x=c(0.9,1.0), legend_y=c(0.45,0.85) )
  }
dev.off()

