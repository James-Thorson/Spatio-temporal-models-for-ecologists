# Data processed from: C:\Users\James.Thorson\Desktop\Work files\AFSC\2022-06 -- Breeding Bird survey
setwd("C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_6")

library(sf)
library(fmesher)
library(TMB)
library(rnaturalearth)

source("../Shared_functions/rmvnorm_prec.R")
source("../Shared_functions/sample_var.R")
source("../Shared_functions/add_legend.R")

# Ozone data
ozone = st_read( "2019_ozone.csv", options=c("X_POSSIBLE_NAMES=X","Y_POSSIBLE_NAMES=Y"), crs=st_crs("+proj=longlat +datum=WGS84") )
  ozone$Daily.Max.8.hour.Ozone.Concentration = as.numeric(ozone$Daily.Max.8.hour.Ozone.Concentration)
ozone <- subset(ozone, Date=="07/01/2019")

#
pop_dens = st_read( "population_density.csv", options=c("X_POSSIBLE_NAMES=X","Y_POSSIBLE_NAMES=Y"), crs=st_crs("+proj=longlat +datum=WGS84") )
  pop_dens$Dens2000 = as.numeric(pop_dens$Dens2000)
  pop_dens$Dens2020 = as.numeric(pop_dens$Dens2020)
pop_dens_stars = stars::st_rasterize(pop_dens)
plot(pop_dens_stars['Dens2000'])

# Define states
state_set = c("Virginia","North Carolina","South Carolina","Georgia","Florida","Maryland","Delaware","District of Columbia","New Jersey","Pennsylvania","New York")
states_sf = ne_states( c("United States of America","Canada"), return="sf")
states_sf = states_sf[pmatch(state_set, states_sf$name),]
domain_sf = st_union( states_sf )

# Get pop-dens for nearest sample
index = st_nearest_feature( ozone, pop_dens )
ozone$pop_dens = pop_dens[index,'Dens2020']
plot(ozone$pop_dens)

#
library(mgcv)
DF = data.frame( "AQI"=ozone$Daily.Max.8.hour.Ozone.Concentration, "Dens"=ozone$pop_dens$Dens2020/1000, st_coordinates(ozone) )
Gam = gam( AQI ~ te(X,Y) + Dens, data=DF )
summary(Gam)
plot(Gam)

# Get pop-dens for nearest point
#plotdata = sp::spsample(as(states_sf,"Spatial"), type="regular", cellsize=0.25)
plotdata = st_make_grid( domain_sf, cellsize=0.25 )
plotdata = st_intersection( plotdata, domain_sf )
index = st_nearest_feature( plotdata, pop_dens )
plotdata = data.frame( "Dens"=pop_dens$Dens2020[index], st_coordinates(st_centroid(plotdata)) )
plotdata$pred = predict(Gam, newdata=plotdata )
plotdata <- st_as_sf(plotdata, coords=c("X","Y"), crs=st_crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#plotdata = stars::st_rasterize(plotdata)
plot(plotdata['pred'], add=TRUE)

###################
# Unequal distance 2D autoregressive
###################

#ozone = ozone[sample(1:nrow(ozone),size=100),]

# Create extrapolation grid
cellsize = 0.25
grid = st_make_grid( states_sf, cellsize=cellsize )
grid = st_intersection( grid, domain_sf )

# create mesh
mesh = fm_mesh_2d( st_coordinates(st_centroid(grid)), refine=TRUE, cutoff=0.2)
# Create matrices in INLA
spde <- fm_fem(mesh, order=2)
# create projection matrix from vertices to samples
A_is = fm_evaluator( mesh, loc=st_coordinates(ozone) )$proj$A
# create projection matrix from vertices to grid
A_gs = fm_evaluator( mesh, loc=st_coordinates(st_centroid(grid)) )$proj$A

# formula
formula = ~ 0 + log(Dens2000)
options(na.action='na.pass')
DF_i = data.frame(pop_dens[st_nearest_feature(ozone, pop_dens),])
X_ij = model.matrix(formula, data= DF_i )
DF_g = data.frame(pop_dens[st_nearest_feature(grid, pop_dens),])
X_gj = model.matrix( formula, data=DF_g )

# Expansion rates
a_g = as.numeric(st_area(grid))
d_g = DF_g$Dens2020

# COmpile
compile( "SPDE.cpp", flags="-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign" )
dyn.load( dynlib("SPDE") )

# Build object
Data = list( "y_i"=ozone$Daily.Max.8.hour.Ozone.Concentration, "A_is"=A_is,
             "A_gs"=A_gs, "M0"=spde$c0, "M1"=spde$g1,
             "M2"=spde$g2, "X_ij"=X_ij, "X_gj"=X_gj, "a_g"=a_g, "d_g"=d_g )
Params = list( "beta0"=0, "ln_tau"=0, "ln_kappa"=0, "ln_sigma"=0,
               "beta_j"=rep(0,ncol(X_ij)), "omega_s"=rnorm(mesh$n) )
Obj = MakeADFun( data=Data, parameters=Params, random="omega_s" )

########## START IN-LINE SNIPPET
# Optimize
Opt = nlminb( start=Obj$par, obj=Obj$fn, gr=Obj$gr )
Opt$SD = sdreport( Obj, bias.correct=TRUE, getJointPrecision=TRUE, bias.correct.control=list(sd=TRUE) )

# Sample standard errors
SE_g = sample_var( obj = Obj,
                   var_name = "yhat_g",
                   mu = Obj$env$last.par.best,
                   prec = Opt$SD$jointPrecision )
########## END IN-LINE SNIPPET

# Plot
report = Obj$report()
grid_sf = st_sf(grid, Ozone=report$yhat_g, "SE"=SE_g, Dens2020=(DF_g$Dens2020) )

png( "Mapped_ozone.png", width=7, height=7, res=200, units="in")
  par( mfrow=c(2,2), mar=c(1,1,1,1) )
  plot( st_geometry(states_sf), main=expression(Mesh ~~ and ~~ data) )
  plot( mesh, add=TRUE )
  plot( ozone[,'Daily.Max.8.hour.Ozone.Concentration'], add=TRUE, pch=20, border=NA )
  add_legend( legend=round(range(ozone$Daily.Max.8.hour.Ozone.Concentration),3) )
  plot( grid_sf['Ozone'], key.pos=NULL, reset=FALSE, main=expression(Predicted ~~ ozone ~~ (ppm)), border=NA, nbreaks=1000 )
  add_legend( legend=round(range(grid_sf$Ozone),3) )
  plot( st_geometry(states_sf), add=TRUE )
  plot( grid_sf['SE'], key.pos=NULL, reset=FALSE, main=expression(Standard ~~ Error ~~ (ppm)), border=NA, nbreaks=1000 )
  add_legend( legend=round(range(grid_sf$SE),3) )
  plot( st_geometry(states_sf), add=TRUE )
  plot( grid_sf['Dens2020'], key.pos=NULL, reset=FALSE, main=expression(Density ~~ (Count/km^2)), border=NA, logz=TRUE, nbreaks=1000 )
  add_legend( legend=round(range(grid_sf$Dens2020),2) )
  plot( st_geometry(states_sf), add=TRUE )
dev.off()

# Calculate averages
library(ggplot2)
Summary = summary(Opt$SD,"report")
ParHat = data.frame("Value"=c("Y1","Y2","Y3"), "Est"=Summary[,1], "SE"=Summary[,2])
  ParHat = rbind(ParHat, data.frame("Value"=c("Z1","Z2","Z3"), "Est"=Summary[,3], "SE"=Summary[,4]) )
ggplot(ParHat, aes(x=Value, y=Est)) +
  geom_point( ) +
  geom_errorbar(aes(ymin=Est-SE, ymax=Est+SE), width=.2)
ggsave( filename="Totals.png", device="png" )


###############
# Simulation design
#  Takes several hours
###############

design_set = c("simple_random", "systematic", "opportunistic", "opport. & fit covariate")
n_reps = 250
Y_drz = array( NA, dim=c(length(design_set),n_reps,2,3), dimnames=list("Design"=design_set,1:n_reps,"Estimator"=c("plugin","epsilon"),c("Est","SE","True")) )
sample_size = 50
metric = "Y2"
ozone_d = vector("list", length=length(design_set))

d = r = 1
for( r in 1:n_reps ){
for( d in seq_along(design_set) ){
  set.seed( 101 + r )
  sim = Obj$simulate()

  if( design_set[d] == "simple_random" ){
    included_samples = sample( 1:length(sim$y_g), size=sample_size, replace=FALSE )
  }else if( design_set[d] == "systematic" ){
    included_samples = seq( 1, length(sim$y_g), length=sample_size )
  }else if( design_set[d] %in% c("opportunistic","opport. & fit covariate") ){
    included_samples = sample( 1:length(sim$y_g), size=sample_size, prob=d_g, replace=FALSE )
  }else{stop()}

  ysim_i = sim$y_g[included_samples]
  ozone_d[[d]] = grid[included_samples]

  # Lowest resolution mesh
  mesh_sim = fm_mesh_2d( st_coordinates(st_centroid(grid)), plot.delay=NULL, refine=TRUE, cutoff=1)
  # Create matrices in INLA
  spde_sim <- fm_fem(mesh_sim, order=2)
  # create projection matrix from vertices to samples
  A_is = fm_evaluator( mesh_sim, loc=st_coordinates(st_centroid(grid[included_samples])) )$proj$A
  # create projection matrix from vertices to grid
  A_gs = fm_evaluator( mesh_sim, loc=st_coordinates(st_centroid(grid)) )$proj$A

  # Remake covariates
  if( design_set[d] == "opport. & fit covariate" ){
    DF_i = data.frame(pop_dens[st_nearest_feature(grid[included_samples], pop_dens),])
    X_ij = model.matrix(formula, data= DF_i )
    DF_g = data.frame(pop_dens[st_nearest_feature(grid, pop_dens),])
    X_gj = model.matrix( formula, data=DF_g )
  }else{
    X_ij = matrix(ncol=0, nrow=sample_size)
    X_gj = matrix(ncol=0, nrow=length(grid))
  }

  # Build object
  Data_sim = list( "y_i"=ysim_i, "A_is"=A_is,
               "A_gs"=A_gs, "M0"=spde_sim$c0, "M1"=spde_sim$g1,
               "M2"=spde_sim$g2, "X_ij"=X_ij, "X_gj"=X_gj, "a_g"=a_g, "d_g"=d_g )
  Params_sim = list( "beta0"=0, "ln_tau"=0, "ln_kappa"=0, "ln_sigma"=0,
                 "beta_j"=rep(0,ncol(X_ij)), "omega_s"=rnorm(mesh_sim$n) )
  Obj_sim = MakeADFun( data=Data_sim, parameters=Params_sim, random="omega_s" )
    Obj_sim$env$beSilent()

  # Optimize
  Opt_sim = nlminb( start=Obj_sim$par, gr=Obj_sim$gr, obj=Obj_sim$fn )
  # Custom check for convergence
  if( all(abs(Opt_sim$par)<5) ){
    Opt_sim$SD = sdreport( Obj_sim, bias.correct=TRUE, bias.correct.control=list(sd=TRUE) )
    Y_drz[d,r,"plugin",1:2] = summary(Opt_sim$SD,"report")[metric,c('Estimate','Std. Error')]
    Y_drz[d,r,"plugin",3] = sim[[metric]]
    Y_drz[d,r,"epsilon",1:2] = summary(Opt_sim$SD,"report")[metric,c('Est. (bias.correct)','Std. (bias.correct)')]
    Y_drz[d,r,"epsilon",3] = sim[[metric]]
  }
}}

# Plot design
png( file="Design.png", width=6, height=4, res=200, units='in' )
  par( mfrow=c(1,3), mar=c(1,1,1,0), mgp=c(2,0.5,0) )
  for(d in seq_along(design_set[1:3])){
    plot( st_geometry(states_sf), main=design_set[d] )
    #plot( plot_grid['Dens2020'], add=TRUE, logz=TRUE )
    plot( ozone_d[[d]], add=TRUE, col="red", cex=3)
  }
dev.off()

# Plot errors
Error_dr = Y_drz[,,,'Est'] - Y_drz[,,,'True']
DF = data.frame( expand.grid(dimnames(Error_dr[1:3,,])), "Error"=as.numeric(Error_dr[1:3,,]) )
ggplot(DF, aes(x=Design,y=Error,colour=Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0) # + geom_jitter(width = 0.2)
ggsave( filename="Design_errors.png", device="png" )

# Preferential sampling mitigation
DF = data.frame( expand.grid(dimnames(Error_dr[4,,,drop=FALSE])), "Error"=as.numeric(Error_dr[4,,,drop=FALSE]) )
ggplot(DF, aes(Design,Error,colour=Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0) # + geom_jitter(width = 0.2)
ggsave( filename="Preferential_sampling.png", device="png" )

# Plot errors
Error_dr = Y_drz[,,,'Est'] - Y_drz[,,,'True']
DF = data.frame( expand.grid(dimnames(Error_dr[,,])), "Error"=as.numeric(Error_dr[,,]) )
ggplot(DF, aes(x=Design,y=Error,colour=Estimator)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0) # + geom_jitter(width = 0.2)
ggsave( filename="Design_errors_combined.png", device="png" )


