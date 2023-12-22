

setwd("C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_1")
source( "../Shared_functions/add_legend.R" )
set.seed(101)

########## START IN-TEXT SNIPPET
library(mvtnorm)
library(stars)

# parameters
SD_km = 100
peak_km = 2 # km
best_elevation = 1
SD_logelev = 0.5
n_indiv = 500

# Simulate elevation and density
get_elevation = function(loc){
  e0 = dmvnorm(c(0,0), mean=c(0,0), sigma=SD_km^2*diag(2))
  peak_km * dmvnorm(loc, mean=c(0,0), sigma=SD_km^2*diag(2)) / e0
}
get_density = function(loc){
  elev = get_elevation(loc)
  exp(-1 * ( log(elev/best_elevation) )^2 / SD_logelev^2)
}

# Get values on grid for plotting
loc_gz = expand.grid( "x"=seq(-200,200,length=100),
                      "y"=seq(-200,200,length=100) )
elev_g = get_elevation(loc_gz)
d_g = get_density(loc_gz)

# Function for rejection sampling for locations of individuals
max_density = optim( fn=get_density, par=c(0,100),
                     control=list(fnscale=-1))$value
simulate_location = function( ... ){
  loc = NULL
  while(is.null(loc)){
    samp = runif(n=2,min=-200,max=200)
    D = get_density(samp)
    rand = runif(n=1,min=0,max=max_density)
    if(rand<D) loc = samp
  }
  return(loc)
}
loc_i = t(sapply( 1:n_indiv, FUN=simulate_location ))
########## END IN-TEXT SNIPPET

# Plot
png( paste0("Spatial distribution.png"), width=6, height=2, res=200, units="in" )
  par( mfrow=c(1,3), mgp=c(2,0.5,0), mar=c(3,3,1,1) )
  out = sf::st_as_sf(cbind("d"=d_g,"elev"=elev_g,loc_gz), coords = c("x","y"))
  out = stars::st_rasterize(out)
  image(out["elev"], axes=TRUE, col=viridisLite::magma(10), reset=FALSE, main="Elevation")
  add_legend( round(range(out[["elev"]]),2), col=viridisLite::magma(10), text_col="white" )
  image(out["d"], axes=TRUE, col=viridisLite::viridis(10), main="Habitat suitability")
  add_legend( round(range(out[["d"]]),2), col=viridisLite::viridis(10), text_col="white" )
  plot( loc_i, xlim=c(-200,200), ylim=c(-200,200), pch=20, cex=0.5, xlab="x", ylab="y", main="Individual locations")
dev.off()

####################
# Format and fit as GLM
####################

########## START IN-TEXT SNIPPET
# make sample-level data frame
samples = data.frame("x"=loc_i[,1], "y"=loc_i[,2] )
samples = st_as_sf( samples, coords=c("x","y") )

# Get count in each grid cell
grid_size = 10
grid = st_make_grid( st_bbox(c(xmin=-200, xmax=200, ymin=-200, ymax=200)),
                     cellsize=grid_size )
grid_i = st_intersects( samples, grid )
N_i = tapply( rep(1,nrow(samples)),
              INDEX = factor(unlist(grid_i),levels=1:length(grid)),
              FUN = sum )

# Convert to a data frame for `glm`
Data = data.frame( st_coordinates(st_centroid(grid)),
                   "N"=ifelse(is.na(N_i),0,N_i) )
Data$elev = get_elevation( Data[,c("X","Y")] )

# Fit with canned GLM software
fit = glm( N ~ log(elev) + I(log(elev)^2), data=Data, family=poisson )
########## END IN-TEXT SNIPPET
Table = round(summary(fit)$coef,3)
colnames(Table) = c("Estimate", "Standard Error", "T value", "Probability")
rownames(Table) = c("Intercept", "log(elevation)", "log(elevation)^2")
write.csv( Table, file=paste0("glm_summary_res",grid_size,".csv") )

# plot using sf
grid_sf = st_sf(grid, N=N_i)
png( paste0("gridded_density_res",grid_size,".png"), width=4, height=3, res=200, units="in")
  breaks = pretty(c(0,N_i),n=11)
  plot( grid_sf, axes=TRUE, reset=FALSE, pal=sf.colors(n=length(breaks)-1, alpha=0.2), breaks=breaks, main="Gridded abundance" )
  plot( samples, add=TRUE, pch=20 )
dev.off()

##################
# Fit in TMB
##################

########## START IN-TEXT SNIPPET
# Compile and load TMB
library(TMB)
compile( "poisson_glm.cpp" )
dyn.load("poisson_glm")

# Make covariate data
formula = ~ log(elev) + I(log(elev)^2)
X_ij = model.matrix( formula, data=Data )

# Built inputs for TMB object
data = list( "y_i"=na.omit(Data)$N, "X_ij"=X_ij )
params = list( "beta_j"=rep(0,ncol(X_ij)) )

# Build object
Obj = MakeADFun( data=data, parameters=params )

# Optimize and get standard errors
Opt = nlminb( objective=Obj$fn, grad=Obj$gr, start=Obj$par )
Opt$SD = sdreport( Obj )
########## END IN-TEXT SNIPPET

# Record stuff
log_mu = Obj$report()$log_mu
Table = round(summary(Opt$SD,"fixed"),3)
colnames(Table) = c("Estimate", "Standard Error")
rownames(Table) = c("Intercept", "log(elevation)", "log(elevation)^2")
write.csv( Table, file=paste0("tmb_summary_res",grid_size,".csv") )

#
grid_sf$est_mu = exp(log_mu)
grid_sf$true_mu = get_density( Data[,c("X","Y")] )
  grid_sf$true_mu = grid_sf$true_mu / sum(grid_sf$true_mu) * n_indiv
png( paste0("estimated_density_res",grid_size,".png"), width=6, height=6, res=200, units="in")
  par( mfrow=c(2,2) )
  plot( grid_sf["true_mu"], axes=TRUE, key.pos = NULL, pal=viridisLite::viridis(10), reset = FALSE )
  add_legend( round(range(grid_sf[["true_mu"]]),2), col=viridisLite::viridis(10), text_col="white" )
  plot( grid_sf["N"], axes=TRUE, key.pos = NULL, pal=viridisLite::viridis(10), reset = FALSE )
  add_legend( round(range(grid_sf[["N"]],na.rm=TRUE),2), col=viridisLite::viridis(10), text_col="white" )
  plot( grid_sf["est_mu"], axes=TRUE, key.pos = NULL, pal=viridisLite::viridis(10), reset = FALSE )
  add_legend( round(range(grid_sf[["est_mu"]]),2), col=viridisLite::viridis(10), text_col="white" )
  xylim = range(c(grid_sf[["true_mu"]],grid_sf[["est_mu"]]))
  plot( x=grid_sf$est_mu, y=grid_sf$true_mu, log="xy", pch=20, cex=2, xlim=xylim, ylim=xylim )
  abline( a=0, b=1, lty="dotted" )
dev.off()

set.seed(101)
########## START IN-TEXT SNIPPET
# Simulate 100 new data sets
Sim = sapply(1:100, FUN=function(i) Obj$simulate()$y_i )

# Make function for plots
plot_ecdf = function( y, ysim, ... ){
  plot(ecdf(ysim), xlim=max(ysim)*c(0,1), ... )
  x_intercept = y
  y_min = mean(ysim<y)
  y_max = mean(ysim<=y)
  y_intercept = runif(1,y_min,y_max)
  arrows( x0=x_intercept, x1=x_intercept, y0=0, y1=y_intercept )
  arrows( x0=x_intercept, x1=0, y0=y_intercept, y1=y_intercept )
}

# Make plot showing two observations
par( mar=c(3,3,1,1), mgp=c(2,0.5,0), mfrow=c(1,2) )
plot_ecdf( y=Data$N[4], ysim=Sim[4,], main="ECDF: sample 4" )
plot_ecdf( y=Data$N[6], ysim=Sim[6,], main="ECDF: sample 6" )
########## END IN-TEXT SNIPPET

png( file=paste0("simulation_residuals_res",grid_size,".png"), width=6, height=3, res=200, units="in")
  par( mar=c(3,3,1,1), mgp=c(2,0.5,0), mfrow=c(1,2) )
  plot_ecdf( y=Data$N[4], ysim=Sim[4,], main="ECDF: sample 4" )
  plot_ecdf( y=Data$N[6], ysim=Sim[6,], main="ECDF: sample 6" )
dev.off()

# Simulation residuals
########## START IN-TEXT SNIPPET
library(DHARMa)

# Create DHARMa object with PIT residuals
dharmaRes = createDHARMa(simulatedResponse = Sim,
  observedResponse = data$y_i,
  fittedPredictedResponse = exp(Obj$report()$log_mu),
  integer = FALSE)

# Plot
plot( dharmaRes )
########## END IN-TEXT SNIPPET

png( file=paste0("DHARMa_residuals_res",grid_size,".png"), width=7, height=4, res=200, units="in")
  plot( dharmaRes )
dev.off()

################
# Sample-based interval
################

########## START IN-TEXT SNIPPET
# Sample new parameters from estimation covariance
beta_rj = mvtnorm::rmvnorm( n=500, mean=Opt$par, sigma=Opt$SD$cov.fixed )
x_i = seq( 0.05, 2, length=1000 )
X_ij = model.matrix( formula, data=data.frame("elev"=x_i) )

# Calculate response curve for each simulated parameter-vector
yhat_ir = scale(X_ij %*% t(beta_rj), scale=FALSE)
ybar_i = scale(X_ij %*% Opt$par, scale=FALSE)
ytrue_i = scale(-1 * (log(x_i) - log(best_elevation))^2 / SD_logelev^2, scale=FALSE)
yci_zi = apply( yhat_ir, MARGIN=1, FUN=quantile, prob=c(0.05,0.95) )

# Plot marginal effect
plot( x=x_i, y=ybar_i, log="x", type="l", lwd=2, col='blue', xlab="elevation", ylab="log-relative density", ylim=c(-30,2) )
lines( x=x_i, y=ytrue_i, lwd=2, col="black")
polygon( x=c(x_i,rev(x_i)), y=c(yci_zi[1,],rev(yci_zi[2,])), col=rgb(0,0,1,0.2) )
########## END IN-TEXT SNIPPET

# Plot
png( file=paste0("response_curve_res",grid_size,".png"), width=4, height=4, res=200, units="in")
  par( xaxs="i" )
  plot( x=x_i, y=ybar_i, log="x", type="l", lwd=2, col='blue', xlab="elevation", ylab="log-relative density", ylim=c(-30,2) )
  lines( x=x_i, y=ytrue_i, lwd=2, col="black")
  polygon( x=c(x_i,rev(x_i)), y=c(yci_zi[1,],rev(yci_zi[2,])), col=rgb(0,0,1,0.2) )
dev.off()

################
# Interval using marginaleffects package
################

########## START IN-TEXT SNIPPET
# Source marginaleffects functions
library(marginaleffects)
library(ggplot2)
source( "../Shared_functions/marginaleffects.R" )

# Define input expected by marginaleffects functions
fit = list( obj = Obj,
            opt = Opt,
            data = Data,
            formula = formula,
            parhat = Obj$env$parList() )
class(fit) = "custom_tmb"

# Get prediction
pred = predictions( fit, param="beta_j", newdata=data.frame("elev"=x_i) )
pred$true = -1 * (log(x_i) - log(best_elevation))^2 / SD_logelev^2
pred$true = pred$true - mean(pred$true)
pred[,c("estimate","conf.low","conf.high")] = pred[,c("estimate","conf.low","conf.high")] - mean(pred[,"estimate"])

# Plot using ggplot
ggplot( as.data.frame(pred), aes(elev, estimate)) +
  geom_line( aes(y=estimate), color="blue", size=3 ) +
  geom_ribbon( aes( x=elev, ymin=conf.low, ymax=conf.high), fill=rgb(0,0,1,0.2) ) +
  geom_line( aes(y=true), color="black", size=3 ) +
  scale_x_continuous(trans='log2') +
  labs(y="Predicted response")
########## END IN-TEXT SNIPPET
ggsave( paste0("marginaleffects_curve_res",grid_size,".png"), device=png, width=4, height=4 )
