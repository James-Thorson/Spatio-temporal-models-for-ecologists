
setwd( R'(C:\Users\James.Thorson\Desktop\Git\Spatio-temporal-models-for-ecologists\Chap_7)' )
source( "../Shared_functions/rmvnorm_prec.R" )

library(sf)
library(Matrix)
library(TMB)
sf_area = rnaturalearth::ne_countries( scale=10, country="Finland", return="sf")

# extract main island
sf_area = st_cast(st_geometry(sf_area),"POLYGON")[[1]]
sf_area = st_sfc( sf_area, crs="+proj=longlat +datum=WGS84" )
plot(sf_area)

# make intro grid
cellsize = 0.5
sf_fullgrid = st_make_grid( sf_area, cellsize=cellsize )
sf_grid = st_make_valid(st_intersection( sf_fullgrid, sf_area ))
sf_grid = sf_grid[ st_area(sf_grid)>(0.2*max(st_area(sf_grid))) ]
#grid = elevatr::get_elev_point( locations=st_point_on_surface(sf_grid), src = "aws" )
#grid$elevation = ifelse( grid$elevation<1, 1, grid$elevation)

# Get adjacency
st_rook = function(a, b = a, ...) st_relate(a, b, pattern = "F***1****", ... )
grid_A = st_rook( sf_grid, sparse=TRUE )
A_ss = as(grid_A,"sparseMatrix")
A_ss = as(A_ss, "TsparseMatrix")
I_ss = Matrix::sparseMatrix( i=1:length(sf_grid), j=1:length(sf_grid), x=rep(1,length(sf_grid)) )
I_ss = as(I_ss, "TsparseMatrix")

# Calculate min and max eigenvalues, range( 1 / eigen(as.matrix(grid_A))$values )
# Using igraph for sparse-matrix calculations
graph = igraph::graph_from_adjacency_matrix( A_ss, weighted=TRUE )
f2 = function(x, extra=NULL) { as.vector(x %*% A_ss) }
rho_min = 1 / Re(igraph::arpack(f2, sym=FALSE, options=list(n=nrow(A_ss), nev=3, ncv=8, which="SR"))$values[1])
rho_max = 1 / Re(igraph::arpack(f2, sym=FALSE, options=list(n=nrow(A_ss), nev=3, ncv=8, which="LR"))$values[1])

# Make ICAR precision
rho = rho_max * 0.75
var_X1 = 0.5^2
var_X2 = 0.5^2
#Q = diag(rowSums(A)) - rho*A
Q = diag(rep(1,nrow(A_ss))) - rho*A_ss

#
X1 = rmvnorm_prec( mu=rep(0,nrow(A_ss)), prec=Q/var_X1, n.sims=1 )
X2 = rmvnorm_prec( mu=rep(0,nrow(A_ss)), prec=Q/var_X2, n.sims=1 )
lambda_s = exp( 1 + X1 + X2 )

#
n_obs = nrow(A_ss)
s_i = sample( 1:nrow(A_ss), size=n_obs, replace=FALSE )
c_i = rpois( n_obs,  lambda_s[s_i] )

#
formula = ~ 0 # + X1 + X2
X_sk = model.matrix( formula, data.frame("X1"=X1,"X2"=X2) )

#
compile( "CAR.cpp" ) # framework='TMBad'
dyn.load( dynlib("CAR") )

# Compile
Params = list( "beta0" = 0,
               "rho_prime" = 0,
               "ln_sigma" = 0,
               "gamma_k" = rep(0,ncol(X_sk)),
               "omega_s" = rep(0,nrow(A_ss)) )
Data = list( "c_i" = c_i,
             "rho_bounds" = c( rho_min, rho_max ),
             "s_i" = s_i-1,
             "X_sk" = X_sk,
             "I_ss" = I_ss,
             "A_ss" = A_ss )

# Build and optimize
obj = MakeADFun( data=Data, parameters=Params, random="omega_s", DLL="CAR", silent=TRUE )
opt = nlminb( objective=obj$fn, grad=obj$gr, start=obj$par )
opt$SD = sdreport( obj, bias.correct=TRUE )
report = obj$report()

# When sampling all locations, sqrt(sum(lambda_s)) = SE(sumlambda)
summary(opt$SD, 'report')
sum(lambda_s)
sqrt(sum(lambda_s))

sink('check_hypotheses.txt')
cat("> summary(opt$SD, 'report')\n")
summary(opt$SD, 'report')
cat("> sum(lambda_s)\n")
sum(lambda_s)
cat("> sqrt(sum(c_i))\n")
sqrt(sum(c_i))
sink()

###########
# Simulate test
###########

n_obs = 100
var_X1 = 0.6 ^ 2   # sigmaX1=0.6, sigmaX2=0.3 works well
var_X2 = 0.3 ^ 2

n_rep = 100
results_crz = array(NA, dim=c(2,n_rep,3), dimnames=list("X1"=c("No","Yes"),"Rep"=1:n_rep,"Param"=c("True","Est","SE")) )
for(r in seq_len(n_rep)){
for(c in 1:2){
  set.seed(r)
  X1 = rmvnorm_prec( mu=rep(0,nrow(A_ss)), prec=Q/var_X1, n.sims=1 )
  X2 = rmvnorm_prec( mu=rep(0,nrow(A_ss)), prec=Q/var_X2, n.sims=1 )
  lambda_s = exp( 1 + X1 + X2 )
  s_i = sample( 1:nrow(A_ss), size=n_obs, replace=FALSE )
  c_i = rpois( n_obs,  lambda_s[s_i] )
  results_crz[c,r,"True"] = sum(lambda_s)
  # Different configurations
  if( c==1 ) formula = ~ 0
  if( c==2 ) formula = ~ 0 + X1
  if( c==3 ) formula = ~ 0 + X1 + X2
  X_sk = model.matrix( formula, data.frame("X1"=X1,"X2"=X2) )
  # Build and optimize
  Params$gamma_k = rep(0,ncol(X_sk))
  Data = list( "c_i" = c_i,
               "rho_bounds" = c( rho_min, rho_max ),
               "s_i" = s_i-1,
               "X_sk" = X_sk,
               "I_ss" = I_ss,
               "A_ss" = A_ss )
  obj = MakeADFun( data=Data, parameters=Params, random="omega_s", DLL="CAR", silent=TRUE )
  opt = nlminb( objective=obj$fn, grad=obj$gr, start=obj$par )
  opt$SD = sdreport( obj, bias.correct=TRUE, bias.correct.control=list(sd=TRUE) )
  results_crz[c,r,c("Est","SE")] = summary(opt$SD, 'report')['sumlambda',c('Est. (bias.correct)','Std. (bias.correct)')]
}}

#
apply((results_crz[,,'Est']-results_crz[,,'True'])/results_crz[,,'True'], MARGIN=1, FUN=mean, na.rm=TRUE)
apply(results_crz[,,'SE'], MARGIN=1, FUN=mean, na.rm=TRUE)

Relative_error = (results_crz[,,'Est']-results_crz[,,'True']) /results_crz[,,'True']
DF1 = cbind( expand.grid(dimnames(Relative_error)), "Metric"="relative error", "Value"=as.vector(Relative_error) )
DF2 = cbind( expand.grid(dimnames(results_crz[,,'SE'])), "Metric"="standard error", "Value"=as.vector(results_crz[,,'SE']) )
DF = rbind(DF1,DF2)

library( ggplot2 )
#ggplot(ParHat, aes(x=Config, y=Value)) +
#  geom_point( ) +
#  geom_errorbar(aes(ymin=Est-SE, ymax=Est+SE), width=.2) +
ggplot(DF, aes(X1, Value, colour=X1) ) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap( vars(Metric), scales="free" )
ggsave( filename="Simulation_results.png", device="png", width=6, height=3 )
