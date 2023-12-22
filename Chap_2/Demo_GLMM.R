


setwd("C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_2")

#data_dir = "C:/Users/James.Thorson/Desktop/Work files/AFSC/2022-05 -- Barro Colorado data/"
#load( paste0(data_dir,"bci.stem.Rdata31Aug2012/bci.stem1.rdata") )
#vismba = subset( bci.stem1, sp %in% "vismba" )
#saveRDS( vismba, "vismba.rds")
vismba = readRDS( "vismba.rds" )

#       
library(sf)
samples = data.frame("x"=vismba$gx, "y"=vismba$gy, "agb"=vismba$agb )
samples = st_as_sf( samples, coords=c("x","y") )

grid = st_make_grid( st_bbox(c(xmin=0, xmax=1000, ymin=0, ymax=500)), cellsize=c(250,125) )
grid_i = st_intersects( samples, grid )
Count_i = tapply( samples$agb, INDEX=factor(unlist(grid_i),levels=1:length(grid)), FUN=length )
Data = data.frame( st_coordinates(st_centroid(grid)), "Count"=ifelse(is.na(Count_i),0,Count_i) )

grid_sf = st_sf(grid, Count=Count_i)
png( "gridded_density.png", width=5, height=2.5, res=200, units="in")
  plot( grid_sf, axes=TRUE, reset=FALSE, pal=sf.colors(n=10, alpha=0.2), breaks=seq(0,max(Data$Count),length=11) )
  plot( samples, add=TRUE, pch=20 )
dev.off()

#############
# Compare with GLMM
#############

########## START IN-TEXT SNIPPET
library(lme4)
Data$Site = factor(1:nrow(Data))
Lme = glmer( Count ~ 1 + (1|Site), data=Data, family=poisson(link = "log") )
########## END IN-TEXT SNIPPET
capture.output(Lme, file="Lme_output.txt" )

#################
# Demo
#################

########## START IN-TEXT SNIPPET
# Joint negative-log-likelihood
joint_nll = function( prandom, pfixed, plist, Data, random ){
  # Read in values
  phat_fixed = relist(pfixed, plist[setdiff(names(plist),random)])
  phat_random = relist(prandom, plist[random])
  phat = c(phat_fixed,phat_random)
  # Data likelihood
  jnll = 0
  for( i in 1:nrow(Data) ){
    jnll = jnll - dpois(Data$Count[i],lambda=exp(phat$eps[i]),log=TRUE)
  }
  # Random effect distribution
  for( i in 1:nrow(Data) ){
    jnll = jnll - dnorm(phat$eps[i],mean=phat$logmu,sd=exp(phat$logsd),log=TRUE)
  }
  return(jnll)
}
########## END IN-TEXT SNIPPET

########## START IN-TEXT SNIPPET
# Outer optimization function
marg_nll = function( pfixed, plist, Data, random, jnll, what="laplace" ){
  # Inner optimizer
  prandom = unlist(plist[random])
  inner = nlminb(start=prandom, objective=jnll, pfixed=pfixed, plist=plist, Data=Data, random=random)

  # Calculate the Laplace approximation
  inner$hessian = optimHess(par=inner$par, fn=jnll, pfixed=pfixed, plist=plist, Data=Data, random=random)
  inner$laplace = inner$objective + 0.5*log(det(inner$hessian)) - length(prandom)/2*log(2*pi)
  if(what=="laplace") return(inner$laplace)
  if(what=="full") return(inner)
}
########## END IN-TEXT SNIPPET

########## START IN-TEXT SNIPPET
# Define inputs
plist = list("logmu"=0, "logsd"=0, "eps"=rep(0,nrow(Data)))
random = "eps"

# Identify MLE using outer optimizer
opt1 = nlminb( objective = marg_nll,
               start = unlist(plist[setdiff(names(plist),random)]),
               jnll = joint_nll,
               plist = plist,
               Data = Data,
               random = random,
               control = list(trace=1) )

# Compute finite-difference approx. to outer Hessian
Hess = optimHess( par = opt1$par,
                  fn = marg_nll,
                  plist = plist,
                  Data = Data,
                  random = random,
                  jnll = joint_nll )
########## END IN-TEXT SNIPPET
parhat = cbind( "Estimate"=opt1$par, "SE"=sqrt(diag(solve(Hess))) )
write.csv( signif(parhat,3), file="opt1_output.csv" )

###############
# Improve R by providing gradient
###############

########## START IN-TEXT SNIPPET
# Function to calculate finite difference gradient
grad = function( pfixed, ... ){
  delta = 0.0001
  gr = rep(0,length(pfixed))
  for(i in seq_along(gr)){
    dvec = rep(0,length(pfixed))
    dvec[i] = delta
    # Calculate central finite difference
    val1 = marg_nll(pfixed=pfixed+dvec, ... )
    val0 = marg_nll(pfixed=pfixed-dvec, ...)
    gr[i] = (val1 - val0) / (2*delta)
  }
  return(gr)
}

# Re-estimate parameters
opt2 = nlminb( objective = marg_nll,
               gradient = grad,
               start = unlist(plist[setdiff(names(plist),random)]),
               jnll = joint_nll,
               plist = plist,
               Data = Data,
               random = random,
               control = list(trace=1) )
Hess = optimHess( par = opt2$par,
                  fn = marg_nll,
                  plist = plist,
                  Data = Data,
                  random = random,
                  jnll = joint_nll )
########## END IN-TEXT SNIPPET
final_grad = grad( opt2$par, plist=plist, Data=Data, random=random, jnll=joint_nll )
parhat = cbind( "Estimate"=opt2$par, "SE"=sqrt(diag(solve(Hess))), "final_grad"=final_grad )
write.csv( signif(parhat,c(3,3,3,3,12)), file="opt2_output.csv" )

##############
# Compare with TMB
##############

########## START IN-TEXT SNIPPET
# Compile and load TMB model
library(TMB)
compile("poisson_glmm.cpp")
dyn.load("poisson_glmm")

# Define inputs
data = list( "y_i"=Data$Count )
params = list( "eps_i"=rep(0,nrow(Data)),
               "ln_mu"=0,
               "ln_sd"=0 )

# Compile TMB object while declaring random effects for inner loop
Obj = MakeADFun( data=data, parameters=params, random="eps_i" )

# Outer optimizer using gradients
Opt = nlminb( start=Obj$par, obj=Obj$fn, grad=Obj$gr )

# Compute standard errors including joint precision
Opt$SD = sdreport( Obj, getJointPrecision=TRUE )
########## END IN-TEXT SNIPPET
parhat = summary(Opt$SD, "fixed")
write.csv( signif(cbind(parhat, Obj$gr(Opt$par)),c(3,3,3,3,12)), file="tmb_glmm.csv")

#############
# Replicate jointPrecision calculation
#############

########## START IN-TEXT SNIPPET
# Extract predicted random effects
fn = function( pfixed ){
  inner = marg_nll( pfixed=pfixed, plist=plist, Data=Data, random=random, jnll=joint_nll, what="full" )
  return(inner$par)
}
# calculate outer Jacobian (gradient) matrix
G = numDeriv::jacobian( func=fn, x=opt2$par )
# Inner and outer Hessian matrices
H1 = marg_nll( pfixed=opt2$par, plist=plist, Data=Data, random=random, jnll=joint_nll, what="full" )$hessian
H2 = Hess
# Assemble joint precision
Q_joint <- rbind(
  cbind2( H1, -H1 %*% G ),
  cbind2( -t(G) %*% H1, as.matrix(t(G) %*% H1 %*% G) + H2 )
)
########## END IN-TEXT SNIPPET

# Check
zapsmall(as.matrix(Q_joint - Opt$SD$jointPrecision))

#############
# Replicate variance of ADREPORTed variables
#############

# Get gradient for ADREPORTed variables
fn = function(pvec) unlist(Obj$report(pvec)[c('yhat_i','yhat_sum')])
G = numDeriv::jacobian( func=fn, x=Obj$env$last.par.best )

# Calculate manually
V_manual = G %*% solve(Q_joint) %*% t(G)

# Check
V_tmb = Opt$SD$cov
zapsmall(as.matrix(V_tmb - V_manual))

#############
# Visualize Hessian
#############

png( "TMB_hessian.png", width=4, height=4, res=200, unit="in")
  inner_hessian = Obj$env$spHess(par=Obj$env$last.par.best,random=TRUE)
  Matrix::image( inner_hessian, main="Inner Hessian in TMB" )
dev.off()

# Visualize graph
library(diagram)
n_groups = 3
names <- c("ln_mu", "ln_SD", paste0("eps",1:n_groups), paste0("Y",1:n_groups) )
M <- array(0, dim=c(length(names),length(names)), dimnames=list(names,names))
  M[grep("eps",names),'ln_mu'] = ""
  M[grep("eps",names),'ln_SD'] = ""
  M[cbind(grep("Y",names),grep("eps",names))] = ""
box.type = c( rep("circle",2), rep("diamond",n_groups), rep("square",n_groups) )
png( "graph.png", width=7, height=3, res=200, un="in")
  par( mar=c(0,0,0,0) )
  plotmat( M, pos = c(2,n_groups,n_groups), curve = 0, name = names, lwd = 3, arr.pos = 0.6,
           box.lwd = 2, box.type = box.type, box.prop = 0.5, cex.txt=0.8 )
dev.off()

##############
# Compare with STAN
##############

########## START IN-TEXT SNIPPET
library(tmbstan)
fit <- tmbstan(Obj, chains=1)
########## END IN-TEXT SNIPPET
plot(fit, pars=setdiff(names(fit),"lp__"))

# Plot TMB and STAN
ln_pred_stan = colMeans(extract(fit,'eps_i')[[1]])
ln_pred_laplace = summary(Opt$SD,"random")[,'Estimate']
grid_sf = st_sf(grid, "STAN"=exp(ln_pred_stan), "Laplace"=exp(ln_pred_laplace) )
png( "gridded_predictions.png", width=6, height=2, res=200, units="in")
  par( mfrow=c(1,2) )
  plot( grid_sf[1], axes=TRUE, reset=FALSE, key.pos=NULL, pal=sf.colors(n=10, alpha=0.2), breaks=seq(0,max(Data$Count),length=11) )
  plot( samples, add=TRUE, pch=20 )
  plot( grid_sf[2], axes=TRUE, reset=FALSE, key.pos=NULL, pal=sf.colors(n=10, alpha=0.2), breaks=seq(0,max(Data$Count),length=11) )
  plot( samples, add=TRUE, pch=20 )
dev.off()

#############
# Profile over variance
#############

eps_i = summary(Opt$SD,"random")[,'Estimate']
LL_z = -1 * c( Obj$report()$jnll, Opt$objective )

# Profile
ln_sd_set = seq(-4,2,length=100)
eps_ji = array(NA, dim=c(length(ln_sd_set),length(eps_i)) )
LL_jz = array(NA, dim=c(length(ln_sd_set),length(LL_z)) )
for( j in seq_along(ln_sd_set) ){
  params$ln_sd = ln_sd_set[j]
  map = list("ln_sd"=factor(NA))
  Obj_j = MakeADFun( data=data, parameters=params, random="eps_i", map=map )
  Opt_j = nlminb( start=Obj_j$par, obj=Obj_j$fn, grad=Obj_j$gr )
  Opt_j$SD = sdreport( Obj_j, getReportCovariance=TRUE )
  eps_ji[j,] = summary(Opt_j$SD,"random")[,'Estimate']
  LL_jz[j,] = -1 * c( Obj_j$report()$jnll, Opt_j$objective )
}
neghalflogdetH_j = LL_jz[,2] - LL_jz[,1] - length(eps_i)/(2*log(2*pi))
LL_jz = cbind(LL_jz, -neghalflogdetH_j )

Glm0 = glm( Count ~ 1, data=Data, family=poisson(link = "log") )
Glm1 = glm( Count ~ 0 + Site, data=Data, family=poisson(link = "log") )
png( file="Shrinkage.png", width=5, height=5, res=200, units="in" )
  par( mfrow=c(2,1), mgp=c(2,0.5,0), mar=c(3,3,1,1), xaxs="i" )
  matplot( x=exp(ln_sd_set), y=LL_jz, type="l", lty="solid", col=c("blue","black","red"), lwd=2, log="x", xlab="", ylab="Values" )
    abline( v=exp(Opt$par['ln_sd']), lty="dotted" )
    legend( "left", fill=c("blue","black","red"), legend=c("joint loglike","marginal loglike","0.5*log(|H|)"), bty="n" )
  matplot( x=exp(ln_sd_set), y=eps_ji, type="l", col="black", lwd=2, lty="solid", log="x", xlab=expression(sigma[epsilon]), ylab=expression(epsilon[i]) )
    abline( v=exp(Opt$par['ln_sd']), lty="dotted" )
    points( x=exp(ln_sd_set[2]), y=Glm0$coef, cex=2, col="lightgreen", pch=20 )
    points( x=rep(exp(rev(ln_sd_set)[2]),length(Glm1$coef)), y=Glm1$coef, cex=2, col="darkgreen", pch=20 )
    legend( "bottomleft", fill=c("lightgreen","darkgreen"), legend=c("Intercept-only GLM","Site-factor GLM"), bty="n" )
dev.off()

#############
# Intentionally cause inner-Hessian crash
#############

# Introduce bug
buggy_params = params
buggy_params$eps_i = rep(0, nrow(Data)+1)
Obj = MakeADFun( data=data, parameters=buggy_params, random="eps_i" )

# Display error
my_log = file( "Inner_hessian_error.txt" )
sink( my_log, append = TRUE, type = "output")
Obj$fn(Obj$par)
sink( file=NULL )

# Visualize inner hessian
png( "TMB_hessian_bug.png", width=4, height=4, res=200, unit="in")
  inner_hessian = Obj$env$spHess(random=TRUE)
  Matrix::image( inner_hessian, main="Inner Hessian in TMB\nwith non-invertible inner Hessian" )
dev.off()

########## START IN-TEXT SNIPPET
# Map off fixed effects and treat random as fixed
FE = setdiff( names(params), "eps_i" )
map <- lapply( FE, function(x) factor(params[[x]]*NA) )
names(map) = FE
# Rebuild and run optimizer
Test = MakeADFun( data=data, parameters=buggy_params, map=map )
Opt_random = nlminb( start=Test$par, obj=Test$fn, grad=Test$gr )
########## END IN-TEXT SNIPPET

########## START IN-TEXT SNIPPET
# Explore inner Hessian
inner_hessian = optimHess( par=Opt_random$par, fn=Test$fn, gr=Test$gr )
# Check eigendecomposition
Eigen = eigen(inner_hessian)
bad_vectors = Eigen$vectors[,which(Eigen$values <= 0)]
data.frame( "Par"=names(Test$par), "MLE"=Test$par, "Status"=ifelse(bad_vectors>0.1,"Bad","Good") )
########## END IN-TEXT SNIPPET

#############
# Intentionally cause outer-Hessian problem
#############

# Introduce bug
buggy_params = params
buggy_params$ln_mu = c(0,0)
Obj = MakeADFun( data=data, parameters=buggy_params, random="eps_i" )

Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
SD = sdreport( Obj )
capture.output( SD, file="sdreport_error.txt")

# Display error
sink( file="Outer_hessian_error.txt", append = FALSE, type = "output")
geterrmessage()
sink( file=NULL )

outer_hessian = optimHess( par=Obj$par, fn=Obj$fn, gr=Obj$gr )
Cov = solve( outer_hessian )
