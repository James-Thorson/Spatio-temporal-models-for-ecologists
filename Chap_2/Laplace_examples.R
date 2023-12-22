
chi2 = function(x, k, order=0, log=FALSE ){
  if(order==0 & log==FALSE) out = 1/(2^(k/2)*gamma(k/2)) * x^(k/2-1) * exp(-x/2)
  if(order==0 & log==TRUE) out = log(1/(2^(k/2)*gamma(k/2)) * x^(k/2-1) * exp(-x/2))
  if(order==1) out = (k/2-1)/x - 0.5
  if(order==2) out = -(k/2-1)/x^2
  return(out)
}
# \int N(x|0,\sigma^2) = dN(x|0,\sigma^2) * sigma\sqrt(2\pi)
#                        Sigma = 4; dnorm(0,mean=0,sd=Sigma) * sqrt(2*pi)*Sigma
# \int MVN(x|0,V) = dMVN(x|0,V) * det(V)\sqrt(2\pi)
#                   mvtnorm::dmvnorm(rep(0,3),sigma=2*diag(3)) * sqrt(det(2*diag(3)))*(2*pi)^(3/2)

# Examples of Laplace approximation
setwd("C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_2/")

############ START IN-TEXT CODE
# Calculate derivatives for log-density of chi-square distribution
log_chi2 = function(x, k, order=0 ){
  # Define log of chi-squared distribution
  Expr = expression( log(1/(2^(k/2)*gamma(k/2)) * x^(k/2-1) * exp(-x/2)) )
  # Calculate symbolic derivatives
  if(order>=1) Expr = D(Expr,"x")
  if(order>=2) Expr = D(Expr,"x")
  out = eval(Expr)
}

# Define function for plotting
plot_laplace = function( k, xlim=c(0,k*3), ... ){
  Opt = optimize( log_chi2, interval=xlim, k=k, order=0, maximum=TRUE )
  H = -1 * log_chi2( x=Opt$maximum, k=k, order=2 )
  x = seq(xlim[1],xlim[2],length=1000)
  ytrue = exp(log_chi2(x,k=k))
  plot( x=x, y=ytrue, type="l", main=paste0("k = ", k), lwd=2, ... )
  ynorm = dnorm( x, mean=Opt$maximum, sd=sqrt(1/H) )
    ynorm = ynorm / max(ynorm) * exp(Opt$objective)
  #Equal to:  ynorm = exp(Opt$objective) * exp( -0.5*(x-Opt$maximum)^2*H )
  lines( x=x, y=ynorm, col="blue", lwd=2 )
  dens_laplace = exp( log(sqrt(2*pi)) + Opt$objective - 0.5*log(H) )
  legend("topright", bty="n", legend=signif(dens_laplace,3), text.col="blue")
}

# Plot for three chi-squared distributions
par(mfrow=c(1,3), mar=c(3,3,3,1), mgp=c(2,0.5,0), xaxs="i")
plot_laplace(k = 5)
plot_laplace(k = 25)
plot_laplace(k = 100)
# END IN-TEXT CODE

# Plot
png( file="Laplace_demo.png", width=7, height=3, res=200, unit="in")
  par(mfrow=c(1,3), mar=c(3,3,3,1), mgp=c(2,0.5,0), xaxs="i")
  plot_laplace(k = 5)
  plot_laplace(k = 25)
  plot_laplace(k = 100)
dev.off()

##############
# Compare with rejection sampling
##############

simulate_chi2 = function(k, xlim, max_density, ...){
  while(TRUE){
    X = runif(n=1, min=xlim[1], max=xlim[2])
    dens = exp(log_chi2( x=X, k=k, order=0))
    rand = runif(n=1,min=0,max=max_density)
    if(rand<dens) return(X)
  }
}

plot_samples = function( k, xlim=c(0,k*3), xrange=c(0,k*5), ... ){
  Opt = optimize( log_chi2, interval=xrange, k=k, order=0, maximum=TRUE )
  H = -1 * log_chi2( x=Opt$maximum, k=k, order=2 )
  x = seq(xlim[1],xlim[2],length=1000)
  ytrue = exp(log_chi2(x,k=k))
  # rejection samples
  samples = sapply( 1:1e3, FUN=simulate_chi2, k=k, xlim=xrange, max_density=exp(Opt$objective) )
  Hist = hist(samples, xlim=xlim, freq=FALSE, main=paste0("k = ", k), ...)
  box()
  lines( x=x, y=ytrue, type="l", main=paste0("k = ", k), lwd=2 )
  legend( "topright", bty="n", legend=signif(mean(samples),3), text.col="black" )
}

# Plot
png( file="Rejection_sampling_demo.png", width=7, height=3, res=200, unit="in")
  par(mfrow=c(1,3), mar=c(3,3,3,1), mgp=c(2,0.5,0), xaxs="i")
  plot_samples(k = 5)
  plot_samples(k = 25)
  plot_samples(k = 100)
dev.off()

