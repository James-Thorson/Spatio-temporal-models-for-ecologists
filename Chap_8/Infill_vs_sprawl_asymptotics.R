
setwd( R'(C:\Users\James.Thorson\Desktop\Git\Spatio-temporal-models-for-ecologists\Chap_8)' )
set.seed(101)

########## START IN-TEXT SNIPPET
library(geoR)
# Simulation settings
area_c = c("original"=1, "infill"=1, "sprawl"=10)
nsamp_c = c("original"=50, "infill"=500, "sprawl"=500)
nrep = 100
cov.pars = c("sigma2"=1, "range"=0.4)
results_rcz = array(NA, dim=c(nrep,length(area_c),length(cov.pars)),
                        dimnames=list(NULL,names(area_c),names(cov.pars)) )

# Run experiment
for( config in seq_along(area_c) ){
for( rep in 1:nrep ){
  loc = matrix( runif(n=2*nsamp_c[config], max=sqrt(area_c[config])), ncol=2)
  # Simulate a spatial field
  field = grf( grid = loc,
               cov.pars = cov.pars,
               cov.model= "exponential" )
  # Construct variogram
  var = variog( field )
  # Fit variogram model
  varhat = variofit( var,
                  cov.model= "exponential",
                  fix.nugget = TRUE,
                  limits = geoR::pars.limits("phi"=c(lower=0.01, upper=Inf)),
                  ini.cov.pars = cov.pars )
  # Compile results
  results_rcz[rep,config,] = varhat$cov.pars
}}
########## END IN-TEXT SNIPPET

# Plot
png( "Infill_and_sprawl_asymptotics.png", width=6, height=3, res=200, units='in')
  par( mfrow=c(1,3), mar=c(2,2,2,0), mgp=c(2,0.5,0), oma=c(2,2,0,1) )
  for( config in 1:3 ){
    xplot = seq(0, sqrt(2*max(area_c[config])), length=1000)
    plot( 0, type="n", xlim=c(0,1.4), ylim=c(0,3*cov.pars[1]), ylab="Predicted semivariance", xlab="Distance", yaxt="n", xaxs="i", yaxs='i' )
    for(i in 1:nrep) lines( x=xplot, y=results_rcz[i,config,1]*(1-exp(-xplot/results_rcz[i,config,2])), col=rgb(0,0,0,0.3) )
    lines( x=xplot, y=cov.pars[1]*(1-exp(-xplot/cov.pars[2])), lwd=3, col="blue" )
    legend( "topleft", bty="n", title="Percent beyond domain:", legend=sum(results_rcz[,config,2]>sqrt(2*max(area_c[1]))) )
    if(config%in%c(1,3)) axis(2)
    title( paste0(nsamp_c[config]," samples in ",area_c[config]," km2"))
  }
  mtext( side=1:2, text=c("Distance (km)","Predicted semivariance"), outer=TRUE )
dev.off()
