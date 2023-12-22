
setwd( R'(C:\Users\James.Thorson\Desktop\Git\Spatio-temporal-models-for-ecologists\Appendix)' )

library(viridisLite)

###############
# Scalar -- Exponential population growth
###############

r = 0.2
b0 = 100
maxT = 3

# True
t = seq(0,maxT,length=100)
b = b0 * exp(r*t)

# Linearize
euler = function( maxT, N, r, b0 ){
  b = rep(NA, N)
  deltaT = maxT / N
  b[1] = b0 + (r*deltaT * b0)
  for( n in 2:N ){
    b[n] = b[n-1] + (r*deltaT * b[n-1])
  }
  return( cbind(t=seq(0,maxT,length=N+1), b=c(b0,b)) )
}

png( file="euler_demo.png", width=4, height=4, res=200, units="in" )
  par( mar=c(3,4,1,1), mgp=c(2,0.5,0), tck=-0.02, xaxs="i", yaxs="i" )
  plot( x=t, y=b, lwd=3, type="l", xlab="Time elapsed", ylab="Abundance under\nexponential growth" )
  lines( euler( maxT=3, N=2, r=r, b0=b0 ), lwd=3, type="b", col=viridis(4)[1] )
  lines( euler( maxT=3, N=4, r=r, b0=b0 ), lwd=3, type="b", col=viridis(4)[2] )
  lines( euler( maxT=3, N=8, r=r, b0=b0 ), lwd=3, type="b", col=viridis(4)[3] )
  lines( euler( maxT=3, N=16, r=r, b0=b0 ), lwd=3, type="b", col=viridis(4)[4] )
  legend( "topleft", ncol=2, bty="n", fill=c("black",viridis(4)), legend=c("true","N=2","N=4","N=8","N=16") )
dev.off()

###############
# Bivariate -- Lotka Volterra
###############

# dxdt = alpha*x + beta*y
# dydt = delta*x + gamma*y

theta = pi / 6
alpha = cos(theta)
beta = -sin(theta)
delta = sin(theta)
gamma = cos(theta)
maxT = 10

expA = matrix( c(alpha,beta,delta,gamma), ncol=2, byrow=TRUE )
A = expm::logm(expA)
z0 = c(1,0.1)

# True
t = seq(0,maxT,length=100)
z = matrix(NA,ncol=2,nrow=length(t))
for( i in seq_along(t)) z[i,] = z0 %*% expm::expm(A*t[i])
matplot(z)

# Linearize
euler = function( maxT, N, A, z0 ){
  z = matrix(NA,ncol=2,nrow=N)
  deltaT = maxT / N
  z[1,] = z0 + (z0 %*% (A*deltaT))
  for( n in 2:N ){
    z[n,] = z[n-1,] + (z[n-1,] %*% (A*deltaT) )
  }
  return( cbind(t=seq(0,maxT,length=N+1), z=rbind(z0,z)) )
}

png( file="euler_demo_2D.png", width=4, height=4, res=200, units="in" )
  par( mar=c(3,4,1,1), mgp=c(2,0.5,0), tck=-0.02, xaxs="i", yaxs="i" )
  matplot( x=t, y=z, lwd=3, type="l", xlab="Time elapsed", ylab="Oscillatory movement\nbetween two states", ylim=range(z)*1.5, lty="solid" )
  out = euler( maxT=maxT, N=100, A=A, z0=z0 )
  matplot( x=out[,1], y=out[,-1], lwd=3, type="b", pch=20, add=TRUE )
  legend( "topleft", ncol=2, bty="n", fill=c("black","red"), legend=c("State-1","State-2") )
dev.off()
