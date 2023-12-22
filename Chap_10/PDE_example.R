set.seed(1234)

################# BEGIN IN-TEST SNIPPET
# Define parameters and randomize locations
n <- 2  # dimensions (2 = spatial coordinates)
var <- runif(1) # diffusion rate for each dimension
Sigma <- var * diag(n)  # full diffusion matrix
s0 <- rnorm(n) # randomized location to evaluate
t0 <- runif(1) # randomized time to evaluate

## Fundamental solution
Density <- function(s,t) mvtnorm::dmvnorm(s, mean=rep(0,n), sigma=Sigma*t)

## diffusion -- time derivative
numDeriv::grad( function(t) Density(s0,t), t0 )
# output:  0.3069884

# diffusion -- trace of Hessian
# NOTE: only works for diagonal covariance
H = numDeriv::hessian( function(s) Density(s,t0), s0 )
var/2 * sum(diag(H)) # trace of matrix
# output:  0.3069884
################# BEGIN IN-TEST SNIPPET

## diffusion equation -- space Laplacian derivative
# NOTE: works for non-diagonal covariance
scaled_grad <- function(s) Sigma %*% numDeriv::grad(function(s) Density(s,t0), s)
H = numDeriv::jacobian( func=scaled_grad, x=s0 )  # gradient of scaled gradient
0.5 * sum(diag(H)) # trace of matrix

