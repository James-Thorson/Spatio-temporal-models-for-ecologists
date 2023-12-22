# settings
N = 20
q = 0.8
n_count = 10

# simulate counts
C = rbinom(n=n_count, size=N, prob=q)

# define negative log-marginal likelihood
jnll = function(q, C, Nmax=200){
  prob_n = NULL
  for( N in max(C):Nmax ){
    prob_n = c( prob_n, prod(dbinom(x=C, size=N, prob=q, log=FALSE)) )
  }
  return( -1*log(sum(prob_n)) )
}

# optimize marginal likelihood
optimize( f=jnll, interval=c(0.001,0.999), C=C )$minimum