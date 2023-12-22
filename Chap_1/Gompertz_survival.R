
setwd("C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_1")
set.seed(101)

########## START IN-TEXT SNIPPET
# Calculate lifetime density function f = d/da( 1 - S(a) )
lifetime_density = function(a, b, k){
  # Write expression for 1 - S(a)
  gompertz_mortality = expression( 1-exp(b/k*(1-exp(k*a))) )
  # Symbolic derivative with respect to age a
  density_function = D(gompertz_mortality, "a")
  # Evaluate derivative using specified values
  eval(Expr)
}
########## END IN-TEXT SNIPPET

########## START IN-TEXT SNIPPET
# Paraeters
b = 0.001 # Initial annual mortality rate
k = 0.1   # Acceleration in mortality rate
birth_rate = 0.2 # animals per year
n_t = 100 # Time domain  to simulate

# Calculate peak of lifetime density for use in rejection sampling
Opt = optimize(f=lifetime_density, interval=c(0,100), maximum=TRUE, b=b, k=k)

# Rejection sampling for age at death
simulate_lifespan = function(b,k){
  while(TRUE){ # Repeat as long as necessary
    A = runif(n=1,min=0,max=150)
    M = lifetime_density(a=A, b=b, k=k)
    rand = runif(n=1,min=0,max=Opt$maximum)
    if(rand<M) return(A) # If accepted, save value
  }
}

# Sample births using sequence of exponential distributions
birth = death = vector()
while( max(0,birth,na.rm=TRUE)<n_t ){
  birth = c(birth, max(0,birth,na.rm=TRUE) + rexp(n=1, rate=birth_rate) )
}

# Sample deaths from lifespace after birth
for( i in seq_along(birth) ){
  death[i] = birth[i] + simulate_lifespan( b=b, k=k )
}
########## END IN-TEXT SNIPPET

# Define
survival_probability = function(a, b, k) exp(b/k*(1-exp(k*a)))

# Plot
png( "Gompertz_survival.png", width=6, height=2, res=200, units="in" )
  par( mfrow=c(1,3), xaxs="i", yaxs="i", mar=c(3,3,1,1), mgp=c(2,0.5,0) )
  x = seq(0,n_t,length=1000)
  Y = cbind( survival_probability(a=x, b=b, k=k),
             lifetime_density(a=x, b=b, k=k) )
  matplot( x=x, y=Y, type="l", lwd=2, xlab="Age", ylab="Proportion", lty="solid" )
  matplot( x=cbind(birth,death), y=cbind(seq_along(birth),seq_along(birth)), type="n",
    xlim=c(0,100), ylim=c(0,sum(birth<n_t)), xlab="Year", ylab="Individual timeline" )
  for( i in seq_along(birth) ){
    lines( x=c(birth[i],death[i]), y=c(i,i), lwd=2)
  }
  y = sapply(x, FUN=function(t){sum(t>birth&t<death)})
  plot(x=x,y=y, type="l", lwd=2, ylim=c(0,1.05*max(y)), xlab="Year", ylab="Population size" )
dev.off()

