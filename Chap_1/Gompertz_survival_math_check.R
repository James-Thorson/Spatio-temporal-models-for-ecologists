
# Gompertz mortality function
m = function(a, b, k) b * exp(k*a)
# Claimed survival-at-age given Gompertz mortality
S = function(a, b, k) exp(b/k*(1-exp(k*a)))
# Symbolic derivative for derivative
dS = function(a, b, k){
  Expr = D(expression(exp(b/k*(1-exp(k*a)))), "a")
  eval(Expr)
}

# dS/da = -m * S
# so dS/da + m * S should be zero for any A
A = runif(1,0,100)
b = exp(runif(1)) # Initial annual mortality rate
k = exp(runif(1))   # Acceleration in mortality rate
dS(a=A, b=b, k=k) + m(a=A, b=b, k=k) * S(a=A, b=b, k=k)

