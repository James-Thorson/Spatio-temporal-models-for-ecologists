#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_IVECTOR( options_z );
   // options_z(0) = 0:  eps ~ Normal( mean=mu, sd=sigma )
   // options_z(0) = 1:  eps ~ Gamma( shape=mu, scale=sigma )
   // options_z(1) = 0:  Y = eps
   // options_z(1) = 1:  Y = sqrt(eps)
   // options_z(1) = 2:  Y = exp(eps)
  DATA_SCALAR( mu );
  DATA_SCALAR( sigma );

  // Parameters
  PARAMETER( epsilon );
  PARAMETER_VECTOR( delta );  // calculate gradient of derived quantity Z

  // Define distribution for random variable epsilon
  Type jnll=0, Z=0;
  if(options_z(0)==0) jnll = -1 * dnorm( epsilon, mu, sigma, true );
  if(options_z(0)==1) jnll = -1 * dgamma( epsilon, mu, sigma, true );

  // Define derived quantity, Y = phi(X)
  if(options_z(1)==0) Z = epsilon;
  if(options_z(1)==1) Z = sqrt(epsilon);
  if(options_z(1)==2) Z = exp(epsilon);

  // Add delta*Z for manual epsilon correction
  if( delta.size() > 0 ){
    jnll += delta(0) * Z;
  }

  // Sample-based estimator for true distribution
  SIMULATE{
    Type Zmean_sampled = 0;
    for( int i=0; i<10000; i++ ){
      if(options_z(0)==0) epsilon = rnorm( mu, sigma );
      if(options_z(0)==1) epsilon = rgamma( mu, sigma );
      if(options_z(1)==0) Zmean_sampled += epsilon / 10000;
      if(options_z(1)==1) Zmean_sampled += sqrt(epsilon) / 10000;
      if(options_z(1)==2) Zmean_sampled += exp(epsilon) / 10000;
    }
    REPORT( Zmean_sampled );
  }

  REPORT( Z );
  ADREPORT( Z ); // Built-in TMB SE and epsilon correction
  return jnll;
}
