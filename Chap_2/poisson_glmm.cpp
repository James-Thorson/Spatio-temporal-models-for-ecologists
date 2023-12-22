#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR( y_i )

  // Parameters
  PARAMETER_VECTOR( eps_i );
  PARAMETER_VECTOR( ln_mu );
  PARAMETER( ln_sd );

  // Global variables
  Type jnll = 0;

  // Probability of data conditional on fixed and random effect values
  vector<Type> yhat_i(y_i.size());
  for( int i=0; i<y_i.size(); i++){
    yhat_i(i) = exp( eps_i(i) );
    jnll -= dpois( y_i(i), yhat_i(i), true );
  }

  // Probability of random effects
  for( int i=0; i<y_i.size(); i++){
    jnll -= dnorm( eps_i(i), ln_mu(0), exp(ln_sd), true );
  }
  Type yhat_sum = yhat_i.sum();

  // Return values to R
  REPORT( jnll );  // Used in shrinkage demo
  REPORT( yhat_i );
  REPORT( yhat_sum );
  ADREPORT( yhat_i );
  ADREPORT( yhat_sum );
  return jnll;
}
