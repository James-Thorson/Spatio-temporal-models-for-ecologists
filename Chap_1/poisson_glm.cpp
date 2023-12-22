#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR( y_i )
  DATA_MATRIX( X_ij );

  // Parameters
  PARAMETER_VECTOR( beta_j );

  // Global variables
  Type jnll = 0;
  int n_i = X_ij.rows();
  int n_j = X_ij.cols();
  vector<Type> log_mu(n_i);
  log_mu.setZero();

  // Probability of data conditional on fixed and random effect values
  for( int i=0; i<n_i; i++){
    for( int j=0; j<n_j; j++){
      log_mu(i) += beta_j(j) * X_ij(i,j);
    }
    jnll -= dpois( y_i(i), exp(log_mu(i)), true );
    // Option to simulate new data given parameters
    SIMULATE{
      y_i(i) = rpois( exp(log_mu(i)) );
    }
  }

  // Return values to R
  SIMULATE{REPORT(y_i);} // If simulating new data, hand the value back to R
  REPORT( log_mu );
  return jnll;
}
