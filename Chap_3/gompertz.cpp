#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR( log_b_t );
  DATA_VECTOR( log_bnew_z );
  DATA_VECTOR( simulate_t );

  // Parameters
  PARAMETER( log_d0 );
  PARAMETER( log_sigmaP );
  PARAMETER( log_sigmaM );
  PARAMETER( alpha );
  PARAMETER( rho );
  PARAMETER_VECTOR( log_d_t );

  // Objective funcction
  Type jnll = 0;

  // Probability of random coefficients
  jnll -= dnorm( log_d_t(0), log_d0, exp(log_sigmaP), true );
  for( int t=1; t<log_b_t.size(); t++){
    if( simulate_t(t) == 1 ){
      SIMULATE{
        log_d_t(t) = rnorm( alpha + rho*log_d_t(t-1), exp(log_sigmaP) );
      }
    }
    jnll -= dnorm( log_d_t(t), alpha + rho*log_d_t(t-1), exp(log_sigmaP), true );
  }

  // Probability of data conditional on fixed and random effect values
  for( int t=0; t<log_b_t.size(); t++){
    if( !R_IsNA(asDouble(log_b_t(t))) ){
      jnll -= dnorm( log_b_t(t), log_d_t(t), exp(log_sigmaM), true );
    }
  }

  // Predicted production function
  vector<Type> log_out_z( log_bnew_z.size() );
  for( int t=0; t<log_bnew_z.size(); t++){
    log_out_z(t) = alpha + rho * log_bnew_z(t);
  }
  ADREPORT( log_out_z );
  SIMULATE{ REPORT( log_d_t ); }

  // Reporting
  return jnll;
}
