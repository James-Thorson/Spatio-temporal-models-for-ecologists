#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;

  // Data
  DATA_VECTOR( c_i );  // counts for observation i
  DATA_IVECTOR( s_i ); //
  DATA_VECTOR( rho_bounds );
  DATA_MATRIX( X_sk );
  DATA_SPARSE_MATRIX(I_ss);
  DATA_SPARSE_MATRIX(A_ss);

  // Parameters and random effects
  PARAMETER( beta0 );
  PARAMETER( rho_prime );
  PARAMETER( ln_sigma );
  PARAMETER_VECTOR( gamma_k );
  PARAMETER_VECTOR( omega_s );

  // Global variables
  Type rho = invlogit(rho_prime)*(rho_bounds(1)-rho_bounds(0)) + rho_bounds(0);
  Type jnll = 0;
  vector<Type> lambda_s( omega_s.size() );
  Eigen::SparseMatrix<Type> Q_ss = (I_ss - rho*A_ss) / exp(2 * ln_sigma);

  // Probability of random effects
  jnll += GMRF(Q_ss)( omega_s );

  // Probability of data conditional on random effects
  lambda_s = exp( beta0 + omega_s + X_sk * gamma_k );
  for( int i=0; i<c_i.size(); i++){
    jnll -= dpois( c_i(i), lambda_s(s_i(i)), true );
  }
  Type sumlambda = lambda_s.sum();

  // Reporting
  REPORT( Q_ss );
  REPORT( rho );
  REPORT( sumlambda );
  ADREPORT( sumlambda );
  return jnll;
}
