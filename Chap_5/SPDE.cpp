
#include <TMB.hpp>

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;

  // Data
  DATA_VECTOR( c_i );  // counts for observation i

  // SPDE objects
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);

  // Projection matrices
  DATA_SPARSE_MATRIX(A_is);
  DATA_SPARSE_MATRIX(A_gs);

  // Parameters
  PARAMETER( beta0 );
  PARAMETER( ln_tau );
  PARAMETER( ln_kappa );

  // Random effects
  PARAMETER_VECTOR( omega_s );

  // Objective funcction
  Type jnll = 0;

  // Derived quantities
  Type Range = sqrt(8) / exp( ln_kappa );
  Type SigmaE = 1 / sqrt(4 * M_PI * exp(2*ln_tau) * exp(2*ln_kappa));

  ///////// START IN-LINE CODE
  // Probability of random effects
  Eigen::SparseMatrix<Type> Q = (exp(4*ln_kappa)*M0 + Type(2.0)*exp(2*ln_kappa)*M1 + M2) * exp(2*ln_tau);
  jnll += GMRF(Q)( omega_s );

  // Project using bilinear interpolation
  vector<Type> omega_i( A_is.rows() );
  omega_i = A_is * omega_s;
  ///////// END IN-LINE CODE

  // Probability of data conditional on random effects
  for( int i=0; i<c_i.size(); i++){
    jnll -= dpois( c_i(i), exp(beta0 + omega_i(i)), true );
  }

  // Extrapolation
  vector<Type> logN_g( A_gs.rows() );
  logN_g = beta0 + A_gs*omega_s;

  // Reporting
  REPORT( Q );
  REPORT( Range );
  REPORT( SigmaE );
  REPORT( logN_g );
  ADREPORT( SigmaE );

  return jnll;
}
