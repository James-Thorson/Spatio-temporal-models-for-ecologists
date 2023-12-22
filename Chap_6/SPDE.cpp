
#include <TMB.hpp>

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;

  // Data
  DATA_VECTOR( y_i );  // measurement for observation i

  // SPDE objects
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);

  // Projection matrices
  DATA_SPARSE_MATRIX(A_is);
  DATA_SPARSE_MATRIX(A_gs);

  // Covariate design matrices
  DATA_MATRIX( X_ij );
  DATA_MATRIX( X_gj );

  // Expansion terms
  DATA_VECTOR( a_g );
  DATA_VECTOR( d_g );

  // Parameters
  PARAMETER( beta0 );
  PARAMETER( ln_tau );
  PARAMETER( ln_kappa );
  PARAMETER( ln_sigma );
  PARAMETER_VECTOR( beta_j );

  // Random effects
  PARAMETER_VECTOR( omega_s );

  // Objective funcction
  Type jnll = 0;

  // Derived quantities
  Type Range = sqrt(8) / exp( ln_kappa );
  Type SigmaE = 1 / sqrt(4 * M_PI * exp(2*ln_tau) * exp(2*ln_kappa));

  // Probability of random effects
  Eigen::SparseMatrix<Type> Q = exp(4*ln_kappa)*M0 + Type(2.0)*exp(2*ln_kappa)*M1 + M2;
  jnll += SCALE( GMRF(Q), 1/exp(ln_tau) )( omega_s );
  SIMULATE{
    SCALE( GMRF(Q), 1/exp(ln_tau) ).simulate(omega_s);
  }

  // Probability of data conditional on random effects
  vector<Type> omega_i( A_is.rows() );
  omega_i = A_is * omega_s;
  vector<Type> phi_i( X_ij.rows() );
  vector<Type> yhat_i( y_i.size() );
  phi_i = X_ij * beta_j;
  for( int i=0; i<y_i.size(); i++){
    // shape = 1/CV^2;   scale = mean*CV^2
    yhat_i(i) = exp( beta0 + omega_i(i) + phi_i(i) );
    jnll -= dgamma( y_i(i), 1/exp(2*ln_sigma), yhat_i(i)*exp(2*ln_sigma), true );
    SIMULATE{
      y_i(i) = rgamma( 1/exp(2*ln_sigma), yhat_i(i)*exp(2*ln_sigma) );
    }
  }

  // Extrapolation
  vector<Type> yhat_g( A_gs.rows() );
  yhat_g = exp( beta0 + A_gs*omega_s + X_gj*beta_j );
  SIMULATE{
    vector<Type> y_g( A_gs.rows() );
    for( int g=0; g<y_g.size(); g++ ){
      y_g(g) = rgamma( 1/exp(2*ln_sigma), yhat_g(g)*exp(2*ln_sigma) );
    }
    REPORT( y_g );
  }

  //
  Type Y1 = (yhat_i).sum() / yhat_i.size();
  Type Y2 = (yhat_g.array() * a_g.array()).sum() / (a_g.array()).sum();
  Type Y3 = (yhat_g.array() * a_g.array() * d_g.array()).sum() / (a_g.array() * d_g.array()).sum();

  // Reporting
  REPORT( yhat_g );
  REPORT( yhat_i );
  REPORT( Y1 );
  REPORT( Y2 );
  REPORT( Y3 );
  ADREPORT( Y1 );
  ADREPORT( Y2 );
  ADREPORT( Y3 );
  SIMULATE{REPORT( y_i );}
  return jnll;
}
