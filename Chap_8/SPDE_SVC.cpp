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

  // Parameters
  PARAMETER( ln_tau );
  PARAMETER( ln_kappa );
  PARAMETER( ln_sigma );
  PARAMETER_VECTOR( beta_j );

  // Random effects
  PARAMETER_ARRAY( xi_sj );

  // Objective funcction
  Type jnll = 0;

  // Derived quantities
  Type Range = sqrt(8) / exp( ln_kappa );
  Type SigmaE = 1 / sqrt(4 * M_PI * exp(2*ln_tau) * exp(2*ln_kappa));

  // Probability of random effects
  Eigen::SparseMatrix<Type> Q = exp(4*ln_kappa)*M0 + Type(2.0)*exp(2*ln_kappa)*M1 + M2;
  for( int j=0; j<xi_sj.cols(); j++ ){
    jnll += SCALE( GMRF(Q), 1/exp(ln_tau) )( xi_sj.col(j) );
  }

  // Probability of data conditional on random effects
  matrix<Type> xi_ij( A_is.rows(), xi_sj.cols() );
  xi_ij = A_is * xi_sj.matrix();
  vector<Type> phi_i( X_ij.rows() );
  phi_i = X_ij * beta_j;
  vector<Type> yhat_i( y_i.size() );
  //array<Type> psi_ij( A_is.rows(), xi_sj.cols() );
  //psi_ij = xi_ij.array() * X_ij.array();
  for( int i=0; i<y_i.size(); i++){
    yhat_i(i) = exp( (xi_ij.row(i).array() * X_ij.row(i).array()).sum() + phi_i(i) );
    //yhat_i(i) = exp( psi_ij.row(i).vector().sum() + phi_i(i) );
    // shape = 1/CV^2;   scale = mean*CV^2
    jnll -= dgamma( y_i(i), 1/exp(2*ln_sigma), yhat_i(i)*exp(2*ln_sigma), true );
  }

  // Extrapolation
  matrix<Type> xi_gj( A_gs.rows(), xi_sj.cols() );
  xi_gj = A_gs * xi_sj.matrix();
  vector<Type> phi_g( X_gj.rows() );
  phi_g = X_gj * beta_j;
  vector<Type> yhat_g( A_gs.rows() );
  for( int g=0; g<A_gs.rows(); g++){
    yhat_g(g) = exp( (xi_gj.row(g).array() * X_gj.row(g).array()).sum() + phi_g(g) );
  }

  // Reporting
  REPORT( yhat_g );
  REPORT( yhat_i );
  return jnll;
}
