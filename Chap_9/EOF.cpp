#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;

  // Data
  DATA_VECTOR( y_i );  
  DATA_IVECTOR( t_i );

  // SPDE objects
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);

  // Projection matrices
  DATA_SPARSE_MATRIX(A_is);
  DATA_SPARSE_MATRIX(A_gs);
  DATA_VECTOR( a_g );

  // Parameters
  PARAMETER( beta0 );
  PARAMETER( ln_tau );
  PARAMETER( ln_kappa );
  PARAMETER( ln_sigma );
  PARAMETER_MATRIX( L_ft );

  // Random effects
  PARAMETER_VECTOR( omega_s );
  PARAMETER_MATRIX( epsilon_sf );

  // Probability of random effects
  Type jnll = 0;
  Eigen::SparseMatrix<Type> Q = exp(2*ln_tau) * (exp(4*ln_kappa)*M0 + Type(2.0)*exp(2*ln_kappa)*M1 + M2);
  jnll += GMRF(Q)( omega_s );
  for( int j=0; j<epsilon_sf.cols(); j++ ){
    jnll += GMRF(Q)( epsilon_sf.col(j) );
  }
  
  // Projection to data
  matrix<Type> epsilon_it( y_i.size(), L_ft.cols() );
  epsilon_it = A_is * epsilon_sf * L_ft;
  vector<Type> omega_i( y_i.size() );
  omega_i = A_is * omega_s.matrix();

  // Probability of data conditional on random effects
  for( int i=0; i<y_i.size(); i++){
    jnll -= dnorm( y_i(i), beta0 + omega_i(i) + epsilon_it(i,t_i(i)), exp(ln_sigma), true );
  }

  // Extrapolation
  matrix<Type> epsilon_gf( A_gs.rows(), epsilon_sf.cols() );
  epsilon_gf = A_gs * epsilon_sf;
  matrix<Type> epsilon_gt( A_gs.rows(), L_ft.cols() );
  epsilon_gt = epsilon_gf * L_ft;
  vector<Type> omega_g( A_gs.rows() );
  omega_g = A_gs * omega_s.matrix();
  
  // Prediction
  matrix<Type> yhat_gt( epsilon_gt.rows(), epsilon_gt.cols() );
  vector<Type> a_t( epsilon_gt.cols() );
  a_t.setZero();
  for( int t=0; t<epsilon_gt.cols(); t++ ){
  for( int g=0; g<epsilon_gt.rows(); g++ ){
    yhat_gt(g,t) = beta0 + omega_g(g) + epsilon_gt(g,t);
    a_t(t) += yhat_gt(g,t) * a_g(g);
  }}

  // Reporting
  REPORT( omega_g );
  REPORT( epsilon_gf );
  REPORT( yhat_gt );
  REPORT( a_t );
  return jnll;
}