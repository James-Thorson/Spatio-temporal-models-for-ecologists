#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;

  // Data
  DATA_VECTOR( n_i );  // counts for observation i
  DATA_IVECTOR( c_i ); // category
  DATA_MATRIX( C_gc ); // Phylogenetic design matrix (g)
  DATA_MATRIX( T_hc ); // Traits (h)
  DATA_MATRIX( X_gk ); // Environment (k)
  DATA_MATRIX( X_ik ); // Environment (k)

  // SPDE objects
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);

  // Projection matrices
  DATA_SPARSE_MATRIX(A_is);
  DATA_SPARSE_MATRIX(A_gs);

  // Parameters and random effects
  PARAMETER( ln_kappa );      // Spatial decorrelation
  PARAMETER( ln_sigmaM );  // log-SD for overdisperion
  PARAMETER_VECTOR( ln_sigmaC );
  PARAMETER_VECTOR( ln_sigmaT_h );
  PARAMETER_VECTOR( alpha_c ); // Intercepts
  PARAMETER_MATRIX( lambda_jc );    // Factors
  PARAMETER_MATRIX( beta_kc ); // X_sk %*% beta_kc
  PARAMETER_MATRIX( gamma_sh ); // gamma_sh %*% T_hc
  PARAMETER_MATRIX( delta_sg ); // delta_sg %*% C_gc
  PARAMETER_MATRIX( omega_sj ); // omega_sj %*% L_jc
  //PARAMETER_VECTOR( eps_i ); // Overdispersion

  // Global variables
  Type jnll = 0;
  Type tau = 1.0 / sqrt(4.0 * M_PI * exp(2*ln_kappa));
  Eigen::SparseMatrix<Type> Q = tau*tau * (exp(4*ln_kappa)*M0 + Type(2.0)*exp(2*ln_kappa)*M1 + M2);

  // Assemble predictions at data
  matrix<Type> beta_ic = X_ik * beta_kc;
  matrix<Type> gamma_ic = A_is * gamma_sh * T_hc;
  matrix<Type> delta_ic = A_is * delta_sg * C_gc;
  matrix<Type> omega_ic = A_is * omega_sj * lambda_jc;
  matrix<Type> p_ic = beta_ic + gamma_ic + delta_ic + omega_ic;

  // Assemble predictions at extrapolation grid
  matrix<Type> beta_gc = X_gk * beta_kc;
  matrix<Type> gamma_gc = A_gs * gamma_sh * T_hc;
  matrix<Type> delta_gc = A_gs * delta_sg * C_gc;
  matrix<Type> omega_gc = A_gs * omega_sj * lambda_jc;
  matrix<Type> p_gc = beta_gc + gamma_gc + delta_gc + omega_gc;

  // Probability of phylogeny
  for( int g=0; g<delta_sg.cols(); g++ ){
    jnll += SCALE( GMRF(Q), exp(ln_sigmaC(0)) )( delta_sg.col(g) );
  }

  // Probability of traits
  for( int h=0; h<gamma_sh.cols(); h++ ){
    jnll += SCALE( GMRF(Q), exp(ln_sigmaT_h(h)) )( gamma_sh.col(h) );
  }

  // Probability of factors
  for( int j=0; j<omega_sj.cols(); j++ ){
    jnll += GMRF(Q)( omega_sj.col(j) );
  }

  // Probability of overdispersion
  //for( int i=0; i<n_i.size(); i++){
  //  jnll -= dnorm( eps_i(i), Type(0.0), exp(ln_sigmaM) );
  //}

  // Probability of data conditional on random effects
  Type mu;
  for( int i=0; i<n_i.size(); i++){
    //jnll -= dpois( n_i(i), exp(alpha_c(c_i(i)) + p_ic(i,c_i(i)) + eps_i(i)), true );
    mu = exp( alpha_c(c_i(i)) + p_ic(i,c_i(i)) );
    jnll -= dnbinom2( n_i(i), mu, mu * (1.0 + exp(2.0*ln_sigmaM)), true );
  }

  // Reporting
  REPORT( Q );  // Used in variance calculations
  REPORT( beta_gc );
  REPORT( gamma_gc );
  REPORT( delta_gc );
  REPORT( omega_gc );
  REPORT( p_gc );
  ADREPORT( beta_kc );  // Get covariance even if using REML
  return jnll;
}
