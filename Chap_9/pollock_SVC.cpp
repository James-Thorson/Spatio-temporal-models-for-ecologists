
#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;

  // Data
  DATA_INTEGER( n_t );
  DATA_VECTOR( a_g );  // Area associated with grid-cell g
  DATA_VECTOR( z_g );  // Location
  DATA_VECTOR( b_i );  // response for observation i
  DATA_VECTOR( a_i );  // area-swept
  DATA_IVECTOR( t_i );  
  DATA_VECTOR( x_t );

  // Projection matrices
  DATA_SPARSE_MATRIX(A_is);
  DATA_SPARSE_MATRIX(A_gs);

  // SPDE objects
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);

  // Parameters
  PARAMETER_VECTOR( beta_t );
  PARAMETER( ln_tauO );
  PARAMETER( ln_tauE );
  PARAMETER( ln_tauX );
  PARAMETER( ln_kappa );
  PARAMETER( ln_phi );
  PARAMETER( logit_rhoE );
  PARAMETER( finv_power );

  // Random effects
  PARAMETER_VECTOR( omega_s );
  PARAMETER_VECTOR( xi_s );
  PARAMETER_MATRIX( epsilon_st );

  // Objective funcction
  vector<Type> jnll_comp(4);
  jnll_comp.setZero();
  int n_i = A_is.rows();
  int n_g = A_gs.rows();
  Type rhoE = invlogit( logit_rhoE );

  // Probability of random effects
  Eigen::SparseMatrix<Type> Q = exp(2*ln_tauE) * (exp(4*ln_kappa)*M0 + Type(2.0)*exp(2*ln_kappa)*M1 + M2);
  jnll_comp(1) += SCALE( GMRF(Q), 1/exp(ln_tauO) )( omega_s );
  jnll_comp(3) += SCALE( GMRF(Q), 1/exp(ln_tauX) )( xi_s );
  for( int t=0; t<n_t; t++){
    if( t==0 ){
      jnll_comp(2) += SCALE( GMRF(Q), 1 / pow( 1.0-pow(rhoE,2), 0.5 ) )( epsilon_st.col(t) );
    }else{
      jnll_comp(2) += GMRF(Q)( epsilon_st.col(t) - rhoE*epsilon_st.col(t-1) );
    }
  }

  // True density and abundance
  vector<Type> omega_g( n_g );
  vector<Type> xi_g( n_g );
  matrix<Type> epsilon_gt( n_g, n_t );
  omega_g = A_gs * omega_s;
  xi_g = A_gs * xi_s;
  epsilon_gt = A_gs * epsilon_st;
  array<Type> ln_d_gt( n_g, n_t );
  for( int t=0; t<n_t; t++){
  for( int g=0; g<n_g; g++){
    ln_d_gt(g,t) = beta_t(t) + omega_g(g) + epsilon_gt(g,t) + xi_g(g)*x_t(t);
  }}

  // Probability of data conditional on random effects
  vector<Type> omega_i( n_i );
  vector<Type> xi_i( n_i );
  matrix<Type> epsilon_it( n_i, n_t );
  omega_i = A_is * omega_s;
  xi_i = A_is * xi_s;
  epsilon_it = A_is * epsilon_st;
  vector<Type> bhat_i( n_i );
  bhat_i.setZero();
  for( int i=0; i<n_i; i++){
    if( !isNA(b_i(i)) ){
      bhat_i(i) = a_i(i) * exp( beta_t(t_i(i)) + omega_i(i) + epsilon_it(i,t_i(i)) + xi_i(i)*x_t(t_i(i)) );
      jnll_comp(0) -= dtweedie( b_i(i), bhat_i(i), exp(ln_phi), Type(1.0)+invlogit(finv_power), true );
    }
  }

  // Objective function
  Type jnll = jnll_comp.sum();

  // Derived quantities
  vector<Type> b_t( n_t );
  b_t.setZero();
  vector<Type> zmean_t( n_t );
  zmean_t.setZero();
  for( int t=0; t<n_t; t++){
    for( int g=0; g<n_g; g++){
      b_t(t) += a_g(g) * exp(ln_d_gt(g,t));
    }
    for( int g=0; g<n_g; g++){
      zmean_t(t) += z_g(g) * a_g(g) * exp(ln_d_gt(g,t)) / b_t(t);
    }
  }

  // Reporting
  REPORT( jnll_comp );
  REPORT( jnll );
  REPORT( ln_d_gt );
  REPORT( bhat_i );
  REPORT( b_t );
  ADREPORT( b_t );
  REPORT( zmean_t );
  REPORT( xi_g );
  REPORT( omega_g );
  ADREPORT( zmean_t );

  return jnll;
}
