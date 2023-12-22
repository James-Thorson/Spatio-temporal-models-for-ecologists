#include <TMB.hpp>
#include "make_M_dense.hpp"
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  //using namespace sparse_matrix_exponential;

  // Configuration
  DATA_INTEGER( CTMC_version );
  DATA_SCALAR( DeltaD );

  // Indexing
  DATA_INTEGER( n_t );
  DATA_INTEGER( n_s );
  DATA_VECTOR( c_i );  // counts for observation i
  DATA_IVECTOR( s_i ); //
  DATA_IVECTOR( t_i ); //

  // Covariates
  DATA_MATRIX( X_sz );

  // Constructing Mrate
  DATA_VECTOR( colsumA_s );
  DATA_IMATRIX( At_zz ); // Indices of A_ss

  // Constructing Q
  DATA_VECTOR( rho_bounds );
  DATA_SPARSE_MATRIX(I_ss);
  DATA_SPARSE_MATRIX(A_ss);

  // Parameters and random effects
  PARAMETER( rho_prime );  // CAR decorrelation rate
  PARAMETER( ln_sigmaO );  // log-SD for CAR
  PARAMETER( ln_sigmaB );  // log-SD for beta
  PARAMETER( ln_D );       // log-diffusion (fixed)
  PARAMETER_VECTOR( gamma_z ); // covariate responses
  PARAMETER_VECTOR( beta_t );  // Changes in abudnance
  PARAMETER_MATRIX( ln_D_st ); // log-density matrix

  // Global variables
  Type rho = invlogit(rho_prime)*(rho_bounds(1)-rho_bounds(0)) + rho_bounds(0);
  Type jnll = 0;
  matrix<Type> ln_Dhat_st( n_s, n_t );

  // Construct CAR precision
  Eigen::SparseMatrix<Type> Q_ss = (I_ss - rho*A_ss) / exp(2.0 * ln_sigmaO);

  // Define preference function
  vector<Type> h_s(n_s);
  h_s = X_sz * gamma_z;

  // Assemble movement
  matrix<Type> Mrate_ss( n_s, n_s );
  Mrate_ss = make_M( CTMC_version, n_s, DeltaD, At_zz, ln_D, h_s, colsumA_s );
  matrix<Type> M_ss( n_s, n_s );
  M_ss = expm(Mrate_ss);

  // Initial density
  vector<Type> tmp_s( n_s );
  for( int s=0; s<n_s; s++ ){
    ln_Dhat_st(s,0) = beta_t(0);
  }
  jnll += GMRF(Q_ss)( ln_D_st.col(0) - ln_Dhat_st.col(0) );

  // Project density forward
  for( int t=1; t<n_t; t++ ){
    for( int s=0; s<n_s; s++ ){
      tmp_s(s) = exp( ln_D_st(s,t-1) );
    }
    tmp_s = tmp_s.matrix().transpose() * M_ss;
    for( int s=0; s<n_s; s++ ){
      ln_Dhat_st(s,t) = beta_t(t) + log(tmp_s(s));
    }
    jnll += GMRF(Q_ss)( ln_D_st.col(t) - ln_Dhat_st.col(t) );
  }

  // Project density without process errors
  matrix<Type> proj_st( n_s, n_t );
  proj_st.col(0) = ln_D_st.col(0);
  for( int t=1; t<n_t; t++ ){
    proj_st.col(t) = proj_st.col(t-1).matrix().transpose() * M_ss;
  }
  
  // Probability of random effects
  for( int t=1; t<n_t; t++ ){
    jnll -= dnorm( beta_t(t), beta_t(t-1), exp(ln_sigmaB), true );
  }

  // Probability of data conditional on random effects
  for( int i=0; i<c_i.size(); i++){
    jnll -= dpois( c_i(i), exp(ln_D_st(s_i(i),t_i(i))), true );   
  }

  // Reporting
  REPORT( rho );
  REPORT( Mrate_ss );
  REPORT( M_ss );
  REPORT( h_s );
  REPORT( ln_Dhat_st );
  REPORT( proj_st );
  return jnll;
}
