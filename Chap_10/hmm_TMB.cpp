#include <TMB.hpp>
#include "make_M.hpp"
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( CTMC_version );
  DATA_INTEGER( expm_version );
  DATA_INTEGER( Nmax );
  DATA_SCALAR( DeltaD );
  DATA_VECTOR( colsumA_g );
  DATA_MATRIX( L_gt );
  DATA_MATRIX( X_gz );
  DATA_IMATRIX( At_zz ); // Indices of A_gg

  // Parameters and random effects
  PARAMETER( ln_D );
  PARAMETER_VECTOR( gamma_z );

  // Global variables
  Type jnll = 0;
  int n_g = L_gt.rows();
  int n_t = L_gt.cols();
  matrix<Type> ForwardProb_gt( n_g, n_t );
  matrix<Type> ForwardPred_gt( n_g, n_t );
  matrix<Type> Forward_gt( n_g, n_t );
  ForwardPred_gt.setZero();
  matrix<Type> BackwardProb_gt( n_g, n_t );
  matrix<Type> BackwardPred_gt( n_g, n_t );
  matrix<Type> Backward_gt( n_g, n_t );
  BackwardPred_gt.setZero();
  
  // Calculate movement matrix
  vector<Type> h_g(n_g);
  h_g = X_gz * gamma_z;
  Eigen::SparseMatrix<Type> Mrate_gg( n_g, n_g );
  Mrate_gg = make_M( CTMC_version, n_g, DeltaD, At_zz, ln_D, h_g, colsumA_g );
  sparse_matrix_exponential::config<Type> myconfig = sparse_matrix_exponential::config<Type>();
  myconfig.Nmax = Nmax;
  sparse_matrix_exponential::expm_generator<Type> M_gg( Mrate_gg, myconfig );  

  // expm_series, see sparse_matrix_exponential.hpp expm_generator code
  // Modify rate matrix to be positive semi-definite for numerical stability
  Type rho = 0; 
  Eigen::SparseMatrix<Type> A_gg(Mrate_gg);
  vector<Type> diag_g = Mrate_gg.diagonal();
  if (diag_g.size() > 0) {
    Type M = diag_g[0];
    for (int i=1; i<diag_g.size(); i++) {
      M = TMBad::min(M, diag_g[i]);
    }
    rho = -M;
  }
  A_gg.diagonal().array() += rho;
  // Forwards algorithm  
  sparse_matrix_exponential::expm_series<Type> Mseries_gg( A_gg, Type(Nmax), myconfig );  
  // Backwards algorithm
  Eigen::SparseMatrix<Type> Aprime_gg( A_gg.transpose() );
  sparse_matrix_exponential::expm_series<Type> Mprimeseries_gg( Aprime_gg, Type(Nmax), myconfig );  
  
  // Project forward
  ForwardProb_gt.col(0) = L_gt.col(0);
  for( int t=1; t<n_t; t++ ){
    // Predict movement (and re-normalize to correct for small numerical issues)
    if( expm_version==0 ){
      ForwardPred_gt.col(t) = M_gg( vector<Type>(ForwardProb_gt.col(t-1)) ).array();
    }else{
      ForwardPred_gt.col(t) = Mseries_gg( vector<Type>(ForwardProb_gt.col(t-1)) ).array();
      ForwardPred_gt.col(t) = exp(-rho) * ForwardPred_gt.col(t);
    }
    ForwardPred_gt.col(t) = ForwardPred_gt.col(t) / ForwardPred_gt.col(t).sum();
    // Calculate probability (and re-normalize to avoid accumulating underflow in likelihood)
    ForwardProb_gt.col(t) = ForwardPred_gt.col(t).array() * L_gt.col(t).array();
    jnll -= log( ForwardProb_gt.col(t).sum() );
    ForwardProb_gt.col(t) = ForwardProb_gt.col(t) / ForwardProb_gt.col(t).sum();
  }
  // Accumulator
  Forward_gt = ForwardProb_gt;

  // Project backwards
  for( int g=0; g<n_g; g++ ){ BackwardPred_gt(g,n_t-1) = 1; }
  for( int t=n_t-2; t>=0; t-- ){
    // Calculate probability (and re-normalize to avoid accumulating underflow in likelihood)
    BackwardProb_gt.col(t+1) = BackwardPred_gt.col(t+1).array() * L_gt.col(t+1).array();
    BackwardProb_gt.col(t+1) = BackwardProb_gt.col(t+1) / BackwardProb_gt.col(t+1).sum();
    // Predict movement (and re-normalize to correct for small numerical issues)
    BackwardPred_gt.col(t) = Mprimeseries_gg( vector<Type>(BackwardProb_gt.col(t+1)) ).array();
    BackwardPred_gt.col(t) = exp(-rho) * BackwardPred_gt.col(t);
    BackwardPred_gt.col(t) = BackwardPred_gt.col(t) / BackwardPred_gt.col(t).sum();
  }
  // Accumulator
  Backward_gt = BackwardPred_gt;

  // Reporting
  REPORT( Mrate_gg );
  REPORT( ForwardProb_gt );
  REPORT( ForwardPred_gt );
  REPORT( Forward_gt );
  REPORT( BackwardProb_gt );
  REPORT( BackwardPred_gt );
  REPORT( Backward_gt );
  REPORT( h_g );
  return jnll;
}
