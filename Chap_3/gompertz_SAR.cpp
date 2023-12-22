#include <TMB.hpp>
template<class Type>
Eigen::SparseMatrix<Type> Q_ar1( Type rho, int n_t ){
  // Compute SparseMatrix precision
  Eigen::SparseMatrix<Type> Q( n_t, n_t );
  for(int t=0; t<n_t; t++){
    Q.coeffRef( t, t ) = 1 + pow(rho,2);
    if(t>=1) Q.coeffRef( t-1, t ) = -rho;
    if(t>=1) Q.coeffRef( t, t-1 ) = -rho;
  }
  Q /= (1 - pow(rho,2));
  return Q;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR( log_b_t );
  DATA_VECTOR( log_bnew_z );

  // Parameters
  PARAMETER( log_delta );
  PARAMETER( log_sigmaP );
  PARAMETER( log_sigmaM );
  PARAMETER( alpha );
  PARAMETER( rho );
  PARAMETER_VECTOR( eps_t );
  
  // Global variables
  Type jnll = 0;
  int n_t = log_b_t.size();
  Eigen::SparseMatrix<Type> Q( n_t, n_t );
  Q = Q_ar1( rho, n_t );

  // Probability of random coefficients
  jnll += density::GMRF(Q)( eps_t );

  // Probability of data conditional on fixed and random effect values
  vector<Type> log_d_t( n_t );
  for( int t=0; t<n_t; t++){
    log_d_t(t) = log_delta*pow(rho, t) + alpha/(1-rho) + exp(log_sigmaP)*eps_t(t);
    jnll -= dnorm( log_b_t(t), log_d_t(t), exp(log_sigmaM), true );
  }
  ADREPORT(log_d_t)

  // Predicted production function
  vector<Type> log_out_z( log_bnew_z.size() );
  for( int t=0; t<log_bnew_z.size(); t++){
    log_out_z(t) = alpha + rho * log_bnew_z(t);
  }
  ADREPORT( log_out_z );

  // Reporting
  return jnll;
}
