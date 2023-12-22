#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;

  // Data
  DATA_VECTOR( c_i );  // counts for observation i
  DATA_IVECTOR( e_i );  // 0: Encounter;  1: Count;  2=Biomass
  DATA_MATRIX( X_ik );  // Habitat covariates for samples
  DATA_MATRIX( X_gk );  // Habitat covariates for integration points
  DATA_MATRIX( Q_ij );  // Detectability covariates for samples
  DATA_SPARSE_MATRIX(M0);  // SPDE matrix-1
  DATA_SPARSE_MATRIX(M1);  // SPDE matrix-2
  DATA_SPARSE_MATRIX(M2);  // SPDE matrix-3
  DATA_SPARSE_MATRIX(A_is); // Project vertices to samples
  DATA_SPARSE_MATRIX(A_gs); // Project vertices to integration points

  // Parameters
  PARAMETER( ln_tau );
  PARAMETER( ln_kappa );
  PARAMETER( ln_phi );
  PARAMETER( finv_power );
  PARAMETER_VECTOR( gamma_k );
  PARAMETER_VECTOR( eta_j );
  PARAMETER_VECTOR( omega_s );

  // Global variables
  Type jnll = 0;
  Type phi = exp(ln_phi);
  Type power = Type(1.0) + invlogit(finv_power);
  vector<Type> logmu_g = A_gs*omega_s + X_gk*gamma_k;
  vector<Type> logmu_i = A_is*omega_s + Q_ij*eta_j + X_ik*gamma_k;

  // Probability of random effects
  Eigen::SparseMatrix<Type> Q = (exp(4*ln_kappa)*M0 + Type(2.0)*exp(2*ln_kappa)*M1 + M2) * exp(2*ln_tau);
  jnll += GMRF(Q)( omega_s );

  // Likelihood of data
  for( int i=0; i<c_i.size(); i++){
    // Bernoulli
    if(e_i(i)==0){
      if( c_i(i) > 0 ){
        jnll -= logspace_sub( Type(log(1.0)), -1*exp(logmu_i(i)) );
      }else{
        jnll -= -1*exp(logmu_i(i));
      }
    }
    // Poisson
    if(e_i(i)==1){
      jnll -= dpois( c_i(i), exp(logmu_i(i)), true );
    }
    // Tweedie
    if(e_i(i)==2){
      jnll -= dtweedie( c_i(i), exp(logmu_i(i)), phi, power, true );
    }
  }

  // Reporting
  REPORT( logmu_g );
  return jnll;
}
