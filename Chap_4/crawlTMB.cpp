#include <TMB.hpp>
#include "make_covariance.hpp"
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_MATRIX( y_iz );       // NA for missing/dropped data
  DATA_VECTOR( DeltaT_i ); // NA for first obs of each track
  DATA_VECTOR( error2_i );
  DATA_MATRIX( X_ij );
  DATA_INTEGER( n_factors );
  DATA_IMATRIX( RAM );

  // Parameters
  PARAMETER_VECTOR( sigma2_z );
  PARAMETER_MATRIX( x_iz );
  PARAMETER_MATRIX( beta_jz );

  // Objective funcction
  int n_i = y_iz.rows();
  int n_z = y_iz.cols();
  Type jnll = 0;
  matrix<Type> I_zz(n_z,n_z);
  I_zz.setIdentity();
  matrix<Type> V_zz(n_z,n_z);
  matrix<Type> S_zz(n_z,n_z);
  matrix<Type> gamma_iz = X_ij * beta_jz;
  matrix<Type> Vi_zz(n_z,n_z);

  // Define process variance
  V_zz = make_covariance( sigma2_z, RAM, sigma2_z, n_z, n_factors );

  // Probability of random coefficients
  for( int i=1; i<n_i; i++){
    if( !R_IsNA(asDouble(DeltaT_i(i))) ){
      Vi_zz = DeltaT_i(i) * V_zz;
      jnll += density::MVNORM(Vi_zz)( x_iz.row(i) - (x_iz.row(i-1)+gamma_iz.row(i-1)) );
    }
  }

  // Probability of data conditional on fixed and random effect values
  for( int i=0; i<n_i; i++){
    S_zz = I_zz * error2_i(i);
    if( !R_IsNA(asDouble(y_iz(i,0))) ){
      jnll += density::MVNORM(S_zz)( y_iz.row(i)-x_iz.row(i) );
    }
  }

  // Cumulative sum of covariates
  matrix<Type> Gsum_iz( n_i, n_z );
  Gsum_iz.row(0) = gamma_iz.row(0) + x_iz.row(0);
  for( int i=1; i<n_i; i++){
    Gsum_iz.row(i) = Gsum_iz.row(i-1) + gamma_iz.row(i);
  }

  // Reporting
  REPORT(gamma_iz);
  ADREPORT(Gsum_iz);
  REPORT( V_zz );
  return jnll;
}
