#include <TMB.hpp>
#include "make_covariance.hpp"
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_MATRIX( y_iz );
  DATA_IMATRIX( RAM );
  DATA_IVECTOR( familycode_z );
  DATA_VECTOR( RAMstart );

  // Parameters
  PARAMETER_MATRIX( x_iz );
  PARAMETER_VECTOR( beta_j );

  // Objective funcction
  int n_i = y_iz.rows();
  int n_z = y_iz.cols();
  matrix<Type> V_zz(n_z,n_z);
  Type jnll = 0;

  // Define process variance
  V_zz = make_covariance( beta_j, RAM, RAMstart, n_z, int(0) );

  // Probability of random coefficients
  for( int i=0; i<n_i; i++){
    jnll += density::MVNORM(V_zz)( x_iz.row(i) );
  }

  // Probability of data conditional on fixed and random effect values
  for( int i=0; i<n_i; i++){
  for( int z=0; z<n_z; z++){
    if(familycode_z(z)==1){
      jnll -= dpois( y_iz(i,z), exp(x_iz(i,z)), true );
    }
  }}

  // Reporting
  REPORT( V_zz );
  return jnll;
}
