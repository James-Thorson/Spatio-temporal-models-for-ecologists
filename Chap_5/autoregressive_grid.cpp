
#include <TMB.hpp>
#include "helper_functions.hpp"
#include "conditional_distribution.hpp"
#include "joint_distribution.hpp"

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Settings
  DATA_VECTOR( Options_vec );
  // Slot 0: method for calculating probability of random effects

  // Data
  DATA_VECTOR( c_i );
  DATA_IMATRIX( xy_i );

  // Parameters
  PARAMETER( beta0 );
  PARAMETER( ln_sigma2 );
  PARAMETER( logit_rho );

  // Random effects
  PARAMETER_ARRAY( epsilon_xy );

  // Objective funcction
  Type jnll = 0;
  Type sigma2 = exp(ln_sigma2);
  Type rho = 1 / (1 + exp(-logit_rho));

  //////  START IN-LINE CODE
  //// Probability of random effects
  // Conditional in dimension-Y, joint in dimension-X
  if( Options_vec(0)==1 ){
    matrix<Type> Q_yy = conditional_distribution( epsilon_xy, rho, sigma2, jnll );
    REPORT( Q_yy );
  }

  // Kroenekcer product of precision in both dimensions
  if( Options_vec(0)==2 ){
    matrix<Type> Q_zz = joint_distribution( epsilon_xy, rho, sigma2, jnll );
    REPORT( Q_zz );
  }

  // Calculate using built-in TMB functions
  using namespace density;
  if( Options_vec(0)==3 ){
    // Include "pow(1-pow(rho,2),0.5)" twice for 2D unit variance
    Type scale = pow(sigma2,0.5) / pow(1-pow(rho,2),0.5) / pow(1-pow(rho,2),0.5);
    jnll += SCALE( SEPARABLE(AR1(rho), AR1(rho)), scale )( epsilon_xy );
  }
  //////  END IN-LINE CODE

  // Predict densities
  matrix<Type> D_xy(epsilon_xy.rows(), epsilon_xy.cols());
  for( int x=0; x<epsilon_xy.rows(); x++ ){
  for( int y=0; y<epsilon_xy.cols(); y++ ){
    D_xy(x,y) = exp(beta0 + epsilon_xy(x,y));
  }}

  // Probability of data conditional on random effects
  for( int i=0; i<c_i.size(); i++){
    jnll -= dpois( c_i(i), D_xy(xy_i(i,0),xy_i(i,1)), true );
  }

  // Reporting
  REPORT( D_xy );
  ADREPORT( D_xy );
  return jnll;
}
