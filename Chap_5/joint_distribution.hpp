template<class Type>
matrix<Type> joint_distribution( array<Type> epsilon_xy,
                                  Type rho,
                                  Type sigma2,
                                  Type &jnll_pointer ){

  int n_x = epsilon_xy.rows();
  int n_y = epsilon_xy.cols();

  // Make precision matrix for Q_xx
  matrix<Type> Q_xx(n_x,n_x);
  Q_xx.setZero();
  for(int x=0; x<n_x; x++) Q_xx(x,x) = (1+pow(rho,2));
  for(int x=1; x<n_x; x++){
    Q_xx(x-1,x) = -rho;
    Q_xx(x,x-1) = -rho;
  }
  // Make precision matrix for Q_yy
  matrix<Type> Q_yy(n_y,n_y);
  Q_yy.setZero();
  for(int y=0; y<n_y; y++) Q_yy(y,y) = (1+pow(rho,2));
  for(int y=1; y<n_y; y++){
    Q_yy(y-1,y) = -rho;
    Q_yy(y,y-1) = -rho;
  }

  // Calculate probability
  int n_z = n_x * n_y;
  matrix<Type> Q_zz(n_z, n_z);
  Q_zz = kronecker( Q_yy, Q_xx );
  Q_zz = Q_zz / sigma2;
  vector<Type> epsilon_z(n_z);
  int Count = 0;
  for( int y=0; y<n_y; y++){
  for( int x=0; x<n_x; x++){
    epsilon_z(Count) = epsilon_xy(x,y);
    Count++;
  }}
  jnll_pointer -= dmvnorm( epsilon_z, Q_zz, true );
  return Q_zz;
}
