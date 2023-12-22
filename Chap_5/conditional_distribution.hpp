template<class Type>
matrix<Type> conditional_distribution( array<Type> epsilon_xy,
                                 Type rho,
                                 Type sigma2,
                                 Type &jnll_pointer ){

  int n_x = epsilon_xy.rows();
  int n_y = epsilon_xy.cols();

  // Make precision matrix for Q_yy
  matrix<Type> Q_yy(n_y, n_y);
  Q_yy.setZero();
  for(int y=0; y<n_y; y++) Q_yy(y,y) = (1+pow(rho,2));
  for(int y=1; y<n_y; y++){
    Q_yy(y-1,y) = -rho;
    Q_yy(y,y-1) = -rho;
  }

  // Calculate probability
  vector<Type> Tmp_y(n_y);
  matrix<Type> Q0_yy(n_y,n_y);
  Q0_yy = Q_yy * ( 1-pow(rho,2) ) / sigma2;
  matrix<Type> Q1_yy(n_y,n_y);
  Q1_yy = Q_yy / sigma2;
  for(int x=0; x<n_x; x++){
    for(int y=0; y<n_y; y++){
      if(x==0) Tmp_y(y) = epsilon_xy(0,y);
      if(x>=1) Tmp_y(y) = epsilon_xy(x,y) - rho*epsilon_xy(x-1,y);
    }
    if(x==0) jnll_pointer -= dmvnorm( Tmp_y, Q0_yy, true );
    if(x>=1) jnll_pointer -= dmvnorm( Tmp_y, Q1_yy, true );
  }
  return(Q1_yy);
}
