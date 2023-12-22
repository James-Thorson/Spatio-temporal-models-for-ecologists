template<class Type>
matrix<Type> make_covariance( vector<Type> s2_z,
                              matrix<int> RAM,
                              vector<Type> RAMstart_z,
                              int n_rows,
                              int n_cols ){

  // Define temporary objects
  matrix<Type> Cov_rr(n_rows, n_rows);
  matrix<Type> I_rr( n_rows, n_rows );
  I_rr.setIdentity();
  if( RAM.rows()>0 ){
    // Assemble Structural Equation Model covariance
    matrix<Type> L_rc(n_rows, n_rows);
    matrix<Type> Rho(n_rows, n_rows);
    matrix<Type> Gamma(n_rows, n_rows);
    Rho.setZero();
    Gamma.setZero();
    for(int zI=0; zI<RAM.rows(); zI++){
      if( RAM(zI,3)>0 ){
        if(RAM(zI,0)==1) Rho(RAM(zI,1)-1,RAM(zI,2)-1) = s2_z(RAM(zI,3)-1);
        if(RAM(zI,0)==2) Gamma(RAM(zI,1)-1, RAM(zI,2)-1) = s2_z(RAM(zI,3)-1);
      }else{
        if(RAM(zI,0)==1) Rho( RAM(zI,1)-1, RAM(zI,2)-1 ) = RAMstart_z(zI);
        if(RAM(zI,0)==2) Gamma( RAM(zI,1)-1, RAM(zI,2)-1 ) = RAMstart_z(zI);
      }
    }
    L_rc = I_rr - Rho;
    L_rc = atomic::matinv( L_rc );
    L_rc = L_rc * Gamma;
    Cov_rr = L_rc * L_rc.transpose();
  }else{
    // Assemble loadings matrix
    matrix<Type> L_rc(n_rows, n_cols);
    L_rc.setZero();
    int Count = 0;
    if( n_cols>0 ){
      for(int r=0; r<n_rows; r++){
      for(int c=0; c<n_cols; c++){
        if(r>=c){
          L_rc(r,c) = s2_z(Count);
          Count++;
        }else{
          L_rc(r,c) = 0.0;
        }
      }}
    }
    // Diagonal and equal coariance matrix
    Cov_rr = I_rr * exp( s2_z(Count) );
    // Add factor-model covariance matrix
    Cov_rr += L_rc * L_rc.transpose();
  }
  return Cov_rr;
}
