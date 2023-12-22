rotate_pca <-
function( L_tf, 
          x_sf,
          order = NULL ){
  
  # Eigen-decomposition
  Cov_tmp = L_tf %*% t(L_tf)
  Cov_tmp = 0.5*Cov_tmp + 0.5*t(Cov_tmp) # Ensure symmetric
  Eigen = eigen(Cov_tmp)
  
  # Desired loadings matrix
  L_tf_rot = (Eigen$vectors%*%diag(sqrt(Eigen$values)))[,1:ncol(L_tf),drop=FALSE]

  # My new factors
  require(corpcor)
  H = pseudoinverse(L_tf_rot) %*% L_tf
  x_sf = t(H %*% t(x_sf))
  
  # Get all loadings matrices to be increasing or decreasing
  if( !is.null(order) ){
    for( f in 1:ncol(L_tf) ){
      Lm = lm( L_tf_rot[,f] ~ 1 + I(1:nrow(L_tf)) )
      Sign = sign(Lm$coef[2]) * ifelse(order=="decreasing", -1, 1)
      L_tf_rot[,f] = L_tf_rot[,f] * Sign
      x_sf[,f] = x_sf[,f] * Sign
    }
  }
  
  # return
  out = list( "L_tf"=L_tf_rot, "x_sf"=x_sf, "H"=H)
  return(out)
}
