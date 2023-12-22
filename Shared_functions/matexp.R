matexp <-
function( Mrate,
          log2steps = 0, # Number of Euler-approximation steps
          zap_small = FALSE ){

  require(Matrix)
  require(expm)
  if( (log2steps <=0 ) || (log2steps > 100) ){
    # Full version ... note that expm::expm is faster that Matrix::expm
    out = expm(Mrate)
    return( Matrix(out) )
  }else{
    # Euler approximation
    Mrate = Diagonal(nrow(Mrate)) + Mrate / (2^log2steps)
    for(stepI in seq(1,log2steps,length=log2steps)){
      Mrate = Mrate %*% Mrate
    }
    if(zap_small) Mrate = zapsmall(Mrate)
    return( Mrate )
  }
}