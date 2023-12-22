do_newtonstep <-
function( par,
          fn,
          gr,
          steps = 1 ){

  ## v^* = v - H^{-1} g
  for(i in seq_len(steps)) {
    g = as.numeric( gr(par) )
    h = optimHess( newpar, fn=fn, gr=gr)
    par = par - solve(h, g)
  }
  return(par)
}
