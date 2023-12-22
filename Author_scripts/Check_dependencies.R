library(renv)

full = NULL
for( i in c(1:11) ){
  dir = R'(C:\Users\James.Thorson\Desktop\Git\Spatio-temporal-models-for-ecologists\Chap_)'
  dep = dependencies( paste0(dir,i) )
  full = rbind( full, cbind(dep,"Chapter"=i) )
}

subset( full, Package=="TMBhelper" )
subset( full, Package=="INLA" )
subset( full, Package=="fmesher" )
