# get_coef:  gets response covariates
get_coef.custom_tmb = function(model, ...){
  #out = as.vector(model$parhat$beta_kc)
  #names(out) = rep("beta_kc", length(out))
  #return(out)
  as.vector(model$parhat$beta_kc)
}
# get_vcov:  gets covariance for response covariates
get_vcov.custom_tmb = function(model, ...){
  #out = model$opt$SD$cov
  #dimnames(out) = list( rep("beta_kc",nrow(out)), rep("beta_kc",nrow(out)) )
  #return(out)
  model$opt$SD$cov
}
# set_coef:  sets new value for response covariates
set_coef.custom_tmb = function(model, newpar, ...){
  model$parhat$beta_kc[] <- newpar
  return(model)
}
# get_predict:  calculates response given new covariates
get_predict.custom_tmb = function(model, newdata, ...){
  # build original model.frame
  frame0 = model.frame( formula=model$formula, data=model$data )
  terms0 = terms( frame0 )
  xlevels = .getXlevels( terms0, frame0 )
  # get new design matrix
  terms1 = delete.response( terms0 )
  frame1 = model.frame( terms1, newdata, xlev=xlevels )
  X_ik = model.matrix( terms1, frame1 )
  # make predictions
  beta_kc = model$parhat$beta_kc
  beta_kc[] = get_coef.custom_tmb(model)
  yhat_ic = X_ik %*% beta_kc
  # return in long-form data frame
  out = expand.grid( rowid=seq_along(yhat_ic[,1]), species=model$species )
  out$estimate = as.vector(yhat_ic)
  return(out)
}
options("marginaleffects_model_classes" = "custom_tmb")
