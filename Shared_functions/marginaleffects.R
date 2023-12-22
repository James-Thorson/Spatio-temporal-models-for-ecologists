# Function to get coefficients for TMB model
get_coef.custom_tmb = function(model, param, ...){
  out = model$parhat[[param]]
  names(out) = rep(param, length(out))
  return(out)
}

# Function to get variance-covariance for TMB model
get_vcov.custom_tmb = function(model, param, ...){
  rows = which( names(model$opt$par) == param )
  array( model$opt$SD$cov.fixed[rows,rows],
    dim = rep(length(rows),2),
    dimnames = list(rep(param,length(rows)),rep(param,length(rows))) )
}

# Function to change coefficients for TMB model
set_coef.custom_tmb = function(model, newpar, param, ...){
  model$parhat[[param]] <- newpar
  return(model)
}

# Function to get predictions when changing coefficients
get_predict.custom_tmb = function(model, newdata, param, center=FALSE, ...){
  # build original model.frame
  frame0 = model.frame( formula=model$formula, data=model$data )
  terms0 = terms( frame0 )
  xlevels = .getXlevels( terms0, frame0 )
  # get new design matrix
  terms1 = delete.response( terms0 )
  frame1 = model.frame( terms1, newdata, xlev=xlevels )
  X_ik = model.matrix( terms1, frame1 )
  gamma_k = get_coef.custom_tmb(model, param)
  # Calculate linear predictor and format output
  yhat_i = X_ik %*% gamma_k
  if(center==TRUE) yhat_i = yhat_i - mean(yhat_i)
  out = data.frame( rowid=seq_along(yhat_i[,1]), estimate=yhat_i )
  return(out)
}

# Change marginaleffects options to define `custom_tmb` class
options("marginaleffects_model_classes" = "custom_tmb")

