# Parameters
n_x = 20
n_y = 20
n_t = 10
corr_s = 0.9
corr_t = 0.9
margSD = 0.5
meanD = 3

#
loc_gz = expand.grid("x"=1:n_x, "y"=1:n_y)
log_d_st = matrix(NA, nrow=n_x*n_y, ncol=n_t )

#
R_ss = corr_s^as.matrix(dist(loc_gz))
condSD = margSD * sqrt(1-corr_t^2)
alpha = meanD*(1-corr_t)

# Loop through years
for( tI in seq_len(n_t) ){
  if(tI==1){
    log_d_st[,1] = mvtnorm::rmvnorm(n=1, mean=rep(meanD,n_x*n_y), sigma=margSD^2*R_ss)
  }else{
    log_d_st[,tI] = alpha + corr_t*mvtnorm::rmvnorm(n=1, mean=log_d_st[,tI-1], sigma=condSD^2*R_ss)
  }
}

# Plot
out = sf::st_as_sf(cbind("d"=log_d_st,loc_gz), coords = c("x","y"))
out = stars::st_rasterize(out)
par(mfrow=c(3,3))
for(tI in 1:9){
  image(out[tI], axes=TRUE, col=viridisLite::viridis(10))
}

