
library(sf)
library(TMB)
library(ape)
library(rnaturalearth)
library(elevatr)
library(phylobase)
library(phylosignal)
library(viridisLite)
library(raster)
library(fmesher)

setwd( R'(C:\Users\James.Thorson\Desktop\Git\Spatio-temporal-models-for-ecologists\Chap_11)' )
source( "../Shared_functions/add_legend.R" )

# Load data
DF = read.csv( file="Top20_Samples.csv" )
trait_set = read.csv( "Top20_traits.csv" )

# Load population density
pop_dens = st_read( "../Chap_6/population_density.csv", options=c("X_POSSIBLE_NAMES=X","Y_POSSIBLE_NAMES=Y"), crs=st_crs("+proj=longlat +datum=WGS84") )
  pop_dens$Dens2020 = as.numeric(pop_dens$Dens2020)

# Load copNDVI (saved from rasterdiv)
copNDVI = raster( "../Chap_10/NDVI.tif" )
#library(rasterdiv)

# Get spatial domain
sf_states = ne_states( c("United States of America"), return="sf")
sf_states = sf_states[pmatch(c("Cal", "Oregon", "Washington", "Idaho", "Montana", "Utah", "New Mex", "Arizona", "Wyoming", "Colorad", "Nevada"), sf_states$name_en),]
sf_states = st_union(sf_states)

# Create data-frame
sf_DF = st_as_sf( DF, coords=c("Longitude","Latitude"), crs="+proj=longlat +datum=WGS84")

#
sf_fullgrid = st_make_grid( sf_DF, cellsize=1, square=FALSE )
sf_grid = st_make_valid(st_intersection( sf_fullgrid, sf_states ))
sf_grid = sf_grid[ st_area(sf_grid)>(0.01*max(st_area(sf_grid))) ]    # or 0.01

# make data frame of covariates
df_grid = st_centroid(sf_grid)
df_grid = get_elev_point( df_grid, src = "aws" )
df_grid$log_elevation_km = log( ifelse(df_grid$elevation<1, 1, df_grid$elevation) / 1000 )
df_grid$NDVI = extract( x=copNDVI, y=as(df_grid,"Spatial") )
df_grid$scale_NDVI = scale( df_grid$NDVI )[,1]
df_grid$pop_dens = pop_dens$Dens2020[ st_nearest_feature( sf_grid, pop_dens ) ]
df_grid$log_pop_dens = log(df_grid$pop_dens)
df_grid = data.frame(df_grid)

#
sf_DF = st_intersection( sf_DF, st_union(sf_grid) )
sf_DF$Genus_species = factor(sf_DF$Genus_species)
temp_DF = get_elev_point( sf_DF, src = "aws" )
sf_DF$log_elevation_km = log( ifelse(temp_DF$elevation<1, 1, temp_DF$elevation) / 1000 )
sf_DF$NDVI = extract( x=copNDVI, y=as(sf_DF,"Spatial") )
sf_DF$scale_NDVI = scale( sf_DF$NDVI )[,1]
sf_DF$pop_dens = pop_dens$Dens2020[ st_nearest_feature( sf_DF, pop_dens ) ]
sf_DF$log_pop_dens = log(sf_DF$pop_dens)

#
taxa = levels( sf_DF$Genus_species )

# create mesh
loc_DF = st_coordinates(sf_DF)
loc_grid = st_coordinates(st_centroid(sf_grid))
mesh = fm_mesh_2d( loc_grid, refine=TRUE, cutoff=1)
# Create matrices in INLA
spde <- fm_fem(mesh, order=2)
# create projection matrix from vertices to samples
A_is = fm_evaluator( mesh, loc=loc_DF )$proj$A
# create projection matrix from vertices to grid
A_gs = fm_evaluator( mesh, loc=loc_grid )$proj$A

############## BEGIN SNIPPET FOR TEXTBOOK
# Read tree
tree = read.tree( "Top20_Tree.tre" )
  tree$node.label = c( "root", seq_len(Nnode(tree)-1) )

# Extract path from species to root
root = Ntip(tree) + 1 # index for root
paths = sapply( match(taxa,tree$tip.label), FUN=nodepath, phy=tree, to=root )
names(paths) = taxa

# Convert to triplet form i=species; j=edge; x=edge_length
i = unlist(lapply(seq_along(paths), FUN=function(i){rep(i,length(paths[[i]]))}))
j = unlist(paths)
edge = match( j, tree$edge[,2] )
x = tree$edge.length[edge]

# remove root and renumber
ijx = na.omit(cbind(i,j,x))
ijx[,1:2] = ifelse( ijx[,1:2]>Ntip(tree), ijx[,1:2]-1, ijx[,1:2] )

# Convert to dense design matrix
PhyloDesign_gc = matrix(0, nrow=Nedge(tree), ncol=length(taxa) )
PhyloDesign_gc[ijx[,2:1]] = ijx[,3]
############### END SNIPPET

##############
# Run alternatives
##############

Version = c( "Residual"=TRUE, "Phylo"=TRUE, "Traits"=TRUE, "Habitat"=TRUE )

# Factors
if( Version["Residual"] ){
  Lform_jc = diag(1:nlevels(sf_DF$Genus_species))
}else{
  Lform_jc = matrix( 0, nrow=0, ncol=length(taxa) )
}

# Phylogenetic Design matrix
if( Version["Phylo"] ){
  C_gc = PhyloDesign_gc / mean(colSums(PhyloDesign_gc))
}else{
  C_gc = matrix(0, nrow=0, ncol=length(taxa) )
}

# Traits
if( Version["Traits"]==TRUE ){
  T_ch = model.matrix( ~ 0 + log(Mass) + factor(Primary.Lifestyle) + log(Hand.Wing.Index), data=trait_set[match(taxa,trait_set$Genus_species),] )
  T_hc = t(scale(T_ch))
}else{
  T_hc = matrix(0, nrow=0, ncol=length(taxa) )
}

#
if( Version["Habitat"]==TRUE ){
  #habitat_formula = ~ 0 + log(elevation) + NDVI + log(pop_dens)
  habitat_formula = ~ 0 + poly(log_elevation_km,2,raw=TRUE) + poly(scale_NDVI,2,raw=TRUE) + poly(log_pop_dens,2,raw=TRUE)
  X_gk = model.matrix( habitat_formula, data=df_grid )
  X_ik = model.matrix( habitat_formula, data=sf_DF )
}else{
  X_gk = model.matrix( ~ 0, data=df_grid )
}

#
compile( "phyloSDM.cpp" ) # framework='TMBad'
dyn.load( dynlib("phyloSDM") )

# Compile
rmatrix = function(nrow,ncol,...) matrix(data=rnorm(n=nrow*ncol,...), nrow=nrow, ncol=ncol )
Params = list( "ln_kappa" = log(1),
               "ln_sigmaM" = log(1),
               "ln_sigmaC" = rep(0, ifelse(nrow(C_gc)>0,1,0) ),
               "ln_sigmaT_h" = rep(0,nrow(T_hc)),
               "alpha_c" = rep(0,length(taxa)),
               "lambda_jc" = ifelse(Lform_jc==0,0,1) * rmatrix(sd=0.1, nrow=nrow(Lform_jc), ncol=ncol(Lform_jc)),
               "beta_kc" = matrix(0, nrow=ncol(X_ik), ncol=length(taxa)),
               "gamma_sh" = matrix(0, nrow=nrow(spde$c0), ncol=nrow(T_hc)),
               "delta_sg" = matrix(0, nrow=nrow(spde$c0), ncol=nrow(C_gc)),
               "omega_sj" = matrix(0, nrow=nrow(spde$c0), ncol=nrow(Lform_jc)) )
Data = list( "n_i" = sf_DF$SpeciesTotal,
             "c_i" = as.numeric(sf_DF$Genus_species)-1,
             "C_gc" = C_gc,
             "T_hc" = T_hc,
             "X_gk" = X_gk,
             "X_ik" = X_ik,
             "M0" = spde$c0,
             "M1" = spde$g1,
             "M2" = spde$g2,
             "A_is" = A_is,
             "A_gs" = A_gs )
Map = list( "lambda_jc" = factor(ifelse(Lform_jc==0,NA,Lform_jc)) )
Random = c("gamma_sh", "delta_sg", "omega_sj")

# Use REML
Random = c( Random, "alpha_c", "beta_kc" )

# Build and optimize
obj = MakeADFun( data=Data, parameters=Params, map=Map, random=Random )
opt = nlminb( start=obj$par, obj=obj$fn, gr=obj$gr )
opt$SD = sdreport( obj, getReportCovariance=TRUE )
report = obj$report()
parhat = obj$env$parList()

# plotting stuff
short_names = sapply( taxa, FUN=function(char){
  tmp = strsplit(char,"_")[[1]]
  paste0( substr(tmp[1],1,1), ". ", tmp[2] )
})
common_names = trait_set[match(taxa,trait_set$Genus_species),'Common_name']
trait_names = c( "log(mass)", "Aerial", "Generalist", "Insessorial", "Terrestrial", "log(Hand Wing Index)")

# Plot tree
tree_tmp = tree
  tree_tmp$tip.label = paste0( tree_tmp$tip.label, ", ", common_names )
png( file=paste0("Phylogeny.png"), width=6, height=6, units="in", res=200 )
  par( mar=c(2,0,0,0), mgp=c(2,0.5,0), oma=c(2,0,0,0) )
  plot( tree_tmp, show.node.label=TRUE, cex=0.8 )
  axisPhylo( 1 )
  mtext( side=1, text="Time before present (million years)", outer=TRUE, adj=0)
dev.off()

# Plot covariates
png( file=paste0("Covariates.png"), width=6, height=2.5, units="in", res=200 )
  par( mfrow=c(1,3) )
  plotgrid = st_sf(sf_grid, "log(elevation)"=df_grid$log_elevation_km, "scale(NDVI)"=df_grid$scale_NDVI, "log(pop density)"=df_grid$log_pop_dens, crs=st_crs(sf_grid) )
  for( z in 1:3 ){
    plot(plotgrid[,z],reset=FALSE, key.pos=NULL, border=NA )
    add_legend( round(range(plotgrid[[z]]),2), legend_y=c(0.6,1), legend_x=c(1,1.05) )
  }
dev.off()

#
colnames(T_ch) = letters[1:ncol(T_ch)]
  tree_tmp$tip.label = rep("",length(taxa)) # common_names[match(tree$tip.label,taxa)]
phylo = phylo4d( tree_tmp,
                   tip.data = T_ch,
                   check.node.labels = "drop",
                   missing.data = "warn" )
png( file=paste0("Phylogeny_and_traits.png"), width=6, height=6, units="in", res=200 )
  par( mar=c(0,0,0,0), oma=c(0,0,0,10) )
  barplot( phylo, tip.labels=rep("",length(taxa)) )
  axis(4, at=seq_along(taxa), labels=common_names[match(tree$tip.label,taxa)], line=0, las=2, lwd=0, lwd.ticks=0 )
dev.off()

# Calculate variances
Vcar = mean(diag( solve(report$Q) ))
  trace = function(mat) sum(diag(mat))
  Vbeta = t(parhat$beta_kc) %*% cov(X_gk) %*% parhat$beta_kc
  Vgamma = Vcar * t(T_hc) %*% diag(exp(2*parhat$ln_sigmaT_h)) %*% T_hc
  Vdelta = Vcar * exp(2*parhat$ln_sigmaC) * t(C_gc) %*% C_gc
  Vomega = Vcar * t(parhat$lambda_jc) %*% parhat$lambda_jc
  Vtotal = Vbeta + Vgamma + Vdelta + Vomega
Height_zc = rbind( diag(Vbeta), diag(Vgamma), diag(Vdelta), diag(Vomega) )
  colnames(Height_zc) = common_names
png( file=paste0("Variance_partitioning.png"), width=6, height=7.5, units="in", res=200 )
  par( mar=c(2,15,1,1) )
  barplot( height=Height_zc, horiz=TRUE, las=1, col=viridisLite::viridis(4) )
  legend( "topright", bty="n", fill=viridisLite::viridis(4), legend=c("habitat","traits","phylogeny","residual"), ncol=2)
dev.off()

# Plot densities
colnames(report$p_gc) = common_names
plotgrid = st_sf(sf_grid, report$p_gc, crs=st_crs(sf_grid) )
png( file=paste0("Densities.png"), width=6, height=7.5, units="in", res=200 )
  par( mfrow=c(5,4) )
  for( p in 1:(ncol(plotgrid)-1) ){
    plot(plotgrid[,p], max.plot=20, border=NA, key.pos=NULL, reset=FALSE, pal=viridis )
    add_legend( round(range(plotgrid[[p]]),2), legend_y=c(0.6,1), legend_x=c(1,1.05), col=viridis(10) )
  }
dev.off()

# Plot traits
if( Version["Traits"] ){
  colnames(report$gamma_gc) = common_names
  plotgrid = st_sf(sf_grid, report$gamma_gc, crs=st_crs(sf_grid) )
  png( file=paste0("Traits.png"), width=6, height=7.5, units="in", res=200 )
    par( mfrow=c(5,4) )
    for( p in 1:(ncol(plotgrid)-1) ){
      plot(plotgrid[,p], max.plot=20, border=NA, key.pos=NULL, reset=FALSE, pal=viridis )
      add_legend( round(range(plotgrid[[p]]),2), legend_y=c(0.6,1), legend_x=c(1,1.05), col=viridis(10) )
    }
  dev.off()

  gamma_gh = as.matrix(A_gs %*% parhat$gamma_sh)
  colnames(gamma_gh) = trait_names
  plotgrid = st_sf(sf_grid, gamma_gh, crs=st_crs(sf_grid) )
  png( file=paste0("Traits-responses.png"), width=6, height=4, units="in", res=200 )
    par( mfrow=c(2,3) )
    for( p in 1:(ncol(plotgrid)-1) ){
      Range = range(parhat$gamma_sh)
      plot( plotgrid[,p], max.plot=20, border=NA, key.pos=NULL, reset=FALSE, pal=viridis, breaks=seq(Range[1],Range[2],len=10) )
      if(p==1) add_legend( round(Range,2), legend_y=c(0.6,1), legend_x=c(1,1.05), col=viridis(10) )
    }
  dev.off()
}

# Plot phylo
if( Version["Phylo"] ){
  colnames(report$delta_gc) = common_names
  plotgrid = st_sf(sf_grid, report$delta_gc, crs=st_crs(sf_grid) )
  png( file=paste0("Phylo.png"), width=6, height=7.5, units="in", res=200 )
    par( mfrow=c(5,4) )
    for( p in 1:(ncol(plotgrid)-1) ){
      plot(plotgrid[,p], max.plot=20, border=NA, key.pos=NULL, reset=FALSE, pal=viridis )
      add_legend( round(range(plotgrid[[p]]),2), legend_y=c(0.6,1), legend_x=c(1,1.05), col=viridis(10) )
    }
  dev.off()

  species = c("Hirundo_rustica", "Petrochelidon_pyrrhonota", "Eremophila_alpestris")   # barn swallow, cliff swallow, horned lark
    paths_to_plot = setdiff( unique(unlist(paths[species])), Ntip(tree)+1 )
    paths_to_plot = sort(ifelse( paths_to_plot>Ntip(tree), paths_to_plot-1, paths_to_plot ))
  delta_gg = as.matrix(A_gs %*% parhat$delta_sg)
    colnames(delta_gg) = c( trait_set[match(tree$tip.label,trait_set$Genus_species),'Common_name'], tree$node.label[-1] )
  plotgrid = st_sf(sf_grid, delta_gg[,paths_to_plot], crs=st_crs(sf_grid) )
  png( file=paste0("Ancestors.png"), width=6, height=3, units="in", res=200 )
    par( mfrow=c(2,4) )
    for( p in 1:(ncol(plotgrid)-1) ){
      Range = range(delta_gg[,paths_to_plot])
      plot(plotgrid[,p], max.plot=20, border=NA, key.pos=NULL, reset=FALSE, cex.main=0.8, pal=viridis, breaks=seq(Range[1],Range[2],len=10) )
      if(p==1) add_legend( round(Range,2), legend_y=c(0.6,1), legend_x=c(1,1.05), col=viridis(10) )
    }
  dev.off()
}

# Plot residuals
if( Version["Residual"] ){
  colnames(report$omega_gc) = common_names
  plotgrid = st_sf(sf_grid, report$omega_gc, crs=st_crs(sf_grid) )
  png( file=paste0("Residuals.png"), width=6, height=7.5, units="in", res=200 )
    par( mfrow=c(5,4) )
    for( p in 1:(ncol(plotgrid)-1) ){
      Range = range(c(plotgrid[[p]],-0.1,0.1))
      plot(plotgrid[,p], max.plot=20, border=NA, key.pos=NULL, reset=FALSE, pal=viridis, breaks=seq(Range[1],Range[2],len=10) )
      add_legend( round(Range,2), legend_y=c(0.6,1), legend_x=c(1,1.05), col=viridis(10) )
    }
  dev.off()
}

# Plot Habitat
if( Version["Habitat"] ){
  colnames(report$beta_gc) = common_names
  plotgrid = st_sf(sf_grid, report$beta_gc, crs=st_crs(sf_grid) )
  png( file=paste0("Habitat.png"), width=6, height=7.5, units="in", res=200 )
    par( mfrow=c(5,4) )
    for( p in 1:(ncol(plotgrid)-1) ){
      plot(plotgrid[,p], max.plot=20, border=NA, key.pos=NULL, reset=FALSE, pal=viridis )
      add_legend( round(range(plotgrid[[p]]),2), legend_y=c(0.6,1), legend_x=c(1,1.05), col=viridis(10) )
    }
  dev.off()

  # Use marginaleffects
  library(marginaleffects)
  library(ggplot2)
  source( "BBS_marginaleffects.R" )
  # Bundle materials
  fit = list( obj = obj,
              opt = opt,
              data = df_grid,
              formula = habitat_formula,
              parhat = obj$env$parList(),
              species = common_names )
  class(fit) = "custom_tmb"
  quant = function(x) seq(min(x),max(x),length=21)
  # get_coef( fit )
  # get_vcov( fit )
  # set_coef( fit, newpar=rep(0,length(get_coef(fit))) )
  # get_predict( fit, newpar=rep(0,length(get_coef(fit))), newdata=newdata )

  # Get prediction for elevation
  newdata = datagrid( newdata=df_grid[c('log_elevation_km','scale_NDVI','log_pop_dens')], log_elevation_km=quant, scale_NDVI=mean, log_pop_dens=mean )
    pred = predictions( fit, newdata=newdata )
  ggplot( pred, aes(log_elevation_km, estimate)) +
    geom_line( aes(y=estimate), color="blue", size=1 ) +
    geom_ribbon( aes( x=log_elevation_km, ymin=conf.low, ymax=conf.high), fill=rgb(0,0,1,0.2) ) +
    facet_wrap(vars(species), scales="free", ncol=4) +
    labs(y="Predicted response")
  ggsave( paste0("Covariates-elevation.png"), width=6, height=7.5 )

  # Get prediction for NDVI
  newdata = datagrid( newdata=df_grid, log_elevation_km=mean, scale_NDVI=quant, log_pop_dens=mean )
    pred = predictions( fit, newdata=newdata )
  ggplot( pred, aes(scale_NDVI, estimate)) +
    geom_line( aes(y=estimate), color="blue", size=1 ) +
    geom_ribbon( aes( x=scale_NDVI, ymin=conf.low, ymax=conf.high), fill=rgb(0,0,1,0.2) ) +
    facet_wrap(vars(species), scales="free", ncol=4) +
    labs(y="Predicted response")
  ggsave( paste0("Covariates-NDVI.png"), width=6, height=7.5 )

  # Get prediction for population density
  newdata = datagrid( newdata=df_grid, log_elevation_km=mean, scale_NDVI=mean, log_pop_dens=quant )
    pred = predictions( fit, newdata=newdata )
  ggplot( pred, aes(log_pop_dens, estimate)) +
    geom_line( aes(y=estimate), color="blue", size=1 ) +
    geom_ribbon( aes( x=log_pop_dens, ymin=conf.low, ymax=conf.high), fill=rgb(0,0,1,0.2) ) +
    facet_wrap(vars(species), scales="free", ncol=4) +
    labs(y="Predicted response")
  ggsave( paste0("Covariates-popdens.png"), width=6, height=7.5 )
}

#############
# Biogeographic clustering
#############

# Cluster
library(fastcluster)
Hclust = hclust.vector( report$p_gc, method="ward" )

# Cut and convert shape
png( file=paste0("Clusters-original.png"), width=6, height=6, units="in", res=200 )
  par( mfrow=c(2,2), mar=c(2,2,1,0) )
  plot( Hclust, labels=FALSE, xlab="" )
  for( k in 2:4 ){
    Class_g = cutree( Hclust, k = k )
    plotgrid = st_sf( sf_grid, Class=Class_g, crs=st_crs(sf_grid) )
    plot( plotgrid, reset=FALSE, key.pos=NULL, border=NA, main=paste0("k = ",k) )
  }
dev.off()

# Clusters and averages
png( file=paste0("Clusters-updated.png"), width=6, height=6, units="in", res=200 )
  par( mar=c(1,3,1,0), oma=c(9,0,0,0), mgp=c(2,0.5,0), tck=-0.02 )   # mfrow=c(3,2),
  layout( mat=matrix(1:6,ncol=2,byrow=TRUE), widths=c(0.5,1) )
  for( k in 2:4 ){
    Class_g = cutree( Hclust, k = k )
    # Maps
    plotgrid = st_sf( sf_grid, Class=Class_g, crs=st_crs(sf_grid) )
    plot( plotgrid, reset=FALSE, key.pos=NULL, border=NA, pal=viridis, main="" )
    mtext( side=2, text=paste0("k = ",k), font=2, line=-1 )
    # Class averages
    Class_gz = model.matrix( ~ 0 + factor(Class_g) )
    pmean_gz = (t(report$p_gc) %*% Class_gz) / nrow(report$p_gc)
    matplot( y=pmean_gz, xlab="", type="l", lwd=3, col=viridis(k), lty="solid", xaxt="n" )
    if(k==2) mtext(3, text="Cluster averages" )
    if(k==4) axis(1, at=1:nrow(pmean_gz), labels=common_names, las=3 )
  }
dev.off()

#############
# Plot beta diversity
#############

library(BAT)
library(FD)
library(stars)

# Convert to RasterBrick
density_gc = exp( report$p_gc + rep(1,nrow(report$p_gc)) %o% parhat$alpha_c )
sf_density = st_sf( sf_grid, density_gc, crs=st_crs(sf_grid) )
raster_grid = raster( sf_density, nrow=9, ncol=12 )
stars_grid = st_as_stars( raster_grid )
stars_density = st_rasterize( sf_density, template=stars_grid )
raster_density = as( stars_density, "Raster" )

# Plot beta diversity
raster_beta = raster.beta( raster_density[[1:20]], abund=TRUE )
plot(raster_beta)
dev.new()
plot(raster_density[[c(1,3)]])

# Convert to sf
stars_beta = st_as_stars(raster_beta, raster_grid)
sf_beta = st_as_sf( stars_beta )
names(sf_beta)[1:3] = c( "Beta_total", "Beta_replacement", "Beta_density" )
png( file=paste0("Beta_diversity.png"), width=6, height=2.5, units="in", res=200 )
  par( mfrow=c(1,3) )
  for( p in 1:(ncol(sf_beta)-1) ){
    plot(sf_beta[,p], pal=viridisLite::viridis, reset=FALSE, key.pos=NULL )
    add_legend( round(range(sf_beta[[p]]),2), legend_y=c(0.6,1), legend_x=c(1,1.05), col=viridis(10) )
  }
dev.off()

# Convert back to sf_grid
#stars_beta = st_as_stars(raster_beta, raster_grid)
#mat_beta = apply( stars_beta$Btotal, MARGIN=3, FUN=function(mat){as.vector(mat[,ncol(mat):1])} )
#colnames(mat_beta) = names(raster_beta)
#sf_beta = st_sf( sf_fullgrid, mat_beta, crs=st_crs(sf_fullgrid) )
#sf_beta = st_intersection( sf_beta, sf_states )
#sf_beta = sf_beta[ st_area(sf_beta)>(0.01*max(st_area(sf_beta))), ]
#plot( sf_beta, pal=viridisLite::viridis )

# Functional diversity
rownames(T_ch) = common_names
out = dbFD( x = as.matrix(T_ch),
      a = density_gc )
out$CWM = density_gc %*% T_ch / (rowSums(density_gc) %o% rep(1,ncol(T_ch)))
colnames(out$CWM) = trait_names

plotgrid = st_sf( sf_grid, as.matrix(out$CWM), crs=st_crs(sf_grid) )
png( file=paste0("community_weighted_traits.png"), width=6, height=4, units="in", res=200 )
  par( mfrow=c(2,3) )
  for( p in 1:ncol(out$CWM) ){
    plot(plotgrid[,p], pal=viridisLite::viridis, border=NA, reset=FALSE, key.pos=NULL)
    add_legend( round(range(plotgrid[[p]]),2), legend_y=c(0.6,1), legend_x=c(1,1.05), col=viridis(10) )
  }
dev.off()

Mat = cbind("functional evenness"=out$FEve, "functional divergence"=out$FDiv, "functional dispersion"=out$FDis)
plotgrid = st_sf( sf_grid, Mat, crs=st_crs(sf_grid) )
png( file=paste0("functional_diversity.png"), width=6, height=2.5, units="in", res=200 )
  par( mfrow=c(1,3) )
  for( p in 1:ncol(Mat) ){
    plot(plotgrid[,p], pal=viridisLite::viridis, border=NA, reset=FALSE, key.pos=NULL)
    add_legend( round(range(plotgrid[[p]]),2), legend_y=c(0.6,1), legend_x=c(1,1.05), col=viridis(10) )
  }
dev.off()

