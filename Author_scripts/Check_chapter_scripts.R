
# Options that do not work
#shell( cmd = paste0("R CMD BATCH --no-environ --no-save ",file.path(root_dir, "Chap_1", "Gompertz_survival.R")) )
#testthat::test_file( file.path(root_dir, "Chap_1", "Gompertz_survival.R") )
#shell.exec( file.path(root_dir, "Chap_1", "Gompertz_survival.R") )

root_dir = R'(C:\Users\James.Thorson\Desktop\Git\Spatio-temporal-models-for-ecologists)'

#
unload_DLLs = function(){
  DLLs = getLoadedDLLs()
  paths = sapply( DLLs, \(x) x[["path"]] )
  isTMB = grep("Spatio-temporal-models-for-ecologists", paths )
  for(i in isTMB ){
    dyn.unload( paths[i] )
  }
}
run_script = function( script ){
  # Source script
  source( script )  # In case of intentional errors
  # Unload linked TMB models
  unload_DLLs()
  # Detach packages
  loaded_packages = names(sessionInfo()$otherPkgs)
  if(length(loaded_packages)>0) lapply(paste0('package:', loaded_packages), detach, character.only=TRUE, unload=TRUE)
}

# Chap-1
scripts = c( "Gompertz_survival.R", "Poisson_point_process.R" )
sapply( scripts, FUN=\(x) run_script(file.path(root_dir, "Chap_1", x)) )

# Chap-2 ... ends with intentional error
scripts = c( "Demo_GLMM.R", "Laplace_examples.R" )
sapply( scripts, FUN=\(x) run_script(file.path(root_dir, "Chap_2", x)) )

# Chap-3
scripts = c( "Gompertz_dynamics.R" )
sapply( scripts, FUN=\(x) run_script(file.path(root_dir, "Chap_3", x)) )

# Chap-4
scripts = c( "Simulate_movement.R" )
sapply( scripts, FUN=\(x) run_script(file.path(root_dir, "Chap_4", x)) )

# Chap-5
scripts = c( "Basis_expansion_demo.R", "Spatial_GLMM.R", "SPDE.R" )
sapply( scripts, FUN=\(x) run_script(file.path(root_dir, "Chap_5", x)) )

# Chap-6
scripts = c( "Ozone_example.R", "Epsilon_demo.R", "Nmixture_example.R" )
sapply( scripts, FUN=\(x) run_script(file.path(root_dir, "Chap_6", x)) )

# Chap-7
scripts = c( "Demo_path_diagram.R", "ICAR_covariate_simulation.R", "integrated_model.R" )
sapply( scripts, FUN=\(x) run_script(file.path(root_dir, "Chap_7", x)) )

# Chap-8
scripts = c( "Infill_vs_sprawl_asymptotics.R", "pollock_index.R" )
sapply( scripts, FUN=\(x) run_script(file.path(root_dir, "Chap_8", x)) )

# Chap-9
scripts = c( "sea_ice.R", "pollock_SVC.R" )
sapply( scripts, FUN=\(x) run_script(file.path(root_dir, "Chap_9", x)) )

# Chap-10
scripts = c( "CTMC_demo.R", "CTMC_pacific_cod.R", "CTMC_eagles.R", "PDE_decisions.R" )
sapply( scripts, FUN=\(x) run_script(file.path(root_dir, "Chap_10", x)) )

# Chap-11
scripts = c( "BBS.R" )
sapply( scripts, FUN=\(x) run_script(file.path(root_dir, "Chap_11", x)) )

