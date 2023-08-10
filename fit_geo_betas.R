library(data.table)
library(dplyr)
library(metapop)
library(metapopnorge)
source("run_functions.R")

N_particles <- 20
N_threads <- 10
n_threads_ext <- 20

scenarios <- fread("parameter_files/scenarios_N.csv")

geo_scenarios <- scenarios[Regional==0 & R_const==1,]

baseline_scenario <- geo_scenarios[Name_geo=="geo_0"]


get_geo_params_file <- function(param_update_file, tag=""){
  update_params <- readRDS(param_update_file)
  print("target")
  target <- get_tot_hosp_scenario(baseline_scenario,1, update_params)
  print("finished target")
  params <-parallel::mclapply(2:nrow(geo_scenarios), function(i) {
    find_factor(geo_scenarios[i], target, update_params)}, mc.cores=n_threads_ext, mc.preschedule = FALSE
    )
  print(params)
  saveRDS(rbindlist(params), glue::glue("geo_params_{tag}.RDS"))
}



get_tot_hosp_scenario <- function(scenario, beta_factor, update_params){
  param_sets <- param_sets_from_scenario(scenario, update_params, NULL)
  for(i in 1:length(param_sets)){
    param_sets[[i]]$beta_day <- param_sets[[i]]$beta_day*beta_factor
  }
  results <- run_param_sets(param_sets[1:10], L=212, N_particles=N_particles, N_threads_internal = 1, N_threads_external=N_threads, silent = TRUE,
                            return_summary_function = function(raw_results, params, scenario){
                              a <- raw_results[params$dust_index$tot_hosp,,212]
                              if(!is.null(dim(a))){
                                return(data.table(tot_hosp=colSums(a), sim=1:dim(a)[2]))
                              }else{
                                return(data.table(tot_hosp=sum(a), sim=1))
                              }
                            })
  return(mean(results$tot_hosp))
}

find_factor <- function(scenario, target, update_params){
  print(scenario$Name)
  o <- optimize(function(x){
    r <- abs(get_tot_hosp_scenario(scenario, x, update_params) - target)
    print(glue::glue("{scenario$Name} - {r}"))
    return(r)}, c(0.9, 1.1), tol=0.002)
  print(glue::glue("{scenario$Name} - {o$minimum}"))
  return(data.table(geo=scenario$Geo, beta_fact=o$minimum))

}


get_geo_params_file("run_params_update_lhs_half.RDS", "half")
get_geo_params_file("run_params_update_lhs_standard.RDS", "standard")
