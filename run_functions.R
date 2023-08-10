generate_param_sets <- function(params, update_params, geo_param){
  paramsets <- list()
  for(i in 1:length(update_params)){
    beta_day <- params$beta_day*update_params[[i]]$beta_day*geo_param
    death_prob <- params$prob_death_hosp * update_params[[i]]$death_prob
    
    update_list <-list(susceptibility_asymp=params$susceptibility_asymp *update_params[[i]]$scale,
                         susceptibility_symp=params$susceptibility_symp *update_params[[i]]$scale,
                         hosp_prob=params$hosp_prob*update_params[[i]]$ihr_scale,
                         incidence_steps_measurmenet=1,
                         icu_prob=array(update_params[[i]]$icu_prob, dim=c(params$n, params$n_vac, params$n_strain)),
                         prob_death_hosp = death_prob,
                         prob_death_icu = death_prob,
                         prob_death_non_hosp = death_prob,
                         beta_day=beta_day)
    paramsets[[i]] <- modifyList(params, update_list)
  }
  return(paramsets)
}

param_sets_from_scenario <- function(scenario, update_params, geo_params, hist_file="parameter_files/vac_history_municip.csv"){
  if(scenario$R_const == -1){
    vac_priority_list <- NULL
  }else if(scenario$R_const == 0){
    vac_priority_list <- list(type="history",
                             filename=hist_file,
                             adherence=paste0("parameter_files/", scenario$Name_adherence, ".csv"),
                             start_date=as.Date("2021-01-01"))
  }else if(scenario$R_const == 1){
    vac_priority_list <- list(type="scenario",
                              amount=scenario$Regional,
                              start_day=scenario$Regional_start,
                              adherence=paste0("parameter_files/", scenario$Name_adherence, ".csv"),
                              file=paste0("parameter_files/", scenario$Name_priority_1, ".csv")
                              )
  }
  params <- get_params(param_file = "parameter_files/parameters_vaccination.xlsx",
                       N_age=18,
                       L=212,
                       n_vac=3,
                       n_strain=1,
                       dose_file = "parameter_files/vaccine_doses_PM_12.csv",
                       vaccine_file = "parameter_files/vaccine_profile_H_35_12.csv",
                                        #"parameter_files/prioritization_1_18p_ABCDE.csv",
                       vac_priority=vac_priority_list,
                       regional_file = glue::glue("parameter_files/geo_{scenario$Geo}.csv"),
                       initial_conditions_file = "parameter_files/epi_scenario_municip_ibm.txt",
                       import_file = "parameter_files/import_M.csv",
                       import_age_dist_file = "parameter_files/import_age_dist.txt",
                       mobility_file = "parameter_files/municip_matrix.RDS",
                       regions = TRUE,
                       add_plus_minus_neutral=TRUE)
  params$vax_type <- 1
  #print(dim(params$S_ini)
  params$N_steps <- 212
  
                                        #params$beta_day[,] <- 1
  ## beta_1 <- fix_beta_large(params, params$S_ini, params$I_ini, 1, beta=params$beta_day[1,], use_eig=TRUE)
  ## params$beta_day <- params$beta_day*beta_1*R_vals*factors
  geo_param <- 1
  if(!is.null(geo_params)){
    geo_param <- geo_params[geo==scenario$Geo, beta_fact]
  }
  param_sets <- generate_param_sets(params, update_params, geo_param)
  return(param_sets)
}
