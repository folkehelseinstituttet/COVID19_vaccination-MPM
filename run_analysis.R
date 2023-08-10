
library(dplyr)
library(data.table)
library(metapop)
library(metapopnorge)
source("run_functions.R")


run_scenario <- function(scenario, N_particles, N_threads, update_params, geo_params, hist_file="parameter_files/vac_history_municip.csv"){

  param_sets <- param_sets_from_scenario(scenario, update_params, geo_params, hist_file=hist_file)
  results <- run_param_sets(param_sets, L=212, N_particles=N_particles, N_threads_internal = 1, N_threads_external=N_threads, silent = TRUE)


  ret <- results[, 
               .(time = t,
                 date=as.Date("2021-01-01")+t,
                 name=scenario$Name,
                 sim,
                 inc_I = incidence,
                 cum_I = tot_infected,
                 cum_I1 = tot_infected_age_1,
                 cum_I2 = tot_infected_age_2,
                 cum_I3 = tot_infected_age_3,
                 cum_I4 = tot_infected_age_4,
                 cum_I5 = tot_infected_age_5,
                 cum_I6 = tot_infected_age_6,
                 cum_I7 = tot_infected_age_7,
                 cum_I8 = tot_infected_age_8,
                 cum_I9 = tot_infected_age_9,
                 cum_D1 = D_age_1,
                 cum_D2 = D_age_2,
                 cum_D3 = D_age_3,
                 cum_D4 = D_age_4,
                 cum_D5 = D_age_5,
                 cum_D6 = D_age_6,
                 cum_D7 = D_age_7,
                 cum_D8 = D_age_8,
                 cum_D9 = D_age_9,
                 cum_H1 = tot_hosp_age_1,
                 cum_H2 = tot_hosp_age_2,
                 cum_H3 = tot_hosp_age_3,
                 cum_H4 = tot_hosp_age_4,
                 cum_H5 = tot_hosp_age_5,
                 cum_H6 = tot_hosp_age_6,
                 cum_H7 = tot_hosp_age_7,
                 cum_H8 = tot_hosp_age_8,
                 cum_H9 = tot_hosp_age_9,
                 cum_V1 = tot_vac_age_1,
                 cum_V2 = tot_vac_age_2,
                 cum_V3 = tot_vac_age_3,
                 cum_V4 = tot_vac_age_4,
                 cum_V5 = tot_vac_age_5,
                 cum_V6 = tot_vac_age_6,
                 cum_V7 = tot_vac_age_7,
                 cum_V8 = tot_vac_age_8,
                 cum_V9 = tot_vac_age_9,
                 cum_VA1 = tot_vac_adm_age_1,
                 cum_VA2 = tot_vac_adm_age_2,
                 cum_VA3 = tot_vac_adm_age_3,
                 cum_VA4 = tot_vac_adm_age_4,
                 cum_VA5 = tot_vac_adm_age_5,
                 cum_VA6 = tot_vac_adm_age_6,
                 cum_VA7 = tot_vac_adm_age_7,
                 cum_VA8 = tot_vac_adm_age_8,
                 cum_VA9 = tot_vac_adm_age_9,
                 cum_resp1 = tot_resp_age_1,
                 cum_resp2 = tot_resp_age_2,
                 cum_resp3 = tot_resp_age_3,
                 cum_resp4 = tot_resp_age_4,
                 cum_resp5 = tot_resp_age_5,
                 cum_resp6 = tot_resp_age_6,
                 cum_resp7 = tot_resp_age_7,
                 cum_resp8 = tot_resp_age_8,
                 cum_resp9 = tot_resp_age_9,
                 cum_H_minus=tot_hosp_minus,
                 cum_H_neutral=tot_hosp_neutral,
                 cum_H_plus=tot_hosp_plus,
                 cum_resp_plus=tot_resp_plus,
                 cum_resp_minus=tot_resp_minus,
                 cum_resp_neutral=tot_resp_neutral,
                 cum_I_minus=tot_infected_minus,
                 cum_I_neutral=tot_infected_neutral,
                 cum_I_plus=tot_infected_plus,
                 cum_D_minus=D_minus,
                 cum_D_neutral=D_neutral,
                 cum_D_plus=D_plus,
                 cum_vac_minus=tot_vac_adm_minus,
                 cum_vac_neutral=tot_vac_adm_neutral,
                 cum_vac_plus=tot_vac_adm_plus,
                 cum_vac_dose_1=tot_vac_vac_2,
                 cum_vac_dose_2=tot_vac_vac_3,
                 prev_I=I,
                 prev_A=A,
                 prev_H = hosp,
                 cum_resp=tot_resp,
                 cum_H = tot_hosp,
                 cum_D = D
                 )]
  print(glue::glue("Finishing {scenario$Name}"))
  return(ret)
}

run_scenarios <- function(param_update_file, geo_file, tag="standard", N_particles=20, N_threads_int=10, N_threads_ext=1, scenario_file="parameter_files/scenarios_N.csv"){
  scenarios <- fread(scenario_file)
  update_params <- readRDS(param_update_file)[1:10]
  geo_params <- readRDS(geo_file)

#  geo_params[, geo:=scenarios[Regional==0, Geo][2:15]]
  geo_params <- rbind(data.table(geo=0, beta_fact=1), geo_params)

  
  list_results <- parallel::mclapply(1:nrow(scenarios), function(i) {
       print(paste(i, " / ", nrow(scenarios)))
       run_scenario(scenarios[i],N_particles,N_threads_int, update_params, geo_params)
   },mc.cores=N_threads_ext, mc.preschedule=FALSE)
  all_results <- rbindlist(list_results)
  print(glue::glue("SAVE {tag}"))
  saveRDS(all_results, glue::glue("all_results-regional_{tag}.RDS"))
}

run_adapted_national_scenario <- function(param_update_file, geo_file, tag="standard", N_particles=10, N_threads_int=10, N_threads_ext=1, scenario_file="parameter_files/scenarios_N.csv"){
  scenarios <- fread(scenario_file)
  N_files <- 100
  update_params <- readRDS(param_update_file)[1:N_particles]
  geo_params <- readRDS(geo_file)

#  geo_params[, geo:=scenarios[Regional==0, Geo][2:15]]
  geo_params <- rbind(data.table(geo=0, beta_fact=1), geo_params)

  hist_files <- paste0("parameter_files/new_version_VaksinertPerKommune-seed_",1:N_files,".csv")
  list_results <- parallel::mclapply(1:length(hist_files), function(i){
  	       print(i)
  	       run_scenario(scenarios[1],1,N_threads_int, update_params, geo_params, hist_file=hist_files[i]) %>% mutate(name=i)
	       },
                                     mc.cores=N_threads_ext, mc.preschedule=FALSE)
  all_results <- rbindlist(list_results)
  all_results[, sim:=paste(name, sim)  ]
  saveRDS(all_results, glue::glue("all_results-adapted-national_{tag}.RDS"))
}

fix_hist_file_format <- function(){
  hist_files <- paste0("parameter_files/VaksinertPerKommune-seed_",1:100,"-municip_code.csv") 
  for(hf in hist_files){
    print(hf)
    df <- fread(hf)
    colnames(df) <- c("V1","date_vax","municip_code", "aldersgruppe", "dosenummer", "vaccine","risikogruppe", "n")
    df[, dosenummer:=as.character(dosenummer)]
    df[dosenummer=="1", dosenummer:="første"]
    df[,risikogruppe:=risikogruppe-1]
    df[, aldersgruppe:=as.integer(factor(df$aldersgruppe))+1]
   fwrite(df, paste0("parameter_files/new_",substr(hf,17,nchar(hf))))
  }
}
#fix_hist_file_format()
#run_scenarios("run_params_update_lhs_half.RDS","geo_params_half.RDS" , N_threads_int=10, N_threads_ext=25, tag="half")
#run_scenarios("run_params_update_lhs_standard.RDS", "geo_params_standard.RDS", N_threads_int=10, N_threads_ext=25)


#run_scenarios("run_params_update_lhs_half.RDS","geo_params_half.RDS" , N_threads_int=10, N_threads_ext=25,N_particles=20,
#               scenario_file = "parameter_files/scenarios_N_deltaT.csv", tag="deltaT")

run_adapted_national_scenario("run_params_update_lhs_half.RDS","geo_params_half.RDS" , N_threads_int=10, N_threads_ext=25, tag="half")

## res <- run_scenario(scenarios[1],4,4)
## data <- readRDS("parameter_files/data_regional_prio.RDS")


check_vax <- function(scenario, ret_params=FALSE){
  print(scenario$Name)
  if(scenario$R_const == -1){
    vac_priority_list <- NULL
  }else if(scenario$R_const == 0){
    vac_priority_list <- list(type="story",
                             filename="parameter_files/vac_history_municip.csv",
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
  if(ret_params) return(params)

  n_list <- seq(1, params$n,by=params$age_groups/2)
  ret <- data.table(name=scenario$Name,
                    over_0=sum(params$vaccinations[,n_list,1]),
                    over_10=sum(params$vaccinations[,n_list+1,1]),
                    over_20=sum(params$vaccinations[,n_list+2,1]),
                    over_30=sum(params$vaccinations[,n_list+3,1]),
                    over_40=sum(params$vaccinations[,n_list+4,1]),
                    over_50=sum(params$vaccinations[,n_list+5,1]),
                    over_60=sum(params$vaccinations[,n_list+6,1]),
                    over_70=sum(params$vaccinations[,n_list+7,1]),
                    over_80=sum(params$vaccinations[,n_list+8,1])
                    )
                    
                    
  return(ret)
}
#list_results <- parallel::mclapply(1:nrow(scenarios), function(i) check_vax(scenarios[i]), mc.cores=10, mc.preschedule=FALSE)

#d <- rbindlist(list_results)


## params1 <- check_vax(scenarios[1], ret_params=T)
## params2 <- check_vax(scenarios[3], ret_params=T)


## sum(params2$vaccinations[40,,1])
## sum(params1$vaccinations[40,,1])

## dim(params2$vaccinations)

## ## cum_h <- data$incidence[type=="hosp"][order(date)]
## ## cum_h[, cum:=cumsum(value)]
## ## library(ggplot2)
## ## ggplot(res[name=="1-18p-65p-18p-ABCDE-H-35-H-12-PM-0p1-10-0-1-50-0-M-0p66"]) + geom_line(aes(x=date, y=cum_H, group=sim)) + geom_point(aes(x=date, y=cum), color="blue", data=cum_h)



## ## plot_inc(res, params, "hosp_incidence") + geom_point(aes(x=date - as.Date("2021-01-01"), y=value, color="red"), data=data$incidence[type=="hosp"])
## ##                                         #colnames(res)
## ##   library(ggplot2)



## ## intervals <- c(1,c(21, 17, 21, 24,19,23, 21, 25,45, 25)/params$dt)

## m <- function(r, b) b/(b+r)
## R_from_r <- function(r, params){
##   c2 <- params$pre_sympt_infect*params$pre_sympt_period / (params$pre_sympt_infect*params$pre_sympt_period + params$infectious_period)
##   c3 <- 1 - c2
##   R <- 1/(c2*m(r, 1/params$latent_period)*m(r, 1/params$pre_sympt_period) + c3*m(r, 1/params$latent_period)*m(r, 1/params$pre_sympt_period)*m(r, 1/params$infectious_period))
##  return(R)
## }
## est_R <- function(inc, params){
##   m <- lm(log(I) ~ t, data=data.frame(t=1:length(inc), I=inc))
##   r <- coef(m)[2]/params$dt
##   return(R_from_r(r, params))
## }


## #intervals <- c(1, rep(10, 24))
## Rs <- list()
## for(i in 2:length(intervals)){
##   t1 <- sum(intervals[1:(i-1)])
##   t2 <- sum(intervals[1:i])
##   inc <- res[, mean(incidence), by=t][t1:t2]
##   Rval <- est_R(inc$V1, params)
##   Rs[[length(Rs) + 1]] <- data.frame(t=c(t1, t2-0.5), R=Rval)
## }
## Rs_df <- rbindlist(Rs)
## #plot(Rs)

## plot_inc <- function(res, params, key){
##   res[, day:=floor(time)]
##   new_res <- res[, sum(get(key)), by=.(day, sim)]
##   q1 <- ggplot(new_res) + geom_line(aes(x=day, y=V1, group=sim))
## }

## q1 <- plot_inc(res, params, "incidence")
##                                         #q1 <- ggplot(res) + geom_line(aes(x=t, y=I, group=sim))
## q1 <- plot_inc(res, params, "hosp_incidence") + geom_point(aes(x=date - as.Date("2021-01-01"), y=value, color="red"), data=data$incidence[type=="hosp"])
## q2 <- ggplot(res) + geom_line(aes(x=t, y=get("Rt"), group=sim)) + geom_line(aes(x=x, y=y), data=data.frame(x=1:max(res$t), y=rep(R_vals, each=4)[1:max(res$t)]), color="red") +  geom_line(aes(x=t, y=R),span=0.1, data=Rs_df, color="blue") + geom_hline(yintercept=1)
## gridExtra::grid.arrange(q1,q2)

## res[, sum(incidence), by=sim]

## res[t==243, tot_infected]

## data <- readRDS("parameter_files/data_regional_prio.RDS")


## params$beta_day[,1]



                                        # Tasks
## 1. beta_day
## 1.1 Regional R-values - DONE
## 1.2 Fix Betas - DONE
## 2. Initial conditions - DONE
## 4. Read scenario file - DONE
## 5. Look at mobility
## 6. Refine to output files - DONE

