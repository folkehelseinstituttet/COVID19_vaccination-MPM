
library(dplyr)
library(data.table)
library(metapop)
library(metapopnorge)

# FIT beta

fit_baseline <- function(ihr_fac=1, tag="standard", n=15000, n_best=50, n_threads=200, run_fit=TRUE){


  vac_priority_list <- list(type="history",
                            filename="parameter_files/vac_history_municip.csv",
                            adherence=paste0("parameter_files/vaccine_adherence_H.csv"),
                            start_date=as.Date("2021-01-01"))
  params <- get_params(param_file = "parameter_files/parameters_vaccination.xlsx",
                       N_age=18,
                       L=212,
                       n_vac=3,
                       n_strain=1,
                       dose_file = "parameter_files/vaccine_doses_PM_12.csv",
                       vaccine_file = "parameter_files/vaccine_profile_H_35_12.csv",
                                        #"parameter_files/prioritization_1_18p_ABCDE.csv",
                       vac_priority=vac_priority_list,
                       regional_file = glue::glue("parameter_files/geo_0.csv"),
                       initial_conditions_file = "parameter_files/epi_scenario_municip_ibm.txt",
                       import_file = "parameter_files/import_M.csv",
                       import_age_dist_file = "parameter_files/import_age_dist.txt",
                       mobility_file = "parameter_files/municip_matrix.RDS",
                       regions = TRUE,
                       add_plus_minus_neutral=TRUE)
  dt <- 1
  params$hosp_prob <- params$hosp_prob * ihr_fac
  params <- change_dt(params, dt)
  params <- fix_beta_mode_params(params)
  
  params$beta_mode <- 1
  params$rand_beta_factors <- params$beta_day[1,]/max(params$beta_day)
  
  beta_1 <- fix_beta_large(params, params$S_ini, params$I_ini, 1, beta=params$beta_day[1,], use_eig=TRUE)
  params$log_beta_ini <- log(beta_1*1.1)
  params$log_beta_sd <- 1
  dust_model <- model$new(pars = params,
                          time = 1,
                          n_particles = 1,
                          n_threads = 1
                          )
  
  data <- readRDS("parameter_files/data_regional_prio.RDS")
  age <- data$age[type=="hosp"][order(age)]
  tot_old <- age[9, value]+age[10,value]
  age[age=="80-89", value:=tot_old]
  age <- age[1:9]
  age$prop <- age$value/sum(age$value)

  age_deaths <- data$age[type=="death"][order(age)]
  tot_old <- age_deaths[7, value]+age_deaths[8,value]
  age_deaths <- rbind(data.table(age=c("0-9", "10-19"), value=c(1e-2, 1e-2), type="death"), age_deaths[age=="80-89", value:=tot_old][1:7], fill=TRUE)

  age_icu <- data$age[type=="ICU"][order(age)]
  tot_old <- age_icu[9, value]+age_icu[10,value]
  age_icu <- age_icu[age=="80-89", value:=tot_old][1:9]

  icu_probs <- age_icu$value / age$value
  
  tot_hosp_index <- dust_model$info()$index$tot_hosp
  tot_hosp_inc_index <- dust_model$info()$index$tot_hosp_inc
  hosp_inc_age <- dust_model$info()$index$hosp_inc
  
  
  hosps <- data$incidence[type=="hosp"][order(date)]
  hosps[,week:=1:nrow(hosps) %/% 7 +1]
  weekly <- hosps[week < 30, .(incidence=sum(value)), by=week]
  
  filter_data <- mcstate::particle_filter_data(data = weekly,
                                               time = "week",
                                               rate = 7 / dt)




  params$rand_beta_days <- 1
  params$incidence_steps_measurement <- 7

  hosp_compare <- function(state, observed, pars = NULL) {
    exp_noise <- 1e6
    n_list <- seq(1, pars$n,by=pars$age_groups/2)
    incidence_observed <- observed$incidence
    incidence_modelled <- state[2, ,drop=TRUE]
    tot_hosp <- state[3:dim(state)[1], ,drop=TRUE]
    dist_ll <- 0
    tot_ll <- rep(0, length(incidence_modelled))
                                        #  if(state[1,1] == 233){
    for(i in 1:9){
      
      tot_by_age <- colSums(tot_hosp[n_list + i -1,])
      lambda <- pars$age[i]$prop*observed$incidence + rexp(n=length(incidence_modelled), exp_noise)
      tot_ll <- tot_ll + dpois(x = tot_by_age, lambda = lambda, log = TRUE) 
    }
                                        # }
    ## lambda <- incidence_modelled +
    ##   rexp(n = length(incidence_modelled), rate = exp_noise)
    
                                        #  dpois(x = incidence_observed, lambda = lambda, log = TRUE) +
    return(tot_ll)
  }
  params$age <- age
  n_list <- seq(1, params$n,by=params$age_groups/2)

  index_list <- list(run=c(1,tot_hosp_inc_index, hosp_inc_age), state=unlist(dust_model$info()$index))
  
  index_func <- purrr::partial(function(x, index_list=index_list) index_list, index_list=!!index_list)
  
  filter <- mcstate::particle_filter$new(data = filter_data,
                                         model = model,
                                         n_particles = 15,
                                         n_threads=1,
                                         index=index_func,
                                         compare = hosp_compare)





  basic_params <- copy(params)
  if(run_fit){
    sus_params <- lhs::randomLHS(n, 12)
    sus_params[,1] <- qunif(sus_params[,1], min=0.3, max=0.8)
    sus_params[,2] <- qunif(sus_params[,2], min=1.5, max=2)
    sus_params[,3] <- qunif(sus_params[,3], min=0.7, max=1.3)
    sus_params[,4] <- qunif(sus_params[,4], min=0.65, max=0.9)
    sus_params[,5] <- qunif(sus_params[,5], min=0.65, max=0.9)
    sus_params[,6] <- qunif(sus_params[,6], min=1.0, max=1.4)
    sus_params[,7] <- qunif(sus_params[,7], min=0.8, max=1.1)
    sus_params[,8] <- qunif(sus_params[,8], min=0.8, max=1.1)
    sus_params[,9] <- qunif(sus_params[,9], min=0.9, max=1.6)
    sus_params[,10] <- qunif(sus_params[,10], min=0.4, max=0.9)
    sus_params[,11] <- qunif(sus_params[,11], min=0.4, max=0.8)
    sus_params[,12] <- qunif(sus_params[,12], min=0.3, max=0.7)
    
    
    runs <- list()
    liks <- c()
    
    
    check_params <- function(p, return_history=FALSE){
      scale <- p[4:12]
      params$beta_day <- basic_params$beta_day*c(rep(p[1], 28), rep(p[2], 42), rep(p[3],212-70))*beta_1
      params$susceptibility_asymp <- basic_params$susceptibility_asymp *scale
      params$susceptibility_symp <- basic_params$susceptibility_symp *scale
      
      run <- filter$run(save_history = return_history, pars = params)
      if(return_history){
        return(filter$history(sample.int(filter$n_particles, 1)))
      }
      return(run)
    }
    
    liks <- parallel::mclapply(1:n, function(i) {
      print(i)
      check_params(sus_params[i,])}, mc.cores=n_threads, mc.preschedule = FALSE
      )
    saveRDS(list(lik=liks,
                 params=sus_params), glue::glue("liklihoods_{tag}.RDS"))
    
  }else{
    par_lik <- readRDS(glue::glue("liklihoods_{tag}.RDS"))
    liks <- par_lik$lik
    sus_params <- par_lik$params
    
  }
  
  args <- doBy::which.maxn(unlist(liks), n_best)
  
  colMeans(sus_params[args,])


  param_sets <- list()
  update_param_list <- list()
  
  for(i in args){
    p <- sus_params[i,]
    new_params <- copy(basic_params)
    scale <- c(p[4:12])
    beta_day <- basic_params$beta_day*c(rep(p[1], 28), rep(p[2], 42), rep(p[3],212-70))*beta_1

    update_params <-list(susceptibility_asymp=basic_params$susceptibility_asymp *scale,
                         susceptibility_symp=basic_params$susceptibility_symp *scale,
                         beta_mode=1,
                         hosp_prob=basic_params$hosp_prob,
                         incidence_steps_measurmenet=1,
                         icu_prob=array(icu_probs, dim=c(params$n, params$n_vac, params$n_strain)),
                         beta_day=beta_day)
    update_param_list[[length(update_param_list) + 1]] <- list(
      scale = scale,
      beta_day = c(rep(p[1], 28), rep(p[2], 42), rep(p[3],212-70))*beta_1,
      icu_prob=icu_probs,
      ihr_scale=ihr_fac
    )
      
    param_sets[[length(param_sets) + 1]] <- modifyList(new_params,update_params)

  }
  results <- run_param_sets(param_sets, L=basic_params$N_steps, N_particles=5, N_threads_internal = 1, N_threads_external=50, silent=FALSE)

  for(i in 1:length(args)){

    mean_deaths <- colMeans(results[name==i & t==212, paste0("D_age_", 1:9)])
    factor <- (age_deaths$value / mean_deaths)
    factor[is.na(factor) | is.infinite(factor)] <- 1
    death_prob <- params$prob_death_hosp * factor
    up_list <- list(
      prob_death_hosp = death_prob,
      prob_death_icu = death_prob,
      prob_death_non_hosp = death_prob
    )
    update_param_list[[i]] <- modifyList(update_param_list[[i]], list(
                                                                 death_prob=factor))
    param_sets[[i]] <- modifyList(param_sets[[i]],up_list)

  }
  saveRDS(update_param_list, glue::glue("run_params_update_lhs_{tag}.RDS"))
  results <- run_param_sets(param_sets, L=basic_params$N_steps, N_particles=5, N_threads_internal = 1, N_threads_external=50, silent=FALSE)
  saveRDS(results, glue::glue("calibration_results_lhs_{tag}.RDS"))
}

fit_baseline(1, "standard", 15000, 50, 200, TRUE)


fit_baseline(0.5, "half", 15000, 50, 200, TRUE)


mat_to_df <- function(values){
  out <- list()
  for(j in 1:dim(values)[1]){
    out[[j]] <- data.table(t=(1:ncol(values)), val=values[j, ],sim=j)
    
  }
  return(rbindlist(out))
}


plot_ages <- function(h, filter_data, x_age){

  sim <- list()
  for(i in 1:9){
    sim[[i]] <- mat_to_df(colSums(h[hosp_inc_age[n_list + i -1],,])) %>% mutate(age=i,
                                                                                data=rep(c(0,filter_data$incidence*x_age[i]$prop), dim(h)[2]))
  }
  all_sims <- rbindlist(sim)
  all_sims <- rbind(all_sims, all_sims[, .(val=sum(val),data=sum(data), age="tot"), by=.(t,sim)])
  all_sims <- rbind(all_sims, data.table(t=1:dim(h)[3], val=as.numeric(t(exp(h[8,,]))), sim=rep(1:dim(h)[2], each=dim(h)[3]), data=NA, age="beta"))
  ggplot(all_sims) + geom_line(aes(x=t, y=val, group=sim)) + geom_point(aes(x=t, y=data), color="red", size=5) + facet_wrap(.~age, scales="free_y")
}

plot_tot_ages <- function(h, filter_data, x_age){

  sim <- list()
  for(i in 1:9){
    sim[[i]] <- data.frame(age=i, sim=1:1:dim(h)[2], data=age[i]$prop*sum(filter_data$incidence), val=rowSums(colSums(h[hosp_inc_age[n_list + i -1],,])))
  }
  all_sims <- rbindlist(sim)
  all_sims <- rbind(all_sims, all_sims[, .(val=sum(val),data=sum(data), age="tot"), by=.(sim)])
  ggplot(all_sims) + geom_col(aes(x=age, y=val, group=sim), position="dodge") + geom_point(aes(x=age, y=data), color="red", size=5) 
}



## hosp_compare <- function(state, observed, pars = NULL) {
##   exp_noise <- 1e6

##   incidence_observed <- observed$incidence
##   inc_hosp <- state[hosp_inc_age, ,drop=TRUE]
##   incidence_modelled <- colSums(inc_hosp)
##   inc_by_age <- array(colSums(inc_hosp[n_list,]), dim=c(1,dim(inc_hosp)[2]))
##   for(i in 2:9){
##     inc_by_age <- rbind(inc_by_age, colSums(inc_hosp[n_list + i -1,]))
##   }
##   dist_ll <- c()
##   for(j in 1:dim(inc_by_age)[2]){
##     dist_ll <- c(dist_ll,dmultinom(inc_by_age[,j], size=incidence_modelled[j], prob=age$prop, log=T))
##   }
  
##   lambda <- incidence_modelled +
##     rexp(n = length(incidence_modelled), rate = exp_noise)

##   dpois(x = incidence_observed, lambda = lambda, log = TRUE) #+ dist_ll
## }

## hosp_compare2 <- function(state, observed, pars = NULL) {
##   exp_noise <- 1e6

##   incidence_observed <- observed$incidence
##   inc_hosp <- state[hosp_inc_age, ,drop=TRUE]
##   incidence_modelled <- colSums(inc_hosp)

##   tot_ll <- rep(0, length(incidence_modelled))
##   for(i in 1:9){
##     inc_by_age <- colSums(inc_hosp[n_list + i -1,])
##     lambda <- incidence_observed*age$prop[i] +
##       rexp(n = length(incidence_modelled), rate = exp_noise)

##     tot_ll <- tot_ll + dpois(x = inc_by_age, lambda = lambda, log = TRUE) 
##   }
##   return(tot_ll)
## }

## new_params <- copy(params)

## new_params$susceptibility_asymp <- params$susceptibility_asymp * c(0.6,0.6,1,0.8,1,1,1,0.5,0.5)
## new_params$susceptibility_symp <- params$susceptibility_symp * c(0.6,0.6,1,1,1,1,1,0.5,0.5)

## 

## filter$run(save_history = TRUE, pars = new_params)


## state <- filter$state()

## tot_hosp <- state[hosp_inc_age, ,drop=TRUE]
## tot_hosp <- state[tot_hosp_index, ,drop=TRUE]
## tot_by_age <- array(colSums(tot_hosp[n_list,]), dim=c(1,dim(tot_hosp)[2]))
## for(i in 2:9){
##   tot_by_age <- rbind(tot_by_age, colSums(tot_hosp[n_list + i -1,]))
## }

## rowMeans(tot_by_age)
