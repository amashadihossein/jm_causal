yspline_get <- function(){
  knots <- c(0.25, 0.5, 0.75)
  theta <-cbind(
    grp0 = c(0.1, 0.9, 0.4, 0.2, 0.1, 0.1, 0.1),
    grp1 = c(0.1, 0.3, 0.5, 0.2, 0.1, 0.1, 0.1))
  
  t_obs <- 1:30
  t_obs.normalize <- (t_obs - min(t_obs))/(max(t_obs) - min(t_obs))
  
  basis <- splines::bs(x = t_obs.normalize, knots = knots, degree = 3, 
                       Boundary.knots = c(0, 1), intercept = TRUE)
  
  
  y_spline <- basis %*% theta
  colnames(y_spline) = c("0", "1")
  rownames(y_spline) = as.character(t_obs)
  return(y_spline)
}

aveProfile<- function(n, time, grp, unitnoise){
  # browser()
  time = as.character(time)
  grp = as.character(grp)
  sapply(1:n, function(ii) y_spline[time[ii], grp[ii]] + unitnoise[ii])
}


get_event_p <- function(tobs, biomarker_tobs, 
                        event,
                        w, x1, betaW, betaX1, 
                        alpha = 3,
                        shape_wb = 1.5,
                        t.end){
  # i_debug <<- i_debug+1
  # j_debug <<- 0
  
  #' p(time=t) = p(time <= t) - p(time < t)
  #'           = 1-S(t) - (1-S(t-1))
  #'           = S(t-1) - S(t)
  #'           = exp(-H(t-1)) - exp(-H(t)) 
  
  # browser()
  
  if(max(tobs)<2)
    return(0)
  
  # kind of unnecessary
  # order_tobs <- order(tobs)
  # tobs <- tobs[order_tobs]
  tobs_n <- tobs/t.end
  y_bar <- biomarker_tobs#[order_tobs]
  
  
  h <- function (t,func, gamma_x, beta_w, scale_wb = 1) {
    # browser()
    f_t <- as.vector(func(t))
    exp(log(scale_wb) +
          log(shape_wb) + (shape_wb - 1) * log(t) +  
          gamma_x + # baseline covariate
          beta_w + # confounder
          f_t * alpha )
  }
  
  #iterate for n sample
  event_p <- sapply(1:nrow(y_bar),FUN = function(nn){
    # j_debug <<- j_debug+1
    
    
    # if(i_debug == 60 & j_debug == 89)
    # browser()
    
    death_nn <- event[nn,]
    y_bar_nn <- y_bar[nn,]
    
    # cat("i_debug:j_debug is", i_debug, ":", j_debug,":")
    # cat("death_nn is", death_nn,"\n")
    
    if(is.na(death_nn))
      return(1)
    
    if(death_nn == 1)
      return(1)
    
    gamma_x <- x1[nn] * betaX1 # baseline covariate
    beta_w <- w[nn] * betaW # confounder
    
    f <- approxfun(x = tobs_n, y = y_bar_nn)
    
    tmp <- try(integrate(h, lower = 1/t.end, upper = rev(tobs_n)[1], func = f,
                         gamma_x = gamma_x, beta_w = beta_w)$value)
    if(class(tmp) == "try-error"){
      
      tmp <- integrate(h, lower = 1/t.end, upper = rev(tobs_n)[1], func = f,
                       gamma_x = gamma_x, beta_w = beta_w,
                       stop.on.error = F)$value
      warning("Integration had failure!!!!")
      
    }
    log_s_t_i <- -tmp
    #---
    tmp <- try(integrate(h, lower = 1/t.end, upper = rev(tobs_n)[2], func = f,
                         gamma_x = gamma_x, beta_w = beta_w)$value)
    if(class(tmp) == "try-error"){
      
      tmp <- integrate(h, lower = 1/t.end, upper = rev(tobs_n)[2], func = f,
                       gamma_x = gamma_x, beta_w = beta_w,
                       stop.on.error = F)$value
      warning("Integration had failure!!!!")
      
    }
    log_s_tm1_i <- -tmp
    
    event_p_i <- exp(log_s_tm1_i) - exp(log_s_t_i) 
    return(event_p_i)
  }) 
  
  # browser()
  
  return(event_p)
  
}


make_lt <- function(fit) {
  
  # arrange the lt data for all rows but the first
  most_rows <-
    tibble(time = fit$time) %>% 
    mutate(time_int = str_c("[", time, ", ", time + 1, ")"), 
           n_risk   = fit$n.risk, 
           n_event  = fit$n.event) %>% 
    mutate(n_censored   = n_risk - n_event - lead(n_risk, default = 0),
           hazard_fun   = n_event / n_risk,
           survivor_fun = fit$surv)
  
  # define the values for t = 2 and t = 1
  time_1 <- fit$time[1]
  time_0 <- time_1 - 1
  
  # define the values for the row for which t = 1
  row_1 <-
    tibble(time         = time_0, 
           time_int     = str_c("[", time_0, ", ", time_1, ")"),
           n_risk       = fit$n.risk[1],
           n_event      = NA,
           n_censored   = NA,
           hazard_fun   = NA, 
           survivor_fun = 1)
  
  # make the full life table
  lt <-
    bind_rows(row_1,
              most_rows)
  
  lt
  
}

emp_hz_get <- function(sim_dat){
  # browser()
  fit11.1 <- 
    survfit(data = sim_dat %>% filter(grps == 0),
            Surv(t, event = event) ~ 1)
  
  fit11.2 <- 
    survfit(data = sim_dat %>% filter(grps == 1),
            Surv(t, event = event) ~ 1)
  
  fit11.3 <- 
    survfit(data = sim_dat,
            Surv(t, event = event) ~ 1)
  
  
  lt <-
    bind_rows(make_lt(fit11.1),
              make_lt(fit11.2),
              make_lt(fit11.3)) %>%
    mutate(grp = factor(rep(c("grp = 0", "grp = 1", "overall"), each = n() / 3))) %>% 
    select(grp, everything())
}
