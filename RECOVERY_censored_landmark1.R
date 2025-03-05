# Set seed for reproducibility
set.seed(456)

library(survival)
library(dplyr)
library(readr)

#------------------------------------------------------------
# 1. Read in the cumulative incidence data
#------------------------------------------------------------
data_int <- read_csv("references/AVATAR_Sx.csv") %>% 
  rename(time = 1, incidence = 2) %>% 
  mutate(across(everything(), as.numeric)) %>% 
  na.omit()

data_clinical <- read_csv("references/AVATAR_Ctl.csv") %>% 
  rename(time = 1, incidence = 2) %>% 
  mutate(across(everything(), as.numeric)) %>% 
  na.omit()

#------------------------------------------------------------
# 1A. Ensure a row at time=0 if missing for interpolation
#------------------------------------------------------------
add_zero_incidence <- function(df) {
  if (!any(df$time == 0)) {
    df <- rbind(data.frame(time = 0, incidence = 0), df)
  }
  df <- df %>% arrange(time)
  return(df)
}

data_int <- add_zero_incidence(data_int)
data_clinical <- add_zero_incidence(data_clinical)

#------------------------------------------------------------
# 2. At-risk numbers & time points
#------------------------------------------------------------
time_points <- c(0, 1, 2, 3, 4, 5)
at_risk_int <- c(113, 97, 76, 65, 51, 18)  # Intervention
at_risk_clin <- c(111, 97, 71, 57, 40, 17) # Clinical

#------------------------------------------------------------
# 3. Simulate IPD with precise event and censoring handling
#------------------------------------------------------------
simulate_ipd_with_risk <- function(df, at_risk, time_points, N, group_label) {
  ipd_list <- list()
  df <- df %>% arrange(time)
  
  for (i in seq_len(length(time_points) - 1)) {
    t_start <- time_points[i]
    t_end <- time_points[i + 1]
    
    ci_start <- approx(df$time, df$incidence, xout = t_start, rule = 2)$y
    ci_end <- approx(df$time, df$incidence, xout = t_end, rule = 2)$y
    
    cum_events_start <- N * ci_start / 100
    cum_events_end <- N * ci_end / 100
    
    events_interval <- min(round(cum_events_end - cum_events_start), at_risk[i] - at_risk[i + 1])
    censor_interval <- (at_risk[i] - at_risk[i + 1]) - events_interval
    
    event_times <- if (events_interval > 0) runif(events_interval, min = t_start, max = t_end) else numeric(0)
    event_df <- data.frame(time = event_times, event = rep(1, length(event_times)))
    
    censor_times <- if (censor_interval > 0) runif(censor_interval, min = t_start, max = t_end) else numeric(0)
    censor_df <- data.frame(time = censor_times, event = rep(0, length(censor_times)))
    
    ipd_list[[i]] <- bind_rows(event_df, censor_df)
  }
  
  t_last <- time_points[length(time_points)]
  remaining <- at_risk[length(time_points)]
  last_df <- data.frame(time = rep(t_last, remaining), event = rep(0, remaining))
  
  ipd <- bind_rows(ipd_list, last_df)
  ipd$group <- group_label
  
  return(ipd)
}

ipd_int <- simulate_ipd_with_risk(data_int, at_risk_int, time_points, at_risk_int[1], "intervention")
ipd_clin <- simulate_ipd_with_risk(data_clinical, at_risk_clin, time_points, at_risk_clin[1], "control")

ipd_all <- bind_rows(ipd_int, ipd_clin)
ipd_all$group <- factor(ipd_all$group, levels = c("control", "intervention"))

#------------------------------------------------------------
# 5. Full follow-up Cox model
#------------------------------------------------------------
cox_model_full <- coxph(Surv(time, event) ~ group, data = ipd_all)
print(summary(cox_model_full))

#------------------------------------------------------------
# 6. Censored at 1 year
#------------------------------------------------------------
ipd_censored <- ipd_all
ipd_censored$time <- pmin(ipd_censored$time, 1)
ipd_censored$event[ipd_censored$time == 1] <- 0
cox_model_censored <- coxph(Surv(time, event) ~ group, data = ipd_censored)
print(summary(cox_model_censored))

#------------------------------------------------------------
# 7. Landmark analysis starting from 1 year
#------------------------------------------------------------
ipd_landmark <- ipd_all %>% filter(time > 1)
ipd_landmark$time <- ipd_landmark$time - 1  # Reset time from 1 year
cox_model_landmark <- coxph(Surv(time, event) ~ group, data = ipd_landmark)
print(summary(cox_model_landmark))
