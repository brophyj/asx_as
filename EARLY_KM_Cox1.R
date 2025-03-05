# Set seed for reproducibility
set.seed(456)

# Load required libraries
library(survival)
library(dplyr)
library(survminer)
library(readr)

#------------------------------------------------------------
# 1. Read in the cumulative incidence data for each group
#------------------------------------------------------------
data_tavr <- read_csv("references/EARLY_TAVR.csv")  # Contains two columns: time and incidence (%)
data_clinical <- read_csv("references/EARLY_Ctl.csv")  # Contains time and incidence (%)

# Rename columns for clarity
data_tavr <- data_tavr %>% rename(time = 1, incidence = 2) %>% mutate(across(everything(), as.numeric)) %>% na.omit()
data_clinical <- data_clinical %>% rename(time = 1, incidence = 2) %>% mutate(across(everything(), as.numeric)) %>% na.omit()

#------------------------------------------------------------
# 2. Define the at-risk numbers and the corresponding time points
#------------------------------------------------------------
# Time points (in years) for which numbers at risk are available:
time_points <- c(0, 1, 2, 3, 4, 5)

# Numbers at risk from the published table:
at_risk_tavr <- c(455, 390, 363, 285, 141, 103)            # TAVR group
at_risk_clin <- c(446, 305, 266, 187, 117, 46)                # Clinical Surveillance group

# The baseline (N) for each group is the number at risk at time 0:
N_tavr <- at_risk_tavr[1]
N_clin <- at_risk_clin[1]

#------------------------------------------------------------
# 3. Define a function to reconstruct IPD using both the cumulative incidence curve and the numbers at risk
#------------------------------------------------------------
simulate_ipd_with_risk <- function(df, at_risk, time_points, N, group_label) {
  # df: data.frame with columns 'time' and 'incidence' (cumulative %)
  # at_risk: vector of numbers at risk at the specified time_points (length = number of time_points)
  # time_points: vector of time points (e.g., c(0,1,2,3,4,5))
  # N: baseline sample size (should equal at_risk[1])
  # group_label: group name (e.g., "TAVR" or "Clinical Surveillance")
  
  # Ensure the cumulative incidence data are sorted
  df <- df %>% arrange(time)
  
  ipd_list <- list()
  
  # Loop over intervals defined by the provided time points
  for (i in 1:(length(time_points) - 1)) {
    t_start <- time_points[i]
    t_end <- time_points[i + 1]
    
    # Interpolate the cumulative incidence (%) at the start and end of the interval
    ci_start <- approx(x = df$time, y = df$incidence, xout = t_start, rule = 2)$y
    ci_end   <- approx(x = df$time, y = df$incidence, xout = t_end,   rule = 2)$y
    
    # The expected number of cumulative events (from baseline) at these times:
    cum_events_start <- N * ci_start / 100
    cum_events_end   <- N * ci_end / 100
    
    # Estimated number of new events in this interval (rounding to the nearest integer)
    events_interval <- round(cum_events_end - cum_events_start)
    
    # The observed drop in numbers at risk (from your published numbers) between these time points:
    risk_drop <- at_risk[i] - at_risk[i + 1]
    
    # We assume that the risk drop is due to events plus censoring.
    # So the estimated number censored in this interval is:
    censor_interval <- max(risk_drop - events_interval, 0)
    
    # Simulate event times uniformly between t_start and t_end for the estimated events
    if (events_interval > 0) {
      event_times <- runif(events_interval, min = t_start, max = t_end)
      event_df <- data.frame(time = event_times, event = 1)
    } else {
      event_df <- data.frame(time = numeric(0), event = numeric(0))
    }
    
    # Simulate censoring times uniformly between t_start and t_end for the estimated censorings
    if (censor_interval > 0) {
      censor_times <- runif(censor_interval, min = t_start, max = t_end)
      censor_df <- data.frame(time = censor_times, event = 0)
    } else {
      censor_df <- data.frame(time = numeric(0), event = numeric(0))
    }
    
    # Combine event and censoring data for this interval
    ipd_list[[i]] <- bind_rows(event_df, censor_df)
  }
  
  # For the final time point, assume that all patients remaining in the risk set are censored at t_last.
  t_last <- time_points[length(time_points)]
  remaining <- at_risk[length(time_points)]
  if (remaining > 0) {
    last_df <- data.frame(time = rep(t_last, remaining), event = 0)
  } else {
    last_df <- data.frame(time = numeric(0), event = numeric(0))
  }
  
  # Combine all intervals
  ipd <- bind_rows(ipd_list, last_df)
  ipd$group <- group_label
  
  return(ipd)
}

#------------------------------------------------------------
# 4. Reconstruct the IPD for each group using the numbers at risk
#------------------------------------------------------------
ipd_tavr <- simulate_ipd_with_risk(data_tavr, at_risk_tavr, time_points, N_tavr, "TAVR")
ipd_clin <- simulate_ipd_with_risk(data_clinical, at_risk_clin, time_points, N_clin, "Clinical Surveillance")

# Combine the reconstructed data from both groups
ipd_all <- bind_rows(ipd_tavr, ipd_clin)
ipd_all$group <- factor(ipd_all$group, levels = c("Clinical Surveillance", "TAVR"))

# Inspect a summary of the reconstructed dataset
summary(ipd_all)

#------------------------------------------------------------
# 5. Fit a Cox Proportional Hazards Model Using the Reconstructed IPD
#------------------------------------------------------------
cox_model <- coxph(Surv(time, event) ~ group, data = ipd_all)
summary(cox_model)

#------------------------------------------------------------
# 6. Plot Kaplan-Meier Survival Curves (Overall)
#------------------------------------------------------------
fit_overall <- survfit(Surv(time, event) ~ group, data = ipd_all)
ggsurvplot(fit_overall, data = ipd_all, risk.table = TRUE, conf.int = TRUE,
           ggtheme = theme_minimal(), legend.title = "Group")

#------------------------------------------------------------
# 7. # Landmark Analysis Starting at 1 Year

landmark_time <- 1.0  # 1 year landmark

# Exclude patients who had an event before the landmark time
ipd_landmark <- ipd_all %>% filter(time > landmark_time)

# Reset the time variable so that time = 0 corresponds to the landmark (i.e. subtract 1 year)
ipd_landmark <- ipd_landmark %>% mutate(time = time - landmark_time)

# Fit the Cox model on the landmark cohort
cox_landmark <- coxph(Surv(time, event) ~ group, data = ipd_landmark)
summary(cox_landmark)

# Plot Kaplan-Meier curves for the landmark analysis
fit_landmark <- survfit(Surv(time, event) ~ group, data = ipd_landmark)
ggsurvplot(fit_landmark, data = ipd_landmark, risk.table = TRUE, conf.int = TRUE,
           ggtheme = theme_minimal(), legend.title = "Group")
