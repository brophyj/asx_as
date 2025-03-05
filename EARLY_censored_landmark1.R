library(survival)
library(dplyr)
library(readr)

# Function to simulate individual patient data (IPD) based on cumulative incidence and at-risk numbers
simulate_ipd_with_risk <- function(df, at_risk, time_points, N, group_label) {
  df <- df %>% arrange(time)
  ipd_list <- list()
  
  for (i in 1:(length(time_points) - 1)) {
    t_start <- time_points[i]
    t_end <- time_points[i + 1]
    
    ci_start <- approx(x = df$time, y = df$incidence, xout = t_start, rule = 2)$y
    ci_end   <- approx(x = df$time, y = df$incidence, xout = t_end, rule = 2)$y
    
    cum_events_start <- N * ci_start / 100
    cum_events_end   <- N * ci_end / 100
    
    events_interval <- round(cum_events_end - cum_events_start)
    risk_drop <- at_risk[i] - at_risk[i + 1]
    
    censor_interval <- max(risk_drop - events_interval, 0)
    
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

perform_cox_analysis <- function(dataset1, dataset2, group1_label, group2_label, time_points, at_risk1, at_risk2) {
  data1 <- read_csv(file = dataset1, col_names = FALSE)
  data2 <- read_csv(file = dataset2, col_names = FALSE)
  
  colnames(data1) <- c("time", "incidence")
  colnames(data2) <- c("time", "incidence")
  
  ipd1 <- simulate_ipd_with_risk(data1, at_risk1, time_points, at_risk1[1], group1_label)
  ipd2 <- simulate_ipd_with_risk(data2, at_risk2, time_points, at_risk2[1], group2_label)
  
  ipd_all <- bind_rows(ipd1, ipd2)
  ipd_all$group <- factor(ipd_all$group, levels = c(group2_label, group1_label))
  
  cox_model_full <- coxph(Surv(time, event) ~ group, data = ipd_all)
  
  # Censored Model at 1 year
  ipd_censored <- ipd_all
  ipd_censored$time <- pmin(ipd_censored$time, 1)  # 1 year as censored point
  ipd_censored$event[ipd_censored$time == 1] <- 0
  cox_model_censored <- coxph(Surv(time, event) ~ group, data = ipd_censored)
  
  # Landmark Model starting from 1 year
  landmark_data <- ipd_all %>% filter(time > 1)
  landmark_data$time <- landmark_data$time - 1  # Reset time from 1 year
  cox_model_landmark <- if (nrow(landmark_data) > 0 && any(landmark_data$event == 1)) {
    coxph(Surv(time, event) ~ group, data = landmark_data)
  } else {
    NULL  # No valid data for analysis
  }
  
  return(list(
    Full_Model = summary(cox_model_full),
    Censored_Model = if (!is.null(cox_model_censored)) summary(cox_model_censored) else "No valid data for censored analysis.",
    Landmark_Model = if (!is.null(cox_model_landmark)) summary(cox_model_landmark) else "No valid data for landmark analysis."
  ))
}

# Example usage
result <- perform_cox_analysis(
  dataset1 = "references/EARLY_TAVR.csv",
  dataset2 = "references/EARLY_Ctl.csv",
  group1_label = "TAVR",
  group2_label = "Clinical Surveillance",
  time_points = c(0, 1, 2, 3, 4, 5),
  at_risk1 = c(455, 390, 363, 285, 141, 103),
  at_risk2 = c(446, 305, 266, 187, 117, 46)
)

# Print the results
print(result$Full_Model)
print(result$Censored_Model)
print(result$Landmark_Model)
