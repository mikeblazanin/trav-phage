#Example growth curve code
set.seed(1)
library("ggplot2")
library("dplyr")

gc_data <- read.csv("Example_dat.csv", stringsAsFactors = F)

#Plot data
plot(gc_data$Time_s, gc_data$OD600)
plot(gc_data$Time_s, gc_data$cfu_ml)

#Smooth data
loess_out <- loess(cfu_ml ~ Time_s,
                   data = gc_data,
                   span = 0.5, #this is the "smoothness" tuning parameter
                   degree = 2)

gc_data$cfu_ml_smoothed <- loess_out$fitted

#Plot data with smoothed
plot(gc_data$Time_s, gc_data$cfu_ml)
lines(gc_data$Time_s, gc_data$cfu_ml_smoothed)

#Fit data with Baranyi growth equations

#Define the function that will return a list of predicted densities
# when provided the parameters of the Baranyi growth equation and
# a vector of timepoints
# (Note that parameters in time units are in /hr for this equation)
baranyi_func <- function(r, k, v, q0, m, d0, t_vals) {
  #Copied from Ram et al 2019
  if (anyNA(c(r, k, v, q0, m, d0, t_vals))) {return(NA)}
  if (q0 < 0) {q0 <- 0}
  t_vals_hrs <- t_vals/3600
  a <- t_vals_hrs + 1/m*log((exp(-m*t_vals_hrs)+q0)/(1+q0))
  d <- k/(1-(1-((k/d0)**v))*exp(-r*v*a))**(1/v)
  return(d)
}

#Example use of baranyi func to predict growth for given set of params
gc_data$ex_pred_cfu <- baranyi_func(r = 1,
                            k = 10**9,
                            v = 1,
                            q0 = 1,
                            m = 1,
                            d0 = 10**8,
                            t_vals = gc_data$Time_s)
plot(gc_data$Time_s, gc_data$ex_pred_cfu)

#Comparing example predicted values to actual data
plot(gc_data$Time_s, gc_data$cfu_ml_smoothed)
lines(gc_data$Time_s, gc_data$ex_pred_cfu)

#Now we'll fit the function to the data
# (because my data has a diauxic shift, and the Baranyi function
# doesn't account for a diauxic shift, we'll just use the data points
# before 35000 seconds)
gc_data <- gc_data[gc_data$Time_s < 35000, ]
plot(gc_data$Time_s, gc_data$cfu_ml_smoothed)
lines(gc_data$Time_s, gc_data$ex_pred_cfu)

#Define the function to hand to optim
# this is the one that simply returns the squared error between
#  the observed densities and the predicted densities for a given set
#  of parameters
baranyi_fit_err <- function(params, t_vals, dens_vals) {
  #params <- c("logk" = ..., "logd0" = ..., "r" = ..., "v" = ...,
  #            "m" = ..., "q0" = ...)
  pred_vals <- baranyi_func(r = params["r"],
                            k = 10**params["logk"],
                            v = params["v"],
                            q0 = params["q0"],
                            m = params["m"],
                            d0 = 10**params["logd0"],
                            t_vals = t_vals)
  pred_vals[pred_vals < 0] <- 0
  err <- sum((log10(pred_vals) - log10(dens_vals))**2)
  if (is.infinite(err) | is.na(err)) {return(2*10**300)} else {return(err)}
}

#Now run optim to find the fitted parameters
optim_results <- optim(
  #Provide a vector of initial guesses for the parameters
  par = c("logk" = log10(10**9),
          "logd0" = log10(10**8),
          "r" = 1,
          "v" = 1,
          "m" = 1,
          "q0" = 1),
  fn = baranyi_fit_err,
  dens_vals = gc_data$cfu_ml_smoothed,
  t_vals = gc_data$Time_s,
  method = "BFGS")
  
#Take a look at the resulting optimized parameters
optim_results$par

#Calculate predicted density using fitted values
gc_data$ex_pred_cfu <- baranyi_func(r = optim_results$par["r"],
                            k = 10**optim_results$par["logk"],
                            v = optim_results$par["v"],
                            q0 = optim_results$par["q0"],
                            m = optim_results$par["m"],
                            d0 = 10**optim_results$par["logd0"],
                            t_vals = gc_data$Time_s)

plot(gc_data$Time_s, gc_data$cfu_ml_smoothed)
lines(gc_data$Time_s, gc_data$ex_pred_cfu)