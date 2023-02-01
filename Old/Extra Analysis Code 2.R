##Isolate growth curves: smooth ----
gc_data$sm_loess_25k <- NA
gc_data$sm_loess_3600 <- NA
gc_data$sm_movmed_3 <- NA
gc_data$sm_movmed_3_loess_3600 <- NA
for (my_well in unique(gc_data$uniq_well)) {
  #(leaving out the first half hour of data)
  my_rows <- which(gc_data$uniq_well == my_well &
                     gc_data$Time_s > 1800)
  #Based on later checking, wells 13 and 201 need slight curation
  if (my_well == "2017-A_7x_Anc_Anc_Anc_1_Orig") {
    my_rows <- which(gc_data$uniq_well == my_well &
                       gc_data$Time_s > 10000)
  }
  if (my_well == "2017-E_7x_C_L_E_1_Orig") {
    my_rows <- which(gc_data$uniq_well == my_well &
                       gc_data$Time_s > 4000)
  }
  
  #Calculate the median timestep (timesteps actually vary slightly in
  # the number of seconds they were recorded as differing by)
  med_timestep <- median(gc_data$Time_s[my_rows[2]:my_rows[length(my_rows)]]-
                           gc_data$Time_s[my_rows[1]:my_rows[(length(my_rows)-1)]])
  
  #Smooth with loess (based on window of ~7 hrs, or span of ~0.4)
  gc_data$sm_loess_25k[my_rows] <- 
    loess(cfu_ml ~ Time_s, 
          data = gc_data[my_rows, ],
          span = ((25000/med_timestep)+1)/length(my_rows),
          degree = 2)$fitted
  #(based on window of 60 mins)
  gc_data$sm_loess_3600[my_rows] <- 
    loess(cfu_ml ~ Time_s, 
          data = gc_data[my_rows, ],
          span = ((3600/med_timestep)+1)/length(my_rows),
          degree = 1)$fitted
  gc_data$sm_movmed_3[my_rows] <- 
    gcplyr::moving_median(cfu_ml ~ Time_s,
                          data = gc_data[my_rows, ],
                          window_width_n = 3)
  gc_data$sm_movmed_3_loess_3600[my_rows] <- 
    c(NA, loess(sm_movmed_3 ~ Time_s, 
                data = gc_data[my_rows, ],
                span = ((3600/med_timestep)+1)/length(my_rows),
                degree = 1)$fitted,
      NA)
}

#Define function that calculates derivatives
calc_deriv <- function(density, percapita = FALSE,
                       subset_by = NULL, time = NULL,
                       time_normalize = NULL) {
  #Note! density values must be sorted sequentially into their unique sets already
  
  #Provided a vector of density values, this function returns (by default) the
  # difference between sequential values
  #if percapita = TRUE, the differences of density are divided by density
  #if subset_by is provided, it should be a vector (same length as density),
  # the unique values of which will separate calculations
  #if time_normalize is specified, time should be provided as a simple 
  # numeric (e.g. number of seconds) in some unit
  #Then the difference will be normalized for the time_normalize value
  #(e.g. if time is provided in seconds and the difference per hour is wanted,
  # time_normalize should = 3600)
  
  #Check inputs
  if (!is.null(time) & !is.numeric(time)) {
    stop("time is not numeric")
  }
  if (!is.null(time_normalize)) {
    if (!is.numeric(time_normalize)) {
      stop("time_normalize is not numeric")
    } else if (is.null(time)) {
      stop("time_normalize is specified, but time is not provided")
    }
  }
  
  #Calc derivative
  ans <- c(density[2:length(density)]-density[1:(length(density)-1)])
  #Percapita (if specified)
  if (percapita) {
    ans <- ans/density[1:(length(density)-1)]
  }
  #Time normalize (if specified)
  if (!is.null(time_normalize)) {
    ans <- ans/
      (c(time[2:length(time)]-time[1:(length(time)-1)])/time_normalize)
  }
  #Subset by (if specified)
  if (!is.null(subset_by)) {
    ans[subset_by[2:length(subset_by)] != subset_by[1:(length(subset_by)-1)]] <- NA
  }
  return(c(ans, NA))
}

#Calculate growth per hour from loess curve
gc_data$deriv_sm_loess_3600 <- calc_deriv(gc_data$sm_loess_3600,
                                          subset_by = gc_data$uniq_well,
                                          time = gc_data$Time_s,
                                          time_normalize = 3600)
gc_data$deriv_sm_loess_25k <- calc_deriv(gc_data$sm_loess_25k,
                                         subset_by = gc_data$uniq_well,
                                         time = gc_data$Time_s,
                                         time_normalize = 3600)
gc_data$deriv_sm_movmed_3_loess_3600 <- calc_deriv(gc_data$sm_movmed_3_loess_3600,
                                                   subset_by = gc_data$uniq_well,
                                                   time = gc_data$Time_s,
                                                   time_normalize = 3600)

#Calculate per capita growth per hour from loess curve
gc_data$percap_deriv_sm_loess_3600 <- calc_deriv(gc_data$sm_loess_3600,
                                                 percapita = TRUE,
                                                 subset_by = gc_data$uniq_well,
                                                 time = gc_data$Time_s,
                                                 time_normalize = 3600)
gc_data$percap_deriv_sm_loess_25k <- calc_deriv(gc_data$sm_loess_25k,
                                                percapita = TRUE,
                                                subset_by = gc_data$uniq_well,
                                                time = gc_data$Time_s,
                                                time_normalize = 3600)
gc_data$percap_deriv_sm_movmed_3_loess_3600 <- calc_deriv(gc_data$sm_movmed_3_loess_3600,
                                                          percapita = TRUE,
                                                          subset_by = gc_data$uniq_well,
                                                          time = gc_data$Time_s,
                                                          time_normalize = 3600)

#View samples of original & smoothed curves
# as well as derivatives (per cap & not) of both orig and smoothed curves
if (make_curveplots) {
  for (my_well in sample(unique(gc_data$uniq_well), 1)) {
    my_rows <- which(gc_data$uniq_well == my_well)
    
    print(cowplot::plot_grid(
      ggplot(data = gc_data[my_rows, ],
             aes(x = Time_s, y = cfu_ml)) +
        geom_line(color = "red", lwd = 1, alpha = 0.5) +
        geom_line(aes(x = Time_s, y = sm_loess_3600),
                  color = "blue", lwd = 1, alpha = 0.5) +
        geom_line(aes(x = Time_s, y = sm_loess_25k),
                  color = "green", lwd = 0.5, alpha = 0.5) +
        ggtitle(gc_data[my_rows[1], "uniq_well"]) +
        #geom_hline(yintercept = 383404890) +
        scale_y_continuous(trans = "log10") +
        NULL,
      ggplot(data = gc_data[my_rows, ],
             aes(x = Time_s, y = deriv_sm_loess_25k)) +
        geom_line(color = "blue") +
        NULL,
      ggplot(data = gc_data[my_rows, ],
             aes(x = Time_s, y = percap_deriv_sm_loess_25k)) +
        geom_line(color = "blue") +
        # geom_line(aes(x = Time_s, y = percap_deriv_cfu),
        #           color = "red") +
        NULL,
      ncol = 1, align = "v"))
  }
}

#Group data by unique wells
gc_data <- group_by(gc_data, Date, Proj, Pop, Treat, Isol, Rep_Well, Media,
                    uniq_well, uniq_well_num)

#For later fitting of functions to the data, we just need to identify:
# at first min
#   density
#   time
# at max percap growth rate
#   rate (r)
#   density
#   time
# at first local minima of non-percap growth rate after max growth rate
#  (density of diauxic shift, aka pseudo carrying capacity)
#   density
#   time

find_local_extrema <- gcplyr::find_local_extrema

#Define the window widths (sensitivity) for local extrema search
first_min_time_window <- 7200
max_percap_time_window <- 21600
max_gr_rate_time_window <- 10800
pseudo_K_time_window <- 7200

#Summarize data (note: message is for grouping of **output**)
gc_summarized <- dplyr::summarize(gc_data,
                                  #Find the first minima in total density
                                  first_min_index = (find_local_extrema(sm_loess_3600,
                                                                        x = Time_s,
                                                                        width_limit = first_min_time_window,
                                                                        return_maxima = FALSE,
                                                                        return_endpoints = TRUE,
                                                                        na.rm = TRUE)[1]),
                                  first_min = sm_loess_3600[first_min_index],
                                  first_min_time = Time_s[first_min_index],
                                  #find max per capita growth rate
                                  max_percap_index = 
                                    #min(which(percap_deriv_sm_movmed_3_loess_3600 >= 0.5)),
                                    ifelse(any(percap_deriv_sm_movmed_3_loess_3600 >= 0.5, na.rm = TRUE),
                                           min(which(percap_deriv_sm_movmed_3_loess_3600 >= 0.5)),
                                           min(which(percap_deriv_sm_movmed_3_loess_3600 >= 0.4))),
                                  max_percap_gr_rate = percap_deriv_sm_movmed_3_loess_3600[max_percap_index],
                                  max_percap_gr_time = Time_s[max_percap_index],
                                  max_percap_gr_dens = sm_movmed_3_loess_3600[max_percap_index],
                                  true_max_percap = max(percap_deriv_sm_movmed_3_loess_3600, na.rm = TRUE),
                                  
                                  #find the local minimas in total grow rate (slope of total density)
                                  #(which is the point when the diauxic shift occurs)
                                  pseudo_K_index = (find_local_extrema(deriv_sm_loess_25k,
                                                                       x = Time_s,
                                                                       width_limit = pseudo_K_time_window,
                                                                       return_maxima = FALSE,
                                                                       return_endpoints = FALSE,
                                                                       na.rm = T)[1]),
                                  pseudo_K = sm_loess_25k[pseudo_K_index],
                                  pseudo_K_time = Time_s[pseudo_K_index],
                                  pseudo_K_deriv = deriv_sm_loess_25k[pseudo_K_index]
)

# #Fix run 185
# temp <- gc_data[gc_data$uniq_well == "2017-E_7x_B_C_E_1_Orig", ]
# temp <- group_by(temp, Date, Proj, Pop, Treat, Isol, Rep_Well, Media,
#                     uniq_well, uniq_well_num)
# gc_summarized[which(gc_summarized$uniq_well == "2017-E_7x_B_C_E_1_Orig"), ] <- 
#   dplyr::summarize(
#     temp,
#     #we know all the added nas are at start, derivs have extra na's at the end
#     # but those don't exist in the sm_loess itself
#     num_nas = sum(is.na(sm_loess_3600)),
#     #Find the first minima in total density
#     first_min_index = (find_local_extrema(sm_loess_3600,
#                                           return_maxima = FALSE,
#                                           width_limit = (first_min_time_window/
#                                                            (Time_s[2]-Time_s[1])) + 1,
#                                           na.rm = T,
#                                           remove_endpoints = FALSE)[1]+num_nas),
#     first_min = sm_loess_3600[first_min_index],
#     first_min_time = Time_s[first_min_index],
#     #find peaks in per capita growth rate
#     max_percap_index =
#       #first find all peaks
#       (find_local_extrema(percap_deriv_sm_loess_25k,
#                           return_minima = FALSE,
#                           width_limit = (max_percap_time_window/
#                                            (Time_s[2]-Time_s[1])) + 1,
#                           na.rm = T,
#                           remove_endpoints = F)[
#                             #But save/use the first one that follows the minimum density
#                             match(TRUE, find_local_extrema(percap_deriv_sm_loess_25k,
#                                                            return_minima = FALSE,
#                                                            width_limit = (max_percap_time_window/
#                                                                             (Time_s[2]-Time_s[1])) + 1,
#                                                            na.rm = T,
#                                                            remove_endpoints = F) >= first_min_index)]+num_nas),
#     max_percap_gr_rate = percap_deriv_sm_loess_25k[max_percap_index],
#     max_percap_gr_time = Time_s[max_percap_index],
#     max_percap_gr_dens = sm_loess_25k[max_percap_index],
#     max_percap_gr_timesincemin = max_percap_gr_time - first_min_time,
#     
#     #find the local minimas in total grow rate (slope of total density)
#     #(which is the point when the diauxic shift occurs)
#     pseudo_K_index = (find_local_extrema(deriv_sm_loess_25k,
#                                          return_maxima = FALSE,
#                                          width_limit = (pseudo_K_time_window/
#                                                           (Time_s[2]-Time_s[1])) + 1,
#                                          na.rm = T,
#                                          remove_endpoints = T)[2]+num_nas),
#     pseudo_K = sm_loess_25k[pseudo_K_index],
#     pseudo_K_time = Time_s[pseudo_K_index],
#     pseudo_K_deriv = deriv_sm_loess_25k[pseudo_K_index]
#     )

#Change to data frame for cleanliness
gc_summarized <- as.data.frame(gc_summarized)

#Make output plots for problematic wells
if (make_curveplots) {
  #New run w/ better first min:
  #13 - first min too early, change to NA all points before 10000s
  #29 - no k, drop it (other rep is good)
  #201 - first min too early, change to NA all points before 4000s
  #484- nok, drop it (other rep is good)
  #496- nok, drop it (other rep is good)
  
  wells_check <- gc_summarized$uniq_well[
    is.na(gc_summarized$max_percap_index) | is.na(gc_summarized$pseudo_K_index)]
  #c("2017-A_7x_Anc_Anc_Anc_1_Orig", "2017-E_7x_C_L_E_1_Orig")
  #This is the old bad wells:
  # wells_check <- c("2017-B_7x_C_L_B_1_Orig",
  #                  "2017-A_7x_Anc_Anc_Anc_1_Orig",
  #                  "2017-A_7x_B_L_A_1_Rich",
  #                  "2017-B_7x_C_L_B_1_Orig",
  #                  "2017-C_7x_B_C_C_1_Rich",
  #                  "2017-C_7x_B_C_C_2_Rich",
  #                  "2017-C_7x_C_C_C_1_Orig",
  #                  "2017-C_7x_C_L_C_2_Rich",
  #                  "2017-C_7x_D_G_C_1_Orig",
  #                  "2017-C_7x_D_G_C_1_Rich",
  #                  "2017-C_7x_D_G_C_2_Orig",
  #                  "2017-C_7x_D_G_C_2_Rich",
  #                  "2017-C_7x_E_G_C_1_Rich",
  #                  "2017-C_7x_E_G_C_2_Rich",
  #                  "2017-E_7x_B_C_E_1_Orig",
  #                  "2017-E_7x_C_C_E_1_Rich",
  #                  "2017-E_7x_C_C_E_2_Rich",
  #                  "2019-09-10_125_B_C_A_1_Orig",
  #                  "2019-09-10_125_B_C_A_2_Orig",
  #                  "2019-09-10_125_B_G_A_1_Orig",
  #                  "2019-09-12_125_B_C_D_1_Orig",
  #                  "2019-09-12_125_B_C_D_2_Orig",
  #                  "2019-09-12_125_B_G_D_1_Orig",
  #                  "2019-09-12_125_D_L_D_1_Rich",
  #                  "2019-09-12_125_D_L_D_2_Rich",
  #                  "2019-09-13_125_B_C_E_1_Rich",
  #                  "2019-09-13_125_B_C_E_2_Orig",
  #                  "2019-09-13_125_C_C_E_1_Rich",
  #                  "2019-09-13_125_D_L_E_1_Rich"
  # )
  for (my_well in wells_check) {
    tiff(filename = paste("./Challenging_growth_curve_plots/", my_well, ".tiff", sep = ""),
         width = 5, height = 10, units = "in", res = 300)
    my_rows <- which(gc_data$uniq_well == my_well)
    print(cowplot::plot_grid(
      ggplot(data = gc_data[my_rows, ],
             aes(x = Time_s, y = cfu_ml)) +
        geom_line(color = "red", lwd = 1, alpha = 0.5) +
        geom_line(aes(y = sm_loess_3600),
                  color = "blue", lwd = 1, alpha = 0.5) +
        geom_line(aes(y = sm_loess_25k),
                  color = "green", lwd = 0.5, alpha = 0.5) +
        geom_line(aes(y = sm_movmed_3_loess_3600),
                  color = "brown", alpha = 0.5) +
        ggtitle(gc_data[my_rows[1], "uniq_well"]) +
        #Add point for first minima
        geom_point(data = gc_summarized[gc_summarized$uniq_well == my_well, ],
                   aes(x = first_min_time, y = first_min),
                   color = "green", size = 3) +
        scale_y_continuous(trans = "log10") +
        NULL,
      ggplot(data = gc_data[my_rows, ],
             aes(x = Time_s, y = deriv_sm_loess_25k)) +
        geom_line(color = "green") +
        geom_line(aes(y = deriv_sm_loess_3600), color = "blue", alpha = 0.5) +
        geom_line(aes(y = deriv_sm_movmed_3_loess_3600), 
                  color = "brown", alpha = 0.5) +
        #Add point for pseudo K
        geom_point(data = gc_summarized[gc_summarized$uniq_well == my_well, ],
                   aes(x = pseudo_K_time, y = pseudo_K_deriv),
                   color = "green", size = 3) +
        NULL,
      ggplot(data = gc_data[my_rows, ],
             aes(x = Time_s, y = percap_deriv_sm_loess_25k)) +
        geom_line(color = "green") +
        geom_line(aes(y = percap_deriv_sm_loess_3600), 
                  color = "blue", alpha = 0.5) +
        geom_line(aes(y = percap_deriv_sm_movmed_3_loess_3600), 
                  color = "brown", alpha = 0.5) +
        #Add point for max growth rate
        geom_point(data = gc_summarized[gc_summarized$uniq_well == my_well, ],
                   aes(x = max_percap_gr_time, y = max_percap_gr_rate),
                   color = "green", size = 3) +
        NULL,
      ncol = 1, align = "v"))
    dev.off()
  }
}

#Make output plots for all wells
dir.create("./Growth_curve_plots/", showWarnings = FALSE)
if (make_curveplots) {
  for (my_well in unique(gc_data$uniq_well)) {
    tiff(filename = 
           paste("./Growth_curve_plots/", 
                 formatC(gc_data$uniq_well_num[gc_data$uniq_well == my_well][1],
                         width = 3, format = "d", flag = "0"),
                 "_", my_well, ".tiff", sep = ""),
         width = 5, height = 10, units = "in", res = 300)
    my_rows <- which(gc_data$uniq_well == my_well)
    print(cowplot::plot_grid(
      ggplot(data = gc_data[my_rows, ],
             aes(x = Time_s, y = cfu_ml)) +
        geom_line(color = "red", lwd = 1, alpha = 0.5) +
        geom_line(aes(y = sm_loess_3600),
                  color = "blue", lwd = 1, alpha = 0.5) +
        geom_line(aes(y = sm_loess_25k),
                  color = "green", lwd = 0.5, alpha = 0.5) +
        geom_line(aes(y = sm_movmed_3_loess_3600),
                  color = "brown", alpha = 0.5) +
        ggtitle(gc_data[my_rows[1], "uniq_well"]) +
        #Add point for first minima
        geom_point(data = gc_summarized[gc_summarized$uniq_well == my_well, ],
                   aes(x = first_min_time, y = first_min),
                   color = "green", size = 3) +
        scale_y_continuous(trans = "log10") +
        NULL,
      ggplot(data = gc_data[my_rows, ],
             aes(x = Time_s, y = deriv_sm_loess_25k)) +
        geom_line(color = "green") +
        geom_line(aes(y = deriv_sm_loess_3600), color = "blue", alpha = 0.5) +
        geom_line(aes(y = deriv_sm_movmed_3_loess_3600), 
                  color = "brown", alpha = 0.5) +
        #Add point for pseudo K
        geom_point(data = gc_summarized[gc_summarized$uniq_well == my_well, ],
                   aes(x = pseudo_K_time, y = pseudo_K_deriv),
                   color = "green", size = 3) +
        NULL,
      ggplot(data = gc_data[my_rows, ],
             aes(x = Time_s, y = percap_deriv_sm_loess_25k)) +
        geom_line(color = "green") +
        geom_line(aes(y = percap_deriv_sm_loess_3600), 
                  color = "blue", alpha = 0.5) +
        geom_line(aes(y = percap_deriv_sm_movmed_3_loess_3600), 
                  color = "brown", alpha = 0.5) +
        #Add point for max growth rate
        geom_point(data = gc_summarized[gc_summarized$uniq_well == my_well, ],
                   aes(x = max_percap_gr_time, y = max_percap_gr_rate),
                   color = "green", size = 3) +
        NULL,
      ncol = 1, align = "v"))
    dev.off()
  }
}

##Isolate growth curves: Fit curves to data ----
#define error for logistic fit (for use with optim())
logis_fit_err <- function(params, t_vals, dens_vals, t_offset) {
  #params <- c("logk" = ..., "d0" = ..., "r" = ..., "delta" = ...)
  t_vals_hrs <- t_vals/3600
  t_offset_hrs <- t_offset/3600
  k <- 10**params["logk"]
  pred_vals <- with(as.list(params),
                    k/(1+(((k-d0)/d0)*exp(-r*(t_vals_hrs-t_offset_hrs)))))
  pred_vals[pred_vals < 0] <- 0
  err <- sum((log10(pred_vals) - log10(dens_vals))**2)
  if (is.infinite(err) | is.na(err)) {return(2*10**300)} else {return(err)}
}

baranyi_func <- function(r, k, v, q0, m, d0, t_vals) {
  #Copied from Ram et al 2019
  if (anyNA(c(r, k, v, q0, m, d0, t_vals))) {return(NA)}
  if (q0 < 0) {q0 <- 0}
  t_vals_hrs <- t_vals/3600
  a <- t_vals_hrs + 1/m*log((exp(-m*t_vals_hrs)+q0)/(1+q0))
  d <- k/(1-(1-((k/d0)**v))*exp(-r*v*a))**(1/v)
  return(d)
}

baranyi_func2 <- function(r, k, v, q0, d0, t_vals) {
  #Copied from Ram et al 2019
  #where m = r
  if (anyNA(c(r, k, v, q0, d0, t_vals))) {return(NA)}
  if (q0 < 0) {q0 <- 0}
  t_vals_hrs <- t_vals/3600
  a <- t_vals_hrs + 1/r*log((exp(-r*t_vals_hrs)+q0)/(1+q0))
  d <- k/(1-(1-((k/d0)**v))*exp(-r*v*a))**(1/v)
  return(d)
}

baranyi_func3 <- function(r, k, v, d0, t_vals) {
  #Modified from Ram et al 2019 where a(t) = 1
  # (reducing to logistic with decel param)
  if (anyNA(c(r, k, v, d0, t_vals))) {return(NA)}
  t_vals_hrs <- t_vals/3600
  d <- k/((1-(1-((k/d0)**v))*exp(-r*v*t_vals_hrs))**(1/v))
  return(d)
}

#Seeing what the Baranyi curve looks like
if (F) {
  temp_t <- seq(from = 0, to = 100000, by = 1)
  temp_d <- baranyi_func(r = 0.5, k = 10**10, v = 1, q0 = 0.05, m = 1.5,
                         d0 = 10**7, t_vals = temp_t)
  ggplot(data = data.frame(time = temp_t, dens = temp_d), aes(x = time, y = dens)) +
    geom_line() + scale_y_continuous(trans = "log10")
}

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

baranyi_fit2_err <- function(params, t_vals, dens_vals) {
  #params <- c("logk" = ..., "logd0" = ..., "r" = ..., "v" = ...,
  #            "q0" = ...)
  #   (m = r for this one)
  pred_vals <- baranyi_func2(r = params["r"],
                             k = 10**params["logk"],
                             v = params["v"],
                             q0 = params["q0"],
                             d0 = 10**params["logd0"],
                             t_vals = t_vals)
  pred_vals[pred_vals < 0] <- 0
  err <- sum((log10(pred_vals) - log10(dens_vals))**2)
  if (is.infinite(err) | is.na(err)) {return(2*10**300)} else {return(err)}
}

baranyi_fit3_err <- function(params, t_vals, dens_vals) {
  #params <- c("logk" = ..., "logd0" = ..., "r" = ..., "v" = ...)
  pred_vals <- baranyi_func3(r = params["r"],
                             k = 10**params["logk"],
                             v = params["v"],
                             d0 = 10**params["logd0"],
                             t_vals = t_vals)
  pred_vals[pred_vals < 0] <- 0
  err <- sum((log10(pred_vals) - log10(dens_vals))**2)
  if (is.infinite(err) | is.na(err)) {return(2*10**300)} else {return(err)}
}

#Do fitting
gc_summarized <- cbind(gc_summarized,
                       data.frame("fit_r" = as.numeric(NA), 
                                  "fit_k" = as.numeric(NA),
                                  "fit_d0" = as.numeric(NA), 
                                  "fit_delta" = as.numeric(NA),
                                  "fit_err" = as.numeric(NA),
                                  "fit2_r" = as.numeric(NA), 
                                  "fit2_k" = as.numeric(NA),
                                  "fit2_v" = as.numeric(NA), 
                                  "fit2_q0" = as.numeric(NA),
                                  "fit2_m" = as.numeric(NA),
                                  "fit2_d0" = as.numeric(NA),
                                  "fit2_err" = as.numeric(NA),
                                  "fit3_r" = as.numeric(NA), 
                                  "fit3_k" = as.numeric(NA),
                                  "fit3_v" = as.numeric(NA), 
                                  "fit3_q0" = as.numeric(NA),
                                  "fit3_d0" = as.numeric(NA),
                                  "fit3_err" = as.numeric(NA),
                                  "fit4_r" = as.numeric(NA), 
                                  "fit4_k" = as.numeric(NA),
                                  "fit4_v" = as.numeric(NA), 
                                  "fit4_d0" = as.numeric(NA),
                                  "fit4_err" = as.numeric(NA)))
for (sum_row in 1:nrow(gc_summarized)) {
  #if (!is.na(gc_summarized$pseudo_K[sum_row])) {
  my_well <- gc_summarized$uniq_well[sum_row]
  if(is.na(gc_summarized$pseudo_K_time[sum_row])) {
    pseudo_K_time <- max(gc_data$Time_s[gc_data$uniq_well == my_well])
  } else {pseudo_K_time <- gc_summarized$pseudo_K_time[sum_row]}
  myrows1 <- which(gc_data$uniq_well == my_well &
                     gc_data$Time_s <= pseudo_K_time &
                     gc_data$Time_s >= gc_summarized$max_percap_gr_time[sum_row])
  myrows2 <- which(gc_data$uniq_well == my_well &
                     gc_data$Time_s <= (pseudo_K_time + 45*60) &
                     gc_data$Time_s >= gc_summarized$first_min_time[sum_row])
  
  #Get initial values
  if(!is.na(gc_summarized$pseudo_K[sum_row])) {
    my_pseudo_K <- gc_summarized$pseudo_K[sum_row]
  } else {my_pseudo_K <- 10**mean(log10(gc_summarized$pseudo_K), na.rm = T)}
  if(!is.na(gc_summarized$first_min[sum_row])) {
    my_d0 <- gc_summarized$first_min[sum_row]
  } else {my_d0 <- 10**mean(log10(gc_summarized$first_min), na.rm = T)}
  
  #If init r estimate is NA or below 0.2 or above 2, 
  # use average of all init estimates
  # (this fixes when the fit was a declining func bc init r estimate was
  # too small, and when fit had r that was too steep)
  #Or if it's one of the wells that's high-ish (generally 
  # max percap rate 1.65-2) and needs to start with a lower fit
  if(is.na(gc_summarized$max_percap_gr_rate[sum_row]) |
     gc_summarized$max_percap_gr_rate[sum_row] < 0.2 |
     gc_summarized$max_percap_gr_rate[sum_row] > 2 |
     gc_summarized$uniq_well[sum_row] %in% 
     c("2017-A_7x_A_C_A_1_Rich", "2017-A_7x_Anc_Anc_Anc_2_Rich", 
       "2017-A_7x_B_L_A_2_Rich", "2017-A_7x_C_C_A_2_Orig", 
       "2019-09-12_125_D_L_D_1_Rich", "2019-09-12_125_D_L_D_2_Rich", 
       "2019-09-13_125_D_L_E_1_Rich")) {
    my_r <- mean(gc_summarized$max_percap_gr_rate, na.rm = T)
  } else {
    my_r <- gc_summarized$max_percap_gr_rate[sum_row]
  }
  
  temp1 <- optim(par = c("logk" = log10(my_pseudo_K),
                         "d0" = my_d0,
                         "r" = my_r,
                         "delta" = 0),
                 fn = logis_fit_err,
                 dens_vals = gc_data$sm_loess_25k[myrows1],
                 t_vals = gc_data$Time_s[myrows1],
                 t_offset = gc_summarized$max_percap_gr_time[sum_row],
                 method = "BFGS")
  temp2 <- optim(par = c("logk" = log10(my_pseudo_K),
                         "logd0" = log10(my_d0),
                         "r" = my_r,
                         "v" = 1, 
                         "m" = my_r,
                         "q0" = 0.05),
                 #We're estimating that q0 ~ f/(1-f) * e^(-mt)
                 # where f is the fraction of cells that are
                 # metabolically active when max percap is achieved
                 # (which I'm just guessing arbitrarily to be 0.9)
                 # (and m is the rate of metabolic activation)
                 # "q0" = 0.5/(1-0.5)*
                 # exp(-gc_summarized$max_percap_gr_rate[sum_row]*
                 #       gc_summarized$max_percap_gr_time[sum_row]/
                 #       3600)),
                 fn = baranyi_fit_err,
                 dens_vals = gc_data$cfu_ml[myrows2],
                 t_vals = gc_data$Time_s[myrows2],
                 method = "L-BFGS-B",
                 #logk, logd0, r, v, m, q0
                 lower = c(5, 4, 0, 0, 0, 0),
                 upper = c(11, 10, 10, 50, 10, 100))
  temp3 <- optim(par = c("logk" = log10(my_pseudo_K),
                         "logd0" = log10(my_d0),
                         "r" = my_r,
                         "v" = 1, 
                         "q0" = 0.05),
                 #We're estimating that q0 ~ f/(1-f) * e^(-mt)
                 # where f is the fraction of cells that are
                 # metabolically active when max percap is achieved
                 # (which I'm just guessing arbitrarily to be 0.9)
                 # (and m is the rate of metabolic activation)
                 # "q0" = 0.5/(1-0.5)*
                 # exp(-gc_summarized$max_percap_gr_rate[sum_row]*
                 #       gc_summarized$max_percap_gr_time[sum_row]/
                 #       3600)),
                 fn = baranyi_fit2_err,
                 dens_vals = gc_data$cfu_ml[myrows2],
                 t_vals = gc_data$Time_s[myrows2],
                 method = "BFGS")
  temp4 <- optim(par = c("logk" = log10(my_pseudo_K),
                         "logd0" = log10(my_d0),
                         "r" = my_r,
                         "v" = 1),
                 fn = baranyi_fit3_err,
                 dens_vals = gc_data$cfu_ml[myrows1],
                 t_vals = gc_data$Time_s[myrows1],
                 method = "L-BFGS-B",
                 #logk, logd0, r, v, m, q0
                 lower = c(5, 4, 0, 0),
                 upper = c(11, 10, 10, 50))
  gc_summarized[sum_row, 
                c("fit_r", "fit_k", "fit_d0", "fit_delta", "fit_err",
                  "fit2_r", "fit2_k", "fit2_v", "fit2_q0",
                  "fit2_m", "fit2_d0", "fit2_err",
                  "fit3_r", "fit3_k", "fit3_v", "fit3_q0",
                  "fit3_d0", "fit3_err",
                  "fit4_r", "fit4_k", "fit4_v", "fit4_d0", "fit4_err")] <-
    data.frame("fit_r" = temp1$par["r"], 
               "fit_k" = 10**temp1$par["logk"],
               "fit_d0" = temp1$par["d0"], 
               "fit_delta" = temp1$par["delta"],
               "fit_err" = temp1$value,
               "fit2_r" = temp2$par["r"], 
               "fit2_k" = 10**temp2$par["logk"],
               "fit2_v" = temp2$par["v"], 
               "fit2_q0" = temp2$par["q0"],
               "fit2_m" = temp2$par["m"],
               "fit2_d0" = 10**temp2$par["logd0"],
               "fit2_err" = temp3$value,
               "fit3_r" = temp3$par["r"], 
               "fit3_k" = 10**temp3$par["logk"],
               "fit3_v" = temp3$par["v"], 
               "fit3_q0" = temp3$par["q0"],
               "fit3_d0" = 10**temp3$par["logd0"],
               "fit3_err" = temp3$value,
               "fit4_r" = temp4$par["r"], 
               "fit4_k" = 10**temp4$par["logk"],
               "fit4_v" = temp4$par["v"], 
               "fit4_d0" = 10**temp4$par["logd0"],
               "fit4_err" = temp4$value)
}

#Make plots of pure logistic fits
dir.create("./Growth_curve_plots_fits", showWarnings = F)
if(make_curveplots) {
  for (sum_row in 1:nrow(gc_summarized)) {
    my_well <- gc_summarized$uniq_well[sum_row]
    t_vals <- gc_data$Time_s[gc_data$uniq_well == my_well]
    t_vals_hrs <- gc_data$Time_s[gc_data$uniq_well == my_well]/3600
    offset_hrs <- gc_summarized$max_percap_gr_time[sum_row]/3600
    pred_vals <- with(as.list(
      gc_summarized[sum_row, c("fit_r", "fit_k", "fit_d0", "fit_delta", "fit_err")]),
      fit_k/(1+(((fit_k-fit_d0)/fit_d0)*
                  exp(-fit_r*(t_vals_hrs-offset_hrs)))))
    pred_vals4 <- baranyi_func3(r = gc_summarized[sum_row, "fit4_r"],
                                k = gc_summarized[sum_row, "fit4_k"],
                                v = gc_summarized[sum_row, "fit4_v"],
                                d0 = gc_summarized[sum_row, "fit4_d0"],
                                t_vals = t_vals)
    
    png(paste("./Growth_curve_plots_fits/", 
              gc_summarized$uniq_well_num[sum_row], "_",
              my_well, ".png", sep = ""),
        width = 4, height = 4, units = "in", res = 300)
    print(ggplot(data = gc_data[gc_data$uniq_well == my_well, ],
                 aes(x = Time_s, y = sm_movmed_3_loess_3600)) +
            geom_line(alpha = 0.5) +
            geom_point(size = 0.5, aes(y = cfu_ml)) +
            # geom_line(data.frame(Time_s = t_vals, pred_dens = pred_vals),
            #           mapping = aes(x = Time_s, y = pred_dens),
            #           color = "red", alpha = 0.5) +
            geom_line(data.frame(Time_s = t_vals, pred_dens = pred_vals4),
                      mapping = aes(x = Time_s, y = pred_dens),
                      color = "blue", alpha = 0.5) +
            scale_y_continuous(trans = "log10") +
            geom_vline(aes(xintercept = gc_summarized$max_percap_gr_time[
              gc_summarized$uniq_well == my_well]), lty = 2) +
            geom_vline(aes(xintercept = gc_summarized$pseudo_K_time[
              gc_summarized$uniq_well == my_well]), lty = 2) +
            NULL)
    dev.off()
  }
}

#Ones that don't work so great w/ 0.5 (or 0.4) cutoff:
# 29 - toss
# 38 - toss
# 89 - toss
# 128 - toss
# 145 - toss
# 181 - toss 
# 182 - toss bc other rep is better
# 185 - toss
# 193 - toss
# 204 - toss
# 261 - toss
# 361 - toss
# 369 - toss
# 398 - toss
# 427 - toss
#In all cases, the other rep well worked well and so these wells
#can just be excluded from future analyses.

rows <- which(gc_summarized$uniq_well_num %in%
                c(29, 38, 89, 128, 145, 181, 182, 185, 193,
                  204, 261, 361, 369, 398, 427))
gc_summarized[rows, grep("fit4_", colnames(gc_summarized))] <- NA

#Make plots of baranyi fits
dir.create("./Growth_curve_plots_fits2", showWarnings = F)
if(make_curveplots) {
  for (sum_row in 1:nrow(gc_summarized)) {
    my_well <- gc_summarized$uniq_well[sum_row]
    t_vals <- gc_data$Time_s[gc_data$uniq_well == my_well]
    pred_vals2 <- baranyi_func(r = gc_summarized[sum_row, "fit2_r"],
                               k = gc_summarized[sum_row, "fit2_k"],
                               v = gc_summarized[sum_row, "fit2_v"],
                               q0 = gc_summarized[sum_row, "fit2_q0"],
                               m = gc_summarized[sum_row, "fit2_m"],
                               d0 = gc_summarized[sum_row, "fit2_d0"],
                               t_vals = t_vals)
    pred_vals3 <- baranyi_func2(r = gc_summarized[sum_row, "fit3_r"],
                                k = gc_summarized[sum_row, "fit3_k"],
                                v = gc_summarized[sum_row, "fit3_v"],
                                q0 = gc_summarized[sum_row, "fit3_q0"],
                                d0 = gc_summarized[sum_row, "fit3_d0"],
                                t_vals = t_vals)
    start_vals2 <- baranyi_func("k" = gc_summarized$pseudo_K[sum_row],
                                "d0" = gc_data$cfu_ml[which(gc_data$uniq_well == my_well)[2]],
                                "r" = gc_summarized$max_percap_gr_rate[sum_row],
                                "v" = 1, 
                                "m" = gc_summarized$max_percap_gr_rate[sum_row],
                                #We're estimating that q0 ~ f/(1-f) * e^(-mt)
                                # where f is the fraction of cells that are
                                # metabolically active when max percap is achieved
                                # (which I'm just guessing arbitrarily to be 0.9)
                                # (and m is the rate of metabolic activation)
                                "q0" = .05,
                                # 0.9/(1-0.9)*
                                #  exp(-gc_summarized$max_percap_gr_rate[sum_row]*
                                #        gc_summarized$max_percap_gr_time[sum_row]),
                                t_vals = t_vals)
    start_vals3 <- baranyi_func2("k" = gc_summarized$pseudo_K[sum_row],
                                 "d0" = gc_data$cfu_ml[which(gc_data$uniq_well == my_well)[2]],
                                 "r" = gc_summarized$max_percap_gr_rate[sum_row],
                                 "v" = 1, 
                                 #"m" = gc_summarized$max_percap_gr_rate[sum_row],
                                 #We're estimating that q0 ~ f/(1-f) * e^(-mt)
                                 # where f is the fraction of cells that are
                                 # metabolically active when max percap is achieved
                                 # (which I'm just guessing arbitrarily to be 0.9)
                                 # (and m is the rate of metabolic activation)
                                 "q0" = .05,
                                 # 0.9/(1-0.9)*
                                 #  exp(-gc_summarized$max_percap_gr_rate[sum_row]*
                                 #        gc_summarized$max_percap_gr_time[sum_row]),
                                 t_vals = t_vals)
    
    png(paste("./Growth_curve_plots_fits2/", 
              formatC(gc_summarized$uniq_well_num[sum_row], width = 3, 
                      format = "d", flag = "0"), "_", my_well, ".png", 
              sep = ""),
        width = 4, height = 4, units = "in", res = 300)
    print(ggplot(data = gc_data[gc_data$uniq_well == my_well, ],
                 aes(x = Time_s, y = cfu_ml)) +
            geom_point(size = 0.5) +
            geom_line(data.frame(Time_s = t_vals, pred_dens = start_vals2),
                      mapping = aes(x = Time_s, y = pred_dens),
                      color = "red", alpha = 0.3, lty = 3) +
            geom_line(data.frame(Time_s = t_vals, pred_dens = pred_vals2),
                      mapping = aes(x = Time_s, y = pred_dens),
                      color = "red", alpha = 0.7) +
            geom_line(data.frame(Time_s = t_vals, pred_dens = start_vals3),
                      mapping = aes(x = Time_s, y = pred_dens),
                      color = "blue", alpha = 0.3, lty = 3) +
            geom_line(data.frame(Time_s = t_vals, pred_dens = pred_vals3),
                      mapping = aes(x = Time_s, y = pred_dens),
                      color = "blue", alpha = 0.7) +
            scale_y_continuous(trans = "log10") +
            geom_vline(aes(xintercept = gc_summarized$first_min_time[
              gc_summarized$uniq_well == my_well]), lty = 2) +
            geom_vline(aes(xintercept = (gc_summarized$pseudo_K_time[
              gc_summarized$uniq_well == my_well]+45*60)), lty = 2) +
            geom_vline(aes(xintercept = (gc_summarized$pseudo_K_time[
              gc_summarized$uniq_well == my_well])), lty = 3, alpha = 0.5) +
            ggtitle(paste("err =",
                          signif(gc_summarized$fit2_err[
                            gc_summarized$uniq_well == my_well], 2))) +
            theme_bw() +
            NULL)
    dev.off()
  }
}

#Red curve (full baranyi function, not where m = r) is better
#bad fits of that curve:
#really bad (declining pred): 13, 61, 78, 89, 97, 98, 99, 109, 110, 111, 112, 167
#                             209, 223, 245, 263, 265, 279, 349, 351, 383, 385,
#                             399
#     FIXED (by forcing init r estimates to be above 0.2)
#other really bad: 185 (pseudo k lims are too close)
#     FIXED (by manually using the 2nd minima for pseudo K)
#high r bad: 14, 26, 38, 40, 46, 180, 186
#     FIXED (by forcing init r estimates to be below 2)
#no fit (as expected bc no pseudo K): 29, 153, 154, 155, 156, 204, 398, 484, 496
#meh: 30
#     FIXED (by forcing init r estimates to be below 2)
#After fixing the above, re-checked and saw problems in these curves:
#     2, 16, 28, 31, 456 (too high r), 458 (too high r), 516 (too high r)
#     FIXED (by forcing init r estimates for these runs only to be lower)

if (F) {
  temp <- c(13, 61, 78, 89, 97, 98, 99, 109, 110, 111, 112, 167,
            209, 223, 245, 263, 265, 279, 349, 351, 383, 385,
            399)
  temp <- c(14, 26, 38, 40, 46, 180, 186)
  temp <- c(153, 154, 155, 156, 204, 398)
  temp <- c(2, 16, 28, 31, 456, 458, 516)
  View(gc_summarized[gc_summarized$uniq_well_num %in% temp, ])
  summary(gc_summarized$max_percap_gr_rate[gc_summarized$uniq_well_num %in% temp])
}

#Check how fit results compare to local extrema results
if (make_statplots) {
  print(ggplot(data = gc_summarized,
               aes(x = max_percap_gr_rate, y = fit2_r)) +
          geom_point() +
          scale_x_continuous(trans = "log10") +
          scale_y_continuous(trans = "log10") +
          geom_abline(slope = 1, intercept = 0, lty = 2))
  print(ggplot(data = gc_summarized,
               aes(x = pseudo_K, y = fit2_k, color = Proj, shape = Media)) +
          geom_point() +
          scale_x_continuous(trans = "log10") +
          scale_y_continuous(trans = "log10") +
          geom_abline(slope = 1, intercept = 0, lty = 2))
  print(ggplot(data = gc_summarized,
               aes(x = fit2_d0, y = first_min)) +
          geom_point() +
          scale_y_continuous(trans = "log10") +
          scale_x_continuous(trans = "log10") +
          geom_abline(slope = 1, intercept = 0, lty = 2))
  hist(gc_summarized$fit_err)
  print(ggplot(data = gc_summarized,
               aes(x = fit3_r, y = fit2_r)) +
          geom_point() +
          scale_x_continuous(trans = "log10") +
          scale_y_continuous(trans = "log10") +
          geom_abline(slope = 1, intercept = 0, lty = 2))
}     

#Calculate lag time
gc_summarized <- as.data.frame(gc_summarized)
gc_summarized$fit2_lagtime_hrs <- NA
for (sum_row in 1:nrow(gc_summarized)) {
  if (all(!is.na(gc_summarized[sum_row,
                               c("fit2_r", "fit2_k", "fit2_v",
                                 "fit2_q0", "fit2_m", "fit2_d0")]))) {
    my_well <- gc_summarized$uniq_well[sum_row]
    t_vals <- seq(from = 0, to = max(gc_data$Time_s[gc_data$uniq_well == my_well]),
                  by = 1)
    pred_vals1 <- baranyi_func(r = gc_summarized[sum_row, "fit2_r"],
                               k = gc_summarized[sum_row, "fit2_k"],
                               v = gc_summarized[sum_row, "fit2_v"],
                               q0 = gc_summarized[sum_row, "fit2_q0"],
                               m = gc_summarized[sum_row, "fit2_m"],
                               d0 = gc_summarized[sum_row, "fit2_d0"],
                               t_vals = t_vals)
    pred_deriv <- calc_deriv(density = pred_vals1)
    exp_grow_index <- which.max(pred_deriv)
    gc_summarized[sum_row, "fit2_lagtime_hrs"] <- 
      (((gc_summarized[sum_row, "fit2_d0"] - pred_vals1[exp_grow_index])/
          pred_deriv[exp_grow_index]) + t_vals[exp_grow_index])/3600
    
    if (F) {
      print(ggplot(data = data.frame(time = t_vals, dens = pred_vals1), 
                   aes(x = t_vals, y = dens)) + geom_line() + 
              #scale_y_continuous(trans = "log10") +
              # geom_vline(xintercept = exp_grow_time) +
              # geom_vline(xintercept = exp_grow_time2, lty = 2) +
              geom_vline(xintercept = 3600*gc_summarized[sum_row, "lagtime_hrs"], 
                         lty = 2) +
              geom_hline(yintercept = gc_summarized[sum_row, "fit2_d0"], lty = 2) +
              geom_abline(slope = pred_deriv[exp_grow_index],
                          intercept = pred_vals1[exp_grow_index] - 
                            pred_deriv[exp_grow_index]*t_vals[exp_grow_index]) +
              NULL)
      print(ggplot(data = data.frame(time = t_vals, 
                                     dNds = calc_deriv(density = pred_vals1)), 
                   aes(x = t_vals, y = dNds)) + geom_line() + 
              scale_y_continuous(trans = "log10") +
              geom_vline(xintercept = exp_grow_time) +
              NULL)
      print(ggplot(data = data.frame(time = t_vals, 
                                     dNdsdN = calc_deriv(density = pred_vals1,
                                                         percapita = T)), 
                   aes(x = t_vals, y = dNdsdN)) + geom_line() + 
              scale_y_continuous(trans = "log10") +
              geom_vline(xintercept = exp_grow_time) +
              NULL)
    }
  }
}

#Make log-transformed v
gc_summarized$fit2_v_log10 <- log10(gc_summarized$fit2_v)

#Isolate growth curves: Make example plots (eg for talks) ----
if (make_curveplots) {
  temp <- gc_data[gc_data$uniq_well == "2017-C_7x_Anc_Anc_Anc_1_Orig", ]
  
  tiff("./Example_curve_plots/gc_plot1.tiff", width = 10, height = 10, units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = cfu_ml)) +
    geom_line(color = "red", lwd = 3) +
    labs(y = "Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30)) 
  dev.off()
  
  tiff("./Example_curve_plots/gc_plot2.tiff", width = 10, height = 10, units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = cfu_ml)) +
    geom_line(color = "red", lwd = 3, alpha = 0.6) +
    geom_line(aes(y = sm_loess_25k), color = "blue", lwd = 3, alpha = 0.6) +
    labs(y = "Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30)) 
  dev.off()
  
  tiff("./Example_curve_plots/gc_plot3.tiff", width = 10, height = 10, units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = deriv_sm_loess_25k)) +
    geom_line(color = "blue", lwd = 3) +
    labs(y = "Change in Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30)) 
  dev.off()
  
  tiff("./Example_curve_plots/gc_plot4.tiff", width = 10, height = 10, units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = percap_deriv_sm_loess_25k)) +
    geom_line(color = "blue", lwd = 3) +
    labs(y = "Per-capita Change in Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30)) 
  dev.off()
  
  #Now add the points
  temp2 <- gc_summarized[gc_summarized$uniq_well == "2017-C_7x_Anc_Anc_Anc_1_Orig", ]
  
  tiff("./Example_curve_plots/gc_plot5.tiff", width = 10, height = 10, units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = sm_loess_25k)) +
    geom_line(color = "blue", lwd = 3) +
    labs(y = "Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30)) +
    geom_point(data = temp2, aes(x = first_min_time/3600, y = first_min),
               color = "black", size = 10)
  dev.off()
  
  tiff("./Example_curve_plots/gc_plot6.tiff", width = 10, height = 10, units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = deriv_sm_loess_25k)) +
    geom_line(color = "blue", lwd = 3) +
    labs(y = "Change in Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30)) +
    geom_point(data = temp2, aes(x = pseudo_K_time/3600, y = pseudo_K_deriv),
               color = "black", size = 10)
  dev.off()
  
  tiff("./Example_curve_plots/gc_plot7.tiff", width = 10, height = 10, units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = sm_loess_25k)) +
    geom_line(color = "blue", lwd = 3) +
    labs(y = "Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30)) +
    geom_point(data = temp2, aes(x = pseudo_K_time/3600, y = pseudo_K),
               color = "black", size = 10)
  dev.off()
  
  tiff("./Example_curve_plots/gc_plot8.tiff", width = 10, height = 10, units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = percap_deriv_sm_loess_25k)) +
    geom_line(color = "blue", lwd = 3) +
    labs(y = "Per-capita Change in Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30)) +
    geom_point(data = temp2, aes(x = max_percap_gr_time/3600, 
                                 y = max_percap_gr_rate),
               color = "black", size = 10)
  dev.off()
  
  tiff("./Example_curve_plots/gc_plot9.tiff", width = 10, height = 10, units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = sm_loess_25k)) +
    geom_line(color = "blue", lwd = 3) +
    labs(y = "Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30)) +
    geom_point(data = temp2, aes(x = pseudo_K_time/3600, y = pseudo_K),
               color = "black", size = 10) +
    geom_point(data = temp2, aes(x = max_percap_gr_time/3600, 
                                 y = max_percap_gr_dens),
               color = "black", size = 10)
  dev.off()
  
  tiff("./Example_curve_plots/gc_plot10.tiff", width = 10, height = 10, units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = sm_loess_25k)) +
    geom_line(color = "blue", lwd = 3) +
    labs(y = "Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30)) +
    geom_point(data = temp2, aes(x = first_min_time/3600, y = first_min),
               color = "black", size = 10)
  dev.off()
  
  tiff("./Example_curve_plots/gc_plot11.tiff", width = 10, height = 10, units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = sm_loess_25k)) +
    geom_line(color = "blue", lwd = 3) +
    labs(y = "Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30)) +
    geom_point(data = temp2, aes(x = pseudo_K_time/3600, y = pseudo_K),
               color = "black", size = 10) +
    geom_point(data = temp2, aes(x = max_percap_gr_time/3600, 
                                 y = max_percap_gr_dens),
               color = "black", size = 7, alpha = 0.7) +
    geom_point(data = temp2, aes(x = first_min_time/3600, y = first_min),
               color = "black", size = 7, alpha = 0.7)
  dev.off()
  
  my_well <- temp$uniq_well[1]
  t_vals <- temp$Time_s
  sum_row <- which(gc_summarized$uniq_well == my_well)
  temp$pred_vals_fit2 <- baranyi_func(r = gc_summarized[sum_row, "fit2_r"],
                                      k = gc_summarized[sum_row, "fit2_k"],
                                      v = gc_summarized[sum_row, "fit2_v"],
                                      q0 = gc_summarized[sum_row, "fit2_q0"],
                                      m = gc_summarized[sum_row, "fit2_m"],
                                      d0 = gc_summarized[sum_row, "fit2_d0"],
                                      t_vals = t_vals)
  
  
  tiff("./Example_curve_plots/gc_plot12.tiff", width = 10, height = 10, 
       units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = cfu_ml)) +
    geom_point(size = 4) +
    labs(y = "Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30)) 
  dev.off()
  
  tiff("./Example_curve_plots/gc_plot13.tiff", width = 10, height = 10, 
       units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = cfu_ml)) +
    geom_point(size = 4) +
    geom_line(aes(y = pred_vals_fit2), color = "red", lwd = 3, alpha = 0.5) +
    labs(y = "Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30)) 
  dev.off()
  
  temp <- gc_data[gc_data$uniq_well == "2017-C_7x_Anc_Anc_Anc_1_Rich", ]
  my_well <- temp$uniq_well[1]
  sum_row <- which(gc_summarized$uniq_well == my_well)
  temp <- temp[temp$Time_s <= gc_summarized$pseudo_K_time[sum_row], ]
  t_vals <- temp$Time_s
  temp$pred_vals_fit2 <- baranyi_func(r = gc_summarized[sum_row, "fit2_r"],
                                      k = gc_summarized[sum_row, "fit2_k"],
                                      v = gc_summarized[sum_row, "fit2_v"],
                                      q0 = gc_summarized[sum_row, "fit2_q0"],
                                      m = gc_summarized[sum_row, "fit2_m"],
                                      d0 = gc_summarized[sum_row, "fit2_d0"],
                                      t_vals = t_vals)
  
  tiff("./Example_curve_plots/gc_plot13.tiff", width = 10, height = 10, 
       units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = cfu_ml)) +
    geom_point(size = 4) +
    geom_line(aes(y = pred_vals_fit2), color = "red", lwd = 3, alpha = 0.5) +
    labs(y = "Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30)) 
  dev.off()
  
}

##Explore rep well variation ----
dir.create("./Growth_curve_rep_plots/", showWarnings = FALSE)
if(make_statplots) {
  for (col_i in grep("fit4_", colnames(gc_summarized))) {
    for (proj in unique(gc_summarized$Proj)) {
      for (media in unique(gc_summarized$Media)) {
        temp <- gc_summarized[gc_summarized$Proj == proj &
                                gc_summarized$Media == media, ]
        tiff(paste(
          "./Growth_curve_rep_plots/",
          proj, "_", media,"_", names(gc_summarized)[col_i], ".tiff", 
          sep = ""), 
          width = 4, height = 5, units = "in", res = 150)
        print(
          ggplot(data = temp,
                 aes_string(x = "Date", 
                            y = names(gc_summarized)[col_i])) +
            geom_point(aes(color = Isol)) +
            facet_grid(Pop ~ Treat) +
            ggtitle(paste(proj, media)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
        dev.off()
      }
    }
  }
}

#Flagging isols for variance in r
#7x A C C Orig
#7x A L E Orig
#7x B C C Orig
#7x C L E Orig
#7x D G B Orig
#7x Anc A Rich
#7x Anc B Rich
#7x C C A Rich
#7x D G B Rich
#7x E L B Rich
#125 B C B Orig
#125 B C E Orig
#125 B L A Orig
#125 B L B Orig
#125 B L D Orig
#125 B C E Rich
#125 C C E Rich
#125 D G C Rich
#




#Isolate growth curves: summarize reps into isols, view variable data & distributions ----

#Summarize replicate wells
gc_summarized <- group_by(gc_summarized, Date, Proj, Pop, Treat,
                          Isol, Media)
gc_sum_isols <- summarize_at(
  gc_summarized,
  .funs = c(avg = function(x) {mean(x, na.rm = TRUE)}, 
            sd = function(x) {sd(x, na.rm = TRUE)}),
  .vars = c(
    "first_min",
    "max_percap_gr_rate",
    "max_percap_gr_dens",
    "max_percap_gr_time",
    "pseudo_K",
    "pseudo_K_time",
    "fit4_r",
    "fit4_k",
    "fit4_v",
    "fit4_d0",
    "fit4_err"))
# "fit2_r", "fit2_k", "fit2_v", "fit2_v_log10", 
# "fit2_q0", "fit2_m", "fit2_d0", "fit2_lagtime_hrs",
# "fit3_r", "fit3_k", "fit3_v", 
# "fit3_q0", "fit3_d0"))
gc_sum_isols <- as.data.frame(gc_sum_isols)

for (proj in unique(gc_sum_isols$Proj)) {
  for (media in unique(gc_sum_isols$Media)) {
    temp <- gc_sum_isols[gc_sum_isols$Proj == proj &
                           gc_sum_isols$Media == media, 
                         c(1:6, 24)]
    temp$grpnm <- paste(temp$Date, temp$Proj, temp$Pop, temp$Treat, temp$Isol, temp$Media,
                        sep = "_")
    temp <- temp[order(-temp$fit4_r_sd), ]
    print(temp[1:10, c(7, 8)])
  }
}

#Taking a look at what to cut out:
#for 7x Orig, would like to cut top 7 (above 0.13 sd)
#for 7x Rich would like to cut top 5 (above 0.12 sd)
#for 125 Orig would like to cut top 6 (above 0.16 sd)
#for 125 Rich would like to cut top 5 (above 0.14 sd)
#Will use 0.12 as 7x cutoff, .14 as 125 cutoff

gc_sum_isols[(gc_sum_isols$Proj == "7x" & 
                !is.na(gc_sum_isols$fit4_r_sd) &
                gc_sum_isols$fit4_r_sd > 0.12),
             grep("fit4_", colnames(gc_sum_isols))] <- NA
gc_sum_isols[gc_sum_isols$Proj == "125" & 
               !is.na(gc_sum_isols$fit4_r_sd) &
               gc_sum_isols$fit4_r_sd > 0.14,
             grep("fit4_", colnames(gc_sum_isols))] <- NA


#Take a look at the standard deviations between replicate wells
# Just raw sd vals (w/ red line for mean avg value)
# if (make_statplots) {
#   for (var in c("first_min_", "first_min_time_", 
#                 "max_percap_gr_rate_", "max_percap_gr_time_", 
#                 "max_percap_gr_dens_", 
#                 "max_percap_gr_timesincemin_",
#                 "pseudo_K_", "pseudo_K_time_", 
#                 "pseudo_K_timesincemin_", 
#                 "pseudo_K_timesince_maxpercap_")) {
#     my_sd <- gc_sum_isols[, paste(var, "sd", sep = "")]
#     my_avg <- mean(gc_sum_isols[, paste(var, "avg", sep = "")])
#     hist(my_sd, main = var, 
#          xlim = c(min(my_sd, my_avg, na.rm = T), max(my_sd, my_avg, na.rm = T)))
#     abline(v = my_avg, col = "red", lwd = 2)
#   }
# }
# Sd vals divided by matching avg value (so 1 is the reference)
#if (make_statplots) {
#   for (var in c("first_min_", "first_min_time_", 
#                 "max_percap_gr_rate_", "max_percap_gr_time_", 
#                 "max_percap_gr_dens_", 
#                 "max_percap_gr_timesincemin_",
#                 "pseudo_K_", "pseudo_K_time_", 
#                 "pseudo_K_timesincemin_", 
#                 "pseudo_K_timesince_maxpercap_")) {
#     my_sd <- gc_sum_isols[, paste(var, "sd", sep = "")]
#     my_avg <- gc_sum_isols[, paste(var, "avg", sep = "")]
#     hist(my_sd/my_avg, main = var, xlim = c(0, max(my_sd/my_avg, 1, na.rm = T)))
#     abline(v = 1, col = "red", lwd = 2)
#   }
# }

#Generally, sd's between reps are small relative to the values themselves

dir.create("./Growth_curve_rep_plots_isols/", showWarnings = FALSE)
if(make_statplots) {
  for (col_i in grep("fit_.*_sd", colnames(gc_sum_isols))) {
    temp <- gc_sum_isols
    tiff(paste(
      "./Growth_curve_rep_plots_isols/", names(gc_sum_isols)[col_i], ".tiff", 
      sep = ""), 
      width = 5, height = 5, units = "in", res = 150)
    print(ggplot(data = temp, 
                 aes_string(x = "Pop", y = names(gc_sum_isols)[col_i], 
                            color = "Isol")) +
            geom_point() +
            facet_grid(Proj*Media ~ Treat, scales = "free_y"))
    dev.off()
  }
}


#View all the isols by variable
my_facet_labels <- c("7x" = "Weak Phage", 
                     "125" = "Strong Phage",
                     "C" = "Control", "G" = "Global", "L" = "Local",
                     "Anc" = "Ancstr",
                     "Rich" = "Rich Media", "Orig" = "Original Media")
dir.create("./Growth_curve_variables_plots/", showWarnings = FALSE)
if (make_statplots) {
  my_vars <- c(
    #"first_min_", 
    #"first_min_time_", 
    #"max_percap_gr_rate_", 
    #"max_percap_gr_time_", 
    #"max_percap_gr_dens_", 
    #"max_percap_gr_timesincemin_",
    #"pseudo_K_", 
    #"pseudo_K_time_", 
    #"pseudo_K_timesincemin_",
    #"pseudo_K_timesince_maxpercap_",
    #"fit_r_", "fit_k_", "fit_d0_",
    # "fit2_r_", "fit2_k_", "fit2_v_", "fit2_v_log10_",
    # "fit2_q0_", "fit2_m_", "fit2_d0_", "fit2_lagtime_hrs_"
    #"fit3_r_", "fit3_k_", "fit3_v_", 
    #"fit3_q0_", "fit3_d0_"
    "fit4_r_", "fit4_k_", "fit4_v_", "fit4_d0_"
  )
  for (i in 1:length(my_vars)) {
    var_root <- my_vars[i]
    var_name <- c(
      #"First minimum density (cfu/mL)",
      #"Maximum per-capita growth rate",
      #"Density at maximum per-capita growth rate",
      #"Time until maximum per-capita growth rate",
      #"Density at diauxic shift (cfu/mL)",
      #"Time until diauxic shift (from min)",
      #"Time until diauxic shift (from max percap)",
      #"Fit r", "Fit carrying capacity", "Fit init density",
      #"Maximum Per Capita Growth Rate (r) (/hr)", 
      #"Density at Diauxic Shift (k) (cfu/mL)", 
      #"Deceleration Parameter (v)", 
      #"Deceleration Parameter (log10(v))", 
      #"Fit 2 q0", "Fit 2 m", "fit 2 d0", 
      #"Lag time (hrs)"
      "r", "k", "v", "d0")[i]
    var <- paste(var_root, "avg", sep = "")
    var_sd <- paste(var_root, "sd", sep = "")
    tiff(paste("./Growth_curve_variables_plots/", var, ".tiff", sep = ""),
         width = 5, height = 5, units = "in", res = 300)
    print(ggplot(data = gc_sum_isols,
                 aes(x = Pop, y = get(var), group = Pop, color = Treat)) +
            geom_point(position = position_dodge(0.6)) +
            facet_nested(Proj ~ Media+Treat, scales = "free_y",
                         labeller = labeller(Proj = my_facet_labels,
                                             Treat = my_facet_labels,
                                             Media = my_facet_labels)) +
            scale_x_discrete(limits = c("Anc", LETTERS[1:5])) +
            scale_color_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                               labels = c("Ancestor", "Control", "Local", "Global"),
                               values = my_cols[c(3, 8, 2, 6)]) +
            geom_errorbar(aes(x = Pop, ymin = get(var)-get(var_sd),
                              ymax = get(var)+get(var_sd)),
                          position = position_dodge(0.6),
                          width = 0.2) +
            labs(y = var_name, x = "Population") +
            theme_bw() +
            theme(legend.position = "none",
                  axis.text.x = element_text(size = 7, angle = 90, 
                                             hjust = 1, vjust = 0.5))
    )
    dev.off()
  }
}

#Note however how some curves have super high sds between
# rep wells

#Also Noted a weird clustering in 7x Rich C of low max percap rates
#Looking at the plots, there's nothing I can see that's wrong
# with those curves. It might be a media batch effect

#Let's normalize by same-plate ancestor
ancestors <- gc_sum_isols[gc_sum_isols$Isol == "Anc", ]
for (var in c("first_min_avg",
              "max_percap_gr_rate_avg",
              "max_percap_gr_dens_avg",
              "max_percap_gr_timesincemin_avg",
              "pseudo_K_avg",
              "pseudo_K_timesincemin_avg",
              "pseudo_K_timesince_maxpercap_avg",
              "fit_r_avg", "fit_k_avg", "fit_d0_avg",
              "fit2_r_avg", "fit2_k_avg", "fit2_v_log10_avg", 
              "fit2_q0_avg", "fit2_m_avg", "fit2_d0_avg", "fit2_lagtime_hrs_avg")) {
  new_var <- paste(var, "_rel", sep = "")
  gc_sum_isols[, new_var] <- gc_sum_isols[, var] -
    ancestors[match(paste(gc_sum_isols$Date, gc_sum_isols$Media), 
                    paste(ancestors$Date, ancestors$Media)), var]
}

#Now view the relative variables
my_facet_labels <- c("7x" = "Weak Phage", 
                     "125" = "Strong Phage",
                     "C" = "Control", "G" = "Global", "L" = "Local",
                     "A" = "WT",
                     "Rich" = "Rich Media", "Orig" = "Original Media")
dir.create("./Growth_curve_variables_plots_relative/", showWarnings = F)
if (make_statplots) {
  my_vars <- c("first_min_", 
               "max_percap_gr_rate_", 
               #"max_percap_gr_dens_", 
               #"max_percap_gr_timesincemin_",
               "pseudo_K_", 
               #"pseudo_K_timesincemin_",
               #"pseudo_K_timesince_maxpercap_",
               #"fit_r_", "fit_k_", "fit_d0_",
               "fit2_r_", "fit2_k_", "fit2_v_log10_", 
               "fit2_q0_", "fit2_m_", "fit2_d0_", "fit2_lagtime_hrs_")
  for (i in 1:length(my_vars)) {
    var_root <- my_vars[i]
    var <- paste(var_root, "avg_rel", sep = "")
    var_name <- c("Relative first minimum density",
                  "Relative maximum per-capita growth rate",
                  #"Relative density at maximum per-capita growth rate",
                  #"Relative time until maximum per-capita growth rate",
                  "Relative density at diauxic shift",
                  #"Relative time until diauxic shift (from min)",
                  #"Relative time until diauxic shift (from max percap)",
                  #"Relative Fit r", "Relative Fit carrying capacity", 
                  #"Relative Fit init density",
                  "Relative Fit 2 r", "Relative Fit 2 k", "Relative Fit 2 v", 
                  "Relative Fit 2 q0", "Relative Fit 2 m", "Relative fit 2 d0",
                  "Relative Fit 2 lag time")[i]
    #Note: if you want to view the sd's between wells of
    # Ancestor-normalized values, you'll have to go back to
    # gc_summarized and calculate the relative values there
    # then re-calculate sd. Have not implemented this
    tiff(paste("./Growth_curve_variables_plots_relative/", var, ".tiff", sep = ""),
         width = 5, height = 5, units = "in", res = 300)
    print(ggplot(data = gc_sum_isols[gc_sum_isols$Pop != "Anc", ],
                 aes(x = Pop, y = get(var), group = Pop,
                     color = Treat)) +
            geom_point(position = position_dodge(0.6)) +
            facet_nested(Proj ~ Media+Treat, scales = "free_y",
                         labeller = labeller(Proj = my_facet_labels,
                                             Treat = my_facet_labels,
                                             Media = my_facet_labels)) +
            scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                               labels = c("Control", "Local", "Global"),
                               values = my_cols[c(8, 2, 6)]) +
            geom_hline(yintercept = 0, lty = 2) +
            labs(y = var_name, x = "Population") +
            theme_bw() +
            theme(legend.position = "none") +
            # geom_errorbar(aes(x = Treat, ymin = get(var)-get(var_sd),
            #                   ymax = get(var)+get(var_sd)),
            #               position = position_dodge(0.6),
            #               width = 0.2)
            NULL
    )
    dev.off()
  }
}

#After looking at those, it definitely improves some of the
# points to normalize by ancestor
# and it doesn't make any others do anything weird
#So we should move forward only with relative variables

gc_sum_isols$PPT <- 
  paste(gc_sum_isols$Proj,gc_sum_isols$Pop, gc_sum_isols$Treat, sep = "_")
gc_sum_isols$PPTI <- 
  paste(gc_sum_isols$Proj, gc_sum_isols$Pop, gc_sum_isols$Treat, 
        gc_sum_isols$Isol, sep = "_")

library(lme4)

gc_isols_merged <- 
  dplyr::left_join(gc_sum_isols, 
                   resis_data_isols[, c("Proj", "Pop", "Treat", 
                                        "Isol", "EOP_avg", "EOP_bd")])
gc_isols_merged$resis_cat <-
  ifelse(gc_isols_merged$EOP_bd, "Resis",
         ifelse(gc_isols_merged$EOP_avg > 0.1,
                "Sens", "Part Resis"))
gc_isols_merged$resis_cat <- factor(gc_isols_merged$resis_cat,
                                    levels = c("Sens", "Part Resis", "Resis"))

mixed_model0_7x_O = lmer(fit4_r_avg ~ (1|PPT) + (1|Date),
                         data = gc_isols_merged[gc_isols_merged$Proj == "7x" &
                                                  gc_isols_merged$Media == "Orig",])
mixed_model0_7x_R = lmer(fit4_r_avg ~ (1|PPT) + (1|Date),
                         data = gc_isols_merged[gc_isols_merged$Proj == "7x" &
                                                  gc_isols_merged$Media == "Rich",])
mixed_model0_125_O = lmer(fit4_r_avg ~ (1|PPT) + (1|Date),
                          data = gc_isols_merged[gc_isols_merged$Proj == "125" &
                                                   gc_isols_merged$Media == "Orig",])
mixed_model0_125_R = lmer(fit4_r_avg ~ (1|PPT) + (1|Date),
                          data = gc_isols_merged[gc_isols_merged$Proj == "125" &
                                                   gc_isols_merged$Media == "Rich",])

mixed_model1_7x_O = lmer(fit4_r_avg ~ Treat + (1|PPT) + (1|Date),
                         data = gc_isols_merged[gc_isols_merged$Proj == "7x" &
                                                  gc_isols_merged$Media == "Orig",])
mixed_model1_7x_R = lmer(fit4_r_avg ~ Treat + (1|PPT) + (1|Date),
                         data = gc_isols_merged[gc_isols_merged$Proj == "7x" &
                                                  gc_isols_merged$Media == "Rich",])
mixed_model1_125_O = lmer(fit4_r_avg ~ Treat + (1|PPT) + (1|Date),
                          data = gc_isols_merged[gc_isols_merged$Proj == "125" &
                                                   gc_isols_merged$Media == "Orig",])
mixed_model1_125_R = lmer(fit4_r_avg ~ Treat + (1|PPT) + (1|Date),
                          data = gc_isols_merged[gc_isols_merged$Proj == "125" &
                                                   gc_isols_merged$Media == "Rich",])

mixed_model2_7x_O = lmer(fit4_r_avg ~ resis_cat + (1|PPT) + (1|Date),
                         data = gc_isols_merged[gc_isols_merged$Proj == "7x" &
                                                  gc_isols_merged$Media == "Orig",])
mixed_model2_7x_R = lmer(fit4_r_avg ~ resis_cat + (1|PPT) + (1|Date),
                         data = gc_isols_merged[gc_isols_merged$Proj == "7x" &
                                                  gc_isols_merged$Media == "Rich",])
mixed_model2_125_O = lmer(fit4_r_avg ~ resis_cat + (1|PPT) + (1|Date),
                          data = gc_isols_merged[gc_isols_merged$Proj == "125" &
                                                   gc_isols_merged$Media == "Orig",])
mixed_model2_125_R = lmer(fit4_r_avg ~ resis_cat + (1|PPT) + (1|Date),
                          data = gc_isols_merged[gc_isols_merged$Proj == "125" &
                                                   gc_isols_merged$Media == "Rich",])

mixed_model3_7x_O = lmer(fit4_r_avg ~ resis_cat + Treat + (1|PPT) + (1|Date),
                         data = gc_isols_merged[gc_isols_merged$Proj == "7x" &
                                                  gc_isols_merged$Media == "Orig",])
mixed_model3_7x_R = lmer(fit4_r_avg ~ resis_cat + Treat + (1|PPT) + (1|Date),
                         data = gc_isols_merged[gc_isols_merged$Proj == "7x" &
                                                  gc_isols_merged$Media == "Rich",])
mixed_model3_125_O = lmer(fit4_r_avg ~ resis_cat + Treat + (1|PPT) + (1|Date),
                          data = gc_isols_merged[gc_isols_merged$Proj == "125" &
                                                   gc_isols_merged$Media == "Orig",])
mixed_model3_125_R = lmer(fit4_r_avg ~ resis_cat + Treat + (1|PPT) + (1|Date),
                          data = gc_isols_merged[gc_isols_merged$Proj == "125" &
                                                   gc_isols_merged$Media == "Rich",])


library(pbkrtest)

KRmodcomp(mixed_model1_7x_O, mixed_model0_7x_O)
KRmodcomp(mixed_model2_7x_O, mixed_model0_7x_O)
KRmodcomp(mixed_model3_7x_O, mixed_model0_7x_O)

KRmodcomp(mixed_model1_7x_R, mixed_model0_7x_R)
KRmodcomp(mixed_model2_7x_R, mixed_model0_7x_R)
KRmodcomp(mixed_model3_7x_R, mixed_model0_7x_R)

KRmodcomp(mixed_model1_125_O, mixed_model0_125_O)
KRmodcomp(mixed_model2_125_O, mixed_model0_125_O)
KRmodcomp(mixed_model3_125_O, mixed_model0_125_O)

KRmodcomp(mixed_model1_125_R, mixed_model0_125_R)
KRmodcomp(mixed_model2_125_R, mixed_model0_125_R)
KRmodcomp(mixed_model3_125_R, mixed_model0_125_R)

##Isolate growth curves: summarize isols into pops, add resis & migr data ----

#Summarize isols into pops
gc_sum_isols <- group_by(gc_sum_isols,
                         Proj, Pop, Treat, Media)
gc_sum_pops <- summarize_at(gc_sum_isols,
                            .funs = c(avg = mean, sd = sd,
                                      med = median),
                            na.rm = TRUE,
                            .vars = c(
                              "first_min_avg_rel",
                              "max_percap_gr_rate_avg_rel",
                              "max_percap_gr_dens_avg_rel",
                              "max_percap_gr_timesincemin_avg_rel",
                              "pseudo_K_avg_rel",
                              "pseudo_K_timesincemin_avg_rel",
                              "pseudo_K_timesince_maxpercap_avg_rel",
                              "fit_r_avg_rel", "fit_k_avg_rel", "fit_d0_avg_rel",
                              "fit2_r_avg_rel", "fit2_k_avg_rel", 
                              "fit2_v_log10_avg_rel",
                              "fit2_q0_avg_rel", "fit2_m_avg_rel", "fit2_d0_avg_rel",
                              "fit2_lagtime_hrs_avg_rel",
                              "fit2_r_avg", "fit2_k_avg", "fit2_v_avg",
                              "fit2_q0_avg", "fit2_m_avg", "fit2_d0_avg"))
gc_sum_pops <- as.data.frame(gc_sum_pops)

#View population-summarized mean data (median below)
my_facet_labels <- c("7x" = "Weak Phage", 
                     "125" = "Strong Phage",
                     "C" = "Control", "G" = "Global", "L" = "Local",
                     "A" = "WT",
                     "Rich" = "Rich Media", "Orig" = "Original Media")
if (make_statplots) {
  var_roots <- c("first_min_avg_rel", 
                 # "first_min_time_avg", 
                 "max_percap_gr_rate_avg_rel", 
                 # "max_percap_gr_time_avg", 
                 "max_percap_gr_dens_avg_rel", 
                 "max_percap_gr_timesincemin_avg_rel",
                 "pseudo_K_avg_rel", 
                 # "pseudo_K_time_avg", 
                 "pseudo_K_timesincemin_avg_rel",
                 "pseudo_K_timesince_maxpercap_avg_rel",
                 "fit2_r_avg_rel", "fit2_k_avg_rel", 
                 "fit2_v_log10_avg_rel", "fit2_q0_avg_rel", 
                 "fit2_m_avg_rel", "fit2_d0_avg_rel",
                 "fit2_lagtime_hrs_avg_rel"
  )
  dir.create("./Growth_curve_variables_plots_pops/", showWarnings = F)
  for (var_root in var_roots) {
    var <- paste(var_root, "_avg", sep = "")
    var_sd <- paste(var_root, "_sd", sep = "")
    tiff(paste("./Growth_curve_variables_plots_pops/", var, ".tiff", sep = ""),
         width = 10, height = 10, units = "in", res = 300)
    print(
      ggplot(data = gc_sum_pops[gc_sum_pops$Pop != "Anc", ],
             aes(x = Treat, y = get(var), group = Pop,
                 color = Treat)) +
        geom_point(size = 5, alpha = 0.65) +
        scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                           labels = c("Control", "Local", "Global"),
                           values = my_cols[c(8, 2, 6)]) +
        facet_grid(Proj ~ Media, scales = "free_y",
                   labeller = labeller(Proj = my_facet_labels,
                                       Media = my_facet_labels)) +
        geom_hline(yintercept = 0, lty = 2) +
        ggtitle(var) +
        theme_bw() +
        # geom_errorbar(aes(x = Treat, ymin = get(var)-get(var_sd),
        #                   ymax = get(var)+get(var_sd)),
        #               position = position_dodge(0.3),
        #               width = 0.2)
        NULL)
    dev.off()
  }
}

#Make combined r-k-lagtime-v plot
if(make_statplots) {
  my_facet_labels <- c("7x" = "Weak Phage", 
                       "125" = "Strong Phage",
                       "C" = "Control", "G" = "Global", "L" = "Local",
                       "A" = "WT",
                       "Rich" = "Rich", "Orig" = "Orig")
  
  temp1 <- gc_sum_pops[gc_sum_pops$Pop != "Anc", ]
  temp1$Title <- "r"
  temp1$Media_title <- "Media"
  rplot <-
    ggplot(data = temp1,
           aes(x = Treat, y = fit2_r_avg_rel_avg, group = Pop,
               color = Treat)) +
    geom_point(size = 5, alpha = 0.65) +
    scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                       labels = c("Control", "Local", "Global"),
                       values = my_cols[c(8, 2, 6)]) +
    scale_x_discrete(breaks = c("C", "L", "G"),
                     labels = c("Control", "Local", "Global")) +
    facet_nested(Proj ~ Media_title*Media, scales = "free_y",
                 labeller = labeller(Proj = my_facet_labels,
                                     Media = my_facet_labels)) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(y = "Evolved Change in Growth Rate (r) (/hr)", x = "Treatment") +
    theme_bw() +
    NULL
  #print(rplot)
  
  temp2 <- gc_sum_pops[gc_sum_pops$Pop != "Anc", ]
  temp2$Title <- "k"
  temp2$Media_title <- "Media"
  kplot <-
    ggplot(data = temp2,
           aes(x = Treat, y = fit2_k_avg_rel_avg/10**8, group = Pop,
               color = Treat)) +
    geom_point(size = 5, alpha = 0.65) +
    scale_x_discrete(breaks = c("C", "L", "G"),
                     labels = c("Control", "Local", "Global")) +
    scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                       labels = c("Control", "Local", "Global"),
                       values = my_cols[c(8, 2, 6)]) +
    facet_nested(Proj ~ Media_title*Media, scales = "free_y",
                 labeller = labeller(Proj = my_facet_labels,
                                     Media = my_facet_labels)) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(y = expression(paste("Evolved Change in Diauxic Shift Density (k) (",
                              10^8, " cfu/mL)")), 
         x = "Treatment") +
    theme_bw() +
    NULL
  #print(kplot)
  
  temp3 <- gc_sum_pops[gc_sum_pops$Pop != "Anc", ]
  temp3$Title <- "lag time"
  temp3$Media_title <- "Media"
  lagplot <-
    ggplot(data = temp3,
           aes(x = Treat, y = fit2_lagtime_hrs_avg_rel_avg, group = Pop,
               color = Treat)) +
    geom_point(size = 5, alpha = 0.65) +
    scale_x_discrete(breaks = c("C", "L", "G"),
                     labels = c("Control", "Local", "Global")) +
    scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                       labels = c("Control", "Local", "Global"),
                       values = my_cols[c(8, 2, 6)]) +
    facet_nested(Proj ~ Media_title*Media, scales = "free_y",
                 labeller = labeller(Proj = my_facet_labels,
                                     Media = my_facet_labels)) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(y = "Evolved Change in Lag Time (hrs)", x = "Treatment") +
    theme_bw() +
    NULL
  #print(lagplot)
  
  temp4 <- gc_sum_pops[gc_sum_pops$Pop != "Anc", ]
  temp4$Title <- "v"
  temp4$Media_title <- "Media"
  vplot <-
    ggplot(data = temp4,
           aes(x = Treat, y = fit2_v_log10_avg_rel_avg, group = Pop,
               color = Treat)) +
    geom_point(size = 5, alpha = 0.65) +
    scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                       labels = c("Control", "Local", "Global"),
                       values = my_cols[c(8, 2, 6)]) +
    scale_x_discrete(breaks = c("C", "L", "G"),
                     labels = c("Control", "Local", "Global")) +
    facet_nested(Proj ~ Media_title*Media, scales = "free_y",
                 labeller = labeller(Proj = my_facet_labels,
                                     Media = my_facet_labels)) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(y = "Evolved Change in Deceleration Parameter (v)", x = "Treatment") +
    theme_bw() +
    NULL
  #print(vplot)
  
  tiff("./Growth_curve_variables_plots_pops/r_k_lag_v_combined.tiff",
       width = 12, height = 8, units = "in", res = 300)
  print(cowplot::plot_grid(
    rplot + theme(legend.position = "none", 
                  strip.background.y = element_blank(), 
                  strip.text.y = element_blank(),
                  strip.text.x = element_text(size = 18),
                  axis.title = element_text(size = 22),
                  axis.text.y = element_text(size = 20),
                  axis.text.x = element_text(angle = 45, size = 18, hjust = 1)),
    kplot + theme(legend.position = "none",
                  strip.background.y = element_blank(), 
                  strip.text.y = element_blank(),
                  strip.text.x = element_text(size = 18),
                  axis.title = element_text(size = 22),
                  axis.text.y = element_text(size = 20),
                  axis.text.x = element_text(angle = 45, size = 18, hjust = 1)),
    lagplot + theme(legend.position = "none",
                    strip.background.y = element_blank(), 
                    strip.text.y = element_blank(),
                    strip.text = element_text(size = 18),
                    axis.title = element_text(size = 22),
                    axis.text.y = element_text(size = 20),
                    axis.text.x = element_text(angle = 45, size = 18, hjust = 1)),
    vplot + theme(legend.position = "none",
                  strip.text = element_text(size = 18),
                  axis.title = element_text(size = 22),
                  axis.text.y = element_text(size = 20),
                  axis.text.x = element_text(angle = 45, size = 18, hjust = 1)),
    ncol = 4, rel_widths = c(1, 1, 1, 1)))
  dev.off()
}

#Make simple plots of r, k, lagtime, v for talks
if(make_statplots) {
  my_facet_labels <- c("7x" = "Weak Phage", 
                       "125" = "Strong Phage",
                       "C" = "Control", "G" = "Global", "L" = "Local",
                       "A" = "WT",
                       "Rich" = "Rich", "Orig" = "Orig")
  
  temp1 <- gc_sum_pops[gc_sum_pops$Pop != "Anc" &
                         gc_sum_pops$Media == "Orig", ]
  temp1$Title <- "r"
  temp1$Media_title <- "Media"
  rplot <-
    ggplot(data = temp1,
           aes(x = Treat, y = fit2_r_avg_rel_avg, group = Pop,
               color = Treat)) +
    geom_point(size = 5, alpha = 0.65) +
    scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                       labels = c("Control", "Local", "Global"),
                       values = my_cols[c(8, 2, 6)]) +
    scale_x_discrete(breaks = c("C", "L", "G"),
                     labels = c("Control", "Local", "Global")) +
    facet_nested(Proj ~ ., scales = "free_y",
                 labeller = labeller(Proj = my_facet_labels)) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(y = "Evolved Change in Growth Rate (/hr)", x = "Treatment") +
    theme_bw() +
    NULL
  #print(rplot)
  
  temp2 <- gc_sum_pops[gc_sum_pops$Pop != "Anc" &
                         gc_sum_pops$Media == "Orig", ]
  temp2$Title <- "k"
  temp2$Media_title <- "Media"
  kplot <-
    ggplot(data = temp2,
           aes(x = Treat, y = fit2_k_avg_rel_avg/10**8, group = Pop,
               color = Treat)) +
    geom_point(size = 5, alpha = 0.65) +
    scale_x_discrete(breaks = c("C", "L", "G"),
                     labels = c("Control", "Local", "Global")) +
    scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                       labels = c("Control", "Local", "Global"),
                       values = my_cols[c(8, 2, 6)]) +
    facet_nested(Proj ~., scales = "free_y",
                 labeller = labeller(Proj = my_facet_labels)) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(y = expression(paste("Evolved Change in Diauxic Shift Density (x ",
                              10^8, " cfu/mL)")), 
         x = "Treatment") +
    theme_bw() +
    NULL
  #print(kplot)
  
  temp3 <- gc_sum_pops[gc_sum_pops$Pop != "Anc" &
                         gc_sum_pops$Media == "Orig", ]
  temp3$Title <- "lag time"
  temp3$Media_title <- "Media"
  lagplot <-
    ggplot(data = temp3,
           aes(x = Treat, y = fit2_lagtime_hrs_avg_rel_avg, group = Pop,
               color = Treat)) +
    geom_point(size = 5, alpha = 0.65) +
    scale_x_discrete(breaks = c("C", "L", "G"),
                     labels = c("Control", "Local", "Global")) +
    scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                       labels = c("Control", "Local", "Global"),
                       values = my_cols[c(8, 2, 6)]) +
    facet_nested(Proj ~., scales = "free_y",
                 labeller = labeller(Proj = my_facet_labels)) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(y = "Evolved Change in Lag Time (hrs)", x = "Treatment") +
    theme_bw() +
    NULL
  #print(lagplot)
  
  temp4 <- gc_sum_pops[gc_sum_pops$Pop != "Anc" &
                         gc_sum_pops$Media == "Orig", ]
  temp4$Title <- "v"
  temp4$Media_title <- "Media"
  vplot <-
    ggplot(data = temp4,
           aes(x = Treat, y = fit2_v_log10_avg_rel_avg, group = Pop,
               color = Treat)) +
    geom_point(size = 5, alpha = 0.65) +
    scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                       labels = c("Control", "Local", "Global"),
                       values = my_cols[c(8, 2, 6)]) +
    scale_x_discrete(breaks = c("C", "L", "G"),
                     labels = c("Control", "Local", "Global")) +
    facet_nested(Proj ~., scales = "free_y",
                 labeller = labeller(Proj = my_facet_labels)) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(y = "Change in Deceleration Parameter", x = "Treatment") +
    theme_bw() +
    NULL
  #print(vplot)
  
  tiff("./Growth_curve_variables_plots_pops/r_only.tiff",
       width = 5, height = 6, units = "in", res = 100)
  print(
    rplot + 
      theme(legend.position = "none", 
            strip.text = element_text(size = 18),
            axis.title = element_text(size = 18),
            axis.text.y = element_text(size = 16),
            axis.text.x = element_text(angle = 45, size = 18, hjust = 1)))
  dev.off()
  
  tiff("./Growth_curve_variables_plots_pops/k_only.tiff",
       width = 5, height = 6, units = "in", res = 100)
  print(
    kplot +
      theme(legend.position = "none", 
            strip.text = element_text(size = 18),
            axis.title = element_text(size = 16),
            axis.text.y = element_text(size = 16),
            axis.text.x = element_text(angle = 45, size = 18, hjust = 1),
            plot.margin = unit(c(0.7, 0.15, 0.15, 0.15), "in")))
  dev.off()
  
  tiff("./Growth_curve_variables_plots_pops/lag_only.tiff",
       width = 5, height = 6, units = "in", res = 100)
  print(
    lagplot + 
      theme(legend.position = "none", 
            strip.text = element_text(size = 18),
            axis.title = element_text(size = 18),
            axis.text.y = element_text(size = 16),
            axis.text.x = element_text(angle = 45, size = 18, hjust = 1)))
  dev.off()
  
  tiff("./Growth_curve_variables_plots_pops/v_only.tiff",
       width = 5, height = 6, units = "in", res = 100)
  print(
    vplot +
      theme(legend.position = "none", 
            strip.text = element_text(size = 18),
            axis.title = element_text(size = 18),
            axis.text.y = element_text(size = 16),
            axis.text.x = element_text(angle = 45, size = 18, hjust = 1)))
  dev.off()
  
  tiff("./Growth_curve_variables_plots_pops/r_k_lag_combined.tiff",
       width = 8, height = 9, units = "in", res = 300)
  print(cowplot::plot_grid(
    lagplot + theme(legend.position = "none",
                    strip.background.y = element_blank(), 
                    strip.text.y = element_blank(),
                    strip.text = element_text(size = 18),
                    axis.title = element_text(size = 22),
                    axis.text.y = element_text(size = 20),
                    axis.text.x = element_text(angle = 45, size = 18, hjust = 1)) +
      xlab(""),
    rplot + theme(legend.position = "none", 
                  strip.background.y = element_blank(), 
                  strip.text.y = element_blank(),
                  strip.text.x = element_text(size = 18),
                  axis.title = element_text(size = 22),
                  axis.text.y = element_text(size = 20),
                  axis.text.x = element_text(angle = 45, size = 18, hjust = 1)),
    kplot + theme(legend.position = "none",
                  #strip.background.y = element_blank(), 
                  strip.text.y = element_text(size = 18),
                  strip.text.x = element_text(size = 18),
                  axis.title = element_text(size = 22),
                  axis.text.y = element_text(size = 20),
                  axis.text.x = element_text(angle = 45, size = 18, hjust = 1)) +
      xlab(""),
    ncol = 3, rel_widths = c(.85, .85, 1)))
  dev.off()
}

#View population mean data r vs k
if (make_statplots) {
  tiff("./Growth_curve_variables_plots_pops/fitr_fitk.tiff",
       width = 10.5, height = 10, units = "in", res = 300)
  print(ggplot(data = gc_sum_pops[gc_sum_pops$Pop != "Anc", ],
               aes(x = fit_r_avg_rel_avg, y = fit_k_avg_rel_avg, 
                   fill = Treat, shape = Pop)) +
          geom_point(size = 8, alpha = 0.8) +
          scale_fill_manual(name = "Treatment", breaks = c("C", "L", "G"),
                            labels = c("Control", "Local", "Global"),
                            values = my_cols[c(8, 2, 6)]) +
          scale_shape_manual(values = 21:26) +
          facet_grid(Proj ~ Media, scales = "free",
                     labeller = labeller(Proj = my_facet_labels,
                                         Media = my_facet_labels)) +
          geom_hline(yintercept = 0, lty = 2) +
          geom_vline(xintercept = 0, lty = 2) +
          theme_bw() +
          guides(fill=guide_legend(override.aes=list(shape=21))) +
          # geom_errorbar(aes(x = Treat, ymin = get(var)-get(var_sd),
          #                   ymax = get(var)+get(var_sd)),
          #               position = position_dodge(0.3),
          #               width = 0.2)
          NULL)
  dev.off()
}

#View population-summarized median data
# if (make_statplots) {
#   dir.create("./Growth_curve_variables_plots_pops/", showWarnings = F)
#   for (var_root in c("first_min_avg_rel", 
#                      "max_percap_gr_rate_avg_rel", 
#                      "max_percap_gr_dens_avg_rel", 
#                      "max_percap_gr_timesincemin_avg_rel",
#                      "pseudo_K_avg_rel", 
# #                     "pseudo_K_timesincemin_avg_rel"
#                      "pseudo_K_timesince_maxpercap_avg_rel"
#                      )) {
#     var <- paste(var_root, "_med", sep = "")
#     var_sd <- paste(var_root, "_sd", sep = "")
#     tiff(paste("./Growth_curve_variables_plots_pops/", var, ".tiff", sep = ""),
#          width = 10, height = 10, units = "in", res = 300)
#     print(ggplot(data = gc_sum_pops[gc_sum_pops$Pop != "Anc", ],
#                  aes(x = Treat, y = get(var), group = Pop,
#                      color = Treat)) +
#             geom_point(position = position_dodge(0.1),
#                        size = 5, alpha = 0.8) +
#             geom_hline(yintercept = 1, lty = 2) +
#             scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
#                                labels = c("Control", "Local", "Global"),
#                                values = my_cols[c(8, 2, 6)]) +
#             facet_grid(Proj ~ Media, scales = "free_y",
#                        labeller = labeller(Proj = my_facet_labels,
#                                            Media = my_facet_labels)) +
#             ggtitle(var) +
#             theme_bw() +
#             # geom_errorbar(aes(x = Treat, ymin = get(var)-get(var_sd),
#             #                   ymax = get(var)+get(var_sd)),
#             #               position = position_dodge(0.3),
#             #               width = 0.2)
#             NULL)
#     dev.off()
#   }
# }

#Plot predicted curves for each pop in ea media
if (make_statplots) {
  t_vals <- seq(from = 0, to = max(gc_data$Time_s), by = 5*60)
  gc_sum_pops_curves <- data.frame(
    Proj = rep(gc_sum_pops$Proj, each = length(t_vals)),
    Pop = rep(gc_sum_pops$Pop, each = length(t_vals)),
    Treat = rep(gc_sum_pops$Treat, each = length(t_vals)),
    Media = rep(gc_sum_pops$Media, each = length(t_vals)),
    Time_s = rep(t_vals, nrow(gc_sum_pops)),
    Dens = rep(0, nrow(gc_sum_pops)*length(t_vals)))
  startrow <- 1
  for (i in 1:nrow(gc_sum_pops)) {
    gc_sum_pops_curves$Dens[startrow:(startrow+length(t_vals)-1)] <- 
      baranyi_func(
        r = gc_sum_pops$fit2_r_avg_avg[i],
        k = gc_sum_pops$fit2_k_avg_avg[i],
        v = gc_sum_pops$fit2_v_avg_avg[i], 
        q0 = gc_sum_pops$fit2_q0_avg_avg[i], 
        d0 = gc_sum_pops$fit2_d0_avg_avg[i], 
        m = gc_sum_pops$fit2_m_avg_avg[i],
        t_vals = t_vals)
    startrow <- startrow + length(t_vals)
  }
  
  my_facet_labels <- c("7x" = "Weak Phage", 
                       "125" = "Strong Phage",
                       "C" = "Control", "G" = "Global", "L" = "Local",
                       "A" = "WT",
                       "Rich" = "Rich Media", "Orig" = "Original Media")
  tiff("./Output_figures/pred_gcs.tiff",
       width = 10.5, height = 10, units = "in", res = 300)
  ggplot(data = gc_sum_pops_curves,
         aes(x = Time_s/3600, y = Dens, color = Treat, group = paste(Pop, Treat))) +
    geom_line(lwd = 1, alpha = 0.8) +
    facet_grid(Proj~Media, scales = "free_y",
               labeller = labeller(Proj = my_facet_labels,
                                   Media = my_facet_labels)) +
    theme_bw() +
    scale_y_continuous(trans = "log10") +
    xlim(0, 11) +
    scale_color_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                       labels = c("Ancestor", "Control", "Local", "Global"),
                       values = my_cols[c(7, 8, 2, 6)]) +
    NULL
  dev.off()
}

#Cast measurements in different medias into different columns
# (using population mean data)
gc_sum_pops_wide <- tidyr::pivot_wider(gc_sum_pops,
                                       values_from = ends_with("avg_rel_avg"),
                                       names_from = Media,
                                       id_cols = c("Proj", "Pop", "Treat"))

#Add in resistance & migration data
isol_data <- full_join(gc_sum_pops_wide, 
                       resis_data_sum)
isol_data <- as.data.frame(full_join(isol_data, isol_migr_sum))

isol_data$EOP_avg <- log10(isol_data$EOP_avg)
isol_data$resis <- -isol_data$EOP_avg

#Rename for brevity
colnames(isol_data) <- gsub("_avg_rel_avg", "", colnames(isol_data))

#Check correlations between variables
# gc_var_cors_7x <- cor(isol_data[isol_data$Proj == "7x", 
#                                        c(4:16, 18)])
# gc_var_cors_125 <- cor(isol_data[isol_data$Proj == "125", 
#                                         c(4:16, 18)])
# write.csv(gc_var_cors_7x, "grow_curve_var_correlations_7x.csv")
# write.csv(gc_var_cors_125, "grow_curve_var_correlations_125.csv")

# all vars are positive w/ ea other
# max percap rate neg w/ first min
# max percap dens pos w/ first min
# max percap dens pos w/ max percap timesincemin

#Make correlation figures
if (make_statplots) {
  #Make base figure
  temp <- isol_data[isol_data$Treat != "Anc" &
                      isol_data$Proj == "7x", ]
  colnames(temp) <- plyr::revalue(
    colnames(temp),
    replace = c("fit2_r_Orig" = "r Orig", 
                "fit2_r_Rich" = "r Rich", 
                "fit2_k_Orig" = "k Orig", 
                "fit2_k_Rich" = "k Rich",
                "fit2_lagtime_hrs_Orig" = "lag Orig",
                "fit2_lagtime_hrs_Rich" = "lag Rich",
                "radius_mm_hr_rel_avg" = "agar growth",
                "fit2_v_log10_Orig" = "v Orig",
                "fit2_v_log10_Rich" = "v Rich"))
  p <- GGally::ggpairs(temp,
                       columns = c("r Orig", "r Rich", "k Orig", "k Rich",
                                   "lag Orig", "lag Rich", "v Orig", "v Rich",
                                   "agar growth", "resis"),
                       lower = list(continuous = "smooth"),
                       upper = list(continuous = "smooth"),
                       ggplot2::aes(color = Treat, group = Proj),
                       title = "Weak Phage") +
    theme(strip.text = element_text(size = 9),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  #Change colors
  for (i in 1:p$nrow) {
    for (j in 1:p$ncol) {
      p[i, j] <- p[i, j] +
        scale_color_manual(breaks = c("C", "L", "G"),
                           values = my_cols[c(8, 2, 6)])
    }
  }
  tiff("./Output_figures/Weakphage_cors.tiff", width = 10, height = 10, 
       units = "in", res = 300)
  print(p)
  dev.off()
  
  #Make base figure
  temp <- isol_data[isol_data$Treat != "Anc" &
                      isol_data$Proj == "125", ]
  colnames(temp) <- plyr::revalue(
    colnames(temp),
    replace = c("fit2_r_Orig" = "r Orig", 
                "fit2_r_Rich" = "r Rich", 
                "fit2_k_Orig" = "k Orig", 
                "fit2_k_Rich" = "k Rich",
                "fit2_lagtime_hrs_Orig" = "lag Orig",
                "fit2_lagtime_hrs_Rich" = "lag Rich",
                "radius_mm_hr_rel_avg" = "agar growth",
                "fit2_v_log10_Orig" = "v Orig",
                "fit2_v_log10_Rich" = "v Rich"))
  p <- GGally::ggpairs(temp,
                       columns = c("r Orig", "r Rich", "k Orig", "k Rich",
                                   "lag Orig", "lag Rich", "v Orig", "v Rich",
                                   "agar growth", "resis"),
                       lower = list(continuous = "smooth"),
                       upper = list(continuous = "smooth"),
                       ggplot2::aes(color = Treat, group = Proj),
                       title = "Strong Phage") +
    theme(strip.text = element_text(size = 9),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  #Change colors
  for (i in 1:p$nrow) {
    for (j in 1:p$ncol) {
      p[i, j] <- p[i, j] +
        scale_color_manual(breaks = c("C", "L", "G"),
                           values = my_cols[c(8, 2, 6)])
    }
  }
  tiff("./Output_figures/Strongphage_cors.tiff", width = 10, height = 10, 
       units = "in", res = 300)
  print(p)
  dev.off()
  
  tiff("./Output_figures/Heatcors_weak.tiff", width = 10, height = 10, units = "in", res = 300)
  GGally::ggcorr(isol_data[isol_data$Treat != "Anc" &
                             isol_data$Proj == "7x", 
                           c("fit2_r_Orig", "fit2_k_Orig", 
                             "fit2_v_log10_Orig", 
                             #"fit2_d0_Orig", 
                             "fit2_lagtime_hrs_Orig",
                             "fit2_r_Rich", "fit2_k_Rich",
                             "fit2_v_log10_Rich",
                             #"fit2_d0_Rich", 
                             "fit2_lagtime_hrs_Rich",
                             #"EOP_avg", 
                             "resis", "radius_mm_hr_rel_avg"
                           )],
                 nbreaks = 5,
                 hjust = 0.8,
                 layout.exp = 1.5) +
    ggplot2::labs(title = "Weak Phage")
  dev.off()
  
  tiff("./Output_figures/Heatcors_strong.tiff", width = 10, height = 10, units = "in", res = 300)
  GGally::ggcorr(isol_data[isol_data$Treat != "Anc" &
                             isol_data$Proj == "125", 
                           c("fit2_r_Orig", "fit2_k_Orig", 
                             "fit2_v_log10_Orig", 
                             #"fit2_d0_Orig", 
                             "fit2_lagtime_hrs_Orig",
                             "fit2_r_Rich", "fit2_k_Rich",
                             "fit2_v_log10_Rich",
                             #"fit2_d0_Rich", 
                             "fit2_lagtime_hrs_Rich",
                             #"EOP_avg", 
                             "resis", "radius_mm_hr_rel_avg"
                           )],
                 nbreaks = 5,
                 hjust = 0.8,
                 layout.exp = 1.5) +
    ggplot2::labs(title = "Strong Phage")
  dev.off()
}

##Isolate growth curves: Check for normality ----

#Check for univariate normality
if (make_statplots) {
  for (var in c("fit2_r_Orig", "fit2_k_Orig", 
                "fit2_v_log10_Orig", 
                #"fit2_d0_Orig", 
                "fit2_lagtime_hrs_Orig",
                "fit2_r_Rich", "fit2_k_Rich",
                "fit2_v_log10_Rich",
                #"fit2_d0_Rich", 
                "fit2_lagtime_hrs_Rich",
                #"EOP_avg", 
                "resis", "radius_mm_hr_rel_avg"
  )) {
    print(ggplot(data = isol_data[isol_data$Pop != "Anc", ],
                 aes(sample = get(var))) +
            geom_qq() +
            geom_qq_line() +
            facet_grid(~Proj) +
            ggtitle(var))
    print(ggplot(data = isol_data[isol_data$Pop != "Anc", ],
                 aes(x = get(var))) +
            geom_histogram(bins = 10) +
            facet_grid(~Proj) +
            ggtitle(var))
  }
}

#Define function to make chi-square quantile plots 
# to test for multivariate normality of data or residuals
# (credit to Jonathan Reuning-Scherer)
CSQPlot<-function(vars,label="Chi-Square Quantile Plot"){
  #usually, vars is xxx$residuals or data from one group and label is for plot
  x<-cov(scale(vars),use="pairwise.complete.obs")
  squares<-sort(diag(as.matrix(scale(vars))%*%solve(x)%*%as.matrix(t(scale(vars)))))
  quantiles<-quantile(squares)
  hspr<-quantiles[4]-quantiles[2]
  cumprob<-c(1:length(vars[,1]))/length(vars[,1])-1/(2*length(vars[,1]))
  degf<-dim(x)[1]
  quants<-qchisq(cumprob,df=degf)
  gval<-(quants**(-1+degf/2))/(exp(quants/2)*gamma(degf/2)*(sqrt(2)**degf))
  scale<-hspr / (qchisq(.75,degf)-qchisq(.25,degf))
  se<-(scale/gval)*sqrt(cumprob*(1-cumprob)/length(squares))
  lower<-quants-2*se
  upper<-quants+2*se
  
  plot(quants,squares,col='red',pch=19,cex=1.2,xlab="Chi-Square Quantiles",
       ylab=label,main=paste("Chi-Square Quantiles for",label),ylim=range(upper,lower, squares) , xlim=range(c(0,quants)))
  lines(c(0,100),c(0,100),col=1)
  lines(quants,upper,col="blue",lty=2,lwd=2)
  lines(quants,lower,col="blue",lty=2,lwd=2)
  legend(0,range(upper,lower)[2]*.9,c("Data","95% Conf Limits"),lty=c(0,2),col=c("red","blue"),lwd=c(2,2),
         pch=c(19,NA))
}

#Make multivariate normality plots
for (proj in unique(isol_data$Proj)) {
  CSQPlot(isol_data[isol_data$Proj == proj, c("fit2_r_Orig", "fit2_k_Orig", 
                                              "fit2_v_log10_Orig", 
                                              #"fit2_d0_Orig", 
                                              "fit2_lagtime_hrs_Orig",
                                              "fit2_r_Rich", "fit2_k_Rich",
                                              "fit2_v_log10_Rich",
                                              #"fit2_d0_Rich", 
                                              "fit2_lagtime_hrs_Rich",
                                              #"EOP_avg", 
                                              "resis", "radius_mm_hr_rel_avg"
  )], label = proj)
}

#Most variables are univariate normal (except for resis which is ~binary
# and k Rich and v Orig which has overdispersed tails)

#Surprisingly, both 7x and 125 are somewhat multivariate normal

##Isolate data: run PCA (gc Orig only) ----
isol_data_pca <- list(
  "7x" = isol_data[isol_data$Proj == "7x", ],
  "125" = isol_data[isol_data$Proj == "125", ])


isol_prcomp <- list()
for (i in 1:length(isol_data_pca)) {
  use_cols <- c("fit2_r_Orig", "fit2_k_Orig", 
                #"fit2_v_Orig", 
                #"fit2_d0_Orig", 
                "fit2_lagtime_hrs_Orig",
                #"fit2_r_Rich", "fit2_k_Rich", 
                #"fit2_v_Rich", 
                #"fit2_d0_Rich", "fit2_lagtime_hrs_Rich",
                #"EOP_avg", 
                "resis", "radius_mm_hr_rel_avg")
  isol_prcomp[[i]] <- prcomp(isol_data_pca[[i]][, use_cols],
                             center = TRUE, scale = TRUE, retx = TRUE)
  isol_prcomp[[i]]$x <- cbind(isol_data_pca[[i]][, c("Proj", "Pop", "Treat")],
                              as.data.frame(isol_prcomp[[i]]$x))
  row.names(isol_prcomp[[i]]$rotation) <- 
    plyr::revalue(x = row.names(isol_prcomp[[i]]$rotation),
                  replace = c("fit2_r_Orig" = "r", "fit2_k_Orig" = "k", 
                              "fit2_lagtime_hrs_Orig" = "lag time", 
                              "resis" = "resistance",
                              "radius_mm_hr_rel_avg" = "agar growth"))
}
names(isol_prcomp) <- names(isol_data_pca)

summary(isol_prcomp[[1]])
summary(isol_prcomp[[2]])

if(make_statplots) {
  arrow_len <- 2 #multiplier for arrow lengths for vis purposes
  
  tiff("./Output_figures/weakphage_PCA.tiff", 
       width = 12, height = 10, units = "in", res = 300)
  weak_pca <- ggplot(isol_prcomp[["7x"]]$x, 
                     aes(x = PC1, y = PC2)) +
    ggtitle("Weak Phage") +
    geom_segment(data = as.data.frame(isol_prcomp[["7x"]]$rotation),
                 aes(x = 0, y = 0, xend = arrow_len*PC1, yend = arrow_len*PC2),
                 arrow = arrow(length = unit(0.02, "npc")),
                 alpha = .8, lwd = 2, color = "gray15") +
    ggrepel::geom_text_repel(data = as.data.frame(isol_prcomp[["7x"]]$rotation),
                             aes(x = arrow_len*PC1, y = arrow_len*PC2,
                                 label = row.names(isol_prcomp[["7x"]]$rotation)),
                             size = 14, alpha = .8, color = "gray0", seed = 8,
                             min.segment.length = unit(1, "native"),
                             nudge_x = c(0, -0.1, -0.4, 0, 0.5), 
                             nudge_y = c(0, -0.1, 0.25, 0.2, 0)) + 
    geom_point(aes(color = Treat), size = 10, alpha = 0.7) +
    theme_bw() +
    labs(x = paste("PC1 (", 
                   round((100*((isol_prcomp[["7x"]]$sdev)**2)/
                            sum((isol_prcomp[["7x"]]$sdev)**2))[1], 1),
                   "%)", sep = ""),
         y = paste("PC2 (", 
                   round((100*((isol_prcomp[["7x"]]$sdev)**2)/
                            sum((isol_prcomp[["7x"]]$sdev)**2))[2], 1),
                   "%)", sep = "")) +
    theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 22),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 32),
          plot.title = element_text(size = 36)) +
    scale_color_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                       labels = c("Ancestor", "Control", "Local", "Global"),
                       values = my_cols[c(7, 8, 2, 6)]) +
    #xlim(-1, NA) +
    NULL
  print(weak_pca)
  dev.off()
  
  tiff("./Output_figures/strongphage_PCA.tiff", 
       width = 12, height = 10, units = "in", res = 300)
  strong_pca <- ggplot(isol_prcomp[["125"]]$x, 
                       aes(x = -PC1, y = PC2)) +
    ggtitle("Strong Phage") +
    geom_segment(data = as.data.frame(isol_prcomp[["125"]]$rotation),
                 aes(x = 0, y = 0, xend = -arrow_len*PC1, yend = arrow_len*PC2),
                 arrow = arrow(length = unit(0.02, "npc")),
                 alpha = .8, lwd = 2, color = "gray15") +
    ggrepel::geom_text_repel(data = as.data.frame(isol_prcomp[["125"]]$rotation),
                             aes(x = -arrow_len*PC1, y = arrow_len*PC2,
                                 label = row.names(isol_prcomp[["125"]]$rotation)),
                             size = 14, alpha = .8, color = "gray0", seed = 1,
                             min.segment.length = unit(1, "native"),
                             nudge_x = c(-0.1, -0.05, 0, 0, 0.05), 
                             nudge_y = c(0, -0.1, 0, 0, -0.05)) +
    geom_point(aes(color = Treat), size = 10, alpha = 0.7) +
    theme_bw() +
    labs(x = paste("PC1 (", 
                   round((100*((isol_prcomp[["125"]]$sdev)**2)/
                            sum((isol_prcomp[["125"]]$sdev)**2))[1], 1),
                   "%)", sep = ""),
         y = paste("PC2 (", 
                   round((100*((isol_prcomp[["125"]]$sdev)**2)/
                            sum((isol_prcomp[["125"]]$sdev)**2))[2], 1),
                   "%)", sep = "")) +
    theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 22),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 32),
          plot.title = element_text(size = 36)) +
    scale_color_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                       labels = c("Ancestor", "Control", "Local", "Global"),
                       values = my_cols[c(7, 8, 2, 6)]) +
    #xlim(-1, NA) +
    NULL
  print(strong_pca)
  dev.off()
  
  tiff("./Output_figures/PCA_combined.tiff", 
       width = 24, height = 10, units = "in", res = 300)
  print(cowplot::plot_grid(
    weak_pca + theme(legend.position = "none"), 
    strong_pca + theme(legend.position = "none"),
    cowplot::get_legend(weak_pca),
    #labels = c("Weak Phage", "Strong Phage"),
    vjust = 0,
    rel_widths = c(1, 1, 0.3),
    nrow = 1, align = "h", axis = "tb"))
  dev.off()
}

##Isolate data: run PCA (gc only) ----
isol_gc_prcomp <- list()
for (i in 1:length(isol_data_pca)) {
  use_cols <- c("fit2_r_Orig", "fit2_k_Orig", 
                #"fit2_v_Orig", 
                #"fit2_d0_Orig", 
                "fit2_lagtime_hrs_Orig",
                "fit2_r_Rich", "fit2_k_Rich",
                #"fit2_v_Rich",
                #"fit2_d0_Rich", 
                "fit2_lagtime_hrs_Rich"
                #"EOP_avg", 
                #"resis", "radius_mm_hr_rel_avg"
  )
  isol_gc_prcomp[[i]] <- prcomp(isol_data_pca[[i]][, use_cols],
                                center = TRUE, scale = TRUE, retx = TRUE)
  isol_gc_prcomp[[i]]$x <- cbind(isol_data_pca[[i]][, c("Proj", "Pop", "Treat")],
                                 as.data.frame(isol_gc_prcomp[[i]]$x))
  row.names(isol_gc_prcomp[[i]]$rotation) <- 
    plyr::revalue(x = row.names(isol_gc_prcomp[[i]]$rotation),
                  replace = c("fit2_r_Orig" = "r Orig", 
                              "fit2_r_Rich" = "r Rich", 
                              "fit2_k_Orig" = "k Orig", 
                              "fit2_k_Rich" = "k Rich",
                              "fit2_lagtime_hrs_Orig" = "lag time Orig",
                              "fit2_lagtime_hrs_Rich" = "lag time Rich"))
}
names(isol_gc_prcomp) <- names(isol_data_pca)

summary(isol_gc_prcomp[[1]])
summary(isol_gc_prcomp[[2]])

if(make_statplots) {
  arrow_len <- 2 #multiplier for arrow lengths for vis purposes
  
  tiff("./Output_figures/weakphage_PCA_gconly.tiff", 
       width = 12, height = 10, units = "in", res = 300)
  weak_pca <- ggplot(isol_gc_prcomp[["7x"]]$x, 
                     aes(x = PC1, y = PC2)) +
    ggtitle("Weak Phage") +
    geom_segment(data = as.data.frame(isol_gc_prcomp[["7x"]]$rotation),
                 aes(x = 0, y = 0, xend = arrow_len*PC1, yend = arrow_len*PC2),
                 arrow = arrow(length = unit(0.02, "npc")),
                 alpha = .8, lwd = 2, color = "gray15") +
    ggrepel::geom_text_repel(data = as.data.frame(isol_gc_prcomp[["7x"]]$rotation),
                             aes(x = arrow_len*PC1, y = arrow_len*PC2,
                                 label = row.names(isol_gc_prcomp[["7x"]]$rotation)),
                             size = 14, alpha = .8, color = "gray0", seed = 8,
                             min.segment.length = unit(1, "native"),
                             nudge_x = c(0.4, -0.8, 1, 0.2, 0.3, -0.9),
                             nudge_y = c(0.2, -0.2, 0.4, 0, 0.2, 0)) +
    geom_point(aes(color = Treat), size = 10, alpha = 0.7) +
    theme_bw() +
    labs(x = paste("PC1 (", 
                   round((100*((isol_gc_prcomp[["7x"]]$sdev)**2)/
                            sum((isol_gc_prcomp[["7x"]]$sdev)**2))[1], 1),
                   "%)", sep = ""),
         y = paste("PC2 (", 
                   round((100*((isol_gc_prcomp[["7x"]]$sdev)**2)/
                            sum((isol_gc_prcomp[["7x"]]$sdev)**2))[2], 1),
                   "%)", sep = "")) +
    theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 22),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 32),
          plot.title = element_text(size = 36)) +
    scale_color_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                       labels = c("Ancestor", "Control", "Local", "Global"),
                       values = my_cols[c(7, 8, 2, 6)]) +
    xlim(NA, 3) +
    NULL
  print(weak_pca)
  dev.off()
  
  tiff("./Output_figures/strongphage_PCA_gconly.tiff", 
       width = 12, height = 10, units = "in", res = 300)
  strong_pca <- ggplot(isol_gc_prcomp[["125"]]$x, 
                       aes(x = PC1, y = PC2)) +
    ggtitle("Strong Phage") +
    geom_segment(data = as.data.frame(isol_gc_prcomp[["125"]]$rotation),
                 aes(x = 0, y = 0, xend = arrow_len*PC1, yend = arrow_len*PC2),
                 arrow = arrow(length = unit(0.02, "npc")),
                 alpha = .8, lwd = 2, color = "gray15") +
    ggrepel::geom_text_repel(data = as.data.frame(isol_gc_prcomp[["125"]]$rotation),
                             aes(x = arrow_len*PC1, y = arrow_len*PC2,
                                 label = row.names(isol_gc_prcomp[["125"]]$rotation)),
                             size = 14, alpha = .8, color = "gray0", seed = 1,
                             min.segment.length = unit(1, "native"),
                             nudge_x = c(0.1, -0.6, 0, -0.55, 0.6, 0), 
                             nudge_y = c(0 ,-0.2, 0.3, 0.2, -0.2, -0.1)) +
    geom_point(aes(color = Treat), size = 10, alpha = 0.7) +
    theme_bw() +
    labs(x = paste("PC1 (", 
                   round((100*((isol_gc_prcomp[["125"]]$sdev)**2)/
                            sum((isol_gc_prcomp[["125"]]$sdev)**2))[1], 1),
                   "%)", sep = ""),
         y = paste("PC2 (", 
                   round((100*((isol_gc_prcomp[["125"]]$sdev)**2)/
                            sum((isol_gc_prcomp[["125"]]$sdev)**2))[2], 1),
                   "%)", sep = "")) +
    theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 22),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 32),
          plot.title = element_text(size = 36)) +
    scale_color_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                       labels = c("Ancestor", "Control", "Local", "Global"),
                       values = my_cols[c(7, 8, 2, 6)]) +
    #xlim(-1, NA) +
    NULL
  print(strong_pca)
  dev.off()
  
  tiff("./Output_figures/PCA_gconly_combined.tiff", 
       width = 24, height = 10, units = "in", res = 300)
  print(cowplot::plot_grid(
    weak_pca + theme(legend.position = "none"), 
    strong_pca + theme(legend.position = "none"),
    cowplot::get_legend(weak_pca),
    #labels = c("Weak Phage", "Strong Phage"),
    vjust = 0,
    rel_widths = c(1, 1, 0.3),
    nrow = 1, align = "h", axis = "tb"))
  dev.off()
}

##Isolate data: run PCA (all vars) ----
isol_all_prcomp <- list()
for (i in 1:length(isol_data_pca)) {
  use_cols <- c("fit2_r_Orig", "fit2_k_Orig", 
                "fit2_v_log10_Orig", 
                #"fit2_d0_Orig", 
                "fit2_lagtime_hrs_Orig",
                "fit2_r_Rich", "fit2_k_Rich",
                "fit2_v_log10_Rich",
                #"fit2_d0_Rich", 
                "fit2_lagtime_hrs_Rich",
                #"EOP_avg", 
                "resis", "radius_mm_hr_rel_avg"
  )
  isol_all_prcomp[[i]] <- prcomp(isol_data_pca[[i]][, use_cols],
                                 center = TRUE, scale = TRUE, retx = TRUE)
  isol_all_prcomp[[i]]$x <- cbind(isol_data_pca[[i]][, c("Proj", "Pop", "Treat")],
                                  as.data.frame(isol_all_prcomp[[i]]$x))
  row.names(isol_all_prcomp[[i]]$rotation) <- 
    plyr::revalue(x = row.names(isol_all_prcomp[[i]]$rotation),
                  replace = c("fit2_r_Orig" = "r Orig", 
                              "fit2_r_Rich" = "r Rich", 
                              "fit2_k_Orig" = "k Orig", 
                              "fit2_k_Rich" = "k Rich",
                              "fit2_lagtime_hrs_Orig" = "lag Orig",
                              "fit2_lagtime_hrs_Rich" = "lag Rich",
                              "radius_mm_hr_rel_avg" = "agar growth",
                              "fit2_v_log10_Orig" = "v Orig",
                              "fit2_v_log10_Rich" = "v Rich"))
}
names(isol_all_prcomp) <- names(isol_data_pca)

summary(isol_all_prcomp[["7x"]])
summary(isol_all_prcomp[["125"]])

isol_all_prcomp[["7x"]]$rotation
isol_all_prcomp[["125"]]$rotation

#Make distance biplots (scores and loading directly plotted)
if(make_statplots) {
  arrow_len <- 4 #multiplier for arrow lengths for vis purposes
  
  isol_all_prcomp[["7x"]]$x$Title <- "Weak Phage"
  
  tiff("./Output_figures/weakphage_PCA12_all.tiff", 
       width = 12, height = 10, units = "in", res = 300)
  weak_pca12 <- ggplot(isol_all_prcomp[["7x"]]$x, 
                       aes(y = PC1, x = PC2)) +
    #ggtitle("Weak Phage") +
    facet_grid(Title~.) +
    geom_segment(data = as.data.frame(isol_all_prcomp[["7x"]]$rotation),
                 aes(x = 0, y = 0, yend = arrow_len*PC1, xend = arrow_len*PC2),
                 arrow = arrow(length = unit(0.02, "npc")),
                 alpha = .8, lwd = 2, color = "gray15") +
    ggrepel::geom_text_repel(data = as.data.frame(isol_all_prcomp[["7x"]]$rotation),
                             aes(y = arrow_len*PC1, x = arrow_len*PC2,
                                 label = row.names(isol_all_prcomp[["7x"]]$rotation)),
                             size = 14, alpha = .8, color = "gray0", seed = 8,
                             min.segment.length = unit(1, "native"),
                             nudge_x = c(0,0,0,0,0,0,0,0),
                             nudge_y = c(0,0,0,0,0,0,0,0)) +
    geom_point(aes(color = Treat), size = 10, alpha = 0.7) +
    theme_bw() +
    labs(y = paste("PC1 (", 
                   round((100*((isol_all_prcomp[["7x"]]$sdev)**2)/
                            sum((isol_all_prcomp[["7x"]]$sdev)**2))[1], 1),
                   "%)", sep = ""),
         x = paste("PC2 (", 
                   round((100*((isol_all_prcomp[["7x"]]$sdev)**2)/
                            sum((isol_all_prcomp[["7x"]]$sdev)**2))[2], 1),
                   "%)", sep = "")) +
    theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 22),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 32),
          strip.text = element_text(size = 36)) +
    scale_color_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                       labels = c("Ancestor", "Control", "Local", "Global"),
                       values = my_cols[c(7, 8, 2, 6)]) +
    #xlim(NA, 3) +
    NULL
  print(weak_pca12)
  dev.off()
  
  tiff("./Output_figures/weakphage_PCA13_all.tiff", 
       width = 12, height = 10, units = "in", res = 300)
  weak_pca13 <- ggplot(isol_all_prcomp[["7x"]]$x, 
                       aes(y = PC1, x = PC3)) +
    #ggtitle("Weak Phage") +
    facet_grid(Title~.) +
    geom_segment(data = as.data.frame(isol_all_prcomp[["7x"]]$rotation),
                 aes(x = 0, y = 0, yend = arrow_len*PC1, xend = arrow_len*PC3),
                 arrow = arrow(length = unit(0.02, "npc")),
                 alpha = .8, lwd = 2, color = "gray15") +
    ggrepel::geom_text_repel(data = as.data.frame(isol_all_prcomp[["7x"]]$rotation),
                             aes(y = arrow_len*PC1, x = arrow_len*PC3,
                                 label = row.names(isol_all_prcomp[["7x"]]$rotation)),
                             size = 14, alpha = .8, color = "gray0", seed = 8,
                             min.segment.length = unit(1, "native"),
                             nudge_x = c(0,0,0,0,0,0,0,0),
                             nudge_y = c(0,0,0,0,0,0,0,0)) +
    geom_point(aes(color = Treat), size = 10, alpha = 0.7) +
    theme_bw() +
    labs(y = paste("PC1 (", 
                   round((100*((isol_all_prcomp[["7x"]]$sdev)**2)/
                            sum((isol_all_prcomp[["7x"]]$sdev)**2))[1], 1),
                   "%)", sep = ""),
         x = paste("PC3 (", 
                   round((100*((isol_all_prcomp[["7x"]]$sdev)**2)/
                            sum((isol_all_prcomp[["7x"]]$sdev)**2))[3], 1),
                   "%)", sep = "")) +
    theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 22),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 32),
          strip.text = element_text(size = 36)) +
    scale_color_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                       labels = c("Ancestor", "Control", "Local", "Global"),
                       values = my_cols[c(7, 8, 2, 6)]) +
    #xlim(NA, 3) +
    NULL
  print(weak_pca13)
  dev.off()
  
  isol_all_prcomp[["125"]]$x$Title <- "Strong Phage"
  
  tiff("./Output_figures/strongphage_PCA12_all.tiff", 
       width = 12, height = 10, units = "in", res = 300)
  strong_pca12 <- ggplot(isol_all_prcomp[["125"]]$x, 
                         aes(y = PC1, x = PC2)) +
    #ggtitle("Strong Phage") +
    facet_grid(Title~.) +
    geom_segment(data = as.data.frame(isol_all_prcomp[["125"]]$rotation),
                 aes(x = 0, y = 0, yend = arrow_len*PC1, xend = arrow_len*PC2),
                 arrow = arrow(length = unit(0.02, "npc")),
                 alpha = .8, lwd = 2, color = "gray15") +
    ggrepel::geom_text_repel(data = as.data.frame(isol_all_prcomp[["125"]]$rotation),
                             aes(y = arrow_len*PC1, x = arrow_len*PC2,
                                 label = row.names(isol_all_prcomp[["125"]]$rotation)),
                             size = 14, alpha = .8, color = "gray0", seed = 1,
                             min.segment.length = unit(1, "native"),
                             nudge_x = c(0,0,0,0,0,0,0,0), 
                             nudge_y = c(0,0,0,0,0,0,0,0)) +
    geom_point(aes(color = Treat), size = 10, alpha = 0.7) +
    theme_bw() +
    labs(y = paste("PC1 (", 
                   round((100*((isol_all_prcomp[["125"]]$sdev)**2)/
                            sum((isol_all_prcomp[["125"]]$sdev)**2))[1], 1),
                   "%)", sep = ""),
         x = paste("PC2 (", 
                   round((100*((isol_all_prcomp[["125"]]$sdev)**2)/
                            sum((isol_all_prcomp[["125"]]$sdev)**2))[2], 1),
                   "%)", sep = "")) +
    theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 22),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 32),
          strip.text = element_text(size = 36)) +
    scale_color_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                       labels = c("Ancestor", "Control", "Local", "Global"),
                       values = my_cols[c(7, 8, 2, 6)]) +
    #xlim(-1, NA) +
    NULL
  print(strong_pca12)
  dev.off()
  
  tiff("./Output_figures/strongphage_PCA13_all.tiff", 
       width = 12, height = 10, units = "in", res = 300)
  strong_pca13 <- ggplot(isol_all_prcomp[["125"]]$x, 
                         aes(y = PC1, x = PC3)) +
    #ggtitle("Strong Phage") +
    facet_grid(Title~.) +
    geom_segment(data = as.data.frame(isol_all_prcomp[["125"]]$rotation),
                 aes(x = 0, y = 0, yend = arrow_len*PC1, xend = arrow_len*PC3),
                 arrow = arrow(length = unit(0.02, "npc")),
                 alpha = .8, lwd = 2, color = "gray15") +
    ggrepel::geom_text_repel(data = as.data.frame(isol_all_prcomp[["125"]]$rotation),
                             aes(y = arrow_len*PC1, x = arrow_len*PC3,
                                 label = row.names(isol_all_prcomp[["125"]]$rotation)),
                             size = 14, alpha = .8, color = "gray0", seed = 1,
                             min.segment.length = unit(1, "native"),
                             nudge_x = c(0,0,0,0,0,0,0,0), 
                             nudge_y = c(0,0,0,0,0,0,0,0)) +
    geom_point(aes(color = Treat), size = 10, alpha = 0.7) +
    theme_bw() +
    labs(y = paste("PC1 (", 
                   round((100*((isol_all_prcomp[["125"]]$sdev)**2)/
                            sum((isol_all_prcomp[["125"]]$sdev)**2))[1], 1),
                   "%)", sep = ""),
         x = paste("PC3 (", 
                   round((100*((isol_all_prcomp[["125"]]$sdev)**2)/
                            sum((isol_all_prcomp[["125"]]$sdev)**2))[3], 1),
                   "%)", sep = "")) +
    theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 22),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 32),
          strip.text = element_text(size = 36)) +
    scale_color_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                       labels = c("Ancestor", "Control", "Local", "Global"),
                       values = my_cols[c(7, 8, 2, 6)]) +
    #xlim(-1, NA) +
    NULL
  print(strong_pca13)
  dev.off()
  
  tiff("./Output_figures/PCA_all_combined.tiff", 
       width = 24, height = 22, units = "in", res = 300)
  print(
    cowplot::plot_grid(
      cowplot::plot_grid(
        weak_pca12 + theme(legend.position = "none",
                           strip.background = element_blank(),
                           strip.text = element_blank()), 
        weak_pca13 + theme(legend.position = "none"),
        strong_pca12 + theme(legend.position = "none",
                             strip.background = element_blank(),
                             strip.text = element_blank()),
        strong_pca13 + theme(legend.position = "none"),
        ncol = 2, align = "hv", axis = "lrtb"),
      cowplot::get_legend(weak_pca12),
      rel_widths = c(2, 0.3), ncol = 2))
  dev.off()
}

#Do normalizations needed for correlation biplot
#In short: pca_res$x %*% diag(1/pca_res$sdev)and pca_res$rotation %*% diag(pca_res$sdev)
isol_all_prcomp_corrbiplot <- list()
for (i in 1:length(isol_all_prcomp)) {
  #adjust rotations all in one go
  isol_all_prcomp_corrbiplot[[i]] <-
    list("rotation" = isol_all_prcomp[[i]]$rotation %*% 
           diag(isol_all_prcomp[[i]]$sdev),
         "x" = isol_all_prcomp[[i]]$x,
         "sdev" = isol_all_prcomp[[i]]$sdev)
  #adjust x column by column
  pc_cntr <- 1
  for (j in 1:ncol(isol_all_prcomp_corrbiplot[[i]]$x)) {
    if (grepl("PC", colnames(isol_all_prcomp_corrbiplot[[i]]$x)[j])) {
      isol_all_prcomp_corrbiplot[[i]]$x[, j] <-
        isol_all_prcomp_corrbiplot[[i]]$x[, j] * 1/isol_all_prcomp[[i]]$sdev[pc_cntr]
      pc_cntr <- pc_cntr + 1
    }
  }
  colnames(isol_all_prcomp_corrbiplot[[i]]$rotation) <-
    colnames(isol_all_prcomp[[i]]$rotation)
}
names(isol_all_prcomp_corrbiplot) <- names(isol_all_prcomp)

#Plot correlation biplots
if(make_statplots) {
  arrow_len <- 1.5 #multiplier for arrow lengths for vis purposes
  
  isol_all_prcomp_corrbiplot[["7x"]]$x$Title <- "Weak Phage"
  
  tiff("./Output_figures/weakphage_PCA12_all_corbiplot.tiff", 
       width = 12, height = 10, units = "in", res = 300)
  weak_pca12 <- ggplot(isol_all_prcomp_corrbiplot[["7x"]]$x, 
                       aes(y = PC1, x = PC2)) +
    #ggtitle("Weak Phage") +
    facet_grid(Title~.) +
    geom_segment(data = as.data.frame(isol_all_prcomp_corrbiplot[["7x"]]$rotation),
                 aes(x = 0, y = 0, yend = arrow_len*PC1, xend = arrow_len*PC2),
                 arrow = arrow(length = unit(0.02, "npc")),
                 alpha = .8, lwd = 2, color = "firebrick4") +
    ggrepel::geom_text_repel(data = as.data.frame(isol_all_prcomp_corrbiplot[["7x"]]$rotation),
                             aes(y = arrow_len*PC1, x = arrow_len*PC2,
                                 label = row.names(isol_all_prcomp_corrbiplot[["7x"]]$rotation)),
                             size = 14, alpha = .8, color = "firebrick4", seed = 8,
                             min.segment.length = unit(1, "native"),
                             nudge_x = c(0,0,0,0,0,0,0,0),
                             nudge_y = c(0,0,0,0,0,0,0,0)) +
    geom_point(aes(color = Treat), size = 10, alpha = 0.7) +
    theme_bw() +
    labs(y = paste("PC1 (", 
                   round((100*((isol_all_prcomp_corrbiplot[["7x"]]$sdev)**2)/
                            sum((isol_all_prcomp_corrbiplot[["7x"]]$sdev)**2))[1], 1),
                   "%)", sep = ""),
         x = paste("PC2 (", 
                   round((100*((isol_all_prcomp_corrbiplot[["7x"]]$sdev)**2)/
                            sum((isol_all_prcomp_corrbiplot[["7x"]]$sdev)**2))[2], 1),
                   "%)", sep = "")) +
    theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 22),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 32),
          strip.text = element_text(size = 36)) +
    scale_color_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                       labels = c("Ancestor", "Control", "Local", "Global"),
                       values = my_cols[c(7, 8, 2, 6)]) +
    #xlim(NA, 3) +
    ylim(-1.6, NA) +
    NULL
  print(weak_pca12)
  dev.off()
  
  tiff("./Output_figures/weakphage_PCA13_all_corbiplot.tiff", 
       width = 12, height = 10, units = "in", res = 300)
  weak_pca13 <- ggplot(isol_all_prcomp_corrbiplot[["7x"]]$x, 
                       aes(y = PC1, x = PC3)) +
    #ggtitle("Weak Phage") +
    facet_grid(Title~.) +
    geom_segment(data = as.data.frame(isol_all_prcomp_corrbiplot[["7x"]]$rotation),
                 aes(x = 0, y = 0, yend = arrow_len*PC1, xend = arrow_len*PC3),
                 arrow = arrow(length = unit(0.02, "npc")),
                 alpha = .8, lwd = 2, color = "firebrick4") +
    ggrepel::geom_text_repel(data = as.data.frame(isol_all_prcomp_corrbiplot[["7x"]]$rotation),
                             aes(y = arrow_len*PC1, x = arrow_len*PC3,
                                 label = row.names(isol_all_prcomp_corrbiplot[["7x"]]$rotation)),
                             size = 14, alpha = .8, color = "firebrick4", seed = 8,
                             min.segment.length = unit(1, "native"),
                             nudge_x = c(0,0,0,0,0,0,0,0),
                             nudge_y = c(0,0,0,0,0,0,0,0)) +
    geom_point(aes(color = Treat), size = 10, alpha = 0.7) +
    theme_bw() +
    labs(y = paste("PC1 (", 
                   round((100*((isol_all_prcomp_corrbiplot[["7x"]]$sdev)**2)/
                            sum((isol_all_prcomp_corrbiplot[["7x"]]$sdev)**2))[1], 1),
                   "%)", sep = ""),
         x = paste("PC3 (", 
                   round((100*((isol_all_prcomp_corrbiplot[["7x"]]$sdev)**2)/
                            sum((isol_all_prcomp_corrbiplot[["7x"]]$sdev)**2))[3], 1),
                   "%)", sep = "")) +
    theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 22),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 32),
          strip.text = element_text(size = 36)) +
    scale_color_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                       labels = c("Ancestor", "Control", "Local", "Global"),
                       values = my_cols[c(7, 8, 2, 6)]) +
    #xlim(NA, 3) +
    NULL
  print(weak_pca13)
  dev.off()
  
  isol_all_prcomp_corrbiplot[["125"]]$x$Title <- "Strong Phage"
  
  tiff("./Output_figures/strongphage_PCA12_all_corbiplot.tiff", 
       width = 12, height = 10, units = "in", res = 300)
  strong_pca12 <- ggplot(isol_all_prcomp_corrbiplot[["125"]]$x, 
                         aes(y = PC1, x = PC2)) +
    #ggtitle("Strong Phage") +
    facet_grid(Title~.) +
    geom_segment(data = as.data.frame(isol_all_prcomp_corrbiplot[["125"]]$rotation),
                 aes(x = 0, y = 0, yend = arrow_len*PC1, xend = arrow_len*PC2),
                 arrow = arrow(length = unit(0.02, "npc")),
                 alpha = .8, lwd = 2, color = "firebrick4") +
    ggrepel::geom_text_repel(data = as.data.frame(isol_all_prcomp_corrbiplot[["125"]]$rotation),
                             aes(y = arrow_len*PC1, x = arrow_len*PC2,
                                 label = row.names(isol_all_prcomp_corrbiplot[["125"]]$rotation)),
                             size = 14, alpha = .8, color = "firebrick4", seed = 1,
                             min.segment.length = unit(1, "native"),
                             nudge_x = c(0,0,0,0,0,0,0,0), 
                             nudge_y = c(0,0,0,0,0,0,0,0)) +
    geom_point(aes(color = Treat), size = 10, alpha = 0.7) +
    theme_bw() +
    labs(y = paste("PC1 (", 
                   round((100*((isol_all_prcomp_corrbiplot[["125"]]$sdev)**2)/
                            sum((isol_all_prcomp_corrbiplot[["125"]]$sdev)**2))[1], 1),
                   "%)", sep = ""),
         x = paste("PC2 (", 
                   round((100*((isol_all_prcomp_corrbiplot[["125"]]$sdev)**2)/
                            sum((isol_all_prcomp_corrbiplot[["125"]]$sdev)**2))[2], 1),
                   "%)", sep = "")) +
    theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 22),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 32),
          strip.text = element_text(size = 36)) +
    scale_color_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                       labels = c("Ancestor", "Control", "Local", "Global"),
                       values = my_cols[c(7, 8, 2, 6)]) +
    #xlim(-1, NA) +
    NULL
  print(strong_pca12)
  dev.off()
  
  tiff("./Output_figures/strongphage_PCA13_all_corbiplot.tiff", 
       width = 12, height = 10, units = "in", res = 300)
  strong_pca13 <- ggplot(isol_all_prcomp_corrbiplot[["125"]]$x, 
                         aes(y = PC1, x = PC3)) +
    #ggtitle("Strong Phage") +
    facet_grid(Title~.) +
    geom_segment(data = as.data.frame(isol_all_prcomp_corrbiplot[["125"]]$rotation),
                 aes(x = 0, y = 0, yend = arrow_len*PC1, xend = arrow_len*PC3),
                 arrow = arrow(length = unit(0.02, "npc")),
                 alpha = .8, lwd = 2, color = "firebrick4") +
    ggrepel::geom_text_repel(data = as.data.frame(isol_all_prcomp_corrbiplot[["125"]]$rotation),
                             aes(y = arrow_len*PC1, x = arrow_len*PC3,
                                 label = row.names(isol_all_prcomp_corrbiplot[["125"]]$rotation)),
                             size = 14, alpha = .8, color = "firebrick4", seed = 1,
                             min.segment.length = unit(1, "native"),
                             nudge_x = c(0,0,0,0,0,0,0,0), 
                             nudge_y = c(0,0,0,0,0,0,0,0)) +
    geom_point(aes(color = Treat), size = 10, alpha = 0.7) +
    theme_bw() +
    labs(y = paste("PC1 (", 
                   round((100*((isol_all_prcomp_corrbiplot[["125"]]$sdev)**2)/
                            sum((isol_all_prcomp_corrbiplot[["125"]]$sdev)**2))[1], 1),
                   "%)", sep = ""),
         x = paste("PC3 (", 
                   round((100*((isol_all_prcomp_corrbiplot[["125"]]$sdev)**2)/
                            sum((isol_all_prcomp_corrbiplot[["125"]]$sdev)**2))[3], 1),
                   "%)", sep = "")) +
    theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 22),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 32),
          strip.text = element_text(size = 36)) +
    scale_color_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                       labels = c("Ancestor", "Control", "Local", "Global"),
                       values = my_cols[c(7, 8, 2, 6)]) +
    #xlim(-1, NA) +
    NULL
  print(strong_pca13)
  dev.off()
  
  tiff("./Output_figures/PCA_all_combined_corbiplot.tiff", 
       width = 24, height = 22, units = "in", res = 300)
  print(
    cowplot::plot_grid(
      cowplot::plot_grid(
        weak_pca12 + theme(legend.position = "none",
                           strip.background = element_blank(),
                           strip.text = element_blank()), 
        weak_pca13 + theme(legend.position = "none"),
        strong_pca12 + theme(legend.position = "none",
                             strip.background = element_blank(),
                             strip.text = element_blank()),
        strong_pca13 + theme(legend.position = "none"),
        ncol = 2, align = "hv", axis = "lrtb"),
      cowplot::get_legend(weak_pca12),
      rel_widths = c(2, 0.3), ncol = 2))
  dev.off()
}

##Isolate data: statistical tests ----

#Manova to test for difference in groups
isol_data_7x <- isol_data[isol_data$Proj == "7x" & isol_data$Pop != "Anc", ]
manova_res_7x <- manova(as.matrix(
  isol_data_7x[, c("fit2_r_Orig", "fit2_k_Orig", 
                   "fit2_v_log10_Orig", 
                   #"fit2_d0_Orig", 
                   "fit2_lagtime_hrs_Orig",
                   "fit2_r_Rich", "fit2_k_Rich",
                   "fit2_v_log10_Rich",
                   #"fit2_d0_Rich", 
                   "fit2_lagtime_hrs_Rich",
                   #"EOP_avg", 
                   "resis", "radius_mm_hr_rel_avg")]) ~ 
    isol_data_7x$Treat)
isol_data_125 <- isol_data[isol_data$Proj == "125" & isol_data$Pop != "Anc", ]
manova_res_125 <- manova(as.matrix(
  isol_data_125[, c("fit2_r_Orig", "fit2_k_Orig", 
                    "fit2_v_log10_Orig", 
                    #"fit2_d0_Orig", 
                    "fit2_lagtime_hrs_Orig",
                    "fit2_r_Rich", "fit2_k_Rich",
                    "fit2_v_log10_Rich",
                    #"fit2_d0_Rich", 
                    "fit2_lagtime_hrs_Rich",
                    #"EOP_avg", 
                    "resis", "radius_mm_hr_rel_avg")]) ~ 
    isol_data_125$Treat)

#7x results
summary(manova_res_7x, test = "Wilks")
summary.aov(manova_res_7x)
p.adjust(summary.aov(manova_res_7x)$` Response fit2_lagtime_hrs_Orig`$`Pr(>F)`,
         n = length(summary.aov(manova_res_7x)) + 
           length(summary.aov(manova_res_125)),
         method = "holm")
p.adjust(summary.aov(manova_res_7x)$` Response fit2_lagtime_hrs_Rich`$`Pr(>F)`,
         n = length(summary.aov(manova_res_7x)) + 
           length(summary.aov(manova_res_125)),
         method = "holm")
p.adjust(summary.aov(manova_res_7x)$` Response resis`$`Pr(>F)`,
         n = length(summary.aov(manova_res_7x)) + 
           length(summary.aov(manova_res_125)),
         method = "holm")
TukeyHSD(aov(lm(resis ~ Treat, data = isol_data_7x)))
p.adjust(c(
  "C" = t.test(isol_data_7x$radius_mm_hr_rel_avg[isol_data_7x$Treat == "C"], 
               mu = 0)$p.value,
  "L" = t.test(isol_data_7x$radius_mm_hr_rel_avg[isol_data_7x$Treat == "L"], 
               mu = 0)$p.value,
  "G" = t.test(isol_data_7x$radius_mm_hr_rel_avg[isol_data_7x$Treat == "G"], 
               mu = 0)$p.value),
  method = "holm", n = 6)

#125 results
summary(manova_res_125, test = "Wilks")
summary.aov(manova_res_125)
p.adjust(summary.aov(manova_res_125)$` Response radius_mm_hr_rel_avg`$`Pr(>F)`,
         n = length(summary.aov(manova_res_7x)) + 
           length(summary.aov(manova_res_125)),
         method = "holm")
p.adjust(summary.aov(manova_res_125)$` Response resis`$`Pr(>F)`,
         n = length(summary.aov(manova_res_7x)) + 
           length(summary.aov(manova_res_125)),
         method = "holm")
TukeyHSD(aov(lm(resis ~ Treat, data = isol_data_125)))
p.adjust(c(
  "C" = t.test(isol_data_125$radius_mm_hr_rel_avg[isol_data_125$Treat == "C"], 
               mu = 0)$p.value,
  "L" = t.test(isol_data_125$radius_mm_hr_rel_avg[isol_data_125$Treat == "L"], 
               mu = 0)$p.value,
  "G" = t.test(isol_data_125$radius_mm_hr_rel_avg[isol_data_125$Treat == "G"], 
               mu = 0)$p.value),
  method = "holm", n = 6)
