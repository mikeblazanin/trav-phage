##TODO: 
##      migration stats
##      growth curve analysis
##      growth curve stats
##      

library("ggplot2")
library("dplyr")
library("data.table")

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")

#Need to split this into 2 scripts: one to take the messy data
# (e.g. isolate information in underscore formats, etc)
#and get it cleaned up, then output those clean data
#And another to take the clean data & analyze it
#Then put the messy data in one folder and the clean data in another folder
#Along with the separate scripts
#So when it's published, can just tell people to re-run the clean data script

#need to make Rep vs Pop consistent
#redo all analysis using dplyr
#redo analysis using vectorization
#include std curve for spec in this code (it may alread be there???)
#Do a scripted analysis of how much to smoothe gc data (not manually checking)
# #to be done later:
# #define function to get means of rates
# mean_rates <- function(sub_by_list, my_data)
# analyze additional plate scans
#change resistance to be by treat then by pop
#Change growth rate for ancestor to dotted line

##Experimental evolution migration ----
exper_evol_migr <- read.csv("./Clean_Data/Experimental_evolution_growth.csv")

#Drop points after T14
exper_evol_migr <- exper_evol_migr[exper_evol_migr$Timepoint <= 14, ]

#Calculate total area
exper_evol_migr$area_cm2 <- pi*exper_evol_migr$Width_cm/2*exper_evol_migr$Height_cm/2

#Make plot with all pops shown
ggplot(data = exper_evol_migr,
       aes(x = Timepoint, y = area_cm2, group = paste(Treat, Pop),
           color = Treat)) +
  geom_line() +
  scale_color_manual(values = my_cols[c(4, 5, 7)]) +
  facet_grid(~Proj)

#Summarize
exper_evol_migr <- group_by(exper_evol_migr, Proj, Treat, Timepoint)
exper_evol_summ <- summarize(exper_evol_migr,
                             area_mean = mean(area_cm2),
                             area_sd = sd(area_cm2),
                             area_n = n())

#Make plot of summarized data
my_facet_labels <- c("7x" = "Weak Phage", "125" = "Strong Phage")

ggplot(data = exper_evol_summ, aes(x = Timepoint, y = area_mean,
                                   color = Treat)) +
  geom_point(position = position_dodge(0.2)) + 
  geom_line(size = 1.2, position = position_dodge(0.2)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 11), axis.text.x = element_text(size = 11),
        legend.text = element_text(size = 16)) +
  facet_grid(~Proj, labeller = labeller(Proj = my_facet_labels)) +
  geom_errorbar(aes(ymax = area_mean+area_sd, ymin = area_mean-area_sd),
                width=1, size = .7, position=position_dodge(0.2)) +
  labs(x = "Transfer", 
       y = expression(paste("Mean Area of Growth ( ", cm^2, ")"))) + 
  scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                  labels = c("Control", "Local", "Global"),
                  values = my_cols[c(8, 2, 6)]) +
  scale_x_continuous(breaks = c(0, 7, 14)) +
  NULL

#Make plot with both summarized and non-summarized data
ggplot(data = exper_evol_migr,
               aes(x = Timepoint, y = area_cm2, group = paste(Treat, Pop),
                   color = Treat)) +
         geom_line(alpha = 0.5) +
         facet_grid(~Proj, labeller = labeller(Proj = my_facet_labels)) +
  theme_bw() +
  geom_line(data = exper_evol_summ,
            aes(x = Timepoint, y = area_mean, color = Treat,
                group = Treat),
            size = 1.2) +
  labs(x = "Transfer", 
       y = expression(paste("Area of Growth ( ", cm^2, ")"))) + 
  scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                     labels = c("Control", "Local", "Global"),
                     values = my_cols[c(8, 2, 6)]) +
  scale_x_continuous(breaks = c(0, 7, 14)) +
  NULL
    

##Isolate migration ----
isol_migration <- read.csv("./Clean_Data/Isolate_migration.csv")

#Calculate total area
isol_migration$area_cm2 <- pi*isol_migration$Width_cm/2*isol_migration$Height_cm/2

#Calculate total area relative to same-day ancestor
ancestors <- isol_migration[isol_migration$Isol == "Anc", ]
isol_migration$relative_area <-
  isol_migration$area_cm2/ancestors$area_cm2[
    match(as.Date(isol_migration$end_timestamp),
          as.Date(ancestors$end_timestamp))]

#Plot data
my_facet_labels <- c("7x" = "Weak Phage", 
                     "125" = "Strong Phage",
                     "C" = "Control", "G" = "Global", "L" = "Local",
                     "A" = "WT")

ggplot(isol_migration[isol_migration$Isol != "Anc", ], 
       aes(x = Treat, y = relative_area, color = Treat,
           group = Pop)) +
  geom_point(position = position_dodge(0.5)) +
  facet_grid(~Proj, labeller = labeller(Proj = my_facet_labels)) +
  theme_bw() + 
  labs(y = "T14 Isolate Area of Growth Relative to Ancestor",
       x = "Treatment") +
  geom_hline(yintercept = 1, lty = 2) +
  scale_x_discrete(limits = c("C", "L", "G"),
                   labels = c("Control", "Local", "Global")) +
  scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                     labels = c("Control", "Local", "Global"),
                     values = my_cols[c(8, 2, 6)]) +
  NULL

## Isolate resistance ----
resis_data <- read.csv("./Clean_Data/Isolate_resistance.csv",
                       stringsAsFactors = F)

#calculate EOP for ea isol
resis_data$EOP <- NA
resis_data$bd <- F
#First Handle below-detection points
my_rows <- which(resis_data$PFU == 0)
resis_data$pfu_ml[my_rows] <- 1*resis_data$dilution[my_rows]
resis_data$bd[my_rows] <- T
#Then calculate all EOPs
for (i in 1:nrow(resis_data)) {
  my_sub <- subset(resis_data, resis_data$Date == resis_data$Date[i])
  resis_data$EOP[i] <- resis_data$pfu_ml[i]/
    mean(my_sub[my_sub$Treat == "Anc",]$pfu_ml)
}

#Add flag for old vs new approach
resis_data$approach <- NA
resis_data$approach[1:70] <- "old"
resis_data$approach[71:nrow(resis_data)] <- "new"

ggplot(resis_data[resis_data$Treat != "Anc", ], 
       aes(x = Treat, y = EOP, color = Pop,
           shape = bd, group = Pop)) +
  facet_grid(Proj~approach) +
  geom_point(position = position_dodge(width = 0.5),
             alpha = 0.7) +
  scale_y_continuous(trans = "log10") +
  theme_bw()

# #For local viewing
# ggplot(resis_data[resis_data$Treat != "A", ], aes(x = Treat, y = 1-EOP)) +
#   scale_x_discrete(name = "Treatment") +
#   ylab("Resistance to Phage") +
#   facet_grid(Proj~Pop, labeller = labeller(Proj = my_facet_labels)) +
#   geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.01, 
#                dotsize = 3) +
#   # geom_jitter(width = 0.2, height = 0, size = 2) +
#   theme_bw() + ggtitle("Population") +
#   theme(plot.title = element_text(size = 12, hjust = 0.5), 
#         axis.title = element_text(size = 12),
#         axis.text.x = element_text(color = "black", size = 12)) +
#   geom_hline(yintercept = 0, linetype = "dotted", size = 1.5)
# 
# #For poster
# png("resis_isols.png", width = 14, height = 9, units = "in", res = 300)
# ggplot(resis_data[resis_data$Treat != "A", ], aes(x = Treat, y = 1-EOP)) +
#   scale_x_discrete(name = "Treatment") +
#   ylab("Resistance to Phage") +
#   facet_grid(Proj~Pop, labeller = labeller(Proj = my_facet_labels)) +
#   geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.02, 
#                dotsize = 2) +
#   theme_bw() + ggtitle("Population") +
#   theme(plot.title = element_text(size = 24, hjust = 0.5), 
#         axis.title = element_text(size = 24),
#         axis.text = element_text(color = "black", size = 20),
#         strip.text = element_text(size = 24)) +
#   geom_hline(yintercept = 0, linetype = "dotted", size = 1.5)
# dev.off()
# 
# #resistance vs growth
# gc_resis_data <- merge(resis_data, gc_mppti)
# gc_resis_mppt <- group_by(gc_resis_data, Media, Proj, Pop, Treat)
# gc_resis_mppt <- summarize(gc_resis_mppt, avg_eop = mean(EOP),
#                            avg_gr = mean(gr_max_avg))
# 
# #Make plot of all isolates
# ggplot(gc_resis_data, aes(x = 1-EOP, y = gr_max_avg)) +
#   geom_point() + 
#   facet_grid(Media~., labeller = labeller(Media = my_facet_labels)) +
#   geom_smooth(method = "lm") +
#   labs(x = "Resistance", y = "Per Capita Growth Rate (/hour)")
# 
# summary(lm(gr_max_avg~Media*EOP, data = gc_resis_data))
# 
# #Make plot of all pops
# #For local viewing
# ggplot(gc_resis_mppt, aes(x = 1-avg_eop, y = avg_gr)) +
#   geom_point(size = 2, aes(pch = Treat)) + 
#   scale_shape_manual(values = c(3, 15, 16, 17), name = "Treatment",
#                      breaks = c("A", "C", "G", "L"),
#                      labels = c("WT", "Control", "Global", "Local")) +
#   facet_grid(Media~Proj, labeller = labeller(Media = my_facet_labels,
#                                              Proj = my_facet_labels)) +
#   geom_smooth(method = "lm") +
#   labs(x = "Resistance", y = "Average Per Capita Growth Rate (/hour)") +
#   theme_bw()
# 
# #For poster
# png("resis_gc_tradeoff.png", width = 10, height = 7, units = "in",
#     res = 300)
# ggplot(gc_resis_mppt, aes(x = 1-avg_eop, y = avg_gr)) +
#   geom_point(size = 4, aes(pch = Treat)) + 
#   scale_shape_manual(values = c(3, 15, 16, 17), name = "Treatment",
#                      breaks = c("A", "C", "G", "L"),
#                      labels = c("WT", "Control", "Global", "Local")) +
#   facet_grid(Proj~Media, labeller = labeller(Media = my_facet_labels,
#                                              Proj = my_facet_labels)) +
#   geom_smooth(method = "lm", se = F) +
#   labs(x = "Resistance", y = "Average Maximum Per Capita Growth Rate (/hr)") +
#   theme_bw() +
#   theme(axis.text = element_text(size = 16),
#         strip.text = element_text(size = 20),
#         axis.title = element_text(size = 20),
#         legend.title = element_text(size = 20),
#         legend.text = element_text(size = 16))
# dev.off()
# 
# summary(lm(avg_gr~avg_eop*Media, data = gc_resis_mppt))
# 


##Isolate growth curves ----

#Read data
gc_data <- read.csv("./Clean_Data/Isolate_growth_curves.csv",
                    header = T, stringsAsFactors = F)

#Make unique well identifiers
gc_data$uniq_well <- paste(gc_data$Date,
                           gc_data$Proj,
                           gc_data$Pop,
                           gc_data$Treat,
                           gc_data$Isol,
                           gc_data$Rep_Well,
                           gc_data$Media,
                           sep = "_")

#reorder
gc_data <- gc_data[order(gc_data$uniq_well, gc_data$Time_s), ]

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
  if (!is.numeric(time)) {
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

#Smooth data
gc_data$sm_loess <- NA
for (my_well in unique(gc_data$uniq_well)) {
  my_rows <- which(gc_data$uniq_well == my_well)
  #Smooth with loess
  gc_data$sm_loess[my_rows] <- loess(cfu_ml ~ Time_s, span = 0.4,
                                             data = gc_data[my_rows, ])$fitted
}

#Calculate growth per hour from loess curve
gc_data$deriv_sm_loess <- calc_deriv(gc_data$sm_loess,
                                     subset_by = gc_data$uniq_well,
                                     time = gc_data$Time_s,
                                     time_normalize = 3600)

#Calculate per capita growth per hour from loess curve
gc_data$percap_deriv_sm_loess <- calc_deriv(gc_data$sm_loess,
                                            percapita = TRUE,
                                            subset_by = gc_data$uniq_well,
                                            time = gc_data$Time_s,
                                            time_normalize = 3600)
  

#Calculate per capita growth per hour from original curve
# (for visual comparison)
# gc_data$percap_deriv_cfu <- calc_deriv(gc_data$cfu_ml,
#                                        percapita = TRUE,
#                                        subset_by = gc_data$uniq_well,
#                                        time = gc_data$Time_s,
#                                        time_normalize = 3600)

#View samples of original & smoothed curves
# as well as derivatives (per cap & not) of both orig and smoothed curves
for (my_well in sample(unique(gc_data$uniq_well), 3)) {
  my_rows <- which(gc_data$uniq_well == my_well)
  
  print(cowplot::plot_grid(
    ggplot(data = gc_data[my_rows, ],
           aes(x = Time_s, y = cfu_ml)) +
      geom_line(color = "red", lwd = 1, alpha = 0.5) +
      geom_line(aes(x = Time_s, y = sm_loess),
                color = "blue", lwd = 1, alpha = 0.5) +
      ggtitle(gc_data[my_rows[1], "uniq_well"]) +
      #geom_hline(yintercept = 383404890) +
      NULL,
    ggplot(data = gc_data[my_rows, ],
           aes(x = Time_s, y = deriv_sm_loess)) +
      geom_line(color = "blue") +
      NULL,
    ggplot(data = gc_data[my_rows, ],
           aes(x = Time_s, y = percap_deriv_sm_loess)) +
      geom_line(color = "blue") +
      # geom_line(aes(x = Time_s, y = percap_deriv_cfu),
      #           color = "red") +
      NULL,
    ncol = 1, align = "v"))
}

find_local_extrema <- function(values, 
                               return_maxima = TRUE,
                               return_minima = TRUE,
                               width_limit = NULL,
                               height_limit = NULL,
                               remove_endpoints = TRUE,
                               na.rm = FALSE) {
  #Takes a vector of values and returns a vector of the indices
  # of all local value extrema (by default, returns both maxima and minima)
  # To only return maxima or minima, change return_maxima/return_minima to FALSE

  #width_limit and/or height_limit must be provided
  #Width is how wide the window will be to look for a maxima/minima
  # Narrower width will be more sensitive to narrow local maxima/minima
  # Wider width will be less sensitive to narrow local maxima/minima
  #Height is how high or low a single step is allowed to take
  # e.g. a maxima-finding function will not pass a valley deeper
  # than height_limit
  #Note that this also limits approaches to extrema, so if set too small
  # function may converge on non-peaks
  #If both width_limit and height_limit are provided, steps are limited
  # conservatively (a single step must meet both criteria)
  
  #This function is designed to be compatible with dplyr::group_by and summarize
  
  #Check inputs
  if (!return_maxima & !return_minima) {
    stop("Both return_maxima and return_minima are FALSE, at least one must be TRUE")
  }
  if (is.null(width_limit) & is.null(height_limit)) {
    stop("Either width_limit or height_limit must be provided")
  }
  if (!is.null(width_limit)) {
    if (width_limit%%2 == 0) {
      warning("width_limit must be odd, will use ", width_limit-1, " as width_limit")
      width_limit <- width_limit - 1
    }
  }
  if (is.null(width_limit) & !is.null(height_limit)) {
    warning("height_limit alone tends to be sensitive to height_limit parameter, use with caution")
  }
  if (na.rm == TRUE & sum(is.na(values)) > 0) {
    if (!all(is.na(values[(1+length(values)-sum(is.na(values))):length(values)]))) {
      warning("Removing NAs found within values vector, returned indices will refer to non-NA values")
      print(values)
    }
    values <- values[!is.na(values)]
  } else if(any(is.na(values))) {
    stop("Some provided values are NA and na.rm = FALSE")
  }
  
  #Define sub-function to find limits of the window
  get_window_limits <- function(cnt_pos,
                                width_limit = NULL,
                                height_limit = NULL,
                                looking_for = c("minima", "maxima"),
                                values = NULL) {
    #Check inputs
    if (length(looking_for) > 1) {stop("looking_for must be specified")}
    if (!is.null(height_limit) & is.null(values)) {
        stop("height_limit is specified, but no values are provided")
    }
    if (is.null(width_limit) & is.null(height_limit)) {
      stop("Either width_limit or height_limit must be provided")
    }
    
    #Define window limits
    window_start <- c(NA, NA)
    if (!is.null(width_limit)) { #using width limit
      window_start[1] <- max(c(1, cnt_pos-floor(width_limit/2)))
    }
    if (!is.null(height_limit)) { #using height limig
      #For startpoint height, we want the latest point that is
      #behind of our current point and
      #either:
      # below current height - height limit
      # or above current height + height limit
      #Then we move one place forward 
      # (so it's the last value w/in height limit)
        window_start[2] <- max(c(1,
                                 1+which(1:length(values) < cnt_pos &
                                         (values >= (values[cnt_pos] + height_limit) |
                                            values <= (values[cnt_pos] - height_limit)))))
        #Make sure we're going at least 1 point backwards
        if(window_start[2] >= cnt_pos) {window_start[2] <- cnt_pos-1}
    }
    window_end <- c(NA, NA)
    if (!is.null(width_limit)) { #using width limit
      window_end[1] <- min(c(length(values), cnt_pos+floor(width_limit/2)))
    }
    if (!is.null(height_limit)) { #using height limit
        #For endpoint height, we want the earliest point that is
        #forward of our current point and
        #either:
        # below current height - height limit
        # or above current height + height limit
        #Then we move one place back 
        # (so it's the last value w/in height limit)
        window_end[2] <- min(c(length(values),
                               -1+which(1:length(values) > cnt_pos & #not backwards
                                    (values <= (values[cnt_pos] - height_limit) |
                                      values >= (values[cnt_pos] + height_limit)))))
        #Make sure we're going at least one point forwards
        if (window_end[2] <= cnt_pos) {window_end[2] <- cnt_pos+1}
    }
    return(c(max(window_start, na.rm = T), min(window_end, na.rm = T)))
  }
  
  find_next_extrema <- function(cnt_pos, values,
                        width_limit = NULL,
                        height_limit = NULL,
                        looking_for = c("minima", "maxima")) {
    if (cnt_pos == length(values)) {best_pos <- cnt_pos-1
    } else {best_pos <- cnt_pos+1}
    
    #Save the starting position so we never go backwards
    start_pos <- cnt_pos
    
    ##Looking for next maxima
    if(looking_for == "maxima") {
      while (cnt_pos != best_pos) {
        #Move the previous best pointer to current pointer location
        best_pos <- cnt_pos
        #Get next window limits
        window_lims <- get_window_limits(cnt_pos = cnt_pos,
                                         width_limit = width_limit,
                                         height_limit = height_limit,
                                         looking_for = "maxima",
                                         values = values)
        #Make sure we're not going backwards
        window_lims <- c(max(start_pos, window_lims[1]),
                         max(start_pos, window_lims[2]))
        #Then move current pointer to highest point within window
        # (making sure not to check non-integer indices, or indices below 1 or
        #  higher than the length of the vector)
        cnt_pos <- window_lims[1]-1+which.max(values[window_lims[1]:window_lims[2]])
      }
    ##Looking for next minima
    } else if (looking_for == "minima") {
      while (cnt_pos != best_pos) {
        #Move the previous best pointer to current pointer location
        best_pos <- cnt_pos
        #Get next window limits
        window_lims <- get_window_limits(cnt_pos = cnt_pos,
                                         width_limit = width_limit,
                                         height_limit = height_limit,
                                         looking_for = "minima",
                                         values = values)
        #Make sure we're not going backwards
        window_lims <- c(max(start_pos, window_lims[1]),
                         max(start_pos, window_lims[2]))
        #Then move current pointer to lowest point within window
        # (making sure not to check non-integer indices, or indices below 1 or
        #  higher than the length of the vector)
        cnt_pos <- window_lims[1]-1+which.min(values[window_lims[1]:window_lims[2]])
      }
    }
    return(best_pos)
  }
  
  cnt_pos <- 1
  ##Find first maxima
  maxima_list <- c(find_next_extrema(cnt_pos, values,
                                     width_limit = width_limit,
                                     height_limit = height_limit,
                                     looking_for = "maxima"))
  ##Find first minima
  minima_list <- c(find_next_extrema(cnt_pos, values,
                                     width_limit = width_limit,
                                     height_limit = height_limit,
                                     looking_for = "minima"))
  
  ##Check for next extrema until...
  while (TRUE) {
    #we're finding repeats
    if (any(duplicated(c(minima_list, maxima_list)))) {break}
    #or we hit the end of the values
    if (length(values) %in% c(maxima_list, minima_list)) {
      break
    }
    #Since maxima & minima must alternate, always start with furthest one 
    # we've found so far
    cnt_pos <- max(c(minima_list, maxima_list))
    #we're looking for a maxima next
    if (cnt_pos %in% minima_list) {
      maxima_list <- c(maxima_list,
                       find_next_extrema(cnt_pos, values,
                                         width_limit = width_limit,
                                         height_limit = height_limit,
                                         looking_for = "maxima"))
      #we're looking for a minima next
    } else if (cnt_pos %in% maxima_list) {
      minima_list <- c(minima_list,
                       find_next_extrema(cnt_pos, values,
                                         width_limit = width_limit,
                                         height_limit = height_limit,
                                         looking_for = "minima"))
    }
  }
  
  #Combine maxima & minima values & remove duplicates
  output <- c()
  if (return_maxima) {output <- c(output, maxima_list)}
  if (return_minima) {output <- c(output, minima_list)}
  #If remove endpoints is true, remove first or last values from return
  if (remove_endpoints) {
    if (1 %in% output) {output <- output[-which(output == 1)]}
    if (length(values) %in% output) {
      output <- output[-which(output == length(values))]}
  }
  #Remove duplicates
  output <- unique(output)
  #Order
  output <- output[order(output)]
  
  return(output)
}

#Group data by unique wells
gc_data <- group_by(gc_data, Date, Proj, Pop, Treat, Isol, Rep_Well, Media,
                    uniq_well)

# at max percap growth rate, find:
#     per capita growth rate (r)
#     density
#     time since some low density (lag time)
# at first local minima of non-percap growth rate after max growth rate 
#   (pseudo carrying capacity)
#   at pseudo K, find:
#     density (carrying capacity)
#     time since max percap growth rate
#     time since low density

#Define the window widths (sensitivity) for local extrema search
first_min_time <- 7200
max_percap_time <- 21600
max_gr_rate_time <- 10800
pseudo_K_time <- 7200

#Summarize data
gc_summarized <- summarize(gc_data,
  #Find the first minima in total density
   first_min_index = find_local_extrema(sm_loess,
                                        return_maxima = FALSE,
                                        width_limit = (first_min_time/
                                                         (Time_s[2]-Time_s[1])) + 1,
                                        na.rm = T,
                                        remove_endpoints = FALSE)[1],
   first_min = sm_loess[first_min_index],
   first_min_time = Time_s[first_min_index],
  #find peaks in per capita growth rate
  # but for one well we need a different criteria (see end)
  max_percap_index = ifelse(all(uniq_well != "2017-E_7x_B_C_E_1_50"),
    #first find all peaks
   find_local_extrema(percap_deriv_sm_loess,
                                        return_minima = FALSE,
                                        width_limit = (max_percap_time/
                                                         (Time_s[2]-Time_s[1])) + 1,
                                        na.rm = T,
                                        remove_endpoints = F)[
     #But save/use the first one that follows the minimum density
     match(TRUE, find_local_extrema(percap_deriv_sm_loess,
                                    return_minima = FALSE,
                                    width_limit = (max_percap_time/
                                                     (Time_s[2]-Time_s[1])) + 1,
                                    na.rm = T,
                                    remove_endpoints = F) >= first_min_index)],
   #if we're in that one well, use the second peak
  find_local_extrema(percap_deriv_sm_loess,
                                        return_minima = FALSE,
                                        width_limit = (max_percap_time/
                                                         (Time_s[2]-Time_s[1])) + 1,
                                        na.rm = T,
                                        remove_endpoints = F)[2]),
  max_percap_gr_rate = percap_deriv_sm_loess[max_percap_index],
  max_percap_gr_rate_time = Time_s[max_percap_index],
  max_percap_gr_rate_dens = sm_loess[max_percap_index],
  max_percap_gr_rate_timesincemin = max_percap_gr_rate_time - first_min_time,
  #find the first peak in total growth rate (slope of total density)
  max_grow_rate_index = find_local_extrema(deriv_sm_loess,
                                           return_minima = FALSE,
                                           width_limit = (max_gr_rate_time/
                                                            (Time_s[2]-Time_s[1])) + 1,
                                           na.rm = T,
                                           remove_endpoints = F)[1],
  #find the local minimas in total grow rate (slope of total density)
  # (which is the point when the diauxic shift occurs
  # but for one well we need a different criteria (see end)
  pseudo_K_index = ifelse(all(uniq_well != "2017-E_7x_B_C_E_1_50"),
    find_local_extrema(deriv_sm_loess,
                                      return_maxima = FALSE,
                                      width_limit = (pseudo_K_time/
                                                       (Time_s[2]-Time_s[1])) + 1,
                                      na.rm = T,
                       remove_endpoints = F)[
      #use the first one that follows the first peak in total gr rate
      # (note endpoints will be included, so when no local minima exists
      # the endpoint will simply be used
    match(TRUE, find_local_extrema(deriv_sm_loess,
                                      return_maxima = FALSE,
                                      width_limit = (pseudo_K_time/
                                                       (Time_s[2]-Time_s[1])) + 1,
                                      na.rm = T,
                                   remove_endpoints = F) > max_grow_rate_index)],
    #if we're in that one well, use the third minima
    find_local_extrema(deriv_sm_loess,
                       return_maxima = FALSE,
                       width_limit = (pseudo_K_time/
                                        (Time_s[2]-Time_s[1])) + 1,
                       na.rm = T,
                       remove_endpoints = F)[3]),
  pseudo_K = sm_loess[pseudo_K_index],
  pseudo_K_time = Time_s[pseudo_K_index],
  pseudo_K_deriv = deriv_sm_loess[pseudo_K_index],
  pseudo_K_timesincemin = pseudo_K_time - first_min_time,
  pseudo_K_timesince_maxpercap = pseudo_K_time - max_percap_gr_rate_time
)

#Change to data frame for cleanliness
gc_summarized <- as.data.frame(gc_summarized)

#Make output plots for problematic wells
if (F) {
  wells_check <- c("2017-B_7x_C_L_B_1_50", #isn't working
                   "2017-A_7x_Anc_Anc_Anc_1_50",
                   "2017-A_7x_B_L_A_1_100",
                   "2017-B_7x_C_L_B_1_50",
                   "2017-C_7x_B_C_C_1_100",
                   "2017-C_7x_B_C_C_2_100",
                   "2017-C_7x_C_C_C_1_50",
                   "2017-C_7x_C_L_C_2_100",
                   "2017-C_7x_D_G_C_1_50",
                   "2017-C_7x_D_G_C_1_100",
                   "2017-C_7x_D_G_C_2_50",
                   "2017-C_7x_D_G_C_2_100",
                   "2017-C_7x_E_G_C_1_100",
                   "2017-C_7x_E_G_C_2_100",
                   "2017-E_7x_B_C_E_1_50",
                   "2017-E_7x_C_C_E_1_100",
                   "2017-E_7x_C_C_E_2_100",
                   "2019-09-10_125_B_C_A_1_25-50",
                   "2019-09-10_125_B_C_A_2_25-50",
                   "2019-09-10_125_B_G_A_1_25-50",
                   "2019-09-12_125_B_C_D_1_25-50",
                   "2019-09-12_125_B_C_D_2_25-50",
                   "2019-09-12_125_B_G_D_1_25-50",
                   "2019-09-12_125_D_L_D_1_50-100",
                   "2019-09-12_125_D_L_D_2_50-100",
                   "2019-09-13_125_B_C_E_1_50-100",
                   "2019-09-13_125_B_C_E_2_25-50",
                   "2019-09-13_125_C_C_E_1_50-100",
                   "2019-09-13_125_D_L_E_1_50-100"
  )
  for (my_well in wells_check) {
    tiff(filename = paste("./Challenging_growth_curve_plots/", my_well, ".tiff", sep = ""),
         width = 5, height = 10, units = "in", res = 300)
    my_rows <- which(gc_data$uniq_well == my_well)
    print(cowplot::plot_grid(
      ggplot(data = gc_data[my_rows, ],
             aes(x = Time_s, y = cfu_ml)) +
        geom_line(color = "red", lwd = 1, alpha = 0.5) +
        geom_line(aes(x = Time_s, y = sm_loess),
                  color = "blue", lwd = 1, alpha = 0.5) +
        ggtitle(gc_data[my_rows[1], "uniq_well"]) +
        #Add point for first minima
        geom_point(data = gc_summarized[gc_summarized$uniq_well == my_well, ],
                   aes(x = first_min_time, y = first_min),
                   color = "green", size = 3) +
        NULL,
      ggplot(data = gc_data[my_rows, ],
             aes(x = Time_s, y = deriv_sm_loess)) +
        geom_line(color = "blue") +
        #Add point for pseudo K
        geom_point(data = gc_summarized[gc_summarized$uniq_well == my_well, ],
                   aes(x = pseudo_K_time, y = pseudo_K_deriv),
                   color = "green", size = 3) +
        #Add point for pseudo K2
        # geom_point(data = gc_summarized[gc_summarized$uniq_well == my_well, ],
        #            aes(x = pseudo_K_time2, y = pseudo_K_deriv2),
        #            color = "dark green", size = 2) +
        NULL,
      ggplot(data = gc_data[my_rows, ],
             aes(x = Time_s, y = percap_deriv_sm_loess)) +
        geom_line(color = "blue") +
        #Add point for max growth rate
        geom_point(data = gc_summarized[gc_summarized$uniq_well == my_well, ],
                   aes(x = max_percap_gr_rate_time, y = max_percap_gr_rate),
                   color = "green", size = 3) +
        NULL,
      ncol = 1, align = "v"))
    dev.off()
  }
}

#Make output plots for all wells
if (FALSE) {
  for (my_well in unique(gc_data$uniq_well)) {
    tiff(filename = paste("./Growth_curve_plots/", my_well, ".tiff", sep = ""),
         width = 5, height = 10, units = "in", res = 300)
    my_rows <- which(gc_data$uniq_well == my_well)
    print(cowplot::plot_grid(
      ggplot(data = gc_data[my_rows, ],
             aes(x = Time_s, y = cfu_ml)) +
        geom_line(color = "red", lwd = 1, alpha = 0.5) +
        geom_line(aes(x = Time_s, y = sm_loess),
                  color = "blue", lwd = 1, alpha = 0.5) +
        ggtitle(gc_data[my_rows[1], "uniq_well"]) +
        #Add point for first minima
        geom_point(data = gc_summarized[gc_summarized$uniq_well == my_well, ],
                   aes(x = first_min_time, y = first_min),
                   color = "green", size = 3) +
        NULL,
      ggplot(data = gc_data[my_rows, ],
             aes(x = Time_s, y = deriv_sm_loess)) +
        geom_line(color = "blue") +
        #Add point for pseudo K
        geom_point(data = gc_summarized[gc_summarized$uniq_well == my_well, ],
                   aes(x = pseudo_K_time, y = pseudo_K_deriv),
                   color = "green", size = 3) +
        NULL,
      ggplot(data = gc_data[my_rows, ],
             aes(x = Time_s, y = percap_deriv_sm_loess)) +
        geom_line(color = "blue") +
        #Add point for max growth rate
        geom_point(data = gc_summarized[gc_summarized$uniq_well == my_well, ],
                   aes(x = max_percap_gr_rate_time, y = max_percap_gr_rate),
                   color = "green", size = 3) +
        NULL,
      ncol = 1, align = "v"))
    dev.off()
  }
}

#Combine replicate wells then run PCA ----

#Summarize replicate wells
gc_summarized <- group_by(gc_summarized, Date, Proj, Pop, Treat,
                          Isol, Media)
gc_sum_noreps <- summarize(gc_summarized,
                           n = n(),
                           first_min_avg = mean(first_min),
                           first_min_sd = sd(first_min),
                           first_min_time_avg = mean(first_min_time),
                           first_min_time_sd = sd(first_min_time),
                           max_percap_gr_rate_avg = mean(max_percap_gr_rate),
                           max_percap_gr_rate_sd = sd(max_percap_gr_rate),
                           max_percap_gr_rate_time_avg = mean(max_percap_gr_rate_time),
                           max_percap_gr_rate_time_sd = sd(max_percap_gr_rate_time),
                           max_percap_gr_rate_dens_avg = mean(max_percap_gr_rate_dens),
                           max_percap_gr_rate_dens_sd = sd(max_percap_gr_rate_dens),
                           max_percap_gr_rate_timesincemin_avg = mean(max_percap_gr_rate_timesincemin),
                           max_percap_gr_rate_timesincemin_sd = sd(max_percap_gr_rate_timesincemin),
                           pseudo_K_avg = mean(pseudo_K),
                           pseudo_K_sd = sd(pseudo_K),
                           pseudo_K_time_avg = mean(pseudo_K_time),
                           pseudo_K_time_sd = sd(pseudo_K_time),
                           pseudo_K_timesincemin_avg = mean(pseudo_K_timesincemin),
                           pseudo_K_timesincemin_sd = sd(pseudo_K_timesincemin),
                           pseudo_K_timesince_maxpercap_avg = mean(pseudo_K_timesince_maxpercap),
                           pseudo_K_timesince_maxpercap_sd = sd(pseudo_K_timesince_maxpercap)
)
gc_sum_noreps <- as.data.frame(gc_sum_noreps)

#Take a look at the standard deviations between replicate wells
if (F) {
  for (var in c("first_min_", "first_min_time_", 
                "max_percap_gr_rate_", "max_percap_gr_rate_time_", 
                "max_percap_gr_rate_dens_", 
                "max_percap_gr_rate_timesincemin_",
                "pseudo_K_", "pseudo_K_time_", 
                "pseudo_K_timesincemin_", 
                "pseudo_K_timesince_maxpercap_")) {
    my_sd <- gc_sum_noreps[, paste(var, "sd", sep = "")]
    my_avg <- mean(gc_sum_noreps[, paste(var, "avg", sep = "")])
    hist(my_sd, main = var, 
         xlim = c(min(my_sd, my_avg, na.rm = T), max(my_sd, my_avg, na.rm = T)))
    abline(v = my_avg, col = "red", lwd = 2)
  }
}
if (F) {
  for (var in c("first_min_", "first_min_time_", 
                "max_percap_gr_rate_", "max_percap_gr_rate_time_", 
                "max_percap_gr_rate_dens_", 
                "max_percap_gr_rate_timesincemin_",
                "pseudo_K_", "pseudo_K_time_", 
                "pseudo_K_timesincemin_", 
                "pseudo_K_timesince_maxpercap_")) {
    my_sd <- gc_sum_noreps[, paste(var, "sd", sep = "")]
    my_avg <- gc_sum_noreps[, paste(var, "avg", sep = "")]
    hist(my_sd/my_avg, main = var, xlim = c(0, max(my_sd/my_avg, 1, na.rm = T)))
    abline(v = 1, col = "red", lwd = 2)
  }
}

#Generally, sd's between reps are small relative to the values themselves

#Cast measurements in different medias
# into different columns
gc_sum_noreps$Media[gc_sum_noreps$Media == "50"] <- "Orig"
gc_sum_noreps$Media[gc_sum_noreps$Media == "100"] <- "Rich"
gc_sum_noreps$Media[gc_sum_noreps$Media == "25-50"] <- "Orig"
gc_sum_noreps$Media[gc_sum_noreps$Media == "50-100"] <- "Rich"

gc_sum_noreps <- as.data.table(gc_sum_noreps)
gc_sum_noreps_wide <- data.table::dcast(gc_sum_noreps,
                           Date+Proj+Pop+Treat+Isol ~ Media,
                           value.var = c("first_min_avg", 
                                         "first_min_time_avg", 
                                         "max_percap_gr_rate_avg", 
                                         "max_percap_gr_rate_time_avg", 
                                         "max_percap_gr_rate_dens_avg", 
                                         "max_percap_gr_rate_timesincemin_avg", 
                                         "pseudo_K_avg", 
                                         "pseudo_K_time_avg", 
                                         "pseudo_K_timesincemin_avg", 
                                         "pseudo_K_timesince_maxpercap_avg"))
gc_sum_noreps_wide <- as.data.frame(gc_sum_noreps_wide)

#Check for univariate normality
if (F) {
  for (proj in unique(gc_sum_noreps_wide$Proj)) {
    for (var in c("first_min_avg_Orig", 
                  #"first_min_time_avg_Orig",
                  "max_percap_gr_rate_avg_Orig",
                  #"max_percap_gr_rate_time_avg_Orig",
                  "max_percap_gr_rate_dens_avg_Orig",
                  "max_percap_gr_rate_timesincemin_avg_Orig",
                  "pseudo_K_avg_Orig",
                  #"pseudo_K_time_avg_Orig",
                  "pseudo_K_timesincemin_avg_Orig",
                  #"pseudo_K_timesince_maxpercap_avg_Orig",
                  "first_min_avg_Rich",
                  #"first_min_time_avg_Rich",
                  "max_percap_gr_rate_avg_Rich",
                  #"max_percap_gr_rate_time_avg_Rich",
                  "max_percap_gr_rate_dens_avg_Rich",
                  "max_percap_gr_rate_timesincemin_avg_Rich",
                  "pseudo_K_avg_Rich",
                  #"pseudo_K_time_avg_Rich",
                  "pseudo_K_timesincemin_avg_Rich"
                  #"pseudo_K_timesince_maxpercap_avg_Rich"
                  )) {
      hist(as.numeric(gc_sum_noreps_wide[gc_sum_noreps_wide$Proj == proj, var]), 
           main = paste(proj, var))
      qqnorm(as.numeric(gc_sum_noreps_wide[gc_sum_noreps_wide$Proj == proj, var]), 
             main = paste(proj, var))
    }
  }
}

#TODO here

gc_sum_noreps$first_min_avg_log10 <- 
  log10(gc_sum_noreps$first_min_avg)
gc_sum_noreps$pseudo_K_avg_log10 <- 
  log10(gc_sum_noreps$pseudo_K_avg)
gc_sum_noreps$max_percap_gr_rate_dens_avg_log10 <-
  log10(gc_sum_noreps$max_percap_gr_rate_dens_avg)

if (F) {
  for (var in c("first_min_avg_log10",
                "pseudo_K_avg_log10",
                "max_percap_gr_rate_dens_avg_log10")) {
    hist(as.numeric(gc_sum_noreps[, var]), main = var)
    qqnorm(as.numeric(gc_sum_noreps[, var]), main = var)
  }
}

#Now they look much better!

#Check for multivariate normality

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


CSQPlot(gc_sum_noreps[, c("first_min_avg_log10", 
                          #"first_min_time_avg", 
                          "max_percap_gr_rate_avg", 
                          #"max_percap_gr_rate_time_avg", 
                          "max_percap_gr_rate_dens_avg_log10", 
                          "max_percap_gr_rate_timesincemin_avg",
                          "pseudo_K_avg_log10", 
                          #"pseudo_K_time_avg", 
                          "pseudo_K_timesincemin_avg"
                          #"pseudo_K_timesince_maxpercap_avg"
                          )])
#Our data is still highly non-multivariate normal

#Get principal components
gc_noreps_prin_comp <- princomp(gc_sum_noreps[, c("first_min_avg_log10", 
                                                  #"first_min_time_avg", 
                                                  "max_percap_gr_rate_avg", 
                                                  #"max_percap_gr_rate_time_avg", 
                                                  "max_percap_gr_rate_dens_avg_log10", 
                                                  "max_percap_gr_rate_timesincemin_avg",
                                                  "pseudo_K_avg_log10", 
                                                  #"pseudo_K_time_avg", 
                                                  "pseudo_K_timesincemin_avg"
                                                  #"pseudo_K_timesince_maxpercap_avg"
                                              )],
                                cor = T,
                                scores = T)

#Print summary
print(summary(gc_noreps_prin_comp), digits = 2)
#PC1 - 40% of variance, PC2 is 30%

#From screeplot, keep PC's 1, 2 & 3?
screeplot(gc_noreps_prin_comp, type = "lines", main = "Movie PCA Scree Plot")

#Check loadings of Principal components
print(gc_noreps_prin_comp$loadings, digits = 2, cutoff = 0)

#PC1 - the 3 densities (first min, dens as max percap, pseudo K)
#PC2 - max percap gr rate vs max percap gr rate time since min
#   (with pseudo K in same direction as max percap gr rate)
#   AKA the faster you grow the higher your pseudo K and less time
#    from min to max gr rate and from min to pseudo K
#PC3 - pseudo K time since min

gc_sum_noreps <- cbind(gc_sum_noreps, gc_noreps_prin_comp$scores)

ggplot(data = gc_sum_noreps,
       aes(x = Comp.1, y = Comp.2, color = Treat, shape = Proj)) +
  geom_point() +
  facet_wrap(~Media)

ggplot(data = gc_sum_noreps,
       aes(x = Comp.1, y = Comp.3, color = Treat, shape = Proj)) +
  geom_point() +
  facet_wrap(~Media)

ggplot(data = gc_sum_noreps,
       aes(x = Comp.2, y = Comp.3, color = Treat, shape = Proj)) +
  geom_point() +
  facet_wrap(~Media)




#Making percap growth rate plots
my_facet_labels <- c("100" = "Rich Environment", 
                     "50" = "Adapted Environment",
                     "C" = "Control", "G" = "Global", "L" = "Local",
                     "A" = "WT", "1" = "Weak Phage", "2" = "Strong Phage")

#plot of all isols
gc_mppti$Media <- factor(gc_mppti$Media, levels = c(50, 100))
ggplot(gc_mppti, aes(x = Treat, y = gr_max_avg)) + 
  geom_jitter(width = 0.1, height = 0, size = 2) + 
  facet_grid(Media ~ Pop, labeller = labeller(Media = my_facet_labels)) +
  labs(x = "Treatment", y = "Per Capita Growth Rate (/hour)") +
  theme(axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11)) +
  theme_bw()

#plot of all pops
gc_mppt$Media <- factor(gc_mppt$Media, levels = c(50, 100))

#For local viewing
ggplot(gc_mppt, aes(x = Treat, y = avg_isols), 
       labeller = labeller(Treat = my_facet_labels)) + 
  geom_point(pch = 1, size = 3) +
  facet_grid(Proj~Media, 
             labeller = labeller(Media = my_facet_labels, Proj = my_facet_labels)) + 
  labs(x = "Treatment", y = "Maximum Per Capita Growth Rate (/hour)") + theme_bw() + 
  scale_x_discrete(labels = c("Ancestor", "Control", "Global", "Local"))

#For poster
png(filename = "growth_rate_pops.png", width = 10, height = 7,
    units = "in", res = 300)
ggplot(gc_mppt, aes(x = Treat, y = avg_isols), 
       labeller = labeller(Treat = my_facet_labels)) + 
  geom_jitter(width = 0.075, height = 0, pch = 16, size = 6) +
  facet_grid(Proj~Media, 
             labeller = labeller(Media = my_facet_labels, Proj = my_facet_labels)) + 
  labs(x = "Treatment", y = "Maximum Per Capita Growth Rate (/hr)") + theme_bw() + 
  scale_x_discrete(labels = c("Ancestor", "Control", "Global", "Local")) +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16),
        strip.text = element_text(size = 20))
dev.off()


## Isolate PCA with all variables ----
