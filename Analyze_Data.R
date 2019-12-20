##TODO: 
##      migration stats
##      growth curve analysis
##      growth curve stats

library("ggplot2")
library("dplyr")

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
  scale_color_hue(name = "Treatment", breaks = c("C", "G", "L"),
                  labels = c("Control", "Global", "Local"))

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
       aes(x = Treat, y = relative_area, color = Pop,
           group = Pop)) +
  geom_point(position = position_dodge(0.4)) +
  facet_grid(~Proj, labeller = labeller(Proj = my_facet_labels)) +
  theme_bw() + 
  labs(y = "Isolate Area of Growth Relative to Ancestor",
       x = "Treatment") +
  geom_hline(yintercept = 1, lty = 2) +
  scale_x_discrete(limits = c("C", "L", "G"),
                   labels = c("Control", "Local", "Global")) +
  NULL

##Isolate growth curves ----


##growth curve analysis
# library("minpack.lm")
# library("tidyr")
# library("lubridate")
# library("dplyr")

#Read data
gc_data <- read.csv("./Clean_Data/Isolate_growth_curves.csv",
                    header = T, stringsAsFactors = F)

#Make unique well identifiers
gc_data$uniq_well <- NA
for (i in 1:nrow(gc_data)) {
  gc_data$uniq_well[i] <- paste(gc_data[i, 1:6], collapse = "_")
}

#reorder
gc_data <- gc_data[order(gc_data$uniq_well, gc_data$Time_s), ]

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
  gc_data$sm_loess[my_rows] <- loess(cfu_ml ~ Time_s,span = 0.4,
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
gc_data$percap_deriv_cfu <- calc_deriv(gc_data$cfu_ml,
                                       percapita = TRUE,
                                       subset_by = gc_data$uniq_well,
                                       time = gc_data$Time_s,
                                       time_normalize = 3600)

#View samples of original & smoothed curves
# as well as derivatives (per cap & not) of both orig and smoothed curves
for (my_well in sample(unique(gc_data$uniq_well), 5)) {
  my_rows <- which(gc_data$uniq_well == my_well)
  
  print(cowplot::plot_grid(
    ggplot(data = gc_data[my_rows, ],
           aes(x = Time_s, y = cfu_ml)) +
      geom_line(color = "red", lwd = 1, alpha = 0.5) +
      geom_line(aes(x = Time_s, y = sm_loess),
                color = "blue", lwd = 1, alpha = 0.5) +
      ggtitle(gc_data[my_rows[1], "uniq_well"]) +
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

#Todo:
# find max growth rate (& density at max growth rate)
# find density at first local minima of non-percap growth rate after max growth rate
#   (pseudo carrying capacity)
# find time from some low density to peak per capita growth rate (pseudo lag time)
# find time from max growth rate to pseudo carrying capacity

find_local_extrema <- function(values, 
                               return_maxima = TRUE,
                               return_minima = FALSE,
                               width_limit = NULL,
                               height_limit = NULL,
                               subset = NULL) {
  #Takes a vector of values and returns a vector of the indices
  # of all local value extrema (by default, returns all extrema)
  # To only return maxima or minima, change return_maxima/return_minima to FALSE

  #Either width_limit or height_limit must be provided
  #Width is how wide the window should be to look for a peak
  # Narrower bandwidth will be more sensitive to narrow local maxima
  # Wider bandwidth will be less sensitive to narrow local maxima
  #Height is how deep of a valley the function will search from each local
  # maxima to find other local maxima

  #Subset can be a Boolean vector same length as density & time
  # with TRUE as values that should be included
  # or it can be a vector of the indices of density & time that should be used
  
  #This function is designed to be run with dplyr::group_by and summarize
  
  #Check inputs
  if (!return_maxima & !return_minima) {
    stop("Both return_maxima and return_minima are FALSE, at least one must be TRUE")
  }
  if (is.null(width_limit) & is.null(height_limit)) {
    stop("Either width_limit or height_limit must be provided")
  } else if (!is.null(width_limit) & !is.null(height_limit)) {
    stop("Both width_limit and height_limit are provided, currently only support for one at a time is implemented")
  }
  
  #Start finding peaks
  maxima_list <- c()
  minima_list <- c()
  
  #use width limit to find peaks
  if (!is.null(width_limit)) { 
    ##Check for first maxima
    cnt_pos <- 1
    best_pos <- cnt_pos+1 #arbitrary index != to cnt_pos
    while (cnt_pos != best_pos) {
      #Move the previous best pointer to current pointer location
      best_pos <- cnt_pos
      #Then move current pointer to highest point within window
      # (making sure not to check non-integer indices, or indices below 1 or
      #  higher than the length of the vector)
      cnt_pos <- which.max(values[max(c(1, cnt_pos-floor(width_limit/2))):
                              min(c(length(values), cnt_pos+floor(width_limit/2)))])
    }
    #add value to maxima_list
    maxima_list <- c(maxima_list, best_pos)
    
    ##Check for first minima
    cnt_pos <- 1
    best_pos <- cnt_pos+1 #arbitrary index != to cnt_pos
    while (cnt_pos != best_pos) {
      #Move the previous best pointer to current pointer location
      best_pos <- cnt_pos
      #Then move current pointer to highest point within window
      # (making sure not to check non-integer indices, or indices below 1 or
      #  higher than the length of the vector)
      cnt_pos <- which.min(values[max(c(1, cnt_pos-floor(width_limit/2))):
                                    min(c(length(values), cnt_pos+floor(width_limit/2)))])
    }
    #add value to minima_list
    minima_list <- c(minima_list, best_pos)
    
    ##Check for next extrema until we no longer find new ones
    while (TRUE) {
      #Since maxima & minima must alternate, always start with furthest one 
      # we've found so far
      cnt_pos <- max(c(minima_list, maxima_list))
      #Assign best_pos to arbitrary valid index
      best_pos <- if (cnt_pos != length(values)) {
        best_pos <- cnt_pos+1
      } else {best_pos <- cnt_pos-1}
      #we're looking for a maxima next
      if (cnt_pos %in% minima_list) { 
        while (best_pos != cnt_pos) {
          #Move the previous best pointer to current pointer location
          best_pos <- cnt_pos
          #Then move current pointer to highest point within window
          # (making sure not to check non-integer indices, or indices below 1 or
          #  higher than the length of the vector)
          cnt_pos <- which.max(values[max(c(1, cnt_pos-floor(width_limit/2))):
                                        min(c(length(values), cnt_pos+floor(width_limit/2)))])
        }
        #add value to maxima_list
        maxima_list <- c(maxima_list, best_pos)
      #we're looking for a minima next
      } else if (cnt_pos %in% maxima_list) { 
        while (best_pos != cnt_pos) {
          #Move the previous best pointer to current pointer location
          best_pos <- cnt_pos
          #Then move current pointer to highest point within window
          # (making sure not to check non-integer indices, or indices below 1 or
          #  higher than the length of the vector)
          cnt_pos <- which.min(values[max(c(1, cnt_pos-floor(width_limit/2))):
                                        min(c(length(values), cnt_pos+floor(width_limit/2)))])
        }
        #add value to minima_list
        minima_list <- c(minima_list, best_pos)
      }
      if (best_pos %in% minima_list | best_pos %in% maxima_list) {break}
    }
  #use height limit to find peaks
  } else if (!is.null(height_limit)) {
    ##Check for first maxima
    cnt_pos <- 1
    best_pos <- cnt_pos+1 #arbitrary index != to cnt_pos
    while (cnt_pos != best_pos) {
      #Move the previous best pointer to current pointer location
      best_pos <- cnt_pos
      #Determine start and end of window
      #Then move current pointer to highest point within window
      # (making sure not to check non-integer indices, or indices below 1 or
      #  higher than the length of the vector)
      cnt_pos <- which.max(values[max(c(1, cnt_pos-floor(width_limit/2))):
                                    min(c(length(values), cnt_pos+floor(width_limit/2)))])
    }
    #add value to maxima_list
    maxima_list <- c(maxima_list, best_pos)
    
    #Check for first minima
    
    #Check for next extrema on loop
  }
  
  
  
  
  
  #Check for first maxima
  while (cnt_pos != best_pos) {
    cnt_pos <- 
  
  
  
  #Walk through points to find first local maxima & first local minima
  #whichever is earlier, start from there and find the next opposite type
  #repeat, alternating searching for maxima or minima, saving indices along the way
  #at end, put the types they requested output into an ordered list
  
  #Check by width:
    #Order all points by height
    #Find highest point
    #Find highest point not within width
    #Find highest point not within widths of #1 or #2
    #Repeat until...
  
  #Check by height:
  # Order all points by height
  # Find highest point
  # Determine closest points to peak that fall below height limit
  # Find highest point outside of that width
  #   repeat finding width essentially where it doesn't fall below height limit
#   Repeat
  
  #Make sure to remove the first and last points, if they're included
  ##and check for duplicates
  
  return(list("density" = ans_density, "time" = ans_time))
}
                               
  
  analyze_curves <- function(od_data, time_data, bandwidth = 10, return) {
    #Takes vectors of the od_data and time_data
    #Bandwidth is how wide the window should be to look for a peak
    #Narrower bandwidth will get you an earlier local maxima
    #Wider bandwidth will get you a later more-global maxima
    #Designed to be run w/ group_by
    prev_max_pos <- 0
    cnt_max_pos <- bandwidth
    while (cnt_max_pos != prev_max_pos) {
      prev_max_pos <- cnt_max_pos
      search_start <- pmax(cnt_max_pos-bandwidth, 1) #start of the search window
      search_end <- pmin(search_start + 2*bandwidth, length(od_data)) #end of the search window
      cnt_max_pos <- which.max(od_data[search_start:search_end]) + search_start - 1
    }
    if (return == "max") {return(od_data[cnt_max_pos])
    } else if (return == "maxtime") {return(time_data[cnt_max_pos])}
  }


#Plot curves & peaks
view_peaks <- function(time_data, dens_data, well_contents, plt_point = TRUE,
                       peak_time = NULL, peak_dens = NULL, peak_well_contents,
                       numplots = 9, lwd = 1) {
  #Inputs:  time_data - vector of time values
  #         dens_data - vector of density values, matched with time_data
  #         well_contents - vector of well contents, matched with time_data
  #         peak_time - vector of peaks' time values
  #         peak_dens - vector of peaks' density values, matched with peak_time
  #         peak_well_contents - vectore of peaks' well contents, matched with peak_time
  require(ggplot2)
  for (start_group in seq(from = 1, to = length(unique(well_contents)),
                          by = numplots)) {
    my_wells <- unique(well_contents)[start_group:(start_group+numplots-1)]
    myplot <- ggplot(data = data.frame(
      "dens" = dens_data[well_contents %in% my_wells],
      "time" = time_data[well_contents %in% my_wells],
      "contents" = well_contents[well_contents %in% my_wells]),
      aes(x = time, y = dens)) +
      geom_line(lwd = lwd) +
      facet_wrap(~contents)
    if (plt_point) {
      myplot <- myplot + 
        geom_point(data = data.frame(
          "dens" = peak_dens[peak_well_contents %in% my_wells],
          "time" = peak_time[peak_well_contents %in% my_wells],
          "contents" = peak_well_contents[peak_well_contents %in% my_wells]),
          aes(x = time, y = dens),
          size = 3, pch = 13)
    }
    print(myplot)
  }
}

view_peaks(tidymerge$Timepoint, tidymerge$OD, tidymerge$Well,
           peak_time = outmerge$maxtime, peak_dens = outmerge$max, 
           peak_well_contents = outmerge$Well)

#Group data by indiv growth curves
grp_data1 <- dplyr::group_by(tidycurves[!is.na(tidycurves$smoothed), ], 
                             Well)

#Get OD peak height & time for each growth curve
out_data1 <- dplyr::summarize(grp_data1, 
                              max = analyze_curves(smoothed, Timepoint, 
                                                   bandwidth = 9, return = "max"),
                              maxtime = analyze_curves(smoothed, Timepoint, 
                                                       bandwidth = 9, 
                                                       return = "maxtime"))



#extract maximum percap growth rates for each uniq well
#using non-overlapping averaged smoothing
gc_sm$Time <- as.character(gc_sm$Time)
gc_mpptir <- group_by(gc_sm, Media, Proj, Pop, Treat, Isol, Rep_Well)
gc_mpptir <- summarize(gc_mpptir, max_pcgr = max(pcgr, na.rm = T))
gc_mppti <- group_by(gc_mpptir, Media, Proj, Pop, Treat, Isol)
gc_mppti <- summarize(gc_mppti, gr_max_avg = mean(max_pcgr),
                      dt_min_sd = sd(max_pcgr))
gc_mppt <- group_by(gc_mppti, Media, Proj, Pop, Treat)
gc_mppt <- summarize(gc_mppt, avg_isols = mean(gr_max_avg),
                     sd_isols = sd(gr_max_avg))

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


## Isolate PCA ----
