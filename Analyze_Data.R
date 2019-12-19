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

##Experimental evolution analysis ----
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

##Isolate migration analysis ----
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

##isolate growth curve analysis ----


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

calc_deriv <- function(CFU, percapita = FALSE,
                          subset_by = NULL, time = NULL,
                          time_normalize = NULL) {
  #Note! CFU values must be sorted sequentially into their unique sets already
  
  #Provided a vector of CFU values, this function returns (by default) the
  # difference between sequential values
  #if percapita = TRUE, the differences of cfu are divided by cfu
  #if subset_by is provided, it should be a vector (same length as cfu),
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
  ans <- c(CFU[2:length(CFU)]-CFU[1:(length(CFU)-1)])
  #Percapita (if specified)
  if (percapita) {
    ans <- ans/CFU[1:(length(CFU)-1)]
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
gc_data$percap_deriv_sm_loess <- NA
gc_data$percap_deriv_cfu <- NA
for (my_well in unique(gc_data$uniq_well)) {
  my_rows <- which(gc_data$uniq_well == my_well)
  #Smooth with loess
  gc_data$sm_loess[my_rows] <- predict(loess(cfu_ml ~ Time_s,
                                             span = 0.4,
                                             data = gc_data[my_rows, ]),
                                       gc_data[my_rows, ])
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
for (my_well in sample(unique(gc_data$uniq_well), 20)) {
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

##Isolate resistance analysis
resis_data_7x <- read.csv("74_75_76_Plaquing.csv", header = T, stringsAsFactors = F)

#split out project info & make unique reps
resis_data_7x <- cbind(resis_data_7x[, 1:2], "Pop" = NA, resis_data_7x[, 3:5])
for (i in 1:nrow(resis_data_7x)) {
  my.proj <- resis_data_7x$Proj[i]
  if (my.proj == "P1.1") {
    resis_data_7x$Pop[i] <- "F"
    resis_data_7x$Treat[i] <- "A" #for Ancestor
  } else if (my.proj == "74") {
    resis_data_7x$Pop[i] <- "A"
  } else {
    resis_data_7x$Proj[i] <- substr(my.proj, 1, 2)
    resis_data_7x$Pop[i] <- substr(my.proj, 3, 3)
  }
  if (resis_data_7x$Proj[i] == "75") {
    if (resis_data_7x$Pop[i] == "A") {resis_data_7x$Pop[i] <- "B"
    } else {resis_data_7x$Pop[i] <- "C"}
  } else if (resis_data_7x$Proj[i] == "76") {
    if (resis_data_7x$Pop[i] == "A") {resis_data_7x$Pop[i] <- "D"
    } else {resis_data_7x$Pop[i] <- "E"}
  }
}
resis_data_7x$Proj <- 1

#this is where merge of resis data will be
resis_data <- resis_data_7x

#calculate EOP for ea isol
resis_data$EOP <- NA
for (i in 1:nrow(resis_data)) {
  my_sub <- subset(resis_data, resis_data$Isol == resis_data$Isol[i])
  resis_data$EOP[i] <- resis_data$PFU[i]/my_sub[my_sub$Treat == "A",]$PFU
}

#For local viewing
ggplot(resis_data[resis_data$Treat != "A", ], aes(x = Treat, y = 1-EOP)) +
  scale_x_discrete(name = "Treatment") +
  ylab("Resistance to Phage") +
  facet_grid(Proj~Pop, labeller = labeller(Proj = my_facet_labels)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.01, 
               dotsize = 3) +
  # geom_jitter(width = 0.2, height = 0, size = 2) +
  theme_bw() + ggtitle("Population") +
  theme(plot.title = element_text(size = 12, hjust = 0.5), 
        axis.title = element_text(size = 12),
        axis.text.x = element_text(color = "black", size = 12)) +
  geom_hline(yintercept = 0, linetype = "dotted", size = 1.5)

#For poster
png("resis_isols.png", width = 14, height = 9, units = "in", res = 300)
ggplot(resis_data[resis_data$Treat != "A", ], aes(x = Treat, y = 1-EOP)) +
  scale_x_discrete(name = "Treatment") +
  ylab("Resistance to Phage") +
  facet_grid(Proj~Pop, labeller = labeller(Proj = my_facet_labels)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.02, 
               dotsize = 2) +
  theme_bw() + ggtitle("Population") +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.title = element_text(size = 24),
        axis.text = element_text(color = "black", size = 20),
        strip.text = element_text(size = 24)) +
  geom_hline(yintercept = 0, linetype = "dotted", size = 1.5)
dev.off()

#resistance vs growth
gc_resis_data <- merge(resis_data, gc_mppti)
gc_resis_mppt <- group_by(gc_resis_data, Media, Proj, Pop, Treat)
gc_resis_mppt <- summarize(gc_resis_mppt, avg_eop = mean(EOP),
                           avg_gr = mean(gr_max_avg))

#Make plot of all isolates
ggplot(gc_resis_data, aes(x = 1-EOP, y = gr_max_avg)) +
  geom_point() + 
  facet_grid(Media~., labeller = labeller(Media = my_facet_labels)) +
  geom_smooth(method = "lm") +
  labs(x = "Resistance", y = "Per Capita Growth Rate (/hour)")

summary(lm(gr_max_avg~Media*EOP, data = gc_resis_data))

#Make plot of all pops
#For local viewing
ggplot(gc_resis_mppt, aes(x = 1-avg_eop, y = avg_gr)) +
  geom_point(size = 2, aes(pch = Treat)) + 
  scale_shape_manual(values = c(3, 15, 16, 17), name = "Treatment",
                     breaks = c("A", "C", "G", "L"),
                     labels = c("WT", "Control", "Global", "Local")) +
  facet_grid(Media~Proj, labeller = labeller(Media = my_facet_labels,
                                             Proj = my_facet_labels)) +
  geom_smooth(method = "lm") +
  labs(x = "Resistance", y = "Average Per Capita Growth Rate (/hour)") +
  theme_bw()

#For poster
png("resis_gc_tradeoff.png", width = 10, height = 7, units = "in",
    res = 300)
ggplot(gc_resis_mppt, aes(x = 1-avg_eop, y = avg_gr)) +
  geom_point(size = 4, aes(pch = Treat)) + 
  scale_shape_manual(values = c(3, 15, 16, 17), name = "Treatment",
                     breaks = c("A", "C", "G", "L"),
                     labels = c("WT", "Control", "Global", "Local")) +
  facet_grid(Proj~Media, labeller = labeller(Media = my_facet_labels,
                                             Proj = my_facet_labels)) +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Resistance", y = "Average Maximum Per Capita Growth Rate (/hr)") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))
dev.off()

summary(lm(avg_gr~avg_eop*Media, data = gc_resis_mppt))


## 125 Migration ----
migration_125 <- read.csv("131_125_isol_migration.csv")
migration_125$area <- pi*migration_125$Width..cm.*migration_125$Height..cm./4
migration_125$relative_area <- NA
migration_125$date <- paste(migration_125$Year, migration_125$Month,
                            migration_125$Day, sep = "_")
for (date in unique(migration_125$date)) {
  migration_125$relative_area[which(migration_125$date == date)] <- 
    migration_125$area[which(migration_125$date == date)]/
    migration_125$area[migration_125$date == date & 
                         migration_125$Population == "Anc"]
}

ggplot(data = migration_125, aes(x = paste(Treatment, Population),
                                 y = relative_area)) +
  geom_point()

## Growth curves ----
#Remove first hour of datapoints
gc_data <- gc_data[hour(gc_data$Time) > 0, ]

#Add variables for getting unique subsets
gc_data$mpptir <- paste(gc_data$Media, gc_data$Proj, gc_data$Pop, gc_data$Treat, 
                        gc_data$Isol, gc_data$Rep_Well, sep = ".")
gc_data$mppti <- paste(gc_data$Media, gc_data$Proj, gc_data$Pop, gc_data$Treat, 
                       gc_data$Isol, sep = ".")

#smooth CFU data
smooth_data <- function(my_data, smooth_over, subset_by) {
  #data must be sorted sequentially before fed into function
  #my_data is a vector of the data to be smoothed
  #smooth over is how many sequential entries to average
  #the unique values of subset_by will be what is iterated over
  out_list <- rep(NA, length(my_data))
  cntr = 1
  for (my_uniq in unique(subset_by)) {
    my_sub <- subset(my_data, subset_by == my_uniq)
    out_list[cntr:(cntr+length(my_sub)-smooth_over)] <- 0
    for (i in 1:smooth_over) {
      out_list[(cntr):(cntr+length(my_sub)-smooth_over)] <-
        out_list[(cntr):(cntr+length(my_sub)-smooth_over)] + 
        my_sub[i:(length(my_sub)-smooth_over+i)]
    }  
    cntr <- cntr+length(my_sub)
  }
  out_list <- out_list/smooth_over
  return(out_list)
}

#smooth by averaging non-overlapping ranges
smooth_data_noverlap <- function(my_data, smooth_over, subset_by) {
  #data must be sorted sequentially before fed into function
  #my_data is a vector of the data to be smoothed
  #smooth over is how many sequential entries to average
  #the unique values of subset_by will be what is iterated over
  out_list <- rep(NA, length(my_data))
  cntr = 0
  for (my_uniq in unique(subset_by)) {
    my_sub <- subset(my_data, subset_by == my_uniq)
    for (i in seq(from = 1, to = length(my_sub), by = smooth_over)) {
      if (i+smooth_over < length(my_sub)) {
        out_list[cntr+i] <- mean(my_sub[i:(i+smooth_over-1)])
      }
    }
    cntr <- cntr+length(my_sub)
  }
  return(out_list)
}

gc_sm <- gc_data

gc_data$Smooth_CFU <- smooth_data(gc_data$CFU, 8, subset_by = gc_data$mpptir)
gc_sm$Smooth_noverlap <- smooth_data_noverlap(gc_data$CFU, 5, 
                                              subset_by = gc_data$mpptir)
gc_sm <- gc_sm[is.na(gc_sm$Smooth_noverlap) == F, ]

per_cap_grate <- function(CFU, time, subset_by) {
  ans <- c((CFU[2:length(CFU)]-CFU[1:(length(CFU)-1)])/
             (as.numeric(difftime(time[2:length(time)], time[1:(length(time)-1)], 
                                  units = "hours")) * 
                0.5 * (CFU[2:length(CFU)]+CFU[1:(length(CFU)-1)])), NA)
  ans[subset_by[2:length(subset_by)] != subset_by[1:(length(subset_by)-1)]] <- NA
  return(ans)
}

gc_data$pcgr <- per_cap_grate(gc_data$Smooth_CFU, gc_data$Time, gc_data$mpptir)
gc_sm$pcgr <- per_cap_grate(gc_sm$Smooth_noverlap, gc_sm$Time, gc_sm$mpptir)

# #Make plot of all smoothed growth rates
# for (i in seq(from = 1, to = length(unique(gc_data$mpptir)), by = 25)) {
#   my_sub <- gc_data[(gc_data$mpptir %in% unique(gc_data$mpptir)[i:(i+24)]), ]
#   print(ggplot(data = my_sub, aes(x = Time, y = pcgr)) + geom_line() +
#           facet_wrap(~mpptir))
# }
# 
# #Make plots of all smoothed nonoverlapping growth rates
# for (i in seq(from = 1, to = length(unique(gc_sm$mpptir)), by = 25)) {
#   my_sub <- gc_sm[(gc_sm$mpptir %in% unique(gc_sm$mpptir)[i:(i+24)]), ]
#   print(ggplot(data = my_sub, aes(x = Time, y = pcgr)) + geom_line() +
#           facet_wrap(~mpptir))
# }

# #Code for looking at an example well
# my_well <- "100.1.D.C.A.1"
# my_sub <- gc_data[gc_data$mpptir == my_well, ]
# ggplot(data = my_sub, aes(x = Time, y = CFU)) +
#   geom_line() + labs(y = "Colony Forming Units (CFU)")
# ggplot(data = my_sub, aes(x = Time, y = Smooth_CFU)) +
#   geom_line() + labs(y = "Colony Forming Units (CFU)")
# ggplot(data = my_sub, aes(x = Time, y = pcgr)) +
#   geom_line()
# 
# my_sub <- gc_sm[gc_sm$mpptir == my_well, ]
# ggplot(data = my_sub, aes(x = Time, y = CFU)) + geom_line()
# ggplot(data = my_sub, aes(x = Time, y = Smooth_noverlap)) +
#   geom_line()
# ggplot(data = my_sub, aes(x = Time, y = pcgr)) + geom_line()
# 
# 
# my_sub$diff <- c(my_sub$CFU[2:nrow(my_sub)]-my_sub$CFU[1:(nrow(my_sub)-1)], NA)
# my_sub$smdiff <- c(my_sub$Smooth_CFU[2:nrow(my_sub)]-my_sub$Smooth_CFU[1:(nrow(my_sub)-1)], NA)
# ggplot(data = my_sub, aes(x = Time, y = smdiff)) +
#   geom_line()

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


## Resistance ----
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

## PCA ----
