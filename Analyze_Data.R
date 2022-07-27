## Load packages and color scale ----
library(ggplot2)
library(dplyr)
library(ggh4x)
library(lme4)

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

#Global options (set these to TRUE & run script to have plots output to files)
make_curveplots <- FALSE
make_statplots <- FALSE

##Experimental evolution migration ----
exper_evol_migr <- read.csv("./Clean_Data/Experimental_evolution_growth.csv")

#Reorder so weak phage is first
exper_evol_migr$Proj <- factor(exper_evol_migr$Proj,
                               levels = c("7x", "125"))

#Drop points after T14
exper_evol_migr <- exper_evol_migr[exper_evol_migr$Timepoint <= 14, ]

#Based on Croze, Ottavio A., et al. "Migration of chemotactic bacteria in 
#soft agar: role of gel concentration." Biophysical journal 101.3 (2011): 525-534.
#radius increases ~linearly with time, so we'll normalize radius by time

#Calculate radius/hr
exper_evol_migr$radius_mm_hr <- 
  10*(exper_evol_migr$Width_cm+exper_evol_migr$Height_cm)/
  (2*exper_evol_migr$time_since_inoc)

#Summarize
exper_evol_migr <- group_by(exper_evol_migr, Proj, Treat, Timepoint)
exper_evol_summ <- summarize(exper_evol_migr,
                             pops_n = n(),
                             radius_mm_hr_mean = mean(radius_mm_hr),
                             radius_mm_hr_sd = sd(radius_mm_hr))

if (make_statplots) {
  my_facet_labels <- c("7x" = "Weak Phage", "125" = "Strong Phage")
  
  #Make plot of summarized data
  tiff("./Output_figures/Exper_evol_migr_nopops.tiff",
       width = 7, height = 4, units = "in", res = 300)
  print(ggplot(data = exper_evol_summ, aes(x = Timepoint, y = radius_mm_hr_mean,
                                           color = Treat)) +
          #geom_point(position = position_dodge(0.2)) + 
          geom_line(size = 2, position = position_dodge(0.2), alpha = 0.8) +
          facet_grid(Proj~., labeller = labeller(Proj = my_facet_labels)) +
          # geom_errorbar(aes(ymax = radius_mm_hr_mean+radius_mm_hr_sd, 
          #                   ymin = radius_mm_hr_mean-radius_mm_hr_sd),
          #               width=1, size = .7, position=position_dodge(0.2)) +
          labs(x = "Transfer", 
               y = expression(paste("Bacterial Growth (mm/hr)"))) + 
          scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                             labels = c("Control", "Local", "Global"),
                             values = my_cols[c(8, 2, 6)]) +
          scale_x_continuous(breaks = c(0, 7, 14)) +
          scale_y_continuous(limits = c(0, NA), breaks = c(0, 0.5, 1)) +
          theme_bw() +
          theme(axis.text.y = element_text(size = 12), 
                axis.text.x = element_text(size = 12),
                axis.title = element_text(size = 15),
                legend.text = element_text(size = 14), 
                legend.title = element_text(size = 15),
                strip.text = element_text(size = 15)) +
          NULL)
  dev.off()

  #Make plot with both summarized and non-summarized data
  tiff("./Output_figures/Exper_evol_migr.tiff",
       width = 7, height = 4, units = "in", res = 300)
  print(ggplot(data = exper_evol_migr,
         aes(x = Timepoint, y = radius_mm_hr, 
             group = paste(Treat, Pop),
             color = Treat)) +
    #  geom_point(size = 0.5, alpha = 0.5) +
    geom_line(alpha = 0.5, lwd = .4) +
    facet_grid(Proj~., labeller = labeller(Proj = my_facet_labels)) +
    geom_line(data = exper_evol_summ,
              aes(x = Timepoint, y = radius_mm_hr_mean, color = Treat,
                  group = Treat),
              size = 1.3) +
    labs(x = "Transfer", y = "Soft Agar Growth (mm/hr)") + 
    scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                       labels = c("Control", "Local", "Global"),
                       values = my_cols[c(8, 2, 6)]) +
    scale_x_continuous(breaks = c(0, 7, 14)) +
    ylim(0, NA) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 11), 
          axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 13), 
          legend.title = element_text(size = 14),
          strip.text = element_text(size = 14)) +
    NULL)
  dev.off()
}

##Isolate migration ----
isol_migration <- read.csv("./Clean_Data/Isolate_migration.csv")

#Reorder so weak phage is first
isol_migration$Proj <- factor(isol_migration$Proj,
                              levels = c("7x", "125"))

#Reorder treatments
isol_migration$Treat <- factor(isol_migration$Treat,
                               levels = c("Anc", "C", "L", "G"))

#Calculate radius/hr
isol_migration$radius_mm_hr <-
  10*(isol_migration$Width_cm+isol_migration$Height_cm)/
  (2*isol_migration$time_since_inoc)

#Calculate normalized "delta" values against same-day ancestor
# (note that first batch of isols had no same-day ancestor to normalize to)
ancestors <- isol_migration[isol_migration$Isol == "Anc", ]

isol_migration$radius_mm_hr_del <-
  isol_migration$radius_mm_hr - ancestors$radius_mm_hr[
    match(as.Date(isol_migration$end_timestamp),
          as.Date(ancestors$end_timestamp))]

#Plot data
my_facet_labels <- c("7x" = "Weak Phage", 
                     "125" = "Strong Phage",
                     "C" = "Control", "G" = "Global", "L" = "Local",
                     "Anc" = "Ancestor")

if (make_statplots) {
  #Raw unnormalized radius
  tiff("./Output_figures/Isol_migration.tiff",
       width = 6, height = 4, units = "in", res = 300)
  print(ggplot(isol_migration, 
         aes(x = Pop, y = radius_mm_hr, 
             color = Treat, fill = Treat)) +
    geom_point(position = position_dodge(0.5), alpha = 0.6,
               size = 2) +
    facet_grid(Proj ~ Treat, scales = "free",
                 labeller = labeller(Proj = my_facet_labels,
                                                  Treat = my_facet_labels)) +
    theme_bw() + 
    labs(y = "Soft Agar Growth (mm/hr)",
         x = "Population") +
    # scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
    #                    labels = c("Control", "Local", "Global"),
    #                    values = my_cols[c(8, 2, 6)]) +
    # scale_fill_manual(name = "Treatment", breaks = c("C", "L", "G"),
    #                   labels = c("Control", "Local", "Global"),
    #                   values = my_cols[c(8, 2, 6)]) +
    theme(legend.position = "none") +
    NULL)
  dev.off()
  
  #Delta radius
  tiff("./Output_figures/Isol_migration_del.tiff",
       width = 5, height = 4, units = "in", res = 300)
  print(ggplot(isol_migration[isol_migration$Isol != "Anc", ], 
         aes(x = Pop, y = radius_mm_hr_del, 
             color = Treat, fill = Treat)) +
    geom_point(position = position_dodge(0.5), alpha = 0.6,
               size = 2) +
    facet_grid(Proj~Treat, labeller = labeller(Proj = my_facet_labels,
                                                  Treat = my_facet_labels)) +
    theme_bw() + 
    labs(y = paste("\u0394", "Soft Agar Growth (mm/hr)", sep = ""),
         x = "Population") +
    geom_hline(yintercept = 0, lty = 2) +
    scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                       labels = c("Control", "Local", "Global"),
                       values = my_cols[c(8, 2, 6)]) +
    scale_fill_manual(name = "Treatment", breaks = c("C", "L", "G"),
                      labels = c("Control", "Local", "Global"),
                      values = my_cols[c(8, 2, 6)]) +
    theme(legend.position = "none") +
    NULL)
  dev.off()
}

#Summarize
isol_migr_sum <- 
  summarize(group_by(isol_migration[!is.na(isol_migration$radius_mm_hr_del), ],
                     Proj, Pop, Treat),
            radius_mm_hr_del_avg = mean(radius_mm_hr_del))
isol_migr_sum_isols <-
  summarize(group_by(isol_migration[!is.na(isol_migration$radius_mm_hr_del), ],
                     Proj, Pop, Treat, Isol),
            radius_mm_hr_del_avg = mean(radius_mm_hr_del))

#Make plot of only pop-level data
if (make_statplots) {
  tiff("./Output_figures/Isol_migration_del_pops.tiff",
       width = 4, height = 4, units = "in", res = 300)
  print(ggplot(isol_migr_sum[isol_migr_sum$Treat != "Anc", ], 
               aes(x = Treat, y = radius_mm_hr_del_avg, 
                   color = Treat, fill = Treat)) +
          geom_point(alpha = 0.6, size = 3) +
          facet_grid(Proj~., labeller = labeller(Proj = my_facet_labels)) +
          labs(y = paste("\u0394", "Soft Agar Growth (mm/hr)", sep = ""),
               x = "Treatment") +
          geom_hline(yintercept = 0, lty = 2) +
          scale_x_discrete(breaks = c("C", "L", "G"),
                           labels = c("Control", "Local", "Global")) +
          scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                             labels = c("Control", "Local", "Global"),
                             values = my_cols[c(8, 2, 6)]) +
          scale_fill_manual(name = "Treatment", breaks = c("C", "L", "G"),
                            labels = c("Control", "Local", "Global"),
                            values = my_cols[c(8, 2, 6)]) +
          theme_bw() + 
          theme(legend.position = "none",
                axis.text.y = element_text(size = 11), 
                axis.text.x = element_text(size = 11),
                axis.title.y = element_text(size = 12.5),
                axis.title.x = element_text(size = 14),
                legend.text = element_text(size = 13), 
                legend.title = element_text(size = 14),
                strip.text = element_text(size = 14)) +
          NULL)
  dev.off()
}


## Isolate resistance ----
resis_data <- read.csv("./Clean_Data/Isolate_resistance.csv",
                       stringsAsFactors = F)

#Remove data collected using older approach (which was superceded)
resis_data <- resis_data[grep("2017-[ABCDE]", resis_data$Date, invert = TRUE), ]

#Reorder projects
resis_data$Proj <- factor(resis_data$Proj,
                          levels = c("7x", "125"))

#Reorder treatments
resis_data$Treat <- factor(resis_data$Treat,
                           levels = c("Anc", "C", "L", "G"))

#calculate EOP for ea isol
resis_data$EOP <- NA
resis_data$bd <- F
#First, handle below-detection points
my_rows <- which(resis_data$PFU == 0)
resis_data$pfu_ml[my_rows] <- 1*resis_data$dilution[my_rows]
resis_data$bd[my_rows] <- T
#Then calculate all EOPs
for (i in 1:nrow(resis_data)) {
  my_sub <- subset(resis_data, resis_data$Date == resis_data$Date[i])
  resis_data$EOP[i] <- resis_data$pfu_ml[i]/
    mean(my_sub[my_sub$Treat == "Anc",]$pfu_ml)
}

#Calculate EOP limit & adjust values below limit
eop_limit <- max(resis_data$EOP[resis_data$bd])
resis_data$EOP[resis_data$EOP < eop_limit] <- eop_limit

#Calculate Ancestral EOP range
Anc_EOP <- data.frame(Proj = rep(c("7x", "125"), each = 3),
                      Treat = rep(c("C", "L", "G"), times = 2),
                      min_EOP = rep(c(
                        min(resis_data$EOP[resis_data$Treat == "Anc" &
                                             resis_data$Proj == "7x"]),
                        min(resis_data$EOP[resis_data$Treat == "Anc" &
                                             resis_data$Proj == "125"])),
                        each = 3),
                      max_EOP = rep(c(
                        max(resis_data$EOP[resis_data$Treat == "Anc" &
                                             resis_data$Proj == "7x"]),
                        max(resis_data$EOP[resis_data$Treat == "Anc" &
                                             resis_data$Proj == "125"])),
                        each = 3))
Anc_EOP$Proj <- factor(Anc_EOP$Proj, levels = c("7x", "125"))
Anc_EOP$Treat <- factor(Anc_EOP$Treat, levels = c("C", "L", "G"))

#Assign facet labels
my_facet_labels <- c("7x" = "Weak Phage", 
                     "125" = "Strong Phage",
                     "C" = "Control", "G" = "Global", "L" = "Local",
                     "A" = "WT")

if (make_statplots) {
  #Nice plot
  tiff("./Output_figures/Isol_resis.tiff", width = 5, height = 4,
       units = "in", res = 300)
  print(ggplot(resis_data[resis_data$Treat != "Anc", ],
         aes(x = Pop, y = EOP, 
             color = Treat, fill = Treat,
             shape = bd)) +
    facet_grid(Proj~Treat, labeller = labeller(Proj = my_facet_labels,
                                                  Treat = my_facet_labels)) +
    geom_rect(data = Anc_EOP,
              aes(x = NULL, y = NULL, color = NULL, fill = NULL, shape = NULL,
                  xmin = 0, xmax = 6,
                  ymin = min_EOP, ymax = max_EOP),
              alpha = 0.3) +
    geom_point(aes(size = bd, alpha = bd)) +
    scale_size_manual(values = c(2, 2.5)) +
    scale_alpha_manual(values = c(0.6, 1)) +
    scale_y_continuous(trans = "log10",
                       breaks = 10**(c(0, -2, -4, -6)),
                       labels = c(1, expression(10^-2), expression(10^-4),
                                  expression(10^-6))) +
    theme_bw() +
    #geom_hline(yintercept = 1, lty = 2) +
    geom_hline(yintercept = eop_limit, lty = 3, lwd = 1) +
    scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                       labels = c("Control", "Local", "Global"),
                       values = my_cols[c(8, 2, 6)]) +
    scale_fill_manual(name = "Treatment", breaks = c("C", "L", "G"),
                      labels = c("Control", "Local", "Global"),
                      values = my_cols[c(8, 2, 6)]) +
    scale_shape_manual(name = "Below Limit", values = c(16, 8)) +
    labs(x = "Population", y = "Efficiency of Plaquing Relative to Ancestor") +
    theme(legend.position = "none") +
    NULL)
  dev.off()
}

#Summarize
resis_data_isols <- summarize(group_by(resis_data, Proj, Pop, Treat, Isol),
                              EOP_avg = 10**mean(log10(EOP)),
                              EOP_bd = any(bd))

resis_data_sum <- summarize(group_by(resis_data, Proj, Pop, Treat),
                            EOP_avg = 10**mean(log10(EOP)),
                            EOP_bd = any(bd))

#Make plot of pop-level data
if (make_statplots) {
  tiff("./Output_figures/Isol_resis_pops.tiff", width = 4, height = 4,
       units = "in", res = 300)
  print(ggplot(resis_data_sum[resis_data_sum$Treat != "Anc", ],
               aes(x = Treat, y = EOP_avg, 
                   color = Treat, fill = Treat)) +
          facet_grid(Proj~., labeller = labeller(Proj = my_facet_labels)) +
          geom_point(alpha = 0.6, size = 3, 
                     position = position_dodge(width = .15),
                     aes(group = Pop)) +
          scale_size_manual(values = c(2, 2.5)) +
          scale_alpha_manual(values = c(0.6, 1)) +
          scale_x_discrete(breaks = c("C", "L", "G"),
                           labels = c("Control", "Local", "Global")) +
          scale_y_continuous(trans = "log10",
                             breaks = 10**(c(0, -2, -4, -6)),
                             labels = c(1, expression(10^-2), expression(10^-4),
                                        expression(10^-6))) +
          geom_hline(yintercept = 1, lty = 2) +
          geom_hline(yintercept = eop_limit, lty = 3, lwd = 1) +
          scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                             labels = c("Control", "Local", "Global"),
                             values = my_cols[c(8, 2, 6)]) +
          scale_fill_manual(name = "Treatment", breaks = c("C", "L", "G"),
                            labels = c("Control", "Local", "Global"),
                            values = my_cols[c(8, 2, 6)]) +
          labs(x = "Treatment", y = "Susceptibility to Phage") +
          theme_bw() +
          theme(legend.position = "none",
                axis.text.y = element_text(size = 11), 
                axis.text.x = element_text(size = 11),
                axis.title.y = element_text(size = 12.5),
                axis.title.x = element_text(size = 14),
                legend.text = element_text(size = 13), 
                legend.title = element_text(size = 14),
                strip.text = element_text(size = 14)) +
          NULL)
  dev.off()
}


##Isolate growth curves: read data & prep ----

#Read data
gc_data <- read.csv("./Clean_Data/Isolate_growth_curves.csv",
                    header = T, stringsAsFactors = F)

#For ease of downstream analysis & visualization, 
# for now we'll recode the media into simply "Original" and "Rich"
# Keeping in mind that "Orig" and "Rich" mean different
# medias depending on the project
gc_data$Media[gc_data$Media == "50"] <- "Orig"
gc_data$Media[gc_data$Media == "100"] <- "Rich"
gc_data$Media[gc_data$Media == "25-50"] <- "Orig"
gc_data$Media[gc_data$Media == "50-100"] <- "Rich"

#Reorder projects
gc_data$Proj <- factor(gc_data$Proj, levels = c("7x", "125"))

#Reorder treatments
gc_data$Treat <- factor(gc_data$Treat, levels = c("Anc", "C", "L", "G"))

#Make unique well identifiers
gc_data$uniq_well <- paste(gc_data$Date,
                           gc_data$Proj,
                           gc_data$Pop,
                           gc_data$Treat,
                           gc_data$Isol,
                           gc_data$Rep_Well,
                           gc_data$Media,
                           sep = "_")
gc_data$uniq_well_num <- match(gc_data$uniq_well, unique(gc_data$uniq_well))

#reorder
gc_data <- gc_data[order(gc_data$uniq_well, gc_data$Time_s), ]

##Isolate growth curves: smooth ----
gc_data$sm_movmed3 <- NA
gc_data$sm_movmed3_loess3600 <- NA
gc_data$sm_movmed3_loess25k <- NA
for (my_well in unique(gc_data$uniq_well)) {
  #(leaving out the first half hour of data)
  my_rows <- which(gc_data$uniq_well == my_well &
                     gc_data$Time_s > 1800)
  
  #Calculate the median timestep (timesteps actually vary slightly in
  # the number of seconds they were recorded as differing by)
  med_timestep <- median(gc_data$Time_s[my_rows[2]:my_rows[length(my_rows)]]-
                   gc_data$Time_s[my_rows[1]:my_rows[(length(my_rows)-1)]])
  
  #Do moving-median smoothing
  gc_data$sm_movmed3[my_rows] <- 
    gcplyr::moving_median(cfu_ml ~ Time_s,
                          data = gc_data[my_rows, ],
                          window_width_n = 3)
  
  #Exclude NA vals of sm_movmed3
  my_rows <- which(gc_data$uniq_well == my_well &
                     gc_data$Time_s > 1800 &
                     !is.na(gc_data$sm_movmed3))
  
  #Do loess smoothing (window of 60 mins)
  gc_data$sm_movmed3_loess3600[my_rows] <- 
    loess(sm_movmed3 ~ Time_s, 
          data = gc_data[my_rows, ],
          span = ((3600/med_timestep)+1)/length(my_rows),
          degree = 1)$fitted
  
  #Do loess smoothing (window of ~7 hrs) 
  gc_data$sm_movmed3_loess25k[my_rows] <- 
    loess(sm_movmed3 ~ Time_s, 
          data = gc_data[my_rows, ],
          span = ((25000/med_timestep)+1)/length(my_rows),
          degree = 2)$fitted
}
#Drop plain moving-median column
gc_data <-  select(gc_data, !sm_movmed3)

#Import function that calculates derivatives
calc_deriv <- gcplyr::calc_deriv

#Calculate growth per hour from 25k loess curve
gc_data$deriv_sm_movmed3_loess25k <- 
  calc_deriv(y = gc_data$sm_movmed3_loess25k,
             x = gc_data$Time_s,
             subset_by = gc_data$uniq_well,
             x_scale = 3600)

#Calculate per capita growth per hour from 3600 loess curve
gc_data$percap_deriv_sm_movmed3_loess3600 <- 
  calc_deriv(y = gc_data$sm_movmed3_loess3600,
             x = gc_data$Time_s,
             subset_by = gc_data$uniq_well,
             x_scale = 3600,
             percapita = TRUE)

##Isolate growth curves: ID transitions ----
#These transition points are used for fitting the data

#Group data by unique wells
gc_data <- group_by(gc_data, Date, Proj, Pop, Treat, Isol, Rep_Well, Media,
                    uniq_well, uniq_well_num)

#Import function
find_local_extrema <- gcplyr::find_local_extrema

#Define the window widths (sensitivity) for local extrema search
diauxie_time_window <- 7200

#Summarize data (note: message is for grouping of **output**)
gc_summarized <- dplyr::summarize(
  gc_data,
  #First point where percap growth rate exceeds 0.5
  # (or 0.4 if it never passes 0.5)
  threshold_percap_gr_time =
    Time_s[ifelse(any(percap_deriv_sm_movmed3_loess3600 >= 0.5, na.rm = TRUE),
                  min(which(percap_deriv_sm_movmed3_loess3600 >= 0.5)),
                  min(which(percap_deriv_sm_movmed3_loess3600 >= 0.4)))],
  #First non-endpoint minima in the slope
  diauxie_time = Time_s[find_local_extrema(deriv_sm_movmed3_loess25k,
                                           x = Time_s,
                                           width_limit = diauxie_time_window,
                                           return_maxima = FALSE,
                                           return_endpoints = FALSE,
                                           na.rm = T)[1]]
)

#Change to data frame for cleanliness
gc_summarized <- as.data.frame(gc_summarized)

##Isolate growth curves: make plots of all wells ----
dir.create("./Growth_curve_plots/", showWarnings = FALSE)
if (make_curveplots) {
  #for (my_well in unique(gc_data$uniq_well)) {
  for(my_well in gc_summarized$uniq_well[is.na(gc_summarized$diauxie_time)]) {
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
        geom_point() +
        geom_line(aes(y = sm_movmed3_loess3600),
                  color = "blue", lwd = 1, alpha = 0.5) +
        geom_line(aes(y = sm_movmed3_loess25k),
                  color = "green", lwd = 0.5, alpha = 0.5) +
        ggtitle(gc_data[my_rows[1], "uniq_well"]) +
        scale_y_continuous(trans = "log10") +
        geom_vline(data = gc_summarized[gc_summarized$uniq_well == my_well, ],
                   aes(xintercept = threshold_percap_gr_time), lty = 2) +
        geom_vline(data = gc_summarized[gc_summarized$uniq_well == my_well, ],
                   aes(xintercept = diauxie_time), lty = 2) +
        NULL,
      ggplot(data = gc_data[my_rows, ],
             aes(x = Time_s, y = deriv_sm_movmed3_loess25k)) +
        geom_line(color = "green") +
        geom_vline(data = gc_summarized[gc_summarized$uniq_well == my_well, ], 
                   aes(xintercept = diauxie_time), lty = 2) +
        NULL,
      ggplot(data = gc_data[my_rows, ],
             aes(x = Time_s, y = percap_deriv_sm_movmed3_loess3600)) +
        geom_line(color = "blue") +
        geom_vline(data = gc_summarized[gc_summarized$uniq_well == my_well, ],
                   aes(xintercept = threshold_percap_gr_time), lty = 2) +
        NULL,
      ncol = 1, align = "v"))
    dev.off()
  }
}

#Isolate growth curves: Make example plots (eg for talks) ----
dir.create("./Example_curve_plots/", showWarnings = FALSE)
if (make_curveplots) {
  temp <- gc_data[gc_data$uniq_well == "2017-C_7x_Anc_Anc_Anc_1_Orig", ]
  
  tiff("./Example_curve_plots/gc_plot1.tiff", width = 10, height = 10, 
       units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = cfu_ml)) +
    geom_point(size = 4) +
    labs(y = "Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30)) 
  dev.off()
  
  tiff("./Example_curve_plots/gc_plot2.tiff", width = 10, height = 10, 
       units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = cfu_ml)) +
    geom_point(size = 4) +
    geom_line(aes(y = sm_movmed3_loess3600), color = "blue", lwd = 3, alpha = 0.6) +
    labs(y = "Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30)) 
  dev.off()
  
  tiff("./Example_curve_plots/gc_plot3.tiff", width = 10, height = 10, 
       units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = deriv_sm_movmed3_loess25k)) +
    geom_line(color = "green", lwd = 3) +
    labs(y = "Change in Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30)) 
  dev.off()
  
  tiff("./Example_curve_plots/gc_plot4.tiff", width = 10, height = 10, 
       units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = percap_deriv_sm_movmed3_loess3600)) +
    geom_line(color = "blue", lwd = 3) +
    labs(y = "Per-capita Change in Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30)) 
  dev.off()
  
  #Now add the points
  temp2 <- gc_summarized[gc_summarized$uniq_well == "2017-C_7x_Anc_Anc_Anc_1_Orig", ]
  
  tiff("./Example_curve_plots/gc_plot5.tiff", width = 10, height = 10, units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = deriv_sm_movmed3_loess25k)) +
    geom_line(color = "blue", lwd = 3) +
    labs(y = "Change in Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30)) +
    geom_vline(data = temp2, aes(xintercept = diauxie_time/3600), lty = 2, lwd = 2)
  dev.off()
  
  tiff("./Example_curve_plots/gc_plot6.tiff", width = 10, height = 10, units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = cfu_ml)) +
    geom_point(size = 4) +
    labs(y = "Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30)) +
    geom_vline(data = temp2, aes(xintercept = diauxie_time/3600), lty = 2, lwd = 2)
  dev.off()
  
  tiff("./Example_curve_plots/gc_plot7.tiff", width = 10, height = 10, 
       units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = percap_deriv_sm_movmed3_loess3600)) +
    geom_line(color = "blue", lwd = 3) +
    labs(y = "Per-capita Change in Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30)) +
    geom_vline(data = temp2, aes(xintercept = threshold_percap_gr_time/3600), 
               lty = 2, lwd = 2)
  dev.off()
  
  tiff("./Example_curve_plots/gc_plot8.tiff", width = 10, height = 10, 
       units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = sm_movmed3_loess3600)) +
    geom_line(color = "blue", lwd = 3) +
    labs(y = "Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30)) +
    geom_vline(data = temp2, aes(xintercept = diauxie_time/3600), lty = 2, lwd = 2) +
    geom_vline(data = temp2, aes(xintercept = threshold_percap_gr_time/3600), 
               lty = 2, lwd = 2)
  dev.off()
  
  my_well <- "2017-C_7x_Anc_Anc_Anc_1_Orig"
  t_vals <- temp$Time_s
  sum_row <- which(gc_summarized$uniq_well == my_well)
  temp$pred_vals_fit <- baranyi_func(r = gc_summarized[sum_row, "fit_r"],
                                      k = gc_summarized[sum_row, "fit_k"],
                                      v = gc_summarized[sum_row, "fit_v"],
                                      d0 = gc_summarized[sum_row, "fit_d0"],
                                      t_vals = t_vals)
  
  
  tiff("./Example_curve_plots/gc_plot9.tiff", width = 10, height = 10, 
       units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = cfu_ml)) +
    geom_point(size = 4) +
    geom_line(aes(y = pred_vals_fit), color = "red", lwd = 3, alpha = 0.5) +
    labs(y = "Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30))  +
    geom_vline(data = temp2, aes(xintercept = diauxie_time/3600), lty = 2, lwd = 2) +
    geom_vline(data = temp2, aes(xintercept = threshold_percap_gr_time/3600), 
               lty = 2, lwd = 2)
  dev.off()
  
  tiff("./Example_curve_plots/gc_plot10.tiff", width = 10, height = 10, 
       units = "in", res = 300)
  ggplot(data = temp, aes(x = Time_s/3600, y = cfu_ml)) +
    geom_point(size = 4) +
    geom_line(aes(y = pred_vals_fit), color = "red", lwd = 3, alpha = 0.5) +
    scale_y_continuous(trans = "log10") +
    labs(y = "Bacterial Density", x = "Time (hr)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 30))  +
    geom_vline(data = temp2, aes(xintercept = diauxie_time/3600), lty = 2, lwd = 2) +
    geom_vline(data = temp2, aes(xintercept = threshold_percap_gr_time/3600), 
               lty = 2, lwd = 2)
  dev.off()
}

##Isolate growth curves: Fit curves to data ----
baranyi_func <- function(r, k, v, d0, t_vals) {
  #Modified from Ram et al 2019 with a(t) = 1
  # (equivalent to logistic with a deceleration param)
  if (anyNA(c(r, k, v, d0, t_vals))) {return(NA)}
  t_vals_hrs <- t_vals/3600
  d <- k/((1-(1-((k/d0)**v))*exp(-r*v*t_vals_hrs))**(1/v))
  return(d)
}

baranyi_fit_err <- function(params, t_vals, dens_vals) {
  #params <- c("logk" = ..., "logd0" = ..., "r" = ..., "v" = ...)
  pred_vals <- baranyi_func(r = params["r"],
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
                                  "fit_v" = as.numeric(NA), 
                                  "fit_d0" = as.numeric(NA),
                                  "fit_err" = as.numeric(NA)))
for (sum_row in 1:nrow(gc_summarized)) {
  my_well <- gc_summarized$uniq_well[sum_row]
  
  start_time <- gc_summarized$threshold_percap_gr_time[sum_row]
  
  if(is.na(gc_summarized$diauxie_time[sum_row])) {
    end_time <- max(gc_data$Time_s[gc_data$uniq_well == my_well])
  } else {end_time <- gc_summarized$diauxie_time[sum_row]}
  
  gc_rows <- which(gc_data$uniq_well == my_well &
                     gc_data$Time_s <= end_time &
                     gc_data$Time_s >= start_time)
  
  #Set initial values
  init_d0 <- min(gc_data$cfu_ml[gc_rows])
  init_K <- max(gc_data$cfu_ml[gc_rows])
  init_r <- 1
  init_v <- 1
  
  #Fit
  temp <- optim(par = c("logk" = log10(init_K),
                        "logd0" = log10(init_d0),
                        "r" = init_r,
                        "v" = init_v),
                fn = baranyi_fit_err,
                dens_vals = gc_data$cfu_ml[gc_rows],
                t_vals = gc_data$Time_s[gc_rows],
                method = "L-BFGS-B",
                #logk, logd0, r, v
                lower = c(5, 4, 0, 0),
                upper = c(11, 10, 10, 50))
  
  #Save fit vals
  gc_summarized[sum_row, 
                c("fit_r", "fit_k", "fit_v", "fit_d0", "fit_err")] <-
    data.frame("fit_r" = temp$par["r"], 
               "fit_k" = 10**temp$par["logk"],
               "fit_v" = temp$par["v"], 
               "fit_d0" = 10**temp$par["logd0"],
               "fit_err" = temp$value)
}

##Isolate growth curves: evaluate fits & exclude wells ----
#Make plots of fits
dir.create("./Growth_curve_plots_fits", showWarnings = F)
if(make_curveplots) {
  for (sum_row in 1:nrow(gc_summarized)) {
    my_well <- gc_summarized$uniq_well[sum_row]
    t_vals <- gc_data$Time_s[gc_data$uniq_well == my_well]
    t_vals_hrs <- gc_data$Time_s[gc_data$uniq_well == my_well]/3600

    pred_vals <- baranyi_func(r = gc_summarized[sum_row, "fit_r"],
                               k = gc_summarized[sum_row, "fit_k"],
                               v = gc_summarized[sum_row, "fit_v"],
                               d0 = gc_summarized[sum_row, "fit_d0"],
                               t_vals = t_vals)
    
    png(paste("./Growth_curve_plots_fits/", 
              gc_summarized$uniq_well_num[sum_row], "_",
              my_well, ".png", sep = ""),
        width = 4, height = 4, units = "in", res = 300)
    print(ggplot(data = gc_data[gc_data$uniq_well == my_well, ],
                 aes(x = Time_s, y = sm_movmed3_loess3600)) +
            geom_line(alpha = 0.5) +
            geom_point(size = 0.5, aes(y = cfu_ml)) +
            geom_line(data.frame(Time_s = t_vals, pred_dens = pred_vals),
                      mapping = aes(x = Time_s, y = pred_dens),
                      color = "blue", alpha = 0.5) +
            scale_y_continuous(trans = "log10") +
            geom_vline(aes(xintercept = gc_summarized$threshold_percap_gr_time[
              gc_summarized$uniq_well == my_well]), lty = 2) +
            geom_vline(aes(xintercept = gc_summarized$diauxie_time[
              gc_summarized$uniq_well == my_well]), lty = 2) +
            NULL)
    dev.off()
  }
}

#Looking at the fits, many of the the top 30 worst will be excluded
# from further analyses (in all excluded cases, the other technical 
# replicate was much better). After that, the errors have dropped off
# substantially
# 215 - toss
# 183 - toss
# 15 - toss
# 115 - toss
# 223 - toss
# 197 - toss
# 19 - keep
# 59 - toss
# 445 - toss
# 210 - toss
# 13 - toss
# 44 - toss
# 53 - toss
# 160 - toss
# 238 - toss
# 243 - toss
# 165 - keep
# 185 - keep
# 159 - keep
# 370 - toss
# 211 - keep
# 354 - toss
# 184 - keep
# 334 - toss
# 216 - keep
# 242 - toss
# 212 - keep
# 294 - keep
# 46 - toss
# 364 - keep

gc_summarized[which(gc_summarized$uniq_well_num %in% 
                      c(215, 183, 15, 115, 223, 197, 59, 445, 210, 13, 
                        44, 53, 160, 238, 243, 370, 354, 334, 242, 46)),
              grep("fit_", colnames(gc_summarized))] <- NA

##Isolate growth curves: summarize reps into isols ----
gc_summarized <- group_by(gc_summarized, Date, Proj, Pop, Treat,
                          Isol, Media)
gc_sum_isols <- summarize_at(
  gc_summarized,
  .funs = c(avg = function(x) {mean(x, na.rm = TRUE)}, 
            sd = function(x) {sd(x, na.rm = TRUE)}),
  .vars = c("threshold_percap_gr_time",
            "diauxie_time",
            "fit_r",
            "fit_k",
            "fit_v",
            "fit_d0",
            "fit_err"))
gc_sum_isols <- as.data.frame(gc_sum_isols)

##Isolate growth curves: evaluate rep well variation & exclude isols ----
dir.create("./Growth_curve_rep_plots/", showWarnings = FALSE)
if(make_statplots) {
  for (col_i in grep("fit_", colnames(gc_summarized))) {
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
  
  for (proj in unique(gc_sum_isols$Proj)) {
    for (media in unique(gc_sum_isols$Media)) {
      temp <- gc_sum_isols[gc_sum_isols$Proj == proj &
                             gc_sum_isols$Media == media, 
                           c(1:6, 16)]
      temp$grpnm <- paste(temp$Date, temp$Proj, temp$Pop, temp$Treat, temp$Isol, temp$Media,
                          sep = "_")
      temp <- temp[order(-temp$fit_r_sd), ]
      print(temp[1:10, c(7, 8)])
    }
  }
}

#Taking a look at what to cut out:
#for 7x Orig, would like to cut top 7 (above 0.08/0.13 sd)
#for 7x Rich would like to cut top 5 (above 0.1/0.12 sd)
#for 125 Orig would like to cut top 6 (above 0.11/0.40 sd)
#for 125 Rich would like to cut top 4 (above 0.13/0.14 sd)
#Will use 0.12 as 7x cutoff, .14 as 125 cutoff

gc_sum_isols[(gc_sum_isols$Proj == "7x" & 
                !is.na(gc_sum_isols$fit_r_sd) &
                gc_sum_isols$fit_r_sd > 0.12),
             grep("fit_", colnames(gc_sum_isols))] <- NA
gc_sum_isols[gc_sum_isols$Proj == "125" & 
               !is.na(gc_sum_isols$fit_r_sd) &
               gc_sum_isols$fit_r_sd > 0.14,
             grep("fit_", colnames(gc_sum_isols))] <- NA

##Isolate growth curves: plot all isols ----
my_facet_labels <- c("7x" = "Weak Phage", 
                     "125" = "Strong Phage",
                     "C" = "Control", "G" = "Global", "L" = "Local",
                     "Anc" = "Ancstr",
                     "Rich" = "Rich Media", "Orig" = "Original Media")
dir.create("./Growth_curve_variables_plots/", showWarnings = FALSE)
if (make_statplots) {
  my_vars <- c("fit_r", "fit_k", "fit_v")
  for (i in 1:length(my_vars)) {
    var_root <- my_vars[i]
    var_name <- c("Maximum Per Capita Growth Rate (r) (/hr)", 
                  "Density at Diauxic Shift (k) (cfu/mL)", 
                  "Deceleration Parameter (v)")[i]
    var <- paste(var_root, "_avg", sep = "")
    var_sd <- paste(var_root, "_sd", sep = "")
    tiff(paste("./Growth_curve_variables_plots/", var_root, ".tiff", sep = ""),
         width = 5, height = 5, units = "in", res = 300)
    print(ggplot(data = gc_sum_isols,
                 aes(x = Pop, y = get(var), group = Pop, color = Treat)) +
            geom_point(position = position_dodge(0.6)) +
            facet_nested(Proj ~ Media+Treat, scales = "free",
                         labeller = labeller(Proj = my_facet_labels,
                                             Treat = my_facet_labels,
                                             Media = my_facet_labels)) +
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

##Isolate analysis: merge dataframes ----
isol_data <- full_join(isol_migr_sum_isols, resis_data_isols)
isol_data <- full_join(isol_data, gc_sum_isols)

isol_data <- select(isol_data, 
                    Proj, Pop, Treat, Isol,
                    Date, Media, threshold_percap_gr_time_avg, diauxie_time_avg,
                    fit_r_avg, fit_k_avg, fit_v_avg, fit_d0_avg,
                    radius_mm_hr_del_avg,
                    EOP_avg, EOP_bd)

##Isolate analysis: mixed effects modeling ----


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


