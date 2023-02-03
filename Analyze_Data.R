## Load packages and color scale ----
library(ggplot2)
library(dplyr)
library(ggh4x)
library(lme4)
library(pbkrtest)
library(emmeans)
library(lmerTest)
library(gcplyr)

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
          facet_grid(Proj~., labeller = labeller(Proj = my_facet_labels),
                     scales = "free_y") +
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
                axis.title = element_text(size = 17),
                legend.text = element_text(size = 16), 
                legend.title = element_text(size = 17),
                strip.text = element_text(size = 11)) +
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
          facet_grid(Proj~., labeller = labeller(Proj = my_facet_labels),
                     scales = "free_y") +
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
          theme(axis.text.y = element_text(size = 12), 
                axis.text.x = element_text(size = 12),
                axis.title = element_text(size = 17),
                legend.text = element_text(size = 16), 
                legend.title = element_text(size = 17),
                strip.text = element_text(size = 11)) +
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
          scale_color_manual(breaks = c("Anc", "C", "L", "G"),
                             values = my_cols[c(3, 8, 2, 6)]) +
          scale_fill_manual(breaks = c("Anc", "C", "L", "G"),
                            values = my_cols[c(3, 8, 2, 6)]) +
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
                                                     Treat = my_facet_labels),
                     scales = "free_y") +
          theme_bw() + 
          labs(y = paste("\u0394", "Dispersal (mm/hr)", sep = ""),
               x = "Population") +
          geom_hline(yintercept = 0, lty = 2) +
          scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                             labels = c("Control", "Local", "Global"),
                             values = my_cols[c(8, 2, 6)]) +
          scale_fill_manual(name = "Treatment", breaks = c("C", "L", "G"),
                            labels = c("Control", "Local", "Global"),
                            values = my_cols[c(8, 2, 6)]) +
          theme(legend.position = "none",
                axis.text = element_text(size = 11), 
                axis.title.y = element_text(size = 15),
                axis.title.x = element_text(size = 14),
                strip.text.x = element_text(size = 15),
                strip.text.y = element_text(size = 13)) +
          
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

#Statistics
isol_migration$PPT <- 
  paste(isol_migration$Proj, isol_migration$Pop, isol_migration$Treat)

#Test for differences from Ancestor
mig7x_lmeranc <- lmer(radius_mm_hr ~ (1|PPT) + (1|start_timestamp) + Treat,
                      data = isol_migration[isol_migration$Proj == "7x", ])
mig125_lmeranc <- lmer(radius_mm_hr ~ (1|PPT) + (1|start_timestamp) + Treat,
                       data = isol_migration[isol_migration$Proj == "125", ])

emmeans(mig7x_lmeranc, ~ Treat, contr = "trt.vs.ctrl", 
        ref = which(levels(isol_migration$Treat) == "Anc"))
emmeans(mig125_lmeranc, ~ Treat, contr = "trt.vs.ctrl", 
        ref = which(levels(isol_migration$Treat) == "Anc"))

#Test for differences bt treats
mig7x_lmer <- lmer(radius_mm_hr_del ~ (1|PPT) + Treat,
                 data = isol_migration[isol_migration$Isol != "Anc" &
                                         isol_migration$Proj == "7x", ])
mig125_lmer <- lmer(radius_mm_hr_del ~ (1|PPT) + Treat,
                   data = isol_migration[isol_migration$Isol != "Anc" &
                                           isol_migration$Proj == "125", ])

mig7x_lmernull <- lmer(radius_mm_hr_del ~ (1|PPT),
                   data = isol_migration[isol_migration$Isol != "Anc" &
                                           isol_migration$Proj == "7x", ])
mig125_lmernull <- lmer(radius_mm_hr_del ~ (1|PPT),
                    data = isol_migration[isol_migration$Isol != "Anc" &
                                            isol_migration$Proj == "125", ])

KRmodcomp(mig7x_lmer, mig7x_lmernull)
KRmodcomp(mig125_lmer, mig125_lmernull)

summary(mig125_lmer)

emmeans(mig125_lmer, ~ Treat, contr = "pairwise", adjust = "tukey")

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
          labs(x = "Population", y = "Susceptibility to Parasite Infection") +
          theme(legend.position = "none",
                axis.text = element_text(size = 11), 
                axis.title.y = element_text(size = 15),
                axis.title.x = element_text(size = 14),
                strip.text.x = element_text(size = 15),
                strip.text.y = element_text(size = 13)) +
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
             percapita = TRUE, blank = 0)

##Isolate growth curves: ID transitions ----
#These transition points are used for fitting the data

#Group data by unique wells
gc_data <- group_by(gc_data, Date, Proj, Pop, Treat, Isol, Rep_Well, Media,
                    uniq_well, uniq_well_num)

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
# 
# 301 - toss (bc of bad diauxic shift fit)
# 330 - toss (bc of bad diauxic shift fit)

gc_summarized[which(gc_summarized$uniq_well_num %in% 
                      c(215, 183, 15, 115, 223, 197, 59, 445, 210, 13, 
                        44, 53, 160, 238, 243, 370, 354, 334, 242, 46, 
                        301, 330)),
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
#(Note the existence of strong batch effects in the data, e.g. max per cap
# growth rate. Statistical testing with random effects for batch should
# counteract this, but it's a good qualifier to keep in mind)

#Make colums of lag time & diauxie time in hrs
gc_sum_isols <-
  mutate(gc_sum_isols,
         across(c("threshold_percap_gr_time_avg", "diauxie_time_avg",
                  "threshold_percap_gr_time_sd", "diauxie_time_sd"),
                .fns = list(hr = function(x) {x/3600})))
colnames(gc_sum_isols)[grep("_hr", colnames(gc_sum_isols))] <-
  c("threshold_percap_gr_time_hr_avg", "diauxie_time_hr_avg",
    "threshold_percap_gr_time_hr_sd", "diauxie_time_hr_sd")
                
my_facet_labels <- c("7x" = "Weak Phage", 
                     "125" = "Strong Phage",
                     "C" = "Control", "G" = "Global", "L" = "Local",
                     "Anc" = "Ancstr",
                     "Rich" = "Rich Media", "Orig" = "Original Media")
dir.create("./Growth_curve_variables_plots/", showWarnings = FALSE)
if (make_statplots) {
  my_vars <- c("fit_r", "fit_k", "fit_v", 
               "threshold_percap_gr_time_hr", "diauxie_time_hr")
  for (i in 1:length(my_vars)) {
    var_root <- my_vars[i]
    var_name <- c("Maximum Per Capita Growth Rate (r) (/hr)", 
                  "Density at Diauxic Shift (k) (cfu/mL)", 
                  "Deceleration Parameter (v)",
                  "Lag time (hr)", "Diauxic Shift Time (kt) (hr)")[i]
    
    var <- paste(var_root, "_avg", sep = "")
    var_sd <- paste(var_root, "_sd", sep = "")
    
    tiff(paste("./Growth_curve_variables_plots/", var_root, ".tiff", sep = ""),
         width = 5, height = 5, units = "in", res = 300)
    print(ggplot(data = gc_sum_isols,
                 aes(x = Pop, y = get(var), group = Pop, color = Treat)) +
            geom_point(position = position_dodge(0.6), alpha = 0.7) +
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
    
    #Make main-text plot
    if(var_root == "fit_r") {
      tiff("./Growth_curve_variables_plots/fit_r_maintext.tiff",
           width = 5, height = 4.5, units = "in", res = 300)
      print(ggplot(data = gc_sum_isols[gc_sum_isols$Media == "Orig", ],
                   aes(x = Pop, y = get(var), group = Pop, color = Treat)) +
              geom_point(position = position_dodge(0.6), alpha = 0.7,
                         size = 2) +
              facet_grid(Proj ~ Treat, scales = "free",
                         labeller = labeller(Proj = my_facet_labels,
                                             Treat = my_facet_labels)) +
              scale_color_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                                 labels = c("Ancestor", "Control", "Local", "Global"),
                                 values = my_cols[c(3, 8, 2, 6)]) +
              geom_errorbar(aes(x = Pop, ymin = get(var)-get(var_sd),
                                ymax = get(var)+get(var_sd)),
                            position = position_dodge(0.6),
                            width = 0.2) +
              labs(y = "Maximum Per Capita Growth Rate (/hr)", x = "Population") +
              theme_bw() +
              theme(legend.position = "none",
                    axis.text = element_text(size = 11),
                    axis.title.y = element_text(size = 15),
                    axis.title.x = element_text(size = 14),
                    strip.text.x = element_text(size = 15),
                    strip.text.y = element_text(size = 13))
      )
      dev.off()
    }
  }
}

##Isolate analysis: merge dataframes ----
isol_data <- full_join(isol_migr_sum_isols, resis_data_isols)
isol_data <- full_join(isol_data, gc_sum_isols)

isol_data <- select(isol_data, 
                    Proj, Pop, Treat, Isol,
                    Date, Media, 
                    threshold_percap_gr_time_hr_avg, diauxie_time_hr_avg,
                    fit_r_avg, fit_k_avg, fit_v_avg, fit_d0_avg,
                    radius_mm_hr_del_avg,
                    EOP_avg, EOP_bd)
isol_data$PPT <- paste(isol_data$Proj, isol_data$Pop, isol_data$Treat, sep = "_")

##Isolate analysis: mixed effects modeling ----
isol_data$resis_cat <-
  ifelse(isol_data$EOP_bd, "Resis",
         ifelse(isol_data$EOP_avg > 0.1, "Sens", "Part Resis"))
isol_data$resis_cat <- factor(isol_data$resis_cat,
                              levels = c("Sens", "Part Resis", "Resis"))

mixed_model0_7x_O = lmer(fit_r_avg ~ (1|PPT) + (1|Date),
                    data = isol_data[isol_data$Proj == "7x" &
                                             isol_data$Media == "Orig",])
mixed_model0_7x_R = lmer(fit_r_avg ~ (1|PPT) + (1|Date),
                    data = isol_data[isol_data$Proj == "7x" &
                                             isol_data$Media == "Rich",])
mixed_model0_125_O = lmer(fit_r_avg ~ (1|PPT) + (1|Date),
                    data = isol_data[isol_data$Proj == "125" &
                                             isol_data$Media == "Orig",])
mixed_model0_125_R = lmer(fit_r_avg ~ (1|PPT) + (1|Date),
                    data = isol_data[isol_data$Proj == "125" &
                                             isol_data$Media == "Rich",])

mixed_model1_7x_O = lmer(fit_r_avg ~ Treat + (1|PPT) + (1|Date),
                         data = isol_data[isol_data$Proj == "7x" &
                                                  isol_data$Media == "Orig",])
mixed_model1_7x_R = lmer(fit_r_avg ~ Treat + (1|PPT) + (1|Date),
                         data = isol_data[isol_data$Proj == "7x" &
                                                  isol_data$Media == "Rich",])
mixed_model1_125_O = lmer(fit_r_avg ~ Treat + (1|PPT) + (1|Date),
                          data = isol_data[isol_data$Proj == "125" &
                                                   isol_data$Media == "Orig",])
mixed_model1_125_R = lmer(fit_r_avg ~ Treat + (1|PPT) + (1|Date),
                          data = isol_data[isol_data$Proj == "125" &
                                                   isol_data$Media == "Rich",])

mixed_model2_7x_O = lmer(fit_r_avg ~ resis_cat + (1|PPT) + (1|Date),
                         data = isol_data[isol_data$Proj == "7x" &
                                                  isol_data$Media == "Orig",])
mixed_model2_7x_R = lmer(fit_r_avg ~ resis_cat + (1|PPT) + (1|Date),
                         data = isol_data[isol_data$Proj == "7x" &
                                                  isol_data$Media == "Rich",])
mixed_model2_125_O = lmer(fit_r_avg ~ resis_cat + (1|PPT) + (1|Date),
                          data = isol_data[isol_data$Proj == "125" &
                                                   isol_data$Media == "Orig",])
mixed_model2_125_R = lmer(fit_r_avg ~ resis_cat + (1|PPT) + (1|Date),
                          data = isol_data[isol_data$Proj == "125" &
                                                   isol_data$Media == "Rich",])

mixed_model3_7x_O = lmer(fit_r_avg ~ resis_cat + Treat + (1|PPT) + (1|Date),
                    data = isol_data[isol_data$Proj == "7x" &
                                             isol_data$Media == "Orig",])
mixed_model3_7x_R = lmer(fit_r_avg ~ resis_cat + Treat + (1|PPT) + (1|Date),
                    data = isol_data[isol_data$Proj == "7x" &
                                             isol_data$Media == "Rich",])
mixed_model3_125_O = lmer(fit_r_avg ~ resis_cat + Treat + (1|PPT) + (1|Date),
                    data = isol_data[isol_data$Proj == "125" &
                                             isol_data$Media == "Orig",])
mixed_model3_125_R = lmer(fit_r_avg ~ resis_cat + Treat + (1|PPT) + (1|Date),
                    data = isol_data[isol_data$Proj == "125" &
                                             isol_data$Media == "Rich",])

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


##Isolate analysis: PCA ----

#Cast measurements in different medias into different columns
isol_data_wide <- tidyr::pivot_wider(
  isol_data[complete.cases(isol_data), ],
  values_from = c("threshold_percap_gr_time_hr_avg", "diauxie_time_hr_avg",            
                  "fit_r_avg", "fit_k_avg", "fit_v_avg", "fit_d0_avg"),
  names_from = Media,
  id_cols = c("Proj", "Pop", "Treat", "Isol", "Date",
              "radius_mm_hr_del_avg", "EOP_avg", "EOP_bd", "resis_cat"),
  values_fill = NA)
                  
#Reduce to complete cases
isol_data_wide <- isol_data_wide[complete.cases(isol_data_wide), ]

#Make resis trait as -log10(EOP)
isol_data_wide$resis <- -log10(isol_data_wide$EOP_avg)

#Make k traits as log10
isol_data_wide$fit_k_log_avg_Orig <- log10(isol_data_wide$fit_k_avg_Orig)
isol_data_wide$fit_k_log_avg_Rich <- log10(isol_data_wide$fit_k_avg_Rich)

##Check univariate normality
pca_cols <- c("radius_mm_hr_del_avg",
              "threshold_percap_gr_time_hr_avg_Orig",
              "threshold_percap_gr_time_hr_avg_Rich",
              "diauxie_time_hr_avg_Orig",
              "diauxie_time_hr_avg_Rich",
              "fit_r_avg_Orig", "fit_r_avg_Rich",                      
              "fit_k_log_avg_Orig", "fit_k_log_avg_Rich",
              "resis")

for (col in pca_cols) {
  print(ggplot(data = isol_data_wide,
               aes_string(x = col, fill = "Treat")) +
          geom_histogram(bins = 10) +
          facet_grid(Proj ~ .) +
          ggtitle(col))
}
#Besides resistance (which is bimodal), 
# within conditions (proj) all traits are ~univariate normal

##Now check for multivariate normality:

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

#Make CSQ plot for multivariate normality
CSQPlot(as.data.frame(isol_data_wide[isol_data_wide$Proj == "7x", pca_cols]))
CSQPlot(as.data.frame(isol_data_wide[isol_data_wide$Proj == "125", pca_cols]))
#Our weak phage data is pretty multivariate normal,
#but our strong phage data is pretty not multivariate normal

##Run PCA
isol_pca_7x <- prcomp(isol_data_wide[isol_data_wide$Proj == "7x",
                                     pca_cols],
                      center = TRUE, scale = TRUE, retx = TRUE)
isol_pca_125 <- prcomp(isol_data_wide[isol_data_wide$Proj == "125",
                                     pca_cols],
                      center = TRUE, scale = TRUE, retx = TRUE)

#Merge pca with data
isol_pca_7x$x <- cbind(as.data.frame(isol_pca_7x$x), 
                       isol_data_wide[isol_data_wide$Proj == "7x", ])
isol_pca_125$x <- cbind(as.data.frame(isol_pca_125$x), 
                       isol_data_wide[isol_data_wide$Proj == "125", ])

#Rename rotations for better plotting
replace_names = c("radius_mm_hr_del_avg" = "disp",
                  "threshold_percap_gr_time_hr_avg_Orig" = "lag_Or",
                  "threshold_percap_gr_time_hr_avg_Rich" = "lag_Ri",
                  "diauxie_time_hr_avg_Orig" = "kt_Or",
                  "diauxie_time_hr_avg_Rich" = "kt_Ri",
                  "fit_r_avg_Orig" = "r_Or",
                  "fit_r_avg_Rich" = "r_Ri",
                  "fit_k_log_avg_Orig" = "k_Or",
                  "fit_k_log_avg_Rich" = "k_Ri")

row.names(isol_pca_7x$rotation) <- plyr::revalue(
  x = row.names(isol_pca_7x$rotation), replace = replace_names)
row.names(isol_pca_125$rotation) <- plyr::revalue(
  x = row.names(isol_pca_125$rotation), replace = replace_names)       

summary(isol_pca_7x)
summary(isol_pca_125)

## Isolate analysis: make distance PCA biplots ----
# (scores and loading directly plotted)
if(make_statplots) {
  arrow_len <- 5 #multiplier for arrow lengths for vis purposes
  
  tiff("./Output_figures/weakphage_PCA.tiff", 
       width = 12, height = 10, units = "in", res = 300)
  weak_pca <- ggplot(isol_pca_7x$x, aes(x = PC1, y = PC2)) +
    ggtitle("Weak Phage") +
    geom_point(aes(color = Treat, fill = Treat, shape = Pop), 
               size = 10, alpha = 0.7) +
    geom_segment(data = as.data.frame(isol_pca_7x$rotation),
                 aes(x = 0, y = 0, xend = arrow_len*PC1, yend = arrow_len*PC2),
                 arrow = arrow(length = unit(0.02, "npc")),
                 alpha = .5, lwd = 2, color = "red4") +
    ggrepel::geom_text_repel(data = as.data.frame(isol_pca_7x$rotation),
                             aes(x = arrow_len*PC1, y = arrow_len*PC2,
                                 label = row.names(isol_pca_7x$rotation)),
              size = 8, alpha = .8, color = "red4", seed = 8,
              min.segment.length = unit(1, "native")) +
    theme_bw() +
    labs(x = paste("PC1 (", 
                   round((100*((isol_pca_7x$sdev)**2)/
                           sum((isol_pca_7x$sdev)**2))[1], 1),
                   "%)", sep = ""),
         y = paste("PC2 (", 
                   round((100*((isol_pca_7x$sdev)**2)/
                            sum((isol_pca_7x$sdev)**2))[2], 1),
                   "%)", sep = "")) +
    theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 22),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 32),
          plot.title = element_text(size = 36)) +
    scale_color_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                       labels = c("Ancestor", "Control", "Local", "Global"),
                       values = my_cols[c(7, 8, 2, 6)]) +
    scale_fill_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                       labels = c("Ancestor", "Control", "Local", "Global"),
                       values = my_cols[c(7, 8, 2, 6)]) +
    scale_shape_manual(name = "Population",
                       breaks = c("Anc", "A", "B", "C", "D", "E"),
                       values = c(21, 21, 22, 23, 24, 25)) +
    NULL
  print(weak_pca)
  dev.off()

  tiff("./Output_figures/strongphage_PCA.tiff", 
       width = 12, height = 10, units = "in", res = 300)
  strong_pca <- ggplot(isol_pca_125$x, aes(x = PC1, y = PC2)) +
    ggtitle("Strong Phage") +
    geom_point(aes(color = Treat, fill = Treat, shape = Pop), 
               size = 10, alpha = 0.7) +
    geom_segment(data = as.data.frame(isol_pca_125$rotation),
                 aes(x = 0, y = 0, xend = arrow_len*PC1, yend = arrow_len*PC2),
                 arrow = arrow(length = unit(0.02, "npc")),
                 alpha = .5, lwd = 2, color = "red4") +
    ggrepel::geom_text_repel(data = as.data.frame(isol_pca_125$rotation),
                             aes(x = arrow_len*PC1, y = arrow_len*PC2,
                                 label = row.names(isol_pca_125$rotation)),
                             size = 8, alpha = .8, color = "red4", seed = 8,
                             min.segment.length = unit(1, "native")) +
    theme_bw() +
    labs(x = paste("PC1 (", 
                   round((100*((isol_pca_125$sdev)**2)/
                            sum((isol_pca_125$sdev)**2))[1], 1),
                   "%)", sep = ""),
         y = paste("PC2 (", 
                   round((100*((isol_pca_125$sdev)**2)/
                            sum((isol_pca_125$sdev)**2))[2], 1),
                   "%)", sep = "")) +
    theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 22),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 32),
          plot.title = element_text(size = 36)) +
    scale_color_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                       labels = c("Ancestor", "Control", "Local", "Global"),
                       values = my_cols[c(7, 8, 2, 6)]) +
    scale_fill_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                      labels = c("Ancestor", "Control", "Local", "Global"),
                      values = my_cols[c(7, 8, 2, 6)]) +
    scale_shape_manual(name = "Population",
                       breaks = c("Anc", "A", "B", "C", "D", "E"),
                       values = c(21, 21, 22, 23, 24, 25)) +
    NULL
  print(strong_pca)
  dev.off()
}

## Isolate analysis: make correlation PCA biplots ----

#Do normalizations needed for correlation biplot
#In short: x = x %*% diag(1/sdev)
#          rotation = rotation %*% diag(sdev)
isol_pca_7x$x_corr <- isol_pca_7x$x
isol_pca_7x$x_corr[, grep("PC", colnames(isol_pca_7x$x_corr))] <-
  as.matrix(isol_pca_7x$x_corr[, grep("PC", colnames(isol_pca_7x$x_corr))]) %*%
  diag(1/isol_pca_7x$sdev)

isol_pca_125$x_corr <- isol_pca_125$x
isol_pca_125$x_corr[, grep("PC", colnames(isol_pca_125$x_corr))] <- 
  as.matrix(isol_pca_125$x_corr[, grep("PC", colnames(isol_pca_125$x_corr))]) %*%
  diag(1/isol_pca_125$sdev)

isol_pca_7x$rotation_corr <- isol_pca_7x$rotation %*% diag(isol_pca_7x$sdev)
isol_pca_125$rotation_corr <- isol_pca_125$rotation %*% diag(isol_pca_125$sdev)
colnames(isol_pca_7x$rotation_corr) <- 
  paste("PC", 1:ncol(isol_pca_7x$rotation_corr), sep = "")
colnames(isol_pca_125$rotation_corr) <- 
  paste("PC", 1:ncol(isol_pca_125$rotation_corr), sep = "")

#Plot correlation biplots
if(make_statplots) {
  arrow_len <- 2 #multiplier for arrow lengths for vis purposes
  
  tiff("./Output_figures/weakphage_PCA_corbiplot.tiff", 
       width = 12, height = 10, units = "in", res = 300)
  weak_pca <- ggplot(isol_pca_7x$x_corr, aes(x = PC1, y = PC2)) +
    ggtitle("Weak Phage") +
    geom_point(aes(color = Treat, fill = Treat, shape = Pop), 
               size = 10, alpha = 0.7) +
    geom_segment(data = as.data.frame(isol_pca_7x$rotation_corr),
                 aes(x = 0, y = 0, xend = arrow_len*PC1, yend = arrow_len*PC2),
                 arrow = arrow(length = unit(0.02, "npc")),
                 alpha = .5, lwd = 2, color = "red4") +
    ggrepel::geom_text_repel(data = as.data.frame(isol_pca_7x$rotation_corr),
                             aes(x = arrow_len*PC1, y = arrow_len*PC2,
                                 label = row.names(isol_pca_7x$rotation_corr)),
                             size = 8, alpha = .8, color = "red4", seed = 8,
                             min.segment.length = unit(1, "native")) +
    theme_bw() +
    labs(x = paste("PC1 (", 
                   round((100*((isol_pca_7x$sdev)**2)/
                            sum((isol_pca_7x$sdev)**2))[1], 1),
                   "%)", sep = ""),
         y = paste("PC2 (", 
                   round((100*((isol_pca_7x$sdev)**2)/
                            sum((isol_pca_7x$sdev)**2))[2], 1),
                   "%)", sep = "")) +
    theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 22),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 32),
          plot.title = element_text(size = 36)) +
    scale_color_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                       labels = c("Ancestor", "Control", "Local", "Global"),
                       values = my_cols[c(7, 8, 2, 6)]) +
    scale_fill_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                      labels = c("Ancestor", "Control", "Local", "Global"),
                      values = my_cols[c(7, 8, 2, 6)]) +
    scale_shape_manual(name = "Population",
                       breaks = c("Anc", "A", "B", "C", "D", "E"),
                       values = c(21, 21, 22, 23, 24, 25)) +
    NULL
  print(weak_pca)
  dev.off()
  
  tiff("./Output_figures/strongphage_PCA_corbiplot.tiff", 
       width = 12, height = 10, units = "in", res = 300)
  strong_pca <- ggplot(isol_pca_125$x_corr, aes(x = PC1, y = PC2)) +
    ggtitle("Strong Phage") +
    geom_point(aes(color = Treat, fill = Treat, shape = Pop), 
               size = 10, alpha = 0.7) +
    geom_segment(data = as.data.frame(isol_pca_125$rotation_corr),
                 aes(x = 0, y = 0, xend = arrow_len*PC1, yend = arrow_len*PC2),
                 arrow = arrow(length = unit(0.02, "npc")),
                 alpha = .5, lwd = 2, color = "red4") +
    ggrepel::geom_text_repel(data = as.data.frame(isol_pca_125$rotation_corr),
                             aes(x = arrow_len*PC1, y = arrow_len*PC2,
                                 label = row.names(isol_pca_125$rotation_corr)),
                             size = 8, alpha = .8, color = "red4", seed = 8,
                             min.segment.length = unit(1, "native")) +
    theme_bw() +
    labs(x = paste("PC1 (", 
                   round((100*((isol_pca_125$sdev)**2)/
                            sum((isol_pca_125$sdev)**2))[1], 1),
                   "%)", sep = ""),
         y = paste("PC2 (", 
                   round((100*((isol_pca_125$sdev)**2)/
                            sum((isol_pca_125$sdev)**2))[2], 1),
                   "%)", sep = "")) +
    theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 22),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 32),
          plot.title = element_text(size = 36)) +
    scale_color_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                       labels = c("Ancestor", "Control", "Local", "Global"),
                       values = my_cols[c(7, 8, 2, 6)]) +
    scale_fill_manual(name = "Treatment", breaks = c("Anc", "C", "L", "G"),
                      labels = c("Ancestor", "Control", "Local", "Global"),
                      values = my_cols[c(7, 8, 2, 6)]) +
    scale_shape_manual(name = "Population",
                       breaks = c("Anc", "A", "B", "C", "D", "E"),
                       values = c(21, 21, 22, 23, 24, 25)) +
    NULL
  print(strong_pca)
  dev.off()
}

##Isolate analysis: summarize into pops ----

isol_data <- group_by(isol_data, Proj, Pop, Treat, Media)
isol_pops <- summarize(isol_data[!is.na(isol_data$Media), ],
                       across(.cols = c("threshold_percap_gr_time_hr_avg",
                                        "diauxie_time_hr_avg",
                                        "fit_r_avg",
                                        "fit_k_avg",
                                        "fit_v_avg",
                                        "fit_d0_avg",
                                        "radius_mm_hr_del_avg",
                                        "EOP_avg"),
                              .fns = list(avg = mean, sd = sd),
                              na.rm = TRUE),
                       "EOP_bd" = any(EOP_bd))
isol_pops <- as.data.frame(isol_pops)