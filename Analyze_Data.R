##TODO: 
##      FIX gc summarize so it's not printing out deriv_sm_loess_25k for ea group
##      Re-run fitting Baryani to data
##      Check & fix bad fits
##      Calculate lag time
##      Plot lag time, r, k, v
##      Make PCA: lag time, r, k, v, resistance, migration
##      Stats
##      normalization of migration data (for time, via log transform?)

## Load packages and color scale ----
library("ggplot2")
library("scales")
library("dplyr")
library("data.table")
library("MASS")
library("ggh4x")
library("npmv")

#Okabe and Ito 2008 colorblind-safe qualitative color scale
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")

#Global options
make_curveplots <- FALSE
make_statplots <- TRUE

##Experimental evolution migration ----
exper_evol_migr <- read.csv("./Clean_Data/Experimental_evolution_growth.csv")

#Reorder so weak phage is first
exper_evol_migr$Proj <- factor(exper_evol_migr$Proj,
                               levels = c("7x", "125"))

#Drop points after T14
exper_evol_migr <- exper_evol_migr[exper_evol_migr$Timepoint <= 14, ]

#Calculate total area
exper_evol_migr$area_cm2 <- pi*exper_evol_migr$Width_cm/2*exper_evol_migr$Height_cm/2

#Let's assume that area increases exponentially
# (e.g. what we might naively expect pop size to do)
# then A(t) = C e^(kt)
#So ln(A/C)/t = k (the growth constant of area
#Well C is just the starting area, which we know has a diameter
# of 2-3mm (see 2017-01-21 in lab notebook)
# (even if the density is lower, since all points are being
#  normalized the same way changes in C would simply move
#  all points up/down proportional to t)
# which equals an area of 0.01*pi cm^2 = 0.0314 cm^2
#So k = ln(A/0.0314)/t

#Calculate k
exper_evol_migr$area_k <- log(exper_evol_migr$area_cm2/(0.01*pi))/
  exper_evol_migr$time_since_inoc

#Make plot of all pops
if (make_statplots) {
  ggplot(data = exper_evol_migr,
         aes(x = Timepoint, y = area_k,
             group = paste(Pop, Treat), color = Treat)) +
    geom_line() +
    facet_grid(~Proj)
}
  
#Summarize
exper_evol_migr <- group_by(exper_evol_migr, Proj, Treat, Timepoint)
exper_evol_summ <- summarize(exper_evol_migr,
                             area_mean = mean(area_cm2),
                             area_sd = sd(area_cm2),
                             area_n = n(),
                             area_k_mean = mean(area_k),
                             area_k_sd = sd(area_k))

#Make plot of summarized data
if (make_statplots) {
  my_facet_labels <- c("7x" = "Weak Phage", "125" = "Strong Phage")
  
  ggplot(data = exper_evol_summ, aes(x = Timepoint, y = area_k_mean,
                                     color = Treat)) +
    geom_point(position = position_dodge(0.2)) + 
    geom_line(size = 1.2, position = position_dodge(0.2)) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 11), axis.text.x = element_text(size = 11),
          legend.text = element_text(size = 16)) +
    facet_grid(~Proj, labeller = labeller(Proj = my_facet_labels)) +
    geom_errorbar(aes(ymax = area_k_mean+area_k_sd, 
                      ymin = area_k_mean-area_k_sd),
                  width=1, size = .7, position=position_dodge(0.2)) +
    labs(x = "Transfer", 
         y = expression(paste("Mean Area of Growth per Hour ( ", cm^2, "/hr)"))) + 
    scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                    labels = c("Control", "Local", "Global"),
                    values = my_cols[c(8, 2, 6)]) +
    scale_x_continuous(breaks = c(0, 7, 14)) +
    NULL

  #Make plot with both summarized and non-summarized data
  tiff("./Output_figures/Exper_evol_migr.tiff",
       width = 7, height = 4, units = "in", res = 300)
  ggplot(data = exper_evol_migr,
                 aes(x = Timepoint, y = area_k, 
                     group = paste(Treat, Pop),
                     color = Treat)) +
  #  geom_point(size = 0.5, alpha = 0.5) +
           geom_line(alpha = 0.5, lwd = .4) +
           facet_grid(~Proj, labeller = labeller(Proj = my_facet_labels)) +
    geom_line(data = exper_evol_summ,
              aes(x = Timepoint, y = area_k_mean, color = Treat,
                  group = Treat),
              size = 1.3) +
    labs(x = "Transfer", y = "Total Growth Parameter") + 
    scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                       labels = c("Control", "Local", "Global"),
                       values = my_cols[c(8, 2, 6)]) +
    scale_x_continuous(breaks = c(0, 7, 14)) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 11), 
          axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 13), 
          legend.title = element_text(size = 14),
          strip.text = element_text(size = 14)) +
    NULL
  dev.off()
  
  tiff("./Output_figures/Exper_evol_migr_tall.tiff",
       width = 7, height = 4, units = "in", res = 300)
  ggplot(data = exper_evol_migr,
         aes(x = Timepoint, y = area_k, 
             group = paste(Treat, Pop),
             color = Treat)) +
    #  geom_point(size = 0.5, alpha = 0.5) +
    geom_line(alpha = 0.5, lwd = .4) +
    facet_grid(Proj~., labeller = labeller(Proj = my_facet_labels)) +
    geom_line(data = exper_evol_summ,
              aes(x = Timepoint, y = area_k_mean, color = Treat,
                  group = Treat),
              size = 1.3) +
    labs(x = "Transfer", y = "Total Growth Parameter") + 
    scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                       labels = c("Control", "Local", "Global"),
                       values = my_cols[c(8, 2, 6)]) +
    scale_x_continuous(breaks = c(0, 7, 14)) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 11), 
          axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 13), 
          legend.title = element_text(size = 14),
          strip.text = element_text(size = 14)) +
    NULL
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

#Calculate total area
isol_migration$area_cm2 <- pi*isol_migration$Width_cm/2*isol_migration$Height_cm/2

#Calculate k
isol_migration$area_k <- log(isol_migration$area_cm2/(0.01*pi))/
  isol_migration$time_since_inoc

#Calculate relative values to same-day ancestor
ancestors <- isol_migration[isol_migration$Isol == "Anc", ]

  #For total area
isol_migration$relative_area <-
  isol_migration$area_cm2/ancestors$area_cm2[
    match(as.Date(isol_migration$end_timestamp),
          as.Date(ancestors$end_timestamp))]

  #For k
isol_migration$relative_k <-
  isol_migration$area_k/ancestors$area_k[
    match(as.Date(isol_migration$end_timestamp),
          as.Date(ancestors$end_timestamp))]

#Plot data
my_facet_labels <- c("7x" = "Weak Phage", 
                     "125" = "Strong Phage",
                     "C" = "Control", "G" = "Global", "L" = "Local",
                     "A" = "WT")

if (make_statplots) {
  #Total area
  ggplot(isol_migration, 
         aes(x = Pop, y = area_cm2, 
             color = Treat, fill = Treat)) +
    geom_point(position = position_dodge(0.5), alpha = 0.6,
               size = 2) +
    facet_nested(~Proj+Treat, labeller = labeller(Proj = my_facet_labels,
                                                  Treat = my_facet_labels)) +
    theme_bw() + 
    labs(y = "Area of Growth (cm^2)", x = "Population") +
    # scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
    #                    labels = c("Control", "Local", "Global"),
    #                    values = my_cols[c(8, 2, 6)]) +
    # scale_fill_manual(name = "Treatment", breaks = c("C", "L", "G"),
    #                   labels = c("Control", "Local", "Global"),
    #                   values = my_cols[c(8, 2, 6)]) +
    theme(legend.position = "none") +
    NULL
  
  #Relative area
  ggplot(isol_migration[isol_migration$Isol != "Anc", ], 
         aes(x = Pop, y = relative_area, 
             color = Treat, fill = Treat)) +
    geom_point(position = position_dodge(0.5), alpha = 0.6,
               size = 2) +
    facet_nested(~Proj+Treat, labeller = labeller(Proj = my_facet_labels,
                                                  Treat = my_facet_labels)) +
    theme_bw() + 
    labs(y = "Total Area Relative to Ancestor",
         x = "Population") +
    geom_hline(yintercept = 1, lty = 2) +
    scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                       labels = c("Control", "Local", "Global"),
                       values = my_cols[c(8, 2, 6)]) +
    scale_fill_manual(name = "Treatment", breaks = c("C", "L", "G"),
                      labels = c("Control", "Local", "Global"),
                      values = my_cols[c(8, 2, 6)]) +
    theme(legend.position = "none") +
    NULL
  
  #K values
  ggplot(isol_migration, 
         aes(x = Pop, y = area_k, 
             color = Treat, fill = Treat)) +
    geom_point(position = position_dodge(0.5), alpha = 0.6,
               size = 2) +
    facet_nested(~Proj+Treat, labeller = labeller(Proj = my_facet_labels,
                                                  Treat = my_facet_labels)) +
    theme_bw() + 
    labs(y = "Total Growth Parameter",
         x = "Population") +
    # scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
    #                    labels = c("Control", "Local", "Global"),
    #                    values = my_cols[c(8, 2, 6)]) +
    # scale_fill_manual(name = "Treatment", breaks = c("C", "L", "G"),
    #                   labels = c("Control", "Local", "Global"),
    #                   values = my_cols[c(8, 2, 6)]) +
    theme(legend.position = "none") +
    NULL
  
  #Relative K
  tiff("./Output_figures/Isol_migration.tiff",
       width = 6, height = 4, units = "in", res = 300)
  ggplot(isol_migration[isol_migration$Isol != "Anc", ], 
         aes(x = Pop, y = relative_k, 
             color = Treat, fill = Treat)) +
    geom_point(position = position_dodge(0.5), alpha = 0.6,
               size = 2) +
    facet_nested(~Proj+Treat, labeller = labeller(Proj = my_facet_labels,
                                                  Treat = my_facet_labels)) +
    theme_bw() + 
    labs(y = "Total Growth Parameter Relative to Ancestor",
         x = "Population") +
    geom_hline(yintercept = 1, lty = 2) +
    scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                       labels = c("Control", "Local", "Global"),
                       values = my_cols[c(8, 2, 6)]) +
    scale_fill_manual(name = "Treatment", breaks = c("C", "L", "G"),
                      labels = c("Control", "Local", "Global"),
                      values = my_cols[c(8, 2, 6)]) +
    theme(legend.position = "none") +
    NULL
  dev.off()
  
  tiff("./Output_figures/Isol_migration_tall.tiff",
       width = 5, height = 4, units = "in", res = 300)
  ggplot(isol_migration[isol_migration$Isol != "Anc", ], 
         aes(x = Pop, y = relative_k, 
             color = Treat, fill = Treat)) +
    geom_point(position = position_dodge(0.5), alpha = 0.6,
               size = 2) +
    facet_grid(Proj~Treat, labeller = labeller(Proj = my_facet_labels,
                                                  Treat = my_facet_labels),
               scales = "free_y") +
    theme_bw() + 
    labs(y = "Total Growth Parameter Relative to Ancestor",
         x = "Population") +
    geom_hline(yintercept = 1, lty = 2) +
    scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                       labels = c("Control", "Local", "Global"),
                       values = my_cols[c(8, 2, 6)]) +
    scale_fill_manual(name = "Treatment", breaks = c("C", "L", "G"),
                      labels = c("Control", "Local", "Global"),
                      values = my_cols[c(8, 2, 6)]) +
    theme(legend.position = "none") +
    NULL
  dev.off()
}

#Summarize for later inclusion w/ gc data
isol_migration_temp <- isol_migration[!is.na(isol_migration$relative_k), ]
isol_migration_temp <- group_by(isol_migration_temp,
                           Proj, Pop, Treat)
isol_migr_sum <- summarize(isol_migration_temp,
                           relative_k_avg = mean(relative_k))

## Isolate resistance ----
resis_data <- read.csv("./Clean_Data/Isolate_resistance.csv",
                       stringsAsFactors = F)

#Reorder projects
resis_data$Proj <- factor(resis_data$Proj,
                          levels = c("7x", "125"))

#Reorder treatments
resis_data$Treat <- factor(resis_data$Treat,
                           levels = c("Anc", "C", "L", "G"))

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

if (make_statplots) {
  #Make plot with both old and new approach
  ggplot(resis_data[resis_data$Treat != "Anc", ], 
         aes(x = Treat, y = EOP, color = Pop,
             shape = bd, group = Pop)) +
    facet_grid(Proj~approach) +
    geom_point(position = position_dodge(width = 0.5),
               alpha = 0.7) +
    scale_y_continuous(trans = "log10") +
    theme_bw()
}

#Calculate EOP limit & adjust values below limit
eop_limit <- max(resis_data$EOP[resis_data$bd &
                                  resis_data$approach == "new"])
resis_data$EOP[resis_data$EOP < eop_limit] <- eop_limit

#Assign facet labels
my_facet_labels <- c("7x" = "Weak Phage", 
                     "125" = "Strong Phage",
                     "C" = "Control", "G" = "Global", "L" = "Local",
                     "A" = "WT")

if (make_statplots) {
  #Nice plot
  tiff("./Output_figures/Isol_resis.tiff", width = 6, height = 4,
       units = "in", res = 300)
  ggplot(resis_data[resis_data$Treat != "Anc" &
                      resis_data$approach == "new", ],
         aes(x = Pop, y = EOP, 
             color = Treat, fill = Treat,
             shape = bd)) +
    facet_nested(~Proj+Treat, labeller = labeller(Proj = my_facet_labels,
                                                  Treat = my_facet_labels)) +
    geom_point(aes(size = bd, alpha = bd)) +
    scale_size_manual(values = c(2, 2.5)) +
    scale_alpha_manual(values = c(0.6, 1)) +
    scale_y_continuous(trans = "log10",
                       breaks = 10**(c(0, -2, -4, -6)),
                       labels = c(1, expression(10^-2), expression(10^-4),
                                  expression(10^-6))) +
    theme_bw() +
    geom_hline(yintercept = 1, lty = 2) +
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
    NULL
  dev.off()
  
  tiff("./Output_figures/Isol_resis_tall.tiff", width = 5, height = 4,
       units = "in", res = 300)
  ggplot(resis_data[resis_data$Treat != "Anc" &
                      resis_data$approach == "new", ],
         aes(x = Pop, y = EOP, 
             color = Treat, fill = Treat,
             shape = bd)) +
    facet_grid(Proj~Treat, labeller = labeller(Proj = my_facet_labels,
                                                  Treat = my_facet_labels)) +
    geom_point(aes(size = bd, alpha = bd)) +
    scale_size_manual(values = c(2, 2.5)) +
    scale_alpha_manual(values = c(0.6, 1)) +
    scale_y_continuous(trans = "log10",
                       breaks = 10**(c(0, -2, -4, -6)),
                       labels = c(1, expression(10^-2), expression(10^-4),
                                  expression(10^-6))) +
    theme_bw() +
    geom_hline(yintercept = 1, lty = 2) +
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
    NULL
  dev.off()
}

#Summarize for later inclusion w/ gc data
resis_data_temp <- resis_data[resis_data$approach == "new", ]
resis_data_temp <- group_by(resis_data_temp,
                            Proj, Pop, Treat)
resis_data_sum <- summarize(resis_data_temp,
                            EOP_avg = mean(EOP),
                            EOP_bd = any(bd))

##Isolate growth curves: read & find peaks ----

#Read data
gc_data <- read.csv("./Clean_Data/Isolate_growth_curves.csv",
                    header = T, stringsAsFactors = F)

#For ease of downstream analysis, for now we'll recode the media
# into simply "Original" and "Rich"
# Keeping in mind that for each project, "Orig" and "Rich" mean different
# medias
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

#reorder
gc_data <- gc_data[order(gc_data$uniq_well, gc_data$Time_s), ]

#Add number (for easy reference)
gc_data$uniq_well_num <- match(gc_data$uniq_well, unique(gc_data$uniq_well))

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
gc_data$sm_loess_25k <- NA
gc_data$sm_loess_3600 <- NA
for (my_well in unique(gc_data$uniq_well)) {
  #(leaving out the first half hour of data)
  my_rows <- which(gc_data$uniq_well == my_well &
                     gc_data$Time_s > 1800)
  
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
}

#Calculate growth per hour from loess curve
gc_data$deriv_sm_loess_25k <- calc_deriv(gc_data$sm_loess_25k,
                                     subset_by = gc_data$uniq_well,
                                     time = gc_data$Time_s,
                                     time_normalize = 3600)

#Calculate per capita growth per hour from loess curve
gc_data$percap_deriv_sm_loess_25k <- calc_deriv(gc_data$sm_loess_25k,
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

if (make_curveplots) {
  #View samples of original & smoothed curves
  # as well as derivatives (per cap & not) of both orig and smoothed curves
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
                    uniq_well, uniq_well_num)

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
first_min_time_window <- 7200
max_percap_time_window <- 21600
max_gr_rate_time_window <- 10800
pseudo_K_time_window <- 7200

#Summarize data (note: message is for grouping of **output**)
gc_summarized <- dplyr::summarize(gc_data,
  num_nas = sum(is.na(sm_loess_3600)), #we know all the nas are at start
  #Find the first minima in total density
   first_min_index = (find_local_extrema(sm_loess_3600,
                                        return_maxima = FALSE,
                                        width_limit = (first_min_time_window/
                                                         (Time_s[2]-Time_s[1])) + 1,
                                        na.rm = T,
                                        remove_endpoints = FALSE)[1]+num_nas),
   first_min = sm_loess_3600[first_min_index],
   first_min_time = Time_s[first_min_index],
  #find peaks in per capita growth rate
  max_percap_index =
    #first find all peaks
   (find_local_extrema(percap_deriv_sm_loess_25k,
                                        return_minima = FALSE,
                                        width_limit = (max_percap_time_window/
                                                         (Time_s[2]-Time_s[1])) + 1,
                                        na.rm = T,
                                        remove_endpoints = F)[
     #But save/use the first one that follows the minimum density
     match(TRUE, find_local_extrema(percap_deriv_sm_loess_25k,
                                    return_minima = FALSE,
                                    width_limit = (max_percap_time_window/
                                                     (Time_s[2]-Time_s[1])) + 1,
                                    na.rm = T,
                                    remove_endpoints = F) >= first_min_index)]+num_nas),
  max_percap_gr_rate = percap_deriv_sm_loess_25k[max_percap_index],
  max_percap_gr_time = Time_s[max_percap_index],
  max_percap_gr_dens = sm_loess_25k[max_percap_index],
  max_percap_gr_timesincemin = max_percap_gr_time - first_min_time,

  #find the local minimas in total grow rate (slope of total density)
  #(which is the point when the diauxic shift occurs)
  pseudo_K_index = (find_local_extrema(deriv_sm_loess_25k,
                                      return_maxima = FALSE,
                                      width_limit = (pseudo_K_time_window/
                                                       (Time_s[2]-Time_s[1])) + 1,
                                      na.rm = T,
                                      remove_endpoints = T)[1]+num_nas),
  pseudo_K = sm_loess_25k[pseudo_K_index],
  pseudo_K_time = Time_s[pseudo_K_index],
  pseudo_K_deriv = deriv_sm_loess_25k[pseudo_K_index],
  pseudo_K_timesincemin = pseudo_K_time - first_min_time,
  pseudo_K_timesince_maxpercap = pseudo_K_time - max_percap_gr_time
)

#Change to data frame for cleanliness
gc_summarized <- as.data.frame(gc_summarized)

#Make output plots for problematic wells
if (make_curveplots) {
  wells_check <- c("2017-B_7x_C_L_B_1_Orig",
                   "2017-A_7x_Anc_Anc_Anc_1_Orig",
                   "2017-A_7x_B_L_A_1_Rich",
                   "2017-B_7x_C_L_B_1_Orig",
                   "2017-C_7x_B_C_C_1_Rich",
                   "2017-C_7x_B_C_C_2_Rich",
                   "2017-C_7x_C_C_C_1_Orig",
                   "2017-C_7x_C_L_C_2_Rich",
                   "2017-C_7x_D_G_C_1_Orig",
                   "2017-C_7x_D_G_C_1_Rich",
                   "2017-C_7x_D_G_C_2_Orig",
                   "2017-C_7x_D_G_C_2_Rich",
                   "2017-C_7x_E_G_C_1_Rich",
                   "2017-C_7x_E_G_C_2_Rich",
                   "2017-E_7x_B_C_E_1_Orig",
                   "2017-E_7x_C_C_E_1_Rich",
                   "2017-E_7x_C_C_E_2_Rich",
                   "2019-09-10_125_B_C_A_1_Orig",
                   "2019-09-10_125_B_C_A_2_Orig",
                   "2019-09-10_125_B_G_A_1_Orig",
                   "2019-09-12_125_B_C_D_1_Orig",
                   "2019-09-12_125_B_C_D_2_Orig",
                   "2019-09-12_125_B_G_D_1_Orig",
                   "2019-09-12_125_D_L_D_1_Rich",
                   "2019-09-12_125_D_L_D_2_Rich",
                   "2019-09-13_125_B_C_E_1_Rich",
                   "2019-09-13_125_B_C_E_2_Orig",
                   "2019-09-13_125_C_C_E_1_Rich",
                   "2019-09-13_125_D_L_E_1_Rich"
  )
  for (my_well in wells_check) {
    tiff(filename = paste("./Challenging_growth_curve_plots/", my_well, ".tiff", sep = ""),
         width = 5, height = 10, units = "in", res = 300)
    my_rows <- which(gc_data$uniq_well == my_well)
    print(cowplot::plot_grid(
      ggplot(data = gc_data[my_rows, ],
             aes(x = Time_s, y = cfu_ml)) +
        geom_line(color = "red", lwd = 1, alpha = 0.5) +
        geom_line(aes(x = Time_s, y = sm_loess_25k),
                  color = "blue", lwd = 1, alpha = 0.5) +
        ggtitle(gc_data[my_rows[1], "uniq_well"]) +
        #Add point for first minima
        geom_point(data = gc_summarized[gc_summarized$uniq_well == my_well, ],
                   aes(x = first_min_time, y = first_min),
                   color = "green", size = 3) +
        NULL,
      ggplot(data = gc_data[my_rows, ],
             aes(x = Time_s, y = deriv_sm_loess_25k)) +
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
             aes(x = Time_s, y = percap_deriv_sm_loess_25k)) +
        geom_line(color = "blue") +
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
    tiff(filename = paste("./Growth_curve_plots/", my_well, ".tiff", sep = ""),
         width = 5, height = 10, units = "in", res = 300)
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
        #Add point for first minima
        geom_point(data = gc_summarized[gc_summarized$uniq_well == my_well, ],
                   aes(x = first_min_time, y = first_min),
                   color = "green", size = 3) +
        scale_y_continuous(trans = "log10") +
        NULL,
      ggplot(data = gc_data[my_rows, ],
             aes(x = Time_s, y = deriv_sm_loess_25k)) +
        geom_line(color = "blue") +
        #Add point for pseudo K
        geom_point(data = gc_summarized[gc_summarized$uniq_well == my_well, ],
                   aes(x = pseudo_K_time, y = pseudo_K_deriv),
                   color = "green", size = 3) +
        NULL,
      ggplot(data = gc_data[my_rows, ],
             aes(x = Time_s, y = percap_deriv_sm_loess_25k)) +
        geom_line(color = "blue") +
        #Add point for max growth rate
        geom_point(data = gc_summarized[gc_summarized$uniq_well == my_well, ],
                   aes(x = max_percap_gr_time, y = max_percap_gr_rate),
                   color = "green", size = 3) +
        NULL,
      ncol = 1, align = "v"))
    dev.off()
  }
}

#Note that orig media has strong diauxic shift
# while rich media has little to none

#Looked through all the plots, some of them no diauxic shift is detected
#But that's okay because for most of them the other tech rep had one
#2017 C 7x clc 2 rich (7x clc 1 rich is good)
#7x dgc 2 orig (7x dgc 1 orig is good)
#125 bce 1 rich (125 bce 2 rich is good)
#125 cce 1 rich (125 cce 2 rich is good)

#One isolate (7x dgc 1&2 rich) both replicate wells had no diauxic shift 
#detected
#Look at plots of 7x dgc 1&2 rich
if(make_curveplots) {
  test <- gc_data[gc_data$uniq_well == "2017-C_7x_D_G_C_2_Rich", ]
  test$deriv <- c(test$OD600[2:nrow(test)]-test$OD600[1:(nrow(test)-1)], NA)
  test$deriv_pc <- test$deriv/test$OD600
  ggplot(test, aes(x = Time_s, y = OD600)) + geom_point()
  ggplot(test, aes(x = Time_s, y = deriv)) + geom_point()
  ggplot(test, aes(x = Time_s, y = deriv_pc)) + geom_point()
}
#Unfortunately 7x dgc 1 & 2 rich diauxic shifts we'll have to leave tossed
# (other isols from that pop will have to suffice, since we fail to detect 
# any diauxic shift in those curves, and looking at the curves shows pretty 
# clearly no sensitivity changes will be sufficient

##Isolate growth curves: Fit logistic curve to data ----
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

#Seeing what the Baranyi curve looks like
if (F) {
  temp_t <- seq(from = 0, to = 100000, by = 1)
  temp_d <- baranyi_func(r = 1.5, k = 10**10, v = 1, q0 = 0.25, m = 2,
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
                                  "fit2_err" = as.numeric(NA)))
for (sum_row in 1:nrow(gc_summarized)) {
  if (!is.na(gc_summarized$pseudo_K[sum_row])) {
    my_well <- gc_summarized$uniq_well[sum_row]
    myrows1 <- which(gc_data$uniq_well == my_well &
                      gc_data$Time_s <= gc_summarized$pseudo_K_time[sum_row] &
                      gc_data$Time_s >= gc_summarized$max_percap_gr_time[sum_row])
    myrows2 <- which(gc_data$uniq_well == my_well &
                       gc_data$Time_s <= (gc_summarized$pseudo_K_time[sum_row] + 45*60) &
                       gc_data$Time_s > 0)
    temp1 <- optim(par = c("logk" = log10(gc_summarized$pseudo_K[sum_row]),
                          "d0" = gc_summarized$max_percap_gr_dens[sum_row],
                          "r" = gc_summarized$max_percap_gr_rate[sum_row],
                          "delta" = 0),
                  fn = logis_fit_err,
                  dens_vals = gc_data$sm_loess_25k[myrows1],
                  t_vals = gc_data$Time_s[myrows1],
                  t_offset = gc_summarized$max_percap_gr_time[sum_row],
                  method = "BFGS")
    temp2 <- optim(par = c("logk" = log10(gc_summarized$pseudo_K[sum_row]),
                           "logd0" = log10(gc_data$cfu_ml[myrows2[2]]),
                           "r" = gc_summarized$max_percap_gr_rate[sum_row],
                           "v" = 1, 
                           #"m" = gc_summarized$max_percap_gr_rate[sum_row],
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
                   method = "CG")
    gc_summarized[sum_row, 
                  c("fit_r", "fit_k", "fit_d0", "fit_delta", "fit_err",
                    "fit2_r", "fit2_k", "fit2_v", "fit2_q0",
                    "fit2_m", "fit2_d0", "fit2_err")] <-
      data.frame("fit_r" = temp1$par["r"], 
                 "fit_k" = 10**temp1$par["logk"],
                 "fit_d0" = temp1$par["d0"], 
                 "fit_delta" = temp1$par["delta"],
                 "fit_err" = temp1$value,
                 "fit2_r" = temp2$par["r"], 
                 "fit2_k" = 10**temp2$par["logk"],
                 "fit2_v" = temp2$par["v"], 
                 "fit2_q0" = temp2$par["q0"],
                 #"fit2_m" = temp2$par["m"],
                 "fit2_m" = NA,
                 "fit2_d0" = 10**temp2$par["logd0"],
                 "fit2_err" = temp2$value)
  }
}

dir.create("./Growth_curve_plots_fits", showWarnings = F)
if(make_curveplots) {
  #Make plots of fits
  for (sum_row in 1:nrow(gc_summarized)) {
    my_well <- gc_summarized$uniq_well[sum_row]
    t_vals <- gc_data$Time_s[gc_data$uniq_well == my_well]
    t_vals_hrs <- gc_data$Time_s[gc_data$uniq_well == my_well]/3600
    offset_hrs <- gc_summarized$max_percap_gr_time[sum_row]/3600
    pred_vals <- with(as.list(
      gc_summarized[sum_row, c("fit_r", "fit_k", "fit_d0", "fit_delta", "fit_err")]),
      fit_k/(1+(((fit_k-fit_d0)/fit_d0)*
                  exp(-fit_r*(t_vals_hrs-offset_hrs)))))
    
    png(paste("./Growth_curve_plots_fits/", my_well, ".png", sep = ""),
        width = 4, height = 4, units = "in", res = 300)
    print(ggplot(data = gc_data[gc_data$uniq_well == my_well, ],
                 aes(x = Time_s, y = sm_loess_25k)) +
            geom_line(lwd = 1.5) +
            geom_point(size = 0.5, aes(y = cfu_ml)) +
            geom_line(data.frame(Time_s = t_vals, pred_dens = pred_vals),
                      mapping = aes(x = Time_s, y = pred_dens),
                      color = "red") +
            scale_y_continuous(trans = "log10") +
            geom_vline(aes(xintercept = gc_summarized$max_percap_gr_time[
              gc_summarized$uniq_well == my_well]), lty = 2) +
            geom_vline(aes(xintercept = gc_summarized$pseudo_K_time[
              gc_summarized$uniq_well == my_well]), lty = 2) +
            NULL)
    dev.off()
  }
}

dir.create("./Growth_curve_plots_fits2", showWarnings = F)
if(make_curveplots) {
  #Make plots of fits
  for (sum_row in 1:nrow(gc_summarized)) {
    my_well <- gc_summarized$uniq_well[sum_row]
    t_vals <- gc_data$Time_s[gc_data$uniq_well == my_well]
    pred_vals2 <- baranyi_func2(r = gc_summarized[sum_row, "fit2_r"],
                               k = gc_summarized[sum_row, "fit2_k"],
                               v = gc_summarized[sum_row, "fit2_v"],
                               q0 = gc_summarized[sum_row, "fit2_q0"],
                               #m = gc_summarized[sum_row, "fit2_m"],
                               d0 = gc_summarized[sum_row, "fit2_d0"],
                               t_vals = t_vals)
    start_vals <- baranyi_func2("k" = gc_summarized$pseudo_K[sum_row],
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
            geom_line(data.frame(Time_s = t_vals, pred_dens = pred_vals2),
                      mapping = aes(x = Time_s, y = pred_dens),
                      color = "red") +
            geom_line(data.frame(Time_s = t_vals, pred_dens = start_vals),
                      mapping = aes(x = Time_s, y = pred_dens),
                      color = "blue", alpha = 0.2) +
            scale_y_continuous(trans = "log10") +
            geom_vline(aes(xintercept = min(gc_data$Time_s[gc_data$Time_s > 0])), 
                       lty = 2) +
            geom_vline(aes(xintercept = (gc_summarized$pseudo_K_time[
              gc_summarized$uniq_well == my_well]+45*60)), lty = 2) +
            geom_vline(aes(xintercept = (gc_summarized$pseudo_K_time[
              gc_summarized$uniq_well == my_well])), lty = 3, alpha = 0.5) +
            theme_bw() +
            NULL)
    dev.off()
  }
}

#New new bad fits (added 45 mins to endtime of data)
# (quite a few need to be cut off at first min but not using loess)
# (because there's an initial drop in density)
#13 (need to drop first 3 timepoints, not just first 1)
#26
#38
#148 - nofit
#154 nofit
#155 nofit
#156 nofit
#185 really bad
#456 really bad
#458 really bad
#484 no fit
#496 nofit
#516 really bad
#518 really bad


#New bad fits
#13
#25
#26
#38
#148 - no fit
#154, 155, 156 - no fit
#157, 158
#185 really bad
#...
#518
#519
#523

#Bad fits (one # is meh, usually k too high; two ## is really bad)
#2019-09-10-125-b-g-a-2-orig
#2019-09-10-125-c-g-a-1-orig
#2019-09-10-125-c-g-a-2-orig
#2019-09-10-125-e-l-a-1-orig
#2019-09-11-125-c-g-b-1-orig
#2019-09-11-125-c-g-b-2-orig
#2019-09-12-125-a-c-d-1-orig
#2019-09-12-125-a-c-d-2-orig
#2019-09-12-125-anc-anc-anc-3-orig
#2019-09-12-125-b-g-d-1-orig
#2019-09-12-125-b-g-d-2-orig
#2019-09-12-125-c-g-d-1-orig
##2019-09-12-125-d-l-d-1-rich
##2019-09-12-125-d-l-d-2-rich
#2019-09-12-125-anc-anc-anc-1-orig
#2019-09-12-125-anc-anc-anc-2-orig
#2019-09-12-125-anc-anc-anc-3-orig
#2019-09-13-125-c-g-e-1-orig
#2019-09-13-125-c-g-e-2-orig
#2019-09-13-125-c-l-e-1-orig
#2019-09-13-125-c-l-e-2-orig
#2019-09-13-125-D-g-e-2-rich
##2019-09-13-125-d-l-e-1-rich
##2019-09-13-125-d-l-e-2-rich
#2019-09-13-125-e-l-e-1-orig


#Check how fit results compare to local extrema results
if (make_statplots) {
  ggplot(data = gc_summarized,
         aes(x = max_percap_gr_rate, y = fit2_r)) +
    geom_point() +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10")
  ggplot(data = gc_summarized,
         aes(x = pseudo_K, y = fit2_k, color = Proj, shape = Media)) +
    geom_point() +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10")
  ggplot(data = gc_summarized,
         aes(x = fit2_d0, y = first_min)) +
    geom_point() +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(trans = "log10")
  hist(gc_summarized$fit_err)
}     

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
  
}

#Isolate growth curves: summarize reps into isols, view variable data & distributions ----

#Summarize replicate wells
gc_summarized <- group_by(gc_summarized, Date, Proj, Pop, Treat,
                          Isol, Media)
gc_sum_isols <- summarize_at(gc_summarized,
                            .funs = c(avg = mean, sd = sd),
                            .vars = c(
                              "first_min",
                              "max_percap_gr_rate",
                              "max_percap_gr_dens",
                              "max_percap_gr_timesincemin",
                              "pseudo_K",
                              "pseudo_K_timesincemin",
                              "pseudo_K_timesince_maxpercap",
                              "fit_r",
                              "fit_k",
                              "fit_d0",
                              "fit2_r", "fit2_k", "fit2_v", 
                              "fit2_q0", "fit2_m", "fit2_d0"))
gc_sum_isols <- as.data.frame(gc_sum_isols)

#Take a look at the standard deviations between replicate wells
# Just raw sd vals (w/ red line for mean avg value)
if (make_statplots) {
  for (var in c("first_min_", "first_min_time_", 
                "max_percap_gr_rate_", "max_percap_gr_time_", 
                "max_percap_gr_dens_", 
                "max_percap_gr_timesincemin_",
                "pseudo_K_", "pseudo_K_time_", 
                "pseudo_K_timesincemin_", 
                "pseudo_K_timesince_maxpercap_")) {
    my_sd <- gc_sum_isols[, paste(var, "sd", sep = "")]
    my_avg <- mean(gc_sum_isols[, paste(var, "avg", sep = "")])
    hist(my_sd, main = var, 
         xlim = c(min(my_sd, my_avg, na.rm = T), max(my_sd, my_avg, na.rm = T)))
    abline(v = my_avg, col = "red", lwd = 2)
  }
}
# Sd vals divided by matching avg value (so 1 is the reference)
if (make_statplots) {
  for (var in c("first_min_", "first_min_time_", 
                "max_percap_gr_rate_", "max_percap_gr_time_", 
                "max_percap_gr_dens_", 
                "max_percap_gr_timesincemin_",
                "pseudo_K_", "pseudo_K_time_", 
                "pseudo_K_timesincemin_", 
                "pseudo_K_timesince_maxpercap_")) {
    my_sd <- gc_sum_isols[, paste(var, "sd", sep = "")]
    my_avg <- gc_sum_isols[, paste(var, "avg", sep = "")]
    hist(my_sd/my_avg, main = var, xlim = c(0, max(my_sd/my_avg, 1, na.rm = T)))
    abline(v = 1, col = "red", lwd = 2)
  }
}

#Generally, sd's between reps are small relative to the values themselves

#View all the isols by variable
my_facet_labels <- c("7x" = "Weak Phage", 
                     "125" = "Strong Phage",
                     "C" = "Control", "G" = "Global", "L" = "Local",
                     "Anc" = "Ancstr",
                     "Rich" = "Rich Media", "Orig" = "Original Media")
if (make_statplots) {
  my_vars <- c("first_min_", 
               #"first_min_time_", 
               "max_percap_gr_rate_", 
               #"max_percap_gr_time_", 
               "max_percap_gr_dens_", 
               "max_percap_gr_timesincemin_",
               "pseudo_K_", 
               #"pseudo_K_time_", 
               "pseudo_K_timesincemin_",
               "pseudo_K_timesince_maxpercap_",
               "fit_r_", "fit_k_", "fit_d0_",
               "fit2_r_", "fit2_k_", "fit2_v_", 
               "fit2_q0_", "fit2_m_", "fit2_d0_")
  for (var_root in my_vars) {
    i <- which(var_root == my_vars)
    var_name <- c("First minimum density (cfu/mL)",
                  "Maximum per-capita growth rate",
                  "Density at maximum per-capita growth rate",
                  "Time until maximum per-capita growth rate",
                  "Density at diauxic shift (cfu/mL)",
                  "Time until diauxic shift (from min)",
                  "Time until diauxic shift (from max percap)",
                  "Fit r", "Fit carrying capacity", "Fit init density",
                  "Fit 2 r", "Fit 2 k", "Fit 2 v", "Fit 2 q0",
                  "Fit 2 m", "fit 2 d0")[i]
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
            theme(legend.position = "none")
    )
    dev.off()
  }
}

#Note however how some curves have super high sds between
# rep wells
#Particularly in vars: first min, percap dens, percap time, K time

#gc_sum_isols[which.max(gc_sum_isols$first_min_sd), ]
#gc_sum_isols[which.max(gc_sum_isols$max_percap_gr_rate_sd), ]
# which.max(gc_sum_isols$pseudo_K_sd)

#After checking out the above cases (and having already
# manually inspected all the curves)
# I'm satisfied that the rare cases where repwells disagree
# strongly are either cases where no algorithm could
# assign differently because of the shape of the curves
# or where the data itself is strangely different
# between the wells

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
              "fit_r_avg", "fit_k_avg", "fit_d0_avg")) {
  new_var <- paste(var, "_rel", sep = "")
  gc_sum_isols[, new_var] <- gc_sum_isols[, var] -
    ancestors[match(paste(gc_sum_isols$Date, gc_sum_isols$Media), 
                    paste(ancestors$Date, ancestors$Media)), var]
}

#Should density measures be relative too?

#Now view the relative variables
my_facet_labels <- c("7x" = "Weak Phage", 
                     "125" = "Strong Phage",
                     "C" = "Control", "G" = "Global", "L" = "Local",
                     "A" = "WT",
                     "Rich" = "Rich Media", "Orig" = "Original Media")
if (make_statplots) {
  for (i in 1:10) {
    var_root <- c("first_min_", 
                  "max_percap_gr_rate_", 
                  "max_percap_gr_dens_", 
                  "max_percap_gr_timesincemin_",
                  "pseudo_K_", 
                  "pseudo_K_timesincemin_",
                  "pseudo_K_timesince_maxpercap_",
                  "fit_r_", "fit_k_", "fit_d0_")[i]
    var <- paste(var_root, "avg_rel", sep = "")
    var_name <- c("Relative first minimum density",
                  "Relative maximum per-capita growth rate",
                  "Relative density at maximum per-capita growth rate",
                  "Relative time until maximum per-capita growth rate",
                  "Relative density at diauxic shift",
                  "Relative time until diauxic shift (from min)",
                  "Relative time until diauxic shift (from max percap)",
                  "Relative Fit r", "Relative Fit carrying capacity", 
                  "Relative Fit init density")[i]
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
                            "fit_r_avg_rel", "fit_k_avg_rel", "fit_d0_avg_rel"))
gc_sum_pops <- as.data.frame(gc_sum_pops)

#View population-summarized mean data (median below)
my_facet_labels <- c("7x" = "Weak Phage", 
                     "125" = "Strong Phage",
                     "C" = "Control", "G" = "Global", "L" = "Local",
                     "A" = "WT",
                     "Rich" = "Rich Media", "Orig" = "Original Media")
if (make_statplots) {
  for (var_root in c("first_min_avg_rel", 
                     # "first_min_time_avg", 
                     "max_percap_gr_rate_avg_rel", 
                     # "max_percap_gr_time_avg", 
                     "max_percap_gr_dens_avg_rel", 
                     "max_percap_gr_timesincemin_avg_rel",
                     "pseudo_K_avg_rel", 
                     # "pseudo_K_time_avg", 
                     "pseudo_K_timesincemin_avg_rel",
                     "pseudo_K_timesince_maxpercap_avg_rel"
  )) {
    var <- paste(var_root, "_avg", sep = "")
    var_sd <- paste(var_root, "_sd", sep = "")
    tiff(paste("./Growth_curve_variables_plots_pops/", var, ".tiff", sep = ""),
         width = 10, height = 10, units = "in", res = 300)
    print(ggplot(data = gc_sum_pops[gc_sum_pops$Pop != "Anc", ],
                 aes(x = Treat, y = get(var), group = Pop,
                     color = Treat)) +
            geom_point(position = position_dodge(0.1),
                       size = 5, alpha = 0.8) +
            scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                               labels = c("Control", "Local", "Global"),
                               values = my_cols[c(8, 2, 6)]) +
            facet_grid(Proj ~ Media, scales = "free_y",
                       labeller = labeller(Proj = my_facet_labels,
                                           Media = my_facet_labels)) +
            geom_hline(yintercept = 1, lty = 2) +
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
if (make_statplots) {
  for (var_root in c("first_min_avg_rel", 
                     "max_percap_gr_rate_avg_rel", 
                     "max_percap_gr_dens_avg_rel", 
                     "max_percap_gr_timesincemin_avg_rel",
                     "pseudo_K_avg_rel", 
#                     "pseudo_K_timesincemin_avg_rel"
                     "pseudo_K_timesince_maxpercap_avg_rel"
                     )) {
    var <- paste(var_root, "_med", sep = "")
    var_sd <- paste(var_root, "_sd", sep = "")
    tiff(paste("./Growth_curve_variables_plots_pops/", var, ".tiff", sep = ""),
         width = 10, height = 10, units = "in", res = 300)
    print(ggplot(data = gc_sum_pops[gc_sum_pops$Pop != "Anc", ],
                 aes(x = Treat, y = get(var), group = Pop,
                     color = Treat)) +
            geom_point(position = position_dodge(0.1),
                       size = 5, alpha = 0.8) +
            geom_hline(yintercept = 1, lty = 2) +
            scale_color_manual(name = "Treatment", breaks = c("C", "L", "G"),
                               labels = c("Control", "Local", "Global"),
                               values = my_cols[c(8, 2, 6)]) +
            facet_grid(Proj ~ Media, scales = "free_y",
                       labeller = labeller(Proj = my_facet_labels,
                                           Media = my_facet_labels)) +
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

#Cast measurements in different medias into different columns
# (using population mean data)
gc_sum_pops <- as.data.table(gc_sum_pops)
gc_sum_pops_wide <- data.table::dcast(gc_sum_pops,
                           Proj+Pop+Treat ~ Media,
                           value.var = c("first_min_avg_rel_avg", 
                                         "max_percap_gr_rate_avg_rel_avg", 
                                         "max_percap_gr_dens_avg_rel_avg", 
                                         "max_percap_gr_timesincemin_avg_rel_avg", 
                                         "pseudo_K_avg_rel_avg", 
                                         "pseudo_K_timesincemin_avg_rel_avg",
                                         "pseudo_K_timesince_maxpercap_avg_rel_avg"))
gc_sum_pops_wide <- as.data.frame(gc_sum_pops_wide)

#Add in resistance & migration data
isol_data <- full_join(gc_sum_pops_wide, 
                              resis_data_sum)
isol_data <- full_join(isol_data, isol_migr_sum)

isol_data$EOP_avg <- log10(isol_data$EOP_avg)

#Rename for brevity
isol_data <- isol_data[, c(1:13, 16:20)]
colnames(isol_data)[4:18] <- c(
  "min_Orig", "min_Rich",
  "pc_rate_Orig", "pc_rate_Rich",
  "pc_dens_Orig", "pc_dens_Rich",
  "pc_time_Orig", "pc_time_Rich",
  "K_Orig", "K_Rich",
  "K_time_Orig", "K_time_Rich",
  "log(EOP)", "EOP_bd", "Agar_grow")

#Check correlations between variables
gc_var_cors_7x <- cor(isol_data[isol_data$Proj == "7x", 
                                       c(4:16, 18)])
gc_var_cors_125 <- cor(isol_data[isol_data$Proj == "125", 
                                        c(4:16, 18)])
write.csv(gc_var_cors_7x, "grow_curve_var_correlations_7x.csv")
write.csv(gc_var_cors_125, "grow_curve_var_correlations_125.csv")

# all vars are positive w/ ea other
# max percap rate neg w/ first min
# max percap dens pos w/ first min
# max percap dens pos w/ max percap timesincemin

#Make correlation figures
if (make_statplots) {
  tiff("./Output_figures/Weakphage_cors.tiff", width = 10, height = 10, units = "in", res = 300)
  #Make base figure
  p <- GGally::ggpairs(isol_data[isol_data$Treat != "Anc" &
                                       isol_data$Proj == "7x", ],
                  columns = c(4:16, 18),
                  lower = list(continuous = "smooth"),
                  upper = list(continuous = "smooth"),
                  ggplot2::aes(color = Treat, group = Proj),
                  title = "Weak Phage") +
    theme(strip.text = element_text(size = 7),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  #Change colors
  for (i in 1:p$nrow) {
    for (j in 1:p$ncol) {
      p[i, j] <- p[i, j] +
        scale_color_manual(breaks = c("C", "L", "G"),
                             values = my_cols[c(8, 2, 6)])
    }
  }
  print(p)
  dev.off()
  
  tiff("./Output_figures/Strongphage_cors.tiff", width = 10, height = 10, units = "in", res = 300)
  #Make base figure
  p <- GGally::ggpairs(isol_data[isol_data$Treat != "Anc" &
                                   isol_data$Proj == "125", ],
                       columns = c(4:16, 18),
                       lower = list(continuous = "smooth"),
                       upper = list(continuous = "smooth"),
                       ggplot2::aes(color = Treat, group = Proj),
                       title = "Strong Phage") +
    theme(strip.text = element_text(size = 7),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  #Change colors
  for (i in 1:p$nrow) {
    for (j in 1:p$ncol) {
      p[i, j] <- p[i, j] +
        scale_color_manual(breaks = c("C", "L", "G"),
                           values = my_cols[c(8, 2, 6)])
    }
  }
  print(p)
  dev.off()
  
  tiff("./Output_figures/Heatcors_weak.tiff", width = 10, height = 10, units = "in", res = 300)
  GGally::ggcorr(isol_data[isol_data$Treat != "Anc" &
                             isol_data$Proj == "7x", c(4:16, 18)],
                 nbreaks = 5,
                 hjust = 0.8,
                 layout.exp = 1.5) +
    ggplot2::labs(title = "Weak Phage")
  dev.off()
  
  tiff("./Output_figures/Heatcors_strong.tiff", width = 10, height = 10, units = "in", res = 300)
  GGally::ggcorr(isol_data[isol_data$Treat != "Anc" &
                             isol_data$Proj == "125", c(4:16, 18)],
                 nbreaks = 5,
                 hjust = 0.8,
                 layout.exp = 1.5) +
                   ggplot2::labs(title = "Strong Phage")
  dev.off()
}

##Isolate growth curves: run PCA (naively) ----
isol_data_pca_7x <- isol_data[isol_data$Proj == "7x", ]
isol_data_pca_125 <- isol_data[isol_data$Proj == "125", ]
#Don't need to scale because all values are relative to ancestor already
# isol_data_pca_7x[, 4:14] <- scale(isol_data_pca_7x[, 4:14])
# isol_data_pca_125[, 4:14] <- scale(isol_data_pca_125[, 4:14])
isol_princomp_7x <- princomp(isol_data_pca_7x[, c(4:16, 18)], cor = T, scores = T)
isol_princomp_125 <- princomp(isol_data_pca_125[, c(4:16, 18)], cor = T, scores = T)

isol_data_pca_7x <- cbind(isol_data_pca_7x, isol_princomp_7x$scores)
isol_data_pca_125 <- cbind(isol_data_pca_125, isol_princomp_125$scores)

loadings_7x <- 6*as.data.frame(isol_princomp_7x$loadings[])
loadings_125 <- 6*as.data.frame(isol_princomp_125$loadings[])
loadings_7x$varnames <- rownames(loadings_7x)
loadings_125$varnames <- rownames(loadings_125)

summary(isol_princomp_7x)
summary(isol_princomp_125)

if(make_statplots) {
  tiff("./Output_figures/weakphage_PCA.tiff", 
       width = 12, height = 10, units = "in", res = 300)
  ggplot(isol_data_pca_7x,
         aes(x = Comp.1, y = Comp.2)) +
    ggtitle("Weak Phage") +
    geom_segment(data = loadings_7x,
                 aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2),
                 arrow = arrow(length = unit(0.02, "npc")),
                 alpha = .8, lwd = 1, color = "gray50") +
    ggrepel::geom_text_repel(data = loadings_7x,
                             aes(x = Comp.1, y = Comp.2, label = varnames),
              size = 7, alpha = .8, color = "gray50") +
    geom_point(aes(color = Treat), size = 10, alpha = 0.7) +
    theme_bw() +
    labs(x = "PC1", y = "PC2") +
    theme(axis.title = element_text(size = 20),
          legend.title = element_text(size = 24),
          legend.text = element_text(size = 20))
  dev.off()

  tiff("./Output_figures/strongphage_PCA.tiff", 
       width = 12, height = 10, units = "in", res = 300)
  ggplot(isol_data_pca_125,
         aes(x = Comp.1, y = Comp.2)) +
    ggtitle("Strong Phage") +
    geom_segment(data = loadings_125,
                 aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2),
                 arrow = arrow(length = unit(0.02, "npc")),
                 alpha = .8, lwd = 1, color = "gray50") +
    ggrepel::geom_text_repel(data = loadings_125,
                             aes(x = Comp.1, y = Comp.2, label = varnames),
                             size = 7, alpha = .8, color = "gray50") +
    geom_point(aes(color = Treat), size = 10, alpha = 0.7) +
    theme_bw() +
    labs(x = "PC1", y = "PC2") +
    theme(axis.title = element_text(size = 20),
          legend.title = element_text(size = 24),
          legend.text = element_text(size = 20))
  dev.off()
}


#print(summary(isol_princomp), digits = 2)
print(isol_princomp_7x$loadings, cutoff = 0)
print(isol_princomp_125$loadings, cutoff = 0)

##Isolate growth curves: Check for normality ----

#Check for univariate normality
if (make_statplots) {
  for (var_root in c("first_min_avg_rel_avg_", 
                     "max_percap_gr_rate_avg_rel_avg_",
                     "max_percap_gr_dens_avg_rel_avg_",
                     "max_percap_gr_timesincemin_avg_rel_avg_",
                     "pseudo_K_avg_rel_avg_",
                     "pseudo_K_timesincemin_avg_rel_avg_"
  )) {
    for (media in c("Orig", "Rich")) {
      for (proj in unique(gc_sum_pops_wide$Proj)) {
        var <- paste(var_root, media, sep = "")
        # hist(as.numeric(gc_sum_pops_wide[gc_sum_pops_wide$Proj == proj, var]), 
        #      main = paste(proj, var))
        qqnorm(as.numeric(gc_sum_pops_wide[gc_sum_pops_wide$Proj == proj &
                                             gc_sum_pops_wide$Pop != "Anc", 
                                           var]), 
               main = paste(proj, var))
        qqline(as.numeric(gc_sum_pops_wide[gc_sum_pops_wide$Proj == proj &
                                             gc_sum_pops_wide$Pop != "Anc", 
                                           var]))
      }
    }
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
for (proj in unique(gc_sum_pops_wide$Proj)) {
    CSQPlot(gc_sum_pops_wide[gc_sum_pops_wide$Proj == proj, 
                             4:15],
            label = proj)
}

#Note that 125 is nearly multivariate normal
# while 7x is so far from multivariate normal no transformations
# will save it

#So we'll have to use non-parametric methods for
# MANOVA/ANOVA

#Luckily, discriminant analysis is not strongly dependent on
# multivariate normality, as long as we're not planning on using
# it to classify future observations

#Split out data into two separate projects
gc_sum_pops_wide_7x <- gc_sum_pops_wide[gc_sum_pops_wide$Proj == "7x", ]
gc_sum_pops_wide_125 <- gc_sum_pops_wide[gc_sum_pops_wide$Proj == "125", ]

##Isolate growth curves: Discriminant Analysis ----



#Run linear discriminant analysis
cols_to_use <- 6:15
gc_lda_7x <- lda(x = scale(gc_sum_pops_wide_7x[, cols_to_use]),
              grouping = gc_sum_pops_wide_7x$Treat)
gc_lda_125 <- lda(x = scale(gc_sum_pops_wide_125[, cols_to_use]),
                 grouping = gc_sum_pops_wide_125$Treat)

#Add discriminant function scores to gc_sum_pops_wide
# (by matrix multiplying the variables with the scaling matrix provided by lda)
gc_sum_pops_wide_7x <- cbind(gc_sum_pops_wide_7x,
                             scale(gc_sum_pops_wide_7x[, cols_to_use])%*%
                               gc_lda_7x$scaling)
gc_sum_pops_wide_125 <- cbind(gc_sum_pops_wide_125,
                              scale(gc_sum_pops_wide_125[, cols_to_use])%*%
                                gc_lda_125$scaling)

if(make_statplots) {
  #Make plots of LD1 and LD2
  ggplot(data = gc_sum_pops_wide_7x,
         aes(x = LD1, y = LD2, color = Treat)) +
    geom_point(size = 3) +
    ggtitle("7x") +
    theme_bw()
  ggplot(data = gc_sum_pops_wide_125,
         aes(x = LD1, y = LD2, color = Treat)) +
    geom_point(size = 3)  +
    ggtitle("125") +
    theme_bw()
  
  ggplot(data = gc_sum_pops_wide_7x,
         aes(x = Treat, y = LD1)) +
    geom_boxplot() +
    ggtitle("7x")
  ggplot(data = gc_sum_pops_wide_125,
         aes(x = Treat, y = LD1)) +
    geom_boxplot() +
    ggtitle("125")
}

#View loadings
gc_lda_7x
gc_lda_125

#7x
# percap rate Orig - left, down
# percap dens Rich - left, down
# percap time Rich - Right
# pseudo K time Orig - Left
# Global evolved higher percap rate and pseudo K time in Orig
# Control evolved lower percap rate & pseudo K time in Orig
# Local evolved lowest percap rate and pseudo K time in Orig
#125
# 

##Isolate growth curves: statistical tests ----

#With the original variables
# Want to test: whether treats are dift from ea other in ea proj
#               whether treats are dift from 1 (Anc value) in ea proj
nonpartest(max_percap_gr_rate_avg_rel_avg_Orig|
             max_percap_gr_rate_avg_rel_avg_Rich|
             max_percap_gr_dens_avg_rel_avg_Orig|
             max_percap_gr_dens_avg_rel_avg_Rich|
             max_percap_gr_timesincemin_avg_rel_avg_Orig|
             max_percap_gr_timesincemin_avg_rel_avg_Rich|
             pseudo_K_avg_rel_avg_Orig|
             pseudo_K_avg_rel_avg_Rich|
             pseudo_K_timesincemin_avg_rel_avg_Orig|
             pseudo_K_timesincemin_avg_rel_avg_Rich~Treat,
           gc_sum_pops_wide_7x[gc_sum_pops_wide_7x$Treat != "Anc", ],
           plots = F)
nonpartest(max_percap_gr_rate_avg_rel_avg_Orig|
             max_percap_gr_rate_avg_rel_avg_Rich|
             max_percap_gr_dens_avg_rel_avg_Orig|
             max_percap_gr_dens_avg_rel_avg_Rich|
             max_percap_gr_timesincemin_avg_rel_avg_Orig|
             max_percap_gr_timesincemin_avg_rel_avg_Rich|
             pseudo_K_avg_rel_avg_Orig|
             pseudo_K_avg_rel_avg_Rich|
             pseudo_K_timesincemin_avg_rel_avg_Orig|
             pseudo_K_timesincemin_avg_rel_avg_Rich~Treat,
           gc_sum_pops_wide_125[gc_sum_pops_wide_125$Treat != "Anc", ],
           plots = F)

#7x
nonpartest(LD1|LD2~Treat, 
           gc_sum_pops_wide_7x[gc_sum_pops_wide_7x$Treat != "Anc", ],
           permreps=1000, plots = F)

#125
nonpartest(LD1|LD2~Treat, 
           gc_sum_pops_wide_125[gc_sum_pops_wide_125$Treat != "Anc", ],
           permreps=1000, plots = F)
ssnonpartest(LD1|LD2~Treat, 
             gc_sum_pops_wide_125[gc_sum_pops_wide_125$Treat != "Anc", ],
             factors.and.variables = T)
            


  
  
  
  
  
  
  