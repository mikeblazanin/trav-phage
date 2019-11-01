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

##experimental evolution analysis
exper_evol_migr <- read.csv("./Clean_Data/Experimental_evolution_growth.csv")

#Drop points after T14
exper_evol_migr <- exper_evol_migr[exper_evol_migr$Timepoint <= 14, ]

#Calculate total area
exper_evol_migr$area_cm <- pi*exper_evol_migr$Width_cm/2*exper_evol_migr$Height_cm/2

#Make plot with all pops shown
ggplot(data = exper_evol_migr,
       aes(x = Timepoint, y = area_cm, group = paste(Treat, Pop),
           color = Treat)) +
  geom_line() +
  facet_grid(~Proj)

#Summarize
exper_evol_migr <- group_by(exper_evol_migr, Proj, Treat, Timepoint)
exper_evol_summ <- summarize(exper_evol_migr,
                             area_mean = mean(area_cm),
                             area_sd = sd(area_cm),
                             area_n = n())

#Make plot of summarized data
my_facet_labels <- c("1" = "Weak Phage", "2" = "Strong Phage")

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

##isolate migration analysis

isol_end_7x <- read.csv("74_75_76_isol_motility.csv", header = T, stringsAsFactors = F)
isol_start_7x <- read.csv("74_75_76_isol_motility_start.csv", header = T, stringsAsFactors = F)
# isol_end_125 <- read.csv("125_isol_motility.csv", header = T, stringsAsFactors = F)
# isol_start_125 <- read.csv("125_isol_motility_start.csv", header = T, stringsAsFactors = F)

#define function to split Strain data
isol_end_strain_split <- function(isol_end) {
  isol_end <- cbind(isol_end, NA, NA, NA, NA, NA)
  colnames(isol_end)[(ncol(isol_end)-4):ncol(isol_end)] <- c("Proj", "Rep", "Treat", "Isol", "Media")
  for (i in 1:nrow(isol_end)) {
    my_split <- strsplit(isol_end$Strain[i], split = "_")[[1]]
    isol_end$Proj[i] <- my_split[[1]]
    if (isol_end$Proj[i] == "P1.1") {
      isol_end$Rep[i] <- "A"
      isol_end$Treat[i] <- "A" #for Ancestor
      isol_end$Isol[i] <- my_split[[2]]
      isol_end$Media[i] <- "+Mg"
    } else if (isol_end$Proj[i] == "74") {
      isol_end$Rep[i] <- "A"
      isol_end$Treat[i] <- my_split[[2]]
      isol_end$Isol[i] <- my_split[[3]]
      if (length(my_split) > 3) {
        isol_end$Media[i] <- my_split[[4]]
      } else {
        isol_end$Media[i] <- "+Mg"
      }
    } else {
      isol_end$Rep[i] <- my_split[[2]]
      isol_end$Treat[i] <- my_split[[3]]
      isol_end$Isol[i] <- my_split[[4]]
      if (length(my_split) > 4) {
        isol_end$Media[i] <- my_split[[5]]
      } else {
        isol_end$Media[i] <- "+Mg"
      }
    }
  }
  return(isol_end)
}

#actually split Strain data
isol_end_7x <- isol_end_strain_split(isol_end_7x)
# isol_end_125 <- isol_end_strain_split(isol_end_125)

#Assign ancestor: time 0, proj 1 or 2, rep F, Treat A, isol same
ancestor_fix_end <- function(isol_end) {
  my_rows <- isol_end$Proj == "P1.1"
  isol_end$Time[my_rows] <- 0
  isol_end$Rep[my_rows] <- "F"
  isol_end$Treat[my_rows] <- "A"
  return(isol_end)
}

#fix ancestor
isol_end_7x <- ancestor_fix_end(isol_end_7x)
# isol_end_125 <- ancestor_fix_end(isol_end_125)
isol_start_7x$Time[isol_start_7x$Proj == "P1.1"] <- 0
# isol_start_125[isol_start_125$Proj == "P1.1"]$Time <- 0

#make Reps unique (74-A, 75A-B, 75B-C, 76A-D, 76B-E)
#& change project to be 1 (74, 75, 76) or 2 (125)
isol_end_7x <- uniq_reps(isol_end_7x)
# isol_end_125 <- uniq_reps(isol_end_125)

#Standardize project numbers
isol_start_7x$Proj <- 1
isol_end_7x$Proj <- 1
# isol_start_125$Project <- 2
# isol_end_125$Proj <- 2

#combine datasets
isol_end <- rbind(isol_end_7x)
isol_start <- rbind(isol_start_7x)

#remove -Mg data
isol_end <- subset(isol_end, isol_end$Media == "+Mg")
isol_end$Media <- NULL

#compile time information into timestamp
isol_end <- make_timestamp(isol_end, type = "end")
isol_start <- make_timestamp(isol_start, type = "isol start")

#convert timestamp to time since inoculation
isol_end <- calc_timediff(isol_start, isol_end, "isol")

my_facet_labels <- c("1" = "Weak Phage", 
                     "2" = "Strong Phage",
                     "C" = "Control", "G" = "Global", "L" = "Local",
                     "A" = "WT")

#plot isolate variation for each pop
isol_end$Treat.Rep <- paste(isol_end$Treat, isol_end$Rep)
ggplot(isol_end, aes(x = Treat.Rep, y = `Rate (cm/hr)`)) + 
  geom_jitter(position = position_jitter(0), size = 2) + 
  facet_grid(.~Proj, labeller = labeller(Proj = my_facet_labels)) +
  theme_bw() + labs(y = "Migration Rate (cm/hr)", x = "") +
  scale_x_discrete(labels = c("WT", LETTERS[1:5], "A", "B",
                              "D", "E", LETTERS[1:3], "E")) +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        plot.margin = unit(c(0, 0, 0.075, 0), "npc"))
grid.text(label = c("Control", "Global", "Local"),
          x = unit(c(.31, .6, 0.86), "npc"),
          y = unit(0.06, "npc"))
grid.polyline(x = unit(c(0.17, 0.45, 0.49, 0.72, 0.75, 0.975), unit = "npc"),
              y = unit(rep(0.09, 6), unit = "npc"),
              id = c(1, 1, 2, 2, 3, 3), gp = gpar(lwd = 3))

#plot isolate variation for each pop
#for poster
isol_end$Treat.Rep <- paste(isol_end$Treat, isol_end$Rep)
png("isol_mig_rate.png", width = 11, height = 7, units = "in",
    res = 300)
ggplot(isol_end, aes(x = Treat.Rep, y = `Rate (cm/hr)`)) + 
  geom_jitter(position = position_jitter(0), size = 4) + 
  facet_grid(Proj~., labeller = labeller(Proj = my_facet_labels)) +
  theme_bw() + labs(y = "Migration Rate (cm/hr)", x = "") +
  scale_x_discrete(labels = c("WT", LETTERS[1:5], "A", "B",
                              "D", "E", LETTERS[1:3], "E")) +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        plot.margin = unit(c(0.01, 0.01, 0.075, 0.01), "npc"))
grid.text(label = c("Control", "Global", "Local"),
          x = unit(c(.3, .58, 0.825), "npc"),
          y = unit(0.06, "npc"),
          gp = gpar(fontsize = 20))
grid.polyline(x = unit(c(0.17, 0.45, 0.48, 0.69, 0.72, 0.94), unit = "npc"),
              y = unit(rep(0.09, 6), unit = "npc"),
              id = c(1, 1, 2, 2, 3, 3), gp = gpar(lwd = 4))
dev.off()

#have to clean up correlation code

# #plot correlation of T14 pop rate & isol rates
# pop_isol_frame <- data.frame(Time=integer(), Rep=character(), Treat=character(), 
#                              Isol=character(), Pop_Rate=double(), Isol_Rate=double(),
#                              stringsAsFactors = F)
# for (i in 1:nrow(isol_end)) {
#   if (isol_end$Proj[i] != "P1.1") {
#     my.sub <- subset(end_data, end_data$Time == isol_end$Time[i] & 
#                        end_data$Rep == isol_end$Rep[i] & 
#                        substr(end_data$Treat, 1, 1) == isol_end$Treat[i])
#     pop_isol_frame <- rbind(pop_isol_frame, 
#                             data.frame(Time=isol_end$Time[i],
#                                        Rep=isol_end$Rep[i],
#                                        Treat=isol_end$Treat[i],
#                                        Isol=isol_end$Isol[i],
#                                        Pop_Rate=my.sub$`Rate (cm/hr)`[1],
#                                        Isol_Rate=isol_end$`Rate (cm/hr)`[i]))
#   }
# }
# # plot(pop_isol_frame$Isol_Rate~pop_isol_frame$Pop_Rate, xlim = c(0.025, 0.065),
# #      ylim = c(0.025, 0.065))
# # lines(x=c(0.02, 0.07), y=c(0.02, 0.07))
# # summary(lm(pop_isol_frame$Isol_Rate~pop_isol_frame$Pop_Rate))
# my_colors <- c("red", "green", "blue")
# my.pch = c(20, 20, 20)
# tiff(filename = "isol_vs_pop.tiff", width = 10, height = 10, units = "in",
#      compression = "none", res = 300)
# par(mar = my.mar + c(0, 1, 0, 0))
# for (i in 1:length(unique(pop_isol_frame$Treat))) {
#   treat = unique(pop_isol_frame$Treat)[i]
#   my.sub <- subset(pop_isol_frame, pop_isol_frame$Treat == treat)
#   if (i == 1) {
#     plot(my.sub$Isol_Rate~my.sub$Pop_Rate, xlim = c(0.027, 0.061), ylim = c(0.027, 0.061),
#          pch = my.pch[i], xlab = "Population Migration Rate (cm/hr)", ylab = "Isolate Migration Rate (cm/hr)",
#          col = my_colors[i], cex = 2, cex.lab = 2, cex.axis = 2)
#   } else {
#     points(my.sub$Isol_Rate~my.sub$Pop_Rate, pch = my.pch[i], col = my_colors[i],
#            cex = 2)
#   }
# }
# lines(x=c(0, 1), y=c(0, 1), lwd = 2)
# legend(x = 0.053, y = 0.035, legend = unique(mean_rates$Treatment), pch = my.pch, 
#        col = my_colors, cex = 2)
# dev.off()

##isolate growth curve analysis

##standard curve to convert OD to CFU (based on 107 data)
stan_plate <- read.csv("107_Std_Curve.csv", stringsAsFactors = F, header = T)
stan_plate_layout <- read.csv("107_Plate_Layout.csv", stringsAsFactors = F, header = F)
stan_spec <- read.csv("107_Trav_Spec.csv", stringsAsFactors = F, header = T)

#define function to melt plate arrays
get_layout <- function(plate_layout, row_start = 1, col_start = 1) {
  #tidies the data: makes list of Well ID & contents from the plate layout
  layout_list <- data.frame("Well" = character((nrow(plate_layout)-(row_start-1))*(ncol(plate_layout)-(col_start-1))), 
                            "Contents" = character((nrow(plate_layout)-(row_start-1))*(ncol(plate_layout)-(col_start-1))),
                            stringsAsFactors = F)
  cntr = 1
  for (i in row_start:nrow(plate_layout)) {
    for (j in col_start:ncol(plate_layout)) {
      layout_list[cntr, ] <- c(paste(LETTERS[i-row_start+1], j-col_start+1, sep = ""), plate_layout[i, j])
      cntr = cntr + 1
    }
  }
  return(layout_list)
}

#melt plate arrays, combine & exclude empty wells
stan_layout_list <- get_layout(stan_plate_layout, 2, 2)
stan_plate_list <- get_layout(stan_plate, 1, 2)
colnames(stan_plate_list)[2] <- "OD600"
stan_plate_list$Dilution <- stan_layout_list$Contents
stan_plate_list <- subset(stan_plate_list, is.na(stan_plate_list$Dilution) == F)

#in data from Trav-lab spec, convert OD to CFU
stan_spec$CFU <- (stan_spec$OD600..on.Trav.lab.spec.+3.318*10^(-3))/(1.252*10^(-9))

#transfer CFU data to plate list
stan_plate_list <- cbind(stan_plate_list, "CFU"=stan_spec$CFU[match(stan_plate_list$Dilution, 
                                                     stan_spec$Dilution)])

#get regression of data
# ggplot(stan_plate_list, aes(x = CFU, y = as.numeric(OD600))) + geom_point()
stan_plate_list$OD600 <- as.numeric(stan_plate_list$OD600)
my.fit <- lm(CFU ~ OD600, data = stan_plate_list)

##growth curve analysis
library("minpack.lm")
library("tidyr")
library("lubridate")
library("dplyr")

options(stringsAsFactors = F)
growth_97 <- read.csv("97B.csv", header = T, stringsAsFactors = F)
growth_98 <- read.csv("98B.csv", header = T, stringsAsFactors = F)
growth_99 <- read.csv("99.csv", header = T, stringsAsFactors = F)
growth_101 <- read.csv("101.csv", header = T, stringsAsFactors = F)
plate_layout_97_101 <- read.csv("97_98_99_100_101_plate_layout.csv", header = F, stringsAsFactors = F)

#make list of well ID & contents from plate layout
plate_layout_list_97_101 <- get_layout(plate_layout_97_101, 2, 2)

#remove columns from growth data that correspond to water
to_rem <- subset(plate_layout_list_97_101, 
                 plate_layout_list_97_101$Contents == "H2O")
growth_97 <- growth_97[, !(colnames(growth_97) %in% to_rem$Well)]
growth_98 <- growth_98[, !(colnames(growth_98) %in% to_rem$Well)]
growth_99 <- growth_99[, !(colnames(growth_99) %in% to_rem$Well)]
growth_101 <- growth_101[, !(colnames(growth_101) %in% to_rem$Well)]

#replace column names with Isols
isol.name.func <- function(growth_array, layout_list) {
  for (i in 1:ncol(growth_array)) {
    j <- match(colnames(growth_array)[i], layout_list$Well)
    if (!is.na(j)) {
      colnames(growth_array)[i] <- as.character(layout_list$Contents[j])
    }
  }
  return(growth_array)
}
growth_97 <- isol.name.func(growth_97, plate_layout_list_97_101)
growth_98 <- isol.name.func(growth_98, plate_layout_list_97_101)
growth_99 <- isol.name.func(growth_99, plate_layout_list_97_101)
growth_101 <- isol.name.func(growth_101, plate_layout_list_97_101)

#Tidy growth data
tidy_97 <- gather(growth_97, "Well", "OD600", 3:ncol(growth_97), na.rm = T)
tidy_98 <- gather(growth_98, "Well", "OD600", 3:ncol(growth_98), na.rm = T)
tidy_99 <- gather(growth_99, "Well", "OD600", 3:ncol(growth_99), na.rm = T)
tidy_101 <- gather(growth_101, "Well", "OD600", 3:ncol(growth_101), na.rm = T)

#Add isolate & project info
tidy_97$Isol <- "A"
tidy_98$Isol <- "B"
tidy_99$Isol <- "C"
tidy_101$Isol <- "E"
tidy_97$Proj <- 1
tidy_98$Proj <- 1
tidy_99$Proj <- 1
tidy_101$Proj <- 1

#Combine data frames
gc_data <- rbind(tidy_97, tidy_98, tidy_99, tidy_101)

#convert OD to CFU
gc_data$CFU <- predict(my.fit, newdata = gc_data)

#separate out ID info
sep.growth.ID <- function(my.array) {
  my.array <- cbind(my.array, "Media"=NA, "Pop"=NA, "Treat"=NA, "Rep_Well"=NA)
  for (i in 1:nrow(my.array)) {
    my.split <- strsplit(as.character(my.array$Well[i]), split = "-")
    my.array$Media[i] <- my.split[[1]][1]
    my.array$Treat[i] <- my.split[[1]][length(my.split[[1]])-1]
    my.array$Rep_Well[i] <- my.split[[1]][length(my.split[[1]])]
    if (length(my.split[[1]])>3) { #it's not an ancestor well
      if (substr(my.split[[1]][2], 1, 2) == "74") {
        my.array$Pop[i] <- "A"
      } else if (substr(my.split[[1]][2], 1, 2) == "75") {
        if (substr(my.split[[1]][2], 3, 3) == "A") {
          my.array$Pop[i] <- "B"
        } else {my.array$Pop[i] <- "C"}
      } else if (substr(my.split[[1]][2], 1, 2) == "76") {
        if (substr(my.split[[1]][2], 3, 3) == "A") {
          my.array$Pop[i] <- "D"
        } else {my.array$Pop[i] <- "E"} 
      }
    } else { #assign ancestor to be Pop F, Treat A
      my.array$Pop[i] <- "F"
    }
  }
  return(my.array)
}

gc_data <- sep.growth.ID(gc_data)

#reformat time column
gc_data$Time <- strptime(gc_data$Time, format = "%H:%M:%S")

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
