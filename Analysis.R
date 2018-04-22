library("ggplot2")
library("reshape")

#need to make Rep vs Pop consistent
#redo all analysis using dplyr
#redo analysis using vectorization
#include std curve for spec in this code (it may alread be there???)
#Do a scripted analysis of how much to smoothe gc data (not manually checking)
# #to be done later:
# #define function to get means of rates
# mean_rates <- function(sub_by_list, my_data)
# analyze additional plate scans

##experimental evolution analysis
end_data_7x <- read.csv("74_75_76_Data_Measurements.csv", header = T, stringsAsFactors = F)
start_data_7x <- read.csv("74_75_76_Start_Data.csv", header = T, stringsAsFactors = F)
fails_7x <- read.csv("74_75_76_fails.csv", header = T, stringsAsFactors = F)
contam_7x <- read.csv("74_75_76_contaminants.csv", header = T, stringsAsFactors = F)
end_data_125 <- read.csv("125_Data_Measurements.csv", header = T, stringsAsFactors = F)
start_data_125 <- read.csv("125_Start_Data.csv", header = T, stringsAsFactors = F)
fails_125 <- read.csv("125_fails.csv", header = T, stringsAsFactors = F)

##fix Sydney's measurements in 74, 75, 76 (they were taken in inches, not cm)
for (i in 1:nrow(end_data_7x)) {
  if (end_data_7x$Entered.By[i] == "Sydney") {
    end_data_7x$Width..cm.[i] <- 2.54*end_data_7x$Width..cm.[i]
    end_data_7x$Height..cm.[i] <- 2.54*end_data_7x$Height..cm.[i]
  }
}

#define function to parse out Strain data into component parts
strain_split <- function(my_data, num_cols) {
  my_data$Proj <- NA
  my_data$Time <- NA
  my_data$Rep <- NA
  if (num_cols == 4) {my_data$Treat <- NA}
  for (i in 1:nrow(my_data)) {
    my_split <- strsplit(my_data$Strain[i], split = "_")[[1]]
    my_data$Proj[i] <- my_split[[1]]
    my_data$Time[i] <- my_split[[2]]
    if (length(my_split) < num_cols) {
      my_data$Rep[i] <- "A"
    } else {
      my_data$Rep[i] <- my_split[[num_cols]]
    }
    if (num_cols == 4) {
      my_data$Treat[i] <- my_split[[3]]
    }
  }
  return(my_data)
}

#parse out strain data
start_data_7x <- strain_split(start_data_7x, 3)
end_data_7x <- strain_split(end_data_7x, 4)
start_data_125 <- strain_split(start_data_125, 3)
end_data_125 <- strain_split(end_data_125, 4)

#turn treatment values into only first letter for end_data_7x
end_data_7x$Treat <- substr(end_data_7x$Treat, 1, 1)
contam_7x$Treat <- substr(contam_7x$Treat, 1, 1)

#combine 2 projects: end data & start data
end_data <- rbind(end_data_7x, end_data_125)
start_data <- rbind(start_data_7x, start_data_125)
fails_data <- rbind(fails_7x, fails_125)

#define functions to remove failed runs
start_remove_fails <- function(start_data, fails) {
  start_data <- cbind(start_data, NA)
  colnames(start_data)[ncol(start_data)] <- "Remove"
  for (i in 1:nrow(start_data)) {
    for (j in 1:nrow(fails)) {
      if (start_data$Proj[i] == fails$Proj[j] & start_data$Time[i] == fails$Time[j] & start_data$Rep[i] == fails$Rep[j]) {
        start_data$Remove[i] <- TRUE
      }
    }
    if (is.na(start_data$Remove[i])) {start_data$Remove[i] <- FALSE}
  }
  start_data <- subset(start_data, start_data$Remove == F)
  start_data <- start_data[, -ncol(start_data)]
  return(start_data)
}
end_remove_fails <- function(end_data, fails) {
  end_data <- cbind(end_data, NA)
  colnames(end_data)[ncol(end_data)] <- "Remove"
  for (i in 1:nrow(end_data)) {
    for (j in 1:nrow(fails)) {
      if (end_data$Proj[i] == fails$Proj[j] & end_data$Time[i] == fails$Time[j] & end_data$Rep[i] == fails$Rep[j]) {
        end_data$Remove[i] <- TRUE
      }
    }
    if (is.na(end_data$Remove[i])) {end_data$Remove[i] <- FALSE}
  }
  end_data <- subset(end_data, end_data$Remove == F)
  end_data <- end_data[, -ncol(end_data)]
  return(end_data)
}

#actually remove failed runs
start_data <- start_remove_fails(start_data, fails_data)
end_data <- end_remove_fails(end_data, fails_data)

#define functions to remove contaminated lines
end_remove_contam <- function(end_data, contam) {
  end_data <- cbind(end_data, NA)
  colnames(end_data)[ncol(end_data)] <- "Remove"
  for (i in 1:nrow(end_data)) {
    for (j in 1:nrow(contam)) {
      if (end_data$Proj[i] == contam$Proj[j] & end_data$Treat[i] == contam$Treat[j] & end_data$Rep[i] == contam$Rep[j]) {
        end_data$Remove[i] <- TRUE
      }
    }
    if (is.na(end_data$Remove[i])) {end_data$Remove[i] <- FALSE}
  }
  end_data <- subset(end_data, end_data$Remove == F)
  end_data <- end_data[, -ncol(end_data)]
  return(end_data)
}

#actually remove contaminated lines: 75 B G, 76 A L
end_data <- end_remove_contam(end_data, contam_7x)

#define function to change the reps so they're unique (74-A, 75A-B, 75B-C, 76A-D, 76B-E)
#& change projects to 1 (for 74, 75, 76) & 2 (for 125) 
uniq_reps <- function(my_data) {
  for (i in 1:nrow(my_data)) {
    if (my_data$Proj[i] == "74") {
      my_data$Rep[i] <- "A"
    } else if (my_data$Proj[i] == "75") {
      if (my_data$Rep[i] == "A") {
        my_data$Rep[i] <- "B"
      } else if (my_data$Rep[i] == "B") {
        my_data$Rep[i] <- "C"
      }
    } else if (my_data$Proj[i] == "76") {
      if (my_data$Rep[i] == "A") {
        my_data$Rep[i] <- "D"
      } else if (my_data$Rep[i] == "B") {
        my_data$Rep[i] <- "E"
      }
    }
    if (my_data$Proj[i]=="74" | my_data$Proj[i]=="75" | my_data$Proj[i]=="76") {
      my_data$Proj[i] <- 1
    } else if (my_data$Proj[i] == "125") {my_data$Proj[i] <- 2}
  }
  return(my_data)
}
  
#make start & end data have unique reps
start_data <- uniq_reps(start_data)
end_data <- uniq_reps(end_data)

#define function to pull out forisol into separate data frames
#returns the main data frame, then the forisol frame, in a list
forisol_rem <- function(my_data) {
  my_data <- cbind(my_data, NA)
  colnames(my_data)[ncol(my_data)] <- "forisol"
  for (i in 1:nrow(my_data)) {
    if (grepl("forisol", my_data$Time[i])) {
      my_data$forisol[i] <- TRUE
    } else {
      my_data$forisol[i] <- FALSE
    }
  }
  my_isol <- subset(my_data, my_data$forisol == T)
  my_data <- subset(my_data, my_data$forisol == F)[, -ncol(my_data)]
  return(list(my_data, my_isol))
}

#actually remove forisol
end_data <- forisol_rem(end_data)[[1]]
start_data <- forisol_rem(start_data)[[1]]

#define function to remove hyphenated 2nd part from tranfer #'s
rem_hyphens <- function(my_data) {
  for (i in 1:nrow(my_data)) {
    my_split <- strsplit(my_data$Time[i], split = "-")
    my_data$Time[i] <- substr(my_split[[1]][1], 2, nchar(my_split[[1]][1]))
  }
  return(my_data)
}

#actually change transfer values so re-dos don't have hyphenated 2nd part
start_data <- rem_hyphens(start_data)
end_data <- rem_hyphens(end_data)

#define function to compile time info into timestamp
#then get rid of all time columns as well as Strain column
make_timestamp <- function(my_data, type) {
  my_data <- cbind(my_data, NA)
  colnames(my_data)[ncol(my_data)] <- "timestamp"
  my_data$timestamp <- paste(paste(my_data$Year, my_data$Month, my_data$Day, sep = "-"),
                                paste(my_data$Hour, my_data$Minute, sep = ":"), sep = " ")
  my_data$timestamp <- strptime(my_data$timestamp, format = "%Y-%m-%d %H:%M")
  if (type=="start") {my_data <- my_data[, -c(1:6)]
  } else if (type=="end") {my_data <- my_data[, -c(2:8)]
  } else if (type=="isol start") {my_data <- my_data[, -c(1:5)]}
  return(my_data)
}

#compile time information into timestamp
start_data <- make_timestamp(start_data, "start")
end_data <- make_timestamp(end_data, "end")

#fix this below to include proj info

#define function to convert timestamp to time since inoculation
#then calculate rate (average radius spread per hour)
calc_timediff <- function(my_start, my_end, type="pop") {
  my_end <- cbind(my_end, NA)
  colnames(my_end)[ncol(my_end)] <- "timediff"
  for (i in 1:nrow(my_end)) {
    if (type == "isol") {
      sub_start <- subset(my_start, my_start$Proj == my_end$Proj[i] &
                            my_start$Time == my_end$Time[i] &
                            my_start$Isol == my_end$Isol[i])
    } else {
      sub_start <- subset(my_start, my_start$Proj == my_end$Proj[i] & 
                            my_start$Time == my_end$Time[i] & 
                            my_start$Rep == my_end$Rep[i])
    }
    my_end$timediff[i] <- difftime(my_end$timestamp[i], sub_start$timestamp[1], units = "hours")
  }
  my_end <- cbind(my_end, NA)
  colnames(my_end)[ncol(my_end)] <- "Rate (cm/hr)"
  my_end$`Rate (cm/hr)` <- (my_end$Width..cm.+my_end$Height..cm.)/(4*my_end$timediff)
  return(my_end)
}

#actually convert timestamp to time since inoculation, then calc rate
end_data <- calc_timediff(start_data, end_data)

#def func for mean & sd of rates for ea treat & timepoint
sum_rates <- function(my_end_data) {
  #returns a data frame summary
  my_mean_rates <- data.frame(Proj=character(), Treat=character(), Time=integer(), 
                           Mean.Rate=double(), SD.Rate=double(), n=integer(),
                           stringsAsFactors = F)
  for (my.proj in unique(my_end_data$Proj)) {
    my.sub.a <- subset(my_end_data, my_end_data$Proj == my.proj)
    for (my.time in unique(my.sub.a$Time)) {
      my.sub.b <- subset(my.sub.a, my.sub.a$Time == my.time)
      for (my.treat in unique(my.sub.b$Treat)) {
        my.sub.c <- subset(my.sub.b, my.sub.b$Treat == my.treat)
        my_mean_rates <- rbind(my_mean_rates, c(as.character(my.proj), as.character(my.treat), 
                                          as.integer(my.time), as.double(mean(my.sub.c$Rate)), 
                                          as.double(sd(my.sub.c$Rate)), nrow(my.sub.c)), 
                            stringsAsFactors = F)
      }
    }
  }
  colnames(my_mean_rates) <- c("Proj", "Treat", "Time", "Mean_Rate", "SD_Rate", "n")
  return(my_mean_rates)
}

#actually get mean rates 
mean_rates <- sum_rates(end_data)

#change classes
end_data$Time <- as.numeric(end_data$Time)
mean_rates$Time <- as.numeric(mean_rates$Time)
mean_rates$Mean_Rate <- as.numeric(mean_rates$Mean_Rate)
mean_rates$SD_Rate <- as.numeric(mean_rates$SD_Rate)

#limit analysis to first 14 transfers
end_data <- subset(end_data, end_data$Time <= 14)
mean_rates <- subset(mean_rates, mean_rates$Time<=14)

my_facet_labels <- c("1" = "Weak Phage", "2" = "Strong Phage")
#make plot by treatment
ggplot(data = mean_rates, 
       aes(x=Time, y=Mean_Rate, group=Treat, colour=Treat)) +
  geom_line(size = 1.2) + geom_point() + 
  theme(axis.text.y = element_text(size = 11), axis.text.x = element_text(size = 11),
        legend.text = element_text(size = 16)) +
  facet_grid(~Proj, labeller = labeller(Proj = my_facet_labels)) + 
  geom_errorbar(aes(ymin=Mean_Rate-SD_Rate, ymax=Mean_Rate+SD_Rate),
                width=1, size = .7, position=position_dodge(0.2)) +
  labs(x = "Transfer", y = "Mean Migration Rate (cm/hr)") + 
  scale_color_hue(name = "Treatment", breaks = c("C", "G", "L"),
                  labels = c("Control", "Global", "Local")) +
  theme_bw()

#make plot of ea treat to check if pops are stable position
my_facet_labels <- c("1" = "Weak Phage", "2" = "Strong Phage", "C" = "Control",
                     "G" = "Global", "L" = "Local")
ggplot(data = end_data, aes(x=Time, y=`Rate (cm/hr)`, group=Rep, colour=Rep)) +
  geom_line(size = 1.1) + geom_point() + 
  facet_grid(Treat~Proj, labeller = labeller(Proj = my_facet_labels, 
                                             Treat = my_facet_labels)) + 
  theme(axis.text.y = element_text(size = 11), axis.text.x = element_text(size = 11)) +
  labs(x = "Transfer", y = "Migration Rate (cm/hr)") + 
  scale_color_hue(name = "Replicate\nPopulation") +
  theme_bw()


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
  theme_bw() + labs(y = "Migration Rate (cm/hr)") +
  scale_x_discrete(labels = c("WT", LETTERS[1:5], "A", "B",
                              "D", "E", LETTERS[1:3], "E")) +
  theme(axis.text.x = element_text(size = 12, color = "black"))

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

doubling_time <- function(CFU, time, subset_by) {
  ans <- c(as.numeric(difftime(time[2:length(time)], time[1:(length(time)-1)],
             units = "mins"))/
    log2(CFU[2:length(CFU)]/CFU[1:(length(CFU)-1)]), NA)
  ans[subset_by[2:length(subset_by)] != subset_by[1:(length(subset_by)-1)]] <- NA
  return(ans)
}

per_cap_grate <- function(CFU, time, subset_by) {
  ans <- c((CFU[2:length(CFU)]-CFU[1:(length(CFU)-1)])/
    (as.numeric(difftime(time[2:length(time)], time[1:(length(time)-1)], 
                         units = "hours")) * 
       0.5 * (CFU[2:length(CFU)]+CFU[1:(length(CFU)-1)])), NA)
  ans[subset_by[2:length(subset_by)] != subset_by[1:(length(subset_by)-1)]] <- NA
  return(ans)
}

gc_data$doubtime <- doubling_time(gc_data$Smooth_CFU, gc_data$Time,
                                  gc_data$mpptir)
gc_sm$doubtime <- doubling_time(gc_sm$Smooth_noverlap, gc_sm$Time, gc_sm$mpptir)
gc_data$pcgr <- per_cap_grate(gc_data$Smooth_CFU, gc_data$Time, gc_data$mpptir)
gc_sm$pcgr <- per_cap_grate(gc_sm$Smooth_noverlap, gc_sm$Time, gc_sm$mpptir)

# #Remove unrealistic values
# gc_data$doubtime[gc_data$doubtime < 0] <- NA
# gc_sm$doubtime[gc_sm$doubtime < 0] <- NA
# gc_data$doubtime[gc_data$doubtime > 1000] <- NA
# gc_sm$doubtime[gc_sm$doubtime > 1000] <- NA

# #Make plot of all smoothed overlapping doubtimes
# for (i in seq(from = 1, to = length(unique(gc_data$mpptir)), by = 25)) {
#   my_sub <- gc_data[(gc_data$mpptir %in% unique(gc_data$mpptir)[i:(i+24)]), ]
#   print(ggplot(data = my_sub, aes(x = Time, y = doubtime)) + geom_line() +
#           facet_wrap(~mpptir))
# }
# 
# #Make plots of all smoothed nonoverlapping doubtimes
# for (i in seq(from = 1, to = length(unique(gc_sm$mpptir)), by = 25)) {
#   my_sub <- gc_sm[(gc_sm$mpptir %in% unique(gc_sm$mpptir)[i:(i+24)]), ]
#   print(ggplot(data = my_sub, aes(x = Time, y = doubtime)) + geom_line() +
#           facet_wrap(~mpptir))
# }

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
# ggplot(data = my_sub, aes(x = Time, y = doubtime)) +
#   geom_line() + labs(y = "Doubling Time (mins)")
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


#extract minimum doubling times for each uniq well
#using overlapping averaged smoothing
gc_data$Time <- as.character(gc_data$Time)
gc_mpptir <- group_by(gc_data, Media, Proj, Pop, Treat, Isol, Rep_Well)
gc_mpptir <- summarize(gc_mpptir, min_dt = min(doubtime, na.rm = T))
gc_mppti <- group_by(gc_mpptir, Media, Proj, Pop, Treat, Isol)
gc_mppti <- summarize(gc_mppti, dt_min_avg = mean(min_dt),
                      dt_min_sd = sd(min_dt))
gc_mppt <- group_by(gc_mppti, Media, Proj, Pop, Treat)
gc_mppt <- summarize(gc_mppt, avg_isols = mean(dt_min_avg),
                     sd_isols = sd(dt_min_avg))

#Making doubtime plots
my_facet_labels <- c("100" = "Rich Environment", 
                     "50" = "Adapted Environment",
                     "C" = "Control", "G" = "Global", "L" = "Local",
                     "A" = "WT")

#plot of all isols
ggplot(gc_mppti, aes(x = Treat, y = dt_min_avg)) + 
  geom_jitter(width = 0.1, height = 0, size = 2) + 
  facet_grid(Media ~ Pop, labeller = labeller(Media = my_facet_labels)) +
  labs(x = "Treatment", y = "Doubling Time (min)") +
  theme(axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11)) +
  theme_bw()

#plot of all pops
ggplot(gc_mppt, aes(x = Treat, y = avg_isols), labeller = labeller(Treat = my_facet_labels)) + geom_point(pch = 1, size = 3) +
  facet_grid(Media ~ ., labeller = labeller(Media = my_facet_labels)) + 
  labs(x = "Treatment", y = "Doubling Time (min)") + theme_bw() + 
  scale_x_discrete(labels = c("Ancestor", "Control", "Global", "Local"))

#extract minimum doubling times for each uniq well
#using non-overlapping averaged smoothing
gc_sm$Time <- as.character(gc_sm$Time)
gc_mpptir <- group_by(gc_sm, Media, Proj, Pop, Treat, Isol, Rep_Well)
gc_mpptir <- summarize(gc_mpptir, min_dt = min(doubtime, na.rm = T))
gc_mppti <- group_by(gc_mpptir, Media, Proj, Pop, Treat, Isol)
gc_mppti <- summarize(gc_mppti, dt_min_avg = mean(min_dt),
                      dt_min_sd = sd(min_dt))
gc_mppt <- group_by(gc_mppti, Media, Proj, Pop, Treat)
gc_mppt <- summarize(gc_mppt, avg_isols = mean(dt_min_avg),
                     sd_isols = sd(dt_min_avg))

#Making doubtime plots
my_facet_labels <- c("100" = "Rich Environment", 
                     "50" = "Adapted Environment",
                     "C" = "Control", "G" = "Global", "L" = "Local",
                     "A" = "WT")

#plot of all isols
gc_mppti$Media <- factor(gc_mppti$Media, levels = c(50, 100))
ggplot(gc_mppti, aes(x = Treat, y = dt_min_avg)) + 
  geom_jitter(width = 0.1, height = 0, size = 2) + 
  facet_grid(Media ~ Pop, labeller = labeller(Media = my_facet_labels)) +
  labs(x = "Treatment", y = "Doubling Time (min)") +
  theme(axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11)) +
  theme_bw()

#plot of all pops
gc_mppt$Media <- factor(gc_mppt$Media, levels = c(50, 100))
ggplot(gc_mppt, aes(x = Treat, y = avg_isols), labeller = labeller(Treat = my_facet_labels)) + geom_point(pch = 1, size = 3) +
  facet_grid(.~Media, labeller = labeller(Media = my_facet_labels)) + 
  labs(x = "Treatment", y = "Doubling Time (min)") + theme_bw() + 
  scale_x_discrete(labels = c("Ancestor", "Control", "Global", "Local"))

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

#Making doubtime plots
my_facet_labels <- c("100" = "Rich Environment", 
                     "50" = "Adapted Environment",
                     "C" = "Control", "G" = "Global", "L" = "Local",
                     "A" = "WT")

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
ggplot(gc_mppt, aes(x = Treat, y = avg_isols), labeller = labeller(Treat = my_facet_labels)) + geom_point(pch = 1, size = 3) +
  facet_grid(.~Media, labeller = labeller(Media = my_facet_labels)) + 
  labs(x = "Treatment", y = "Per Capita Growth Rate (/hour)") + theme_bw() + 
  scale_x_discrete(labels = c("Ancestor", "Control", "Global", "Local"))


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

# geom_point(position = position_dodge(width = 0), size = 2) +
  

#resistance vs growth
gc_resis_data <- merge(resis_data, gc_mppti)
gc_resis_mppt <- group_by(gc_resis_data, Media, Proj, Pop, Treat)
gc_resis_mppt <- summarize(gc_resis_mppt, avg_eop = mean(EOP),
                           avg_dt = mean(dt_min_avg))

#Make plot of all isolates
ggplot(gc_resis_data, aes(x = 1-EOP, y = dt_min_avg)) +
  geom_point() + 
  facet_grid(Media~., labeller = labeller(Media = my_facet_labels)) +
  geom_smooth(method = "lm") +
  labs(x = "Resistance", y = "Doubling Time (min)")

summary(lm(dt_min_avg~Media*EOP, data = gc_resis_data))

#Make plot of all pops
ggplot(gc_resis_mppt, aes(x = 1-avg_eop, y = avg_dt)) +
  geom_point() + facet_grid(Media~.) +
  geom_smooth(method = "lm")

summary(lm(avg_dt~avg_eop*Media, data = gc_resis_mppt))