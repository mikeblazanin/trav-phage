library("ggplot2")

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
    } else {my_data$Proj[i] <- 2}
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
calc_timediff <- function(my_start, my_end) {
  my_end <- cbind(my_end, NA)
  colnames(my_end)[ncol(my_end)] <- "timediff"
  for (i in 1:nrow(my_end)) {
    sub_start <- subset(my_start, my_start$Proj == my_end$Proj[i] & 
                          my_start$Time == my_end$Time[i] & 
                          my_start$Rep == my_end$Rep[i])
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

#make plot by treatment
ggplot(data = mean_rates, 
       aes(x=Time, y=Mean_Rate, group=Treat, colour=Treat)) +
  geom_line() + geom_point() + theme(axis.text.y = element_blank()) +
  facet_grid(~Proj) + 
  geom_errorbar(aes(ymin=Mean_Rate-SD_Rate, ymax=Mean_Rate+SD_Rate),
                width=0.7, position=position_dodge(0.2))

#make plot of ea treat to check if pops are stable position
ggplot(data = end_data, aes(x=Time, y=`Rate (cm/hr)`, group=Rep, colour=Rep)) +
  geom_line() + geom_point() + facet_grid(Treat~Proj) +
  theme(axis.text.y = element_blank())

#isolate 50% KB migration analysis

isol_end_7x <- read.csv("74_75_76_isol_motility.csv", header = T, stringsAsFactors = F)
isol_start_7x <- read.csv("74_75_76_isol_motility_start.csv", header = T, stringsAsFactors = F)

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

#make Reps unique (74-A, 75A-B, 75B-C, 76A-D, 76B-E)
#& change project to be 1 (74, 75, 76) or 2 (125)
isol_end_7x <- uniq_reps(isol_end_7x)

#problem: makes ancestor into proj 2
#soln: should label ancestor as Time 0, Proj 1, Rep F, Treat A, isol matching
#soln: in start data, need to split out proj's to reps (so then P1.1 can be F)

#compile time information into timestamp
isol_end_7x <- make_timestamp(isol_end_7x, type = "end")
isol_start_7x <- make_timestamp(isol_start_7x, type = "isol start")

#convert timestamp to time since inoculation
isol_end_7x <- calc_timediff(isol_start_7x, isol_end_7x)


isol_end <- cbind(isol_end, NA)
colnames(isol_end)[ncol(isol_end)] <- "timediff"
for (i in 1:nrow(isol_end)) {
  if (isol_end$Proj[i] == "P1.1") {
    my_start <- subset(isol_start, isol_start$Project == isol_end$Proj[i] & isol_start$Isol == isol_end$Isol[i])
  } else {
    my_start <- subset(isol_start, isol_start$Project == isol_end$Proj[i] & isol_start$Isol == isol_end$Isol[i] & isol_start$Time == isol_end$Time[i])
  }
  isol_end$timediff[i] <- difftime(isol_end$timestamp[i], my_start$timestamp[1], units = "hours")
}

#calculate rate (average radius spread per hour)
isol_end <- cbind(isol_end, NA)
colnames(isol_end)[ncol(isol_end)] <- "Rate (cm/hr)"
isol_end$`Rate (cm/hr)` <- (isol_end$Width..cm.+isol_end$Height..cm.)/(4*isol_end$timediff)

#remove -Mg data
no_mg_isol_end <- subset(isol_end, isol_end$Media == "-Mg")
isol_end <- subset(isol_end, isol_end$Media == "+Mg")

#remove >24 hr data
non_24_hr_isol_end <- subset(isol_end, isol_end$timediff > 36)
isol_end <- subset(isol_end, isol_end$timediff < 36)

#plot isolate variation for each pop
tiff(filename = "isol_rates.tiff", width = 15, height = 7, units = "in",
     compression = "none", res = 300)
par(mar = my.mar + c(1, 3, 0, 0), xpd = T)
isol_end <- cbind(isol_end, "Treat.Rep" = paste(isol_end$Treat, isol_end$Rep, sep = "."))
stripchart(isol_end$`Rate (cm/hr)`~isol_end$Treat.Rep, vert = T, pch = 19,
           cex = 2, cex.axis = 2, ylab = "",
           group.names = c("WT", "A", "B", "C", "D", "E", "A", "B", "C", "D",
                           "A", "B", "C", "D"))
mtext("Isolate Migration Rate (cm/hr)", side = 2, cex = 2.5, line = 4)
#xlab = "Population.Treatment",
text(c("Control", "Global", "Local"), x = c(4, 8.5, 12.5), y = 0.022, cex = 2.5)
my.lwd <- 3
lines(x = c(1.8, 6.2), y = c(0.024, 0.024), lwd = my.lwd)
lines(x = c(6.8, 10.2), y = c(0.024, 0.024), lwd = my.lwd)
lines(x = c(10.8, 14.2), y = c(0.024, 0.024), lwd = my.lwd)
dev.off()

#plot correlation of T14 pop rate & isol rates
pop_isol_frame <- data.frame(Time=integer(), Rep=character(), Treat=character(), 
                             Isol=character(), Pop_Rate=double(), Isol_Rate=double(),
                             stringsAsFactors = F)
for (i in 1:nrow(isol_end)) {
  if (isol_end$Proj[i] != "P1.1") {
    my.sub <- subset(end_data, end_data$Time == isol_end$Time[i] & 
                       end_data$Rep == isol_end$Rep[i] & 
                       substr(end_data$Treat, 1, 1) == isol_end$Treat[i])
    pop_isol_frame <- rbind(pop_isol_frame, 
                            data.frame(Time=isol_end$Time[i],
                                       Rep=isol_end$Rep[i],
                                       Treat=isol_end$Treat[i],
                                       Isol=isol_end$Isol[i],
                                       Pop_Rate=my.sub$`Rate (cm/hr)`[1],
                                       Isol_Rate=isol_end$`Rate (cm/hr)`[i]))
  }
}
# plot(pop_isol_frame$Isol_Rate~pop_isol_frame$Pop_Rate, xlim = c(0.025, 0.065),
#      ylim = c(0.025, 0.065))
# lines(x=c(0.02, 0.07), y=c(0.02, 0.07))
# summary(lm(pop_isol_frame$Isol_Rate~pop_isol_frame$Pop_Rate))
my_colors <- c("red", "green", "blue")
my.pch = c(20, 20, 20)
tiff(filename = "isol_vs_pop.tiff", width = 10, height = 10, units = "in",
     compression = "none", res = 300)
par(mar = my.mar + c(0, 1, 0, 0))
for (i in 1:length(unique(pop_isol_frame$Treat))) {
  treat = unique(pop_isol_frame$Treat)[i]
  my.sub <- subset(pop_isol_frame, pop_isol_frame$Treat == treat)
  if (i == 1) {
    plot(my.sub$Isol_Rate~my.sub$Pop_Rate, xlim = c(0.027, 0.061), ylim = c(0.027, 0.061),
         pch = my.pch[i], xlab = "Population Migration Rate (cm/hr)", ylab = "Isolate Migration Rate (cm/hr)",
         col = my_colors[i], cex = 2, cex.lab = 2, cex.axis = 2)
  } else {
    points(my.sub$Isol_Rate~my.sub$Pop_Rate, pch = my.pch[i], col = my_colors[i],
           cex = 2)
  }
}
lines(x=c(0, 1), y=c(0, 1), lwd = 2)
legend(x = 0.053, y = 0.035, legend = unique(mean_rates$Treatment), pch = my.pch, 
       col = my_colors, cex = 2)
dev.off()

##isolate growth curve analysis
library("minpack.lm")
options(stringsAsFactors = F)
setwd("C:/Users/mikeb/Google Drive/Phage-Bacteria Project/Data/97-101_Growth_Curves")
growth_97 <- read.csv("97B.csv", header = T, stringsAsFactors = F)
growth_98 <- read.csv("98B.csv", header = T, stringsAsFactors = F)
growth_99 <- read.csv("99.csv", header = T, stringsAsFactors = F)
plate_layout <- read.csv("Plate_layout.csv", header = F, stringsAsFactors = F)

#make list of well ID & contents from plate layout
layout_list <- data.frame("Location" = character(), "Isol" = character())
for (i in 2:nrow(plate_layout)) {
  for (j in 2:ncol(plate_layout)) {
    layout_list <- rbind(layout_list, data.frame("Location" = paste(plate_layout[i, 1], plate_layout[1, j], sep = ""),
                                                 "Isol" = plate_layout[i, j]))
  }
}

#remove columns from growth data that correspond to water
to_rem <- subset(layout_list, layout_list$Isol == "H2O")
growth_97 <- growth_97[, !(colnames(growth_97) %in% to_rem$Location)]
growth_98 <- growth_98[, !(colnames(growth_98) %in% to_rem$Location)]
growth_99 <- growth_99[, !(colnames(growth_99) %in% to_rem$Location)]

#replace column names with Isols
isol.name.func <- function(my.array) {
  for (i in 1:ncol(my.array)) {
    j <- match(colnames(my.array)[i], layout_list$Location)
    if (!is.na(j)) {
      colnames(my.array)[i] <- as.character(layout_list$Isol[[j]])
    }
  }
  return(my.array)
}
growth_97 <- isol.name.func(growth_97)
growth_98 <- isol.name.func(growth_98)
growth_99 <- isol.name.func(growth_99)

#turn into tidy data
tidy.up <- function(my.array) {
  tidy.array <- data.frame("Isol"=character(), "Time"=character(),
                           "Temp"=character(), "OD600"=double())
  for (i in 3:ncol(my.array)) {
    for (j in 1:nrow(my.array)) {
      tidy.array <- rbind(tidy.array, data.frame("Isol" = colnames(my.array)[i],
                                           "Time" = my.array[j, 1],
                                           "Temp" = my.array[j, 2],
                                           "OD600" = my.array[j, i]))
    }
  }
  return(tidy.array)
}

tidy_97 <- tidy.up(growth_97)
tidy_98 <- tidy.up(growth_98)
tidy_99 <- tidy.up(growth_99)

#convert OD to CFU
tidy_97 <- cbind(tidy_97, "CFU"=(tidy_97$OD600-0.08560402)/(2.917486*10^(-10)))
tidy_98 <- cbind(tidy_98, "CFU"=(tidy_98$OD600-0.08560402)/(2.917486*10^(-10)))
tidy_99 <- cbind(tidy_99, "CFU"=(tidy_99$OD600-0.08560402)/(2.917486*10^(-10)))

#separate out ID info (add isol info)
sep.ID <- function(my.array, my.isol) {
  #my.isol is the isolate letter included in the run (assumes all the same) (eg. "A")
  my.array <- cbind("Media"=NA, "Proj"=NA, "Pop"=NA, "Treat"=NA,
                   "Rep"=NA, my.array)
  for (i in 1:nrow(my.array)) {
    my.split <- strsplit(as.character(my.array$Isol[i]), split = "-")
    my.array$Media[i] <- my.split[[1]][1]
    my.array$Proj[i] <- my.split[[1]][2]
    my.array$Rep[i] <- my.split[[1]][3]
    my.array$Isol[i] <- my.isol
  }
  my.array$Treat <- substr(my.array$Proj, nchar(my.array$Proj),
                          nchar(my.array$Proj))
  my.array$Proj <- substr(my.array$Proj, 1, nchar(my.array$Proj)-1)
  for (i in 1:nrow(my.array)) {
    if (my.array$Proj[i] == "") {
      my.array$Proj[i] <- "A"
      my.array$Pop[i] <- "Z"
    } else { #must be non-ancestor
      if (substr(my.array$Proj[i], 1, 2) == "74") {
        my.array$Pop[i] <- "A"
      } else if (substr(my.array$Proj[i], 1, 2) == "75") {
        if (substr(my.array$Proj[i], 3, 3) == "A") {
          my.array$Pop[i] <- "B"
        } else {
          my.array$Pop[i] <- "C"
        }
      } else { #must be project 76
        if (substr(my.array$Proj[i], 3, 3) == "A") {
          my.array$Pop[i] <- "D"
        } else {
          my.array$Pop[i] <- "E"
        } 
      }
      my.array$Proj[i] <- substr(my.array$Proj[i], 1, 2)
    }
  }
  return(my.array)
}

tidy_97 <- sep.ID(tidy_97, "A")
tidy_98 <- sep.ID(tidy_98, "B")
tidy_99 <- sep.ID(tidy_99, "C")

#combine all runs
gc_data <- rbind(tidy_97, tidy_98, tidy_99)
#reorganize columns
gc_data <- gc_data[c(1:4, 6, 5, 7, 8, 10)]

#prep for calculation
gc_data$Time <- strptime(gc_data$Time, format = "%H:%M:%S")
gc_data <- cbind(gc_data, "M.P.P.T.I.R"=paste(gc_data$Media, 
                                              gc_data$Proj,
                                              gc_data$Pop,
                                              gc_data$Treat,
                                              gc_data$Isol,
                                              gc_data$Rep,
                                              sep = "."),
                 "Smooth_CFU"=NA, "dCFU/hour"=NA)

#non-linear least-squares fit
lambda = 0.01 #error threshold for starting values
frac.change = 1
#starting values for fit: 50 start vals listed first
start.vals <- list(list(a=1059393503, b=1.778213, c=5.056632, 
                   d=2759201448, e=0.3280987, f=10.26474), 
                   list(a=1059393503, b=1.778213, c=5.056632, 
                       d=2759201448, e=0.3280987, f=10.26474))
while (frac.change > lambda) {
  #make frame for fit coeffs
  my.coef <- data.frame(Media = character(), Proj = character(), Pop = character(), 
                        Treat = character(), Isol = character(), Rep = character(),
                        a = double(), b = double(), c = double(), d = double(), e = double(),
                        f = double(), cor = double())
  #make list to save nls models
  nls.list <- list()
  #calculate coeffs for ea growth curve (using current starting values)
  for (my.run in unique(gc_data$M.P.P.T.I.R)) {
    my.sub <- subset(gc_data, gc_data$M.P.P.T.I.R == my.run)
    my.sub <- my.sub[order(my.sub$Time), ]
    my.data <- data.frame(x = as.numeric(format(my.sub$Time, "%H"))+
                            as.numeric(format(my.sub$Time, "%M"))/60, y = my.sub$CFU)
    my.data <- my.data[-(1:3), ]
    if (my.sub$Media[1] == 50) {
      my.nls <- nlsLM(y ~ (a/(1+exp(-b*(x-c))))+(d/(1+exp(-e*(x-f)))),
                      data = my.data, start = start.vals[[1]])
      
    } else if (my.sub$Media[1] == 100) {
      my.nls <- nlsLM(y ~ (a/(1+exp(-b*(x-c))))+(d/(1+exp(-e*(x-f)))),
                      data = my.data, start = start.vals[[2]])
    }
    # plot(my.data$x, my.data$y)
    # lines(my.data$x, predict(my.nls))
    #store coeffs
    my.coef <- rbind(my.coef, data.frame(Media = my.sub$Media[1], Proj = my.sub$Proj[1], 
                                         Pop = my.sub$Pop[1], Treat = my.sub$Treat[1], 
                                         Isol = my.sub$Isol[1], Rep = my.sub$Rep[1],
                                         a = my.nls$m$getPars()[1], b = my.nls$m$getPars()[2], 
                                         c = my.nls$m$getPars()[3], d = my.nls$m$getPars()[4], 
                                         e = my.nls$m$getPars()[5], f = my.nls$m$getPars()[6], 
                                         cor = cor(my.data$y, predict(my.nls))))
    nls.list <- c(nls.list, list(my.nls))
  }
  #change the starting value for ea coeff to be the cnt mean of each coeff
  past.vals <- start.vals
  change <- c()
  for (i in 1:2) {
    my.media <- c(50, 100)[i]
    my.sub = subset(my.coef, my.coef$Media == my.media)
    start.vals[[i]] <- list(a = median(my.sub$a), b = median(my.sub$b),
                            c = median(my.sub$c), d = median(my.sub$d),
                            e = median(my.sub$e), f = median(my.sub$f))
    change <-  c(change, 
      abs((as.numeric(start.vals[[i]])-as.numeric(past.vals[[i]]))/
            as.numeric(past.vals[[i]])))
  }
  #check if the change is small enough (reached local maxima)
  frac.change <- mean(change)
}


# stripchart(my.coef$f ~ paste(my.coef$Media, my.coef$Pop, my.coef$Treat,
#                            my.coef$Isol), vert = T)

#average values across replicate growth curves
my.coef.mean <- my.coef[0, -c(2, 6, 13)]
my.coef$M.P.T.I <- paste(my.coef$Media, my.coef$Pop, my.coef$Treat, my.coef$Isol, sep=".")
for (my.run in unique(my.coef$M.P.T.I)) {
  my.sub <- subset(my.coef, my.coef$M.P.T.I == my.run)
  my.coef.mean <- rbind(my.coef.mean, data.frame("Media" = my.sub$Media[1],
                                     "Pop" = my.sub$Pop[1],
                                     "Treat" = my.sub$Treat[1],
                                     "Isol" = my.sub$Isol[1],
                                     "a" = mean(my.sub$a),
                                     "b" = mean(my.sub$b),
                                     "c" = mean(my.sub$c),
                                     "d" = mean(my.sub$d),
                                     "e" = mean(my.sub$e),
                                     "f" = mean(my.sub$f),
                                     "M.P.T.I" = my.sub$M.P.T.I[1]))
}

#plot averaged values for each pop & treat for each coeff
for (i in 1:length(unique(my.coef.mean$Media))) {
  my.media <- unique(my.coef.mean$Media)[i]
  my.sub <- subset(my.coef.mean, my.coef.mean$Media == my.media)
  my.sub <- my.sub[order(my.sub$Treat, my.sub$Pop), ]
  for (j in 1:6) {
    tiff(filename = paste(my.media, letters[j], "v2_isol_rates.tiff", sep = "_"),
         width = 15, height = 7, units = "in",
         compression = "none", res = 300)
    stripchart(my.sub[, j+4] ~ paste(my.sub$Treat, my.sub$Pop), vert = T, pch = 19,
               cex = 2, cex.axis = 2, ylab = "",
               group.names = c("WT", "CA", "B", "C", "D", "E", "GA", "B", "C", "D",
                               "LA", "B", "C", "D"))
    dev.off()
  }
  # par(mar = my.mar + c(1, 3, 0, 0), xpd = T)
  # 
  # 
  # mtext(expression(paste(mu[max], " (hr"^"-1", ")")), side = 2, cex = 2.5, line = 4)
  # my.y <- c(0.05, 0.32)
  # text(c("Control", "Global", "Local"), x = c(4, 8.5, 12.5), y = my.y[i], cex = 2.5)
  # my.lwd <- 3
  # my.y <- c(0.175, 0.42)
  # lines(x = c(1.8, 6.2), y = c(my.y[i], my.y[i]), lwd = my.lwd)
  # lines(x = c(6.8, 10.2), y = c(my.y[i], my.y[i]), lwd = my.lwd)
  # lines(x = c(10.8, 14.2), y = c(my.y[i], my.y[i]), lwd = my.lwd)
  # dev.off()
}

##Isolate resistance analysis
setwd("C:/Users/mikeb/Google Drive/Phage-Bacteria Project/Data/74_75_76_Analysis")
resis_data <- read.csv("74_75_76_Plaquing.csv", header = T, stringsAsFactors = F)

#split out project info & mke unique reps
resis_data <- cbind(resis_data[, 1:2], "Pop" = NA, resis_data[, 3:5])
for (i in 1:nrow(resis_data)) {
  my.proj <- resis_data$Proj[i]
  if (my.proj == "P1.1") {
    resis_data$Pop[i] <- "A"
    resis_data$Treat[i] <- "A" #for Ancestor
  } else if (my.proj == "74") {
    resis_data$Pop[i] <- "A"
  } else {
    resis_data$Proj[i] <- substr(my.proj, 1, 2)
    resis_data$Pop[i] <- substr(my.proj, 3, 3)
  }
  if (resis_data$Proj[i] == "75") {
    if (resis_data$Pop[i] == "A") {resis_data$Pop[i] <- "B"
    } else {resis_data$Pop[i] <- "C"}
  } else if (resis_data$Proj[i] == "76") {
    if (resis_data$Pop[i] == "A") {resis_data$Pop[i] <- "D"
    } else {resis_data$Pop[i] <- "E"}
  }
}

#plot PFU for each pop
tiff(filename = "resist.tiff", width = 15, height = 7, units = "in",
     compression = "none", res = 300)
par(mar = my.mar + c(1, 3, 0, 0), xpd = T)
stripchart(resis_data$PFU ~ paste(resis_data$Treat, resis_data$Pop), vert = T, pch = 19,
           cex = 2, cex.axis = 2, ylab = "",
           group.names = c("WT", "A", "B", "C", "D", "E", "A", "B", "C", "D",
                           "A", "B", "C", "D"))
mtext("Ancestral PFU Observed", side = 2, cex = 2.5, line = 4)
my.y <- -34
text(c("Control", "Global", "Local"), x = c(4, 8.5, 12.5), y = my.y, cex = 2.5)
my.lwd <- 3
my.y <- -25
lines(x = c(1.8, 6.2), y = c(my.y, my.y), lwd = my.lwd)
lines(x = c(6.8, 10.2), y = c(my.y, my.y), lwd = my.lwd)
lines(x = c(10.8, 14.2), y = c(my.y, my.y), lwd = my.lwd)
dev.off()
