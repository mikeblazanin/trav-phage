library("reshape")

##Experimental evolution cleanup ----
end_data_7x <- read.csv("./Raw_Data/74_75_76_Data_Measurements.csv", header = T, stringsAsFactors = F)
start_data_7x <- read.csv("./Raw_Data/74_75_76_Start_Data.csv", header = T, stringsAsFactors = F)
fails_7x <- read.csv("./Raw_Data/74_75_76_fails.csv", header = T, stringsAsFactors = F)
contam_7x <- read.csv("./Raw_Data/74_75_76_contaminants.csv", header = T, stringsAsFactors = F)
end_data_125 <- read.csv("./Raw_Data/125_Data_Measurements.csv", header = T, stringsAsFactors = F)
start_data_125 <- read.csv("./Raw_Data/125_Start_Data.csv", header = T, stringsAsFactors = F)
fails_125 <- read.csv("./Raw_Data/125_fails.csv", header = T, stringsAsFactors = F)

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
  my_data$Pop <- NA
  if (num_cols == 4) {my_data$Treat <- NA}
  for (i in 1:nrow(my_data)) {
    my_split <- strsplit(my_data$Strain[i], split = "_")[[1]]
    my_data$Proj[i] <- my_split[[1]]
    my_data$Time[i] <- my_split[[2]]
    if (length(my_split) < num_cols) {
      my_data$Pop[i] <- "A"
    } else {
      my_data$Pop[i] <- my_split[[num_cols]]
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

#Change Rep to Pop
colnames(fails_data)[3] <- "Pop"
colnames(contam_7x)[3] <- "Pop"

#define functions to remove failed runs
start_remove_fails <- function(start_data, fails) {
  start_data <- cbind(start_data, NA)
  colnames(start_data)[ncol(start_data)] <- "Remove"
  for (i in 1:nrow(start_data)) {
    for (j in 1:nrow(fails)) {
      if (start_data$Proj[i] == fails$Proj[j] & 
          start_data$Time[i] == fails$Time[j] & 
          start_data$Pop[i] == fails$Pop[j]) {
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
      if (end_data$Proj[i] == fails$Proj[j] & 
          end_data$Time[i] == fails$Time[j] & 
          end_data$Pop[i] == fails$Pop[j]) {
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
      if (end_data$Proj[i] == contam$Proj[j] &
          end_data$Treat[i] == contam$Treat[j] &
          end_data$Pop[i] == contam$Pop[j]) {
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

#define function to change the pops so they're unique (74-A, 75A-B, 75B-C, 76A-D, 76B-E)
#& change projects to 1 (for 74, 75, 76) & 2 (for 125) 
uniq_pops <- function(my_data) {
  for (i in 1:nrow(my_data)) {
    if (my_data$Proj[i] == "74") {
      my_data$Pop[i] <- "A"
    } else if (my_data$Proj[i] == "75") {
      if (my_data$Pop[i] == "A") {
        my_data$Pop[i] <- "B"
      } else if (my_data$Pop[i] == "B") {
        my_data$Pop[i] <- "C"
      }
    } else if (my_data$Proj[i] == "76") {
      if (my_data$Pop[i] == "A") {
        my_data$Pop[i] <- "D"
      } else if (my_data$Pop[i] == "B") {
        my_data$Pop[i] <- "E"
      }
    }
    if (my_data$Proj[i]=="74" | my_data$Proj[i]=="75" | my_data$Proj[i]=="76") {
      my_data$Proj[i] <- "7x"
    } else if (my_data$Proj[i] == "125") {my_data$Proj[i] <- "125"}
  }
  return(my_data)
}
  
#make start & end data have unique pops
start_data <- uniq_pops(start_data)
end_data <- uniq_pops(end_data)

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
make_timestamp <- function(my_data, type) {
  my_data <- cbind(my_data, NA)
  colnames(my_data)[ncol(my_data)] <- "timestamp"
  my_data$timestamp <- paste(paste(my_data$Year, my_data$Month, my_data$Day, sep = "-"),
                                paste(my_data$Hour, my_data$Minute, sep = ":"), sep = " ")
  my_data$timestamp <- strptime(my_data$timestamp, format = "%Y-%m-%d %H:%M")
  # if (type=="start") {my_data <- my_data[, -c(1:6)]
  # } else if (type=="end") {my_data <- my_data[, -c(2:8)]
  # } else if (type=="isol start") {my_data <- my_data[, -c(1:5)]}
  return(my_data)
}

#compile time information into timestamp
start_data <- make_timestamp(start_data, "start")
end_data <- make_timestamp(end_data, "end")

#Put both timestamps in end_data
colnames(end_data)[ncol(end_data)] <- "end_timestamp"
start_data$Strain_clean <- paste(start_data$Proj,
                                 start_data$Pop,
                                 start_data$Time, sep = "_")
end_data$Strain_clean <- paste(end_data$Proj,
                               end_data$Pop,
                               end_data$Time, sep = "_")
end_data$start_timestamp <- start_data$timestamp[match(end_data$Strain_clean,
                                                       start_data$Strain_clean)]
end_data$time_since_inoc <- NA
for (i in 1:nrow(end_data)) {
  end_data$time_since_inoc <- difftime(end_data$end_timestamp,
                                end_data$start_timestamp,
                                units = "hours")
}

#Tidy up & drop unneccesary columns
end_data <- end_data[, c("Proj", "Pop", "Treat", "Time",
                         "start_timestamp", "end_timestamp", "time_since_inoc",
                         "Width..cm.", "Height..cm.")]
end_data$start_timestamp <- as.character(end_data$start_timestamp)
end_data$end_timestamp <- as.character(end_data$end_timestamp)
colnames(end_data)[c(4, 8, 9)] <- c("Timepoint", "Width_cm", "Height_cm")

write.csv(end_data, "./Clean_Data/Experimental_evolution_growth.csv",
          row.names = FALSE)

##isolate migration cleanup ----

isol_end_7x <- read.csv("./Raw_Data/74_75_76_isol_motility.csv", header = T, stringsAsFactors = F)
isol_start_7x <- read.csv("./Raw_Data/74_75_76_isol_motility_start.csv", header = T, stringsAsFactors = F)
isol_end_125 <- read.csv("./Raw_Data/131_125_isol_migration.csv", header = T, stringsAsFactors = F)
isol_start_125 <- read.csv("./Raw_Data/131_125_isol_migration_start.csv", header = T, stringsAsFactors = F)

#define function to split Strain data
isol_end_strain_split <- function(isol_end) {
  isol_end <- cbind(isol_end, NA, NA, NA, NA, NA)
  colnames(isol_end)[(ncol(isol_end)-4):ncol(isol_end)] <- c("Proj", "Pop", "Treat", "Isol", "Media")
  for (i in 1:nrow(isol_end)) {
    my_split <- strsplit(isol_end$Strain[i], split = "_")[[1]]
    isol_end$Proj[i] <- my_split[[1]]
    if (isol_end$Proj[i] == "P1.1") {
      isol_end$Pop[i] <- "A"
      isol_end$Treat[i] <- "A" #for Ancestor
      isol_end$Isol[i] <- my_split[[2]]
      isol_end$Media[i] <- "+Mg"
    } else if (isol_end$Proj[i] == "74") {
      isol_end$Pop[i] <- "A"
      isol_end$Treat[i] <- my_split[[2]]
      isol_end$Isol[i] <- my_split[[3]]
      if (length(my_split) > 3) {
        isol_end$Media[i] <- my_split[[4]]
      } else {
        isol_end$Media[i] <- "+Mg"
      }
    } else {
      isol_end$Pop[i] <- my_split[[2]]
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

#Fix Ancestor
my_rows <- which(isol_end_7x$Proj == "P1.1")
isol_end_7x$Pop[my_rows] <- "Anc"
isol_end_7x$Treat[my_rows] <- "Anc"
isol_end_7x$Isol[my_rows] <- "Anc"
isol_end_7x$Time[my_rows] <- "0"
my_rows <- which(isol_end_125$Population == "Anc")
isol_end_125$Treatment[my_rows] <- "Anc"
isol_end_125$Isolate[my_rows] <- "Anc"
isol_end_125$Time[my_rows] <- 0

#make Reps unique (74-A, 75A-B, 75B-C, 76A-D, 76B-E)
#& change project to be 1 (74, 75, 76) or 2 (125)
isol_end_7x <- uniq_pops(isol_end_7x)

#Add media info to 125 end
isol_end_125$Media <- "+Mg"

#Standardize columns/format in start dataframes
isol_start_7x$Isol[isol_start_7x$Proj == "P1.1"] <- "Anc"
isol_start_125$Isol[isol_start_125$Pop == "Anc"] <- "Anc"

isol_start_7x$Time[isol_start_7x$Isol == "Anc"] <- 0
isol_start_125$Time[isol_start_125$Isol == "Anc"] <- 0

isol_start_7x$Proj <- "7x"
isol_start_125$Proj <- "125"

isol_start_7x <- isol_start_7x[, c("Year", "Month", "Day", "Hour", "Minute",
                                   "Proj", "Isol", "Time")]
colnames(isol_start_7x) <- c("Year", "Month", "Day", "Hour", "Minute",
                             "Proj", "Isol", "Timepoint")
isol_start_125 <- isol_start_125[, c("Year", "Month", "Day", "Hour", "Minute",
                                    "Proj", "Isol", "Time")]
colnames(isol_start_125) <- c("Year", "Month", "Day", "Hour", "Minute",
                             "Proj", "Isol", "Timepoint")

isol_start <- rbind(isol_start_7x, isol_start_125)
isol_start <- unique(isol_start)
isol_start <- isol_start[order(isol_start$Proj,
                               isol_start$Isol), ]

#Standardize columns/format in end dataframes
isol_end_7x <- isol_end_7x[, c("Year", "Month", "Day", "Hour", "Minute", "Second",
                               "Proj", "Pop", "Treat", "Time", "Isol",
                               "Width..cm.", "Height..cm.", "Media")]
colnames(isol_end_7x) <- c("Year", "Month", "Day", "Hour", "Minute", "Second",
                           "Proj", "Pop", "Treat", "Timepoint", "Isol",
                           "Width_cm", "Height_cm", "Media")
isol_end_125 <- isol_end_125[, c("Year", "Month", "Day", "Hour", "Minute", "Second",
                               "Project", "Population", "Treatment", "Time", "Isolate",
                               "Width..cm.", "Height..cm.", "Media")]
colnames(isol_end_125) <- c("Year", "Month", "Day", "Hour", "Minute", "Second",
                           "Proj", "Pop", "Treat", "Timepoint", "Isol",
                           "Width_cm", "Height_cm", "Media")

isol_end_7x$Proj <- "7x"
isol_end_125$Proj <- "125"

isol_end <- rbind(isol_end_7x, isol_end_125)
                                   
#compile time information into timestamp
isol_end <- make_timestamp(isol_end, type = "end")
isol_start <- make_timestamp(isol_start, type = "isol start")

#Put both timestamps in isol_end
colnames(isol_end)[ncol(isol_end)] <- "end_timestamp"
isol_start$start_date <- paste(as.Date(isol_start$timestamp), isol_start$Isol,
                               sep = "_")
isol_end$start_date <- paste(as.Date(isol_end$end_timestamp-24*60*60),
                             isol_end$Isol, sep = "_")
isol_end$start_timestamp <- isol_start$timestamp[match(isol_end$start_date,
                                                       isol_start$start_date)]
isol_end$time_since_inoc <- NA
for (i in 1:nrow(isol_end)) {
  isol_end$time_since_inoc <- difftime(isol_end$end_timestamp,
                                       isol_end$start_timestamp,
                                       units = "hours")
}

#Drop rows with -Mg Media
isol_end <- isol_end[isol_end$Media == "+Mg", ]

#Tidy up & drop unneccesary columns
isol_end <- isol_end[, c("Proj", "Pop", "Treat", "Timepoint", "Isol",
                         "start_timestamp", "end_timestamp", "time_since_inoc",
                         "Width_cm", "Height_cm")]
write.csv(isol_end, "./Clean_Data/Isolate_migration.csv",
          row.names = FALSE)

##isolate growth curve cleanup ----

#in 7x, ancestors are originally coded as pop F, treat A, 
# and isol (according to block they were run in)
#in 125, ancestors are originally coded as pop Anc, treat Anc, 
# and isol (according to block they were run in)
#After cleanup, they're all pop Anc, treat Anc, isol Anc
# (with Date denoting the different experimental blocks)

#Standard curves ----

#Travisano lab standard curve:
# (note: I didn't actually plate out from the plate standard curve, instead
# I just measured OD in the plate and in the spec and since I knew the spec
# standard curve I just connected the two)

trav_stan_spec <- read.csv("./Raw_Data/107_Trav_Spec.csv", stringsAsFactors = F, header = T)
trav_stan_plate <- read.csv("./Raw_Data/107_Std_Curve.csv", stringsAsFactors = F, header = T)
trav_stan_plate_layout <- read.csv("./Raw_Data/107_Plate_Layout.csv", stringsAsFactors = F, header = F)

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

#melt trav std curve plate arrays, combine & exclude empty wells
trav_stan_layout_list <- get_layout(trav_stan_plate_layout, 2, 2)
trav_stan_plate_list <- get_layout(trav_stan_plate, 1, 2)
colnames(trav_stan_plate_list)[2] <- "OD600"
trav_stan_plate_list$Dilution <- trav_stan_layout_list$Contents
trav_stan_plate_list <- subset(trav_stan_plate_list, 
                               is.na(trav_stan_plate_list$Dilution) == F)

# #in data from Trav-lab spec, convert OD to CFU
trav_stan_spec$cfu_ml <- (trav_stan_spec$OD600..on.Trav.lab.spec.+3.318*10^(-3))/(1.252*10^(-9))

#transfer CFU data to plate list
trav_stan_plate_list <- cbind(trav_stan_plate_list, 
                              "cfu_ml"=trav_stan_spec$cfu_ml[match(trav_stan_plate_list$Dilution, 
                                                     trav_stan_spec$Dilution)])

#get regression of data
# ggplot(stan_plate_list, aes(x = cfu_ml, y = as.numeric(OD600))) + geom_point()
trav_stan_plate_list$OD600 <- as.numeric(trav_stan_plate_list$OD600)
trav_std_curve <- lm(cfu_ml ~ OD600, data = trav_stan_plate_list)

#Turner lab standard curve:
turn_std_cv_data <- data.frame("Dil" = 1*1/(2**(0:8)),
                               "Abs_cuv" = c(2.109, 1.098, 0.450, 0.224, 0.108, 0.050,
                                             0.025, 0.010, 0.004),
                               "Abs_plt" = c(0.8319, 0.4711, 0.2591, 0.176, 0.1263,
                                             0.1009, 0.0905, 0.085, 0.0825))
#Change colname to OD600 because that's colname in growthcurve data
colnames(turn_std_cv_data)[3] <- "OD600"
turn_std_cv_data$cfu_ml <- turn_std_cv_data$Dil * 1305600000

turn_std_curve <- lm(cfu_ml ~ OD600, data = turn_std_cv_data)

##growth curve analysis
#library("minpack.lm")
library("tidyr")
#library("lubridate")
library("dplyr")

##7x growth curve cleanup ----

growth_97 <- read.csv("./Raw_Data/7xT14isol_growthcurves/97B.csv", 
                      header = T, stringsAsFactors = F)
growth_98 <- read.csv("./Raw_Data/7xT14isol_growthcurves/98B.csv", 
                      header = T, stringsAsFactors = F)
growth_99 <- read.csv("./Raw_Data/7xT14isol_growthcurves/99.csv", 
                      header = T, stringsAsFactors = F)
growth_101 <- read.csv("./Raw_Data/7xT14isol_growthcurves/101.csv", 
                       header = T, stringsAsFactors = F)
plate_layout_97_101 <- read.csv("./Raw_Data/7xT14isol_growthcurves/97_98_99_100_101_plate_layout.csv", 
                                header = F, stringsAsFactors = F)

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
tidy_97$Proj <- "7x"
tidy_98$Proj <- "7x"
tidy_99$Proj <- "7x"
tidy_101$Proj <- "7x"

#Combine data frames
gc_data_7x <- rbind(tidy_97, tidy_98, tidy_99, tidy_101)

#convert OD to cfu/mL
gc_data_7x$cfu_ml <- predict(trav_std_curve, newdata = gc_data_7x)

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

gc_data_7x <- sep.growth.ID(gc_data_7x)

#reformat time column
colnames(gc_data_7x)[which(colnames(gc_data_7x) == "Time")] <- "Time_s"
gc_data_7x$Time_s <- strptime(gc_data_7x$Time_s, format = "%H:%M:%S")
gc_data_7x$Time_s <- difftime(gc_data_7x$Time_s,
                            strptime("0:00:00", format = "%H:%M:%S"), 
                            units = "secs")

#Rename temperature column
colnames(gc_data_7x)[2] <- "Temp_C"

#Add Date column for experimental blocks
gc_data_7x$Date <- paste("2017-", gc_data_7x$Isol, sep = "")

#Adjust Ancestor wells
my_rows <- which(gc_data_7x$Treat == "A")
gc_data_7x$Pop[my_rows] <- "Anc"
gc_data_7x$Treat[my_rows] <- "Anc"
gc_data_7x$Isol[my_rows] <- "Anc"

#125 growth curve cleanup ----
growth_125 <- list(rep(NA, 5))
i <- 1
for (fil in list.files("./Raw_Data/129_125T14isol_growthcurves/")) {
  growth_125[[i]] <- read.csv(paste("./Raw_Data/129_125T14isol_growthcurves/", fil, sep = ""),
                              header = T, stringsAsFactors = F)
  growth_125[[i]] <- cbind(data.frame("Isol" = substr(fil, nchar(fil)-4, nchar(fil)-4),
                                      "Date" = strsplit(fil, split = "_")[[1]][2]),
                           growth_125[[i]])
  i <- i + 1
}

for (i in 1:length(growth_125)) {
  #Drop last row ("Date of measurement")
  growth_125[[i]] <- growth_125[[i]][-c(nrow(growth_125[[i]])), ]
  #Rename columns
  colnames(growth_125[[i]])[3:4] <- c("Time_s", "Temp_C")
  #strip "s" off of time
  growth_125[[i]][, "Time_s"] <- gsub("s", "", growth_125[[i]][, "Time_s"])
  #strip degrees C off of Temp_C
  growth_125[[i]][, "Temp_C"] <- gsub(" °C", "", growth_125[[i]][, "Temp_C"])
}

#Get layout information
isol_layout_125 <- read.csv("./Raw_data/129_isolate_layout.csv",
                            header = F, stringsAsFactors = F)
media_layout_125 <- read.csv("./Raw_data/129_media_layout.csv",
                             header = F, stringsAsFactors = F)

isol_layout_125 <- get_layout(isol_layout_125, 2, 2)
media_layout_125 <- get_layout(media_layout_125, 2, 2)

isol_layout_125$Contents[isol_layout_125$Contents == ""] <- NA
isol_layout_125 <- cbind(isol_layout_125, data.frame("Pop" = NA, "Treat" = NA, "Rep_Well" = NA))
my_split <- strsplit(isol_layout_125$Contents, "-")
for (i in 1:length(my_split)) {
  if (length(my_split[[i]]) > 1) {
    isol_layout_125$Pop[i] <- my_split[[i]][1]
    isol_layout_125$Treat[i] <- my_split[[i]][2]
    isol_layout_125$Rep_Well[i] <- my_split[[i]][3]
  }
}
isol_layout_125 <- isol_layout_125[, c("Well", "Pop", "Treat", "Rep_Well")]

media_layout_125$Contents[media_layout_125$Contents == "H2O"] <- NA
media_layout_125$Contents[media_layout_125$Contents == "Orig"] <- "25-50"
media_layout_125$Contents[media_layout_125$Contents == "Rich"] <- "50-100"
colnames(media_layout_125)[2] <- "Media"

#Melt density data, then join with media & isolate information
#Then drop wells with missing info (e.g. empty wells)
for (i in 1:length(growth_125)) {
  growth_125[[i]] <- pivot_longer(growth_125[[i]],
                                  cols = -c("Date", "Isol", "Time_s", "Temp_C"),
                                  names_to = "Well",
                                  values_to = "OD600")
  growth_125[[i]] <- left_join(growth_125[[i]],
                               isol_layout_125,
                               by = "Well")
  growth_125[[i]] <- left_join(growth_125[[i]],
                               media_layout_125,
                               by = "Well")
  growth_125[[i]] <- growth_125[[i]][complete.cases(growth_125[[i]]), ]
}

#Merge 125 data into one dataframe
gc_data_125 <- growth_125[[1]]
for (i in 2:length(growth_125)) {
  gc_data_125 <- rbind(gc_data_125, growth_125[[i]])
}

#Adjust Ancestor Isols
gc_data_125$Isol <- as.character(gc_data_125$Isol)
gc_data_125$Isol[gc_data_125$Pop == "Anc"] <- "Anc"

#Add project
gc_data_125$Proj <- "125"

#Add cfu/ml
gc_data_125$cfu_ml <- predict(turn_std_curve, newdata = gc_data_125)

#Reorder columns
gc_data_125 <- gc_data_125[, c("Date", "Proj", "Pop", "Treat", "Isol", "Rep_Well", 
                               "Media", "Time_s", "Temp_C", "OD600", "cfu_ml")]
gc_data_7x <- gc_data_7x[, c("Date", "Proj", "Pop", "Treat", "Isol", "Rep_Well", 
                             "Media", "Time_s", "Temp_C", "OD600", "cfu_ml")]

#Merge
gc_data_all <- rbind(gc_data_7x, gc_data_125)

#Output
write.csv(gc_data_all, "./Clean_Data/Isolate_growth_curves.csv",
          row.names = F)


##Isolate resistance cleanup ----
resis_7x_old <- read.csv("./Raw_Data/74_75_76_Plaquing.csv", header = T, stringsAsFactors = F)
resis_new <- read.csv("./Raw_Data/130_new_resis_assays.csv", header = T, stringsAsFactors = F)

#split out project info & make unique reps
resis_7x_old <- cbind(resis_7x_old[, 1:2], "Pop" = NA, resis_7x_old[, 3:5])
for (i in 1:nrow(resis_7x_old)) {
  my.proj <- resis_7x_old$Proj[i]
  if (my.proj == "P1.1") {
    resis_7x_old$Pop[i] <- "F"
    resis_7x_old$Treat[i] <- "A" #for Ancestor
  } else if (my.proj == "74") {
    resis_7x_old$Pop[i] <- "A"
  } else {
    resis_7x_old$Proj[i] <- substr(my.proj, 1, 2)
    resis_7x_old$Pop[i] <- substr(my.proj, 3, 3)
  }
  if (resis_7x_old$Proj[i] == "75") {
    if (resis_7x_old$Pop[i] == "A") {resis_7x_old$Pop[i] <- "B"
    } else {resis_7x_old$Pop[i] <- "C"}
  } else if (resis_7x_old$Proj[i] == "76") {
    if (resis_7x_old$Pop[i] == "A") {resis_7x_old$Pop[i] <- "D"
    } else {resis_7x_old$Pop[i] <- "E"}
  }
}
resis_7x_old$Proj <- "7x"

###Make two dataframes compatible for merge

##With old data
resis_7x_old$Date <- paste("2017-", resis_7x_old$Isol, sep = "")
resis_7x_old$Timepoint <- 14

#Reformat Ancestor rows
resis_7x_old[resis_7x_old$Treat == "A", c("Pop", "Treat", "Isol")] <- "Anc"
resis_7x_old[resis_7x_old$Treat == "Anc", c("Timepoint")] <- 0

#We used 45.5 uL of -7 dilution for old resistance assays, so PFU/mL needs
# to be scaled accordingly
resis_7x_old$dilution <- 10**7 * 1000/45.5
resis_7x_old$pfu_ml <- resis_7x_old$PFU * resis_7x_old$dilution

#Drop unnecessary columns
resis_7x_old <- resis_7x_old[, c("Date", "Proj", "Pop", "Treat", "Timepoint",
                                 "Isol", "PFU", "dilution", "pfu_ml")]

##With new data
for (i in 8:16) {
  resis_new[resis_new[, i] == "clr" | resis_new[, i] == "...", i] <- "clear"
  resis_new[resis_new[, i] == "too", i] <- "too many"
}

#Calculate pfu & dilution
resis_new$PFU <- NA
resis_new$dilution <- NA
resis_new$pfu_ml <- NA

for (i in 1:nrow(resis_new)) {
  for (j in 8:16) {
    if (!is.na(suppressWarnings(as.numeric(resis_new[i, j])))) {
      resis_new$PFU[i] <- resis_new[i, j]
      #100 fold for 10 uL, times the dilution of the tube
      resis_new$dilution[i] <- 100*10**as.numeric( 
        substr(colnames(resis_new)[j],
               nchar(colnames(resis_new)[j]),
               nchar(colnames(resis_new)[j])))
      break
    }
  }
}
resis_new$pfu_ml <- as.numeric(resis_new$PFU) * resis_new$dilution

#Drop unneeded columns
resis_new <- resis_new[, c("Date", "Proj", "Pop", "Treat", "Timepoint",
                                 "Isol", "PFU", "dilution", "pfu_ml")]

#Merge
resis_data <- rbind(resis_7x_old, resis_new)

#Write
write.csv(resis_data,
          "./Clean_Data/Isolate_resistance.csv",
          row.names = FALSE)
