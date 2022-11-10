library(ggplot2)

#Read in modeling results data
data <- lapply(
  X = grep(pattern = "Analysis.*/.*_tot.csv", value = TRUE,
           x = list.files("./Modeling/", recursive = TRUE, full.names = TRUE)),
  FUN = read.csv)
#Add names
names(data) = 
  gsub(pattern = "\\./Modeling//", replacement = "",
       x = grep(pattern = "Modeling//(Analysis.*/.*_tot.csv)", value = TRUE,
                x = list.files("./Modeling/", 
                               recursive = TRUE, full.names = TRUE)))

#Calculate "relative resistance" col
data <- lapply(X = data,
               FUN = function(x) {
                 cbind(x, data.frame("relativeR" = 1/x[, "relativeI"]))})

#Add cols for phage distribution and variables manipulated
for (i in 1:length(data)) {
  if(strsplit(names(data)[i], split = "/")[[1]][1] == "Analysis") {
    data[[i]] <- cbind(data.frame(distrib = "global"), data[[i]])
  } else if (
    strsplit(names(data)[i], split = "/")[[1]][1] == "Analysis_GaussPhage") {
    data[[i]] <- cbind(data.frame(distrib = "local"), data[[i]])
  } else if (
    strsplit(names(data)[i], split = "/")[[1]][1] == "Analysis_GaussPhage_Wide") {
    data[[i]] <- cbind(data.frame(distrib = "global_gauss"), data[[i]])
  } else if (
    strsplit(names(data)[i], split = "/")[[1]][1] == "Analysis_NoPhage") {
    data[[i]] <- cbind(data.frame(distrib = "no_paras"), data[[i]])
  }
  
  data[[i]] <- cbind(
    data.frame(
      vars_manip_1 = strsplit(strsplit(names(data)[i], split = "/")[[1]][2],
                              split = "_")[[1]][1],
      vars_manip_2 = strsplit(strsplit(names(data)[i], split = "/")[[1]][2],
                              split = "_")[[1]][3],
      vars_manip = paste(strsplit(strsplit(names(data)[i], split = "/")[[1]][2],
                                  split = "_")[[1]][1:3], collapse = "_")),
    
    data[[i]])
}

#Merge data
library(magrittr)
library(dplyr)
library(purrr)
data_mrg <- data %>% reduce(full_join)
data_mrg$Cell2_Cell1 <- data_mrg$Cell2_population/data_mrg$Cell_population

#Make plots
for (vars_manip in unique(data_mrg$vars_manip)) {
  myrows <- which(data_mrg$vars_manip == vars_manip)
  my_data <- data_mrg[myrows, ]
  
  #Get limits for log10(Cell2/Cell) across the dift initial distributions
  mybreaks <- pretty(x = log10(my_data$Cell2_Cell1), n = 6)
  
  #axis labels labeller for the non-resistance axis
  mylabeller <-
    as_labeller(c(relativecA = "Relative Attractant Consumption",
                  relativeChi = "Relative Chemotactic Sensitivity",
                  relativecR = "Relative Growth Rate",
                  relativeY = "Relative Growth Yield"))
  
  #facet labels labeller
  mylabller <-
    as_labeller(c(global = "Global Parasite Distribution",
                  local = "Local Parasite Distribution",
                  global_gauss = "Global Parasite Distribution (Gaussian)",
                  no_paras = "No Parasites"))
  
  #Make file
  png(
    paste(sep = "", "./Model_plots/", vars_manip, ".png"),
    width = 6, height = 5, units = "in", res = 150)
  #Make plot
  print(
    ggplot(my_data,
           aes_string(x = "relativeR", 
                      y = paste("relative", my_data$vars_manip_2[1], sep = ""))) +
      geom_contour_filled(aes(z = log10(Cell2_population/Cell_population)),
                          breaks = mybreaks) +
      geom_point(aes(color = log10(Cell2_population/Cell_population)),
                 size = 3) +
      facet_wrap(~distrib) +
      scale_x_continuous(trans = "log10") +
      scale_y_continuous(trans = "log2") +
      scale_fill_viridis_d(drop = FALSE) +
      scale_color_continuous(name = "Fitness", type = "viridis",
                             limits = c(min(mybreaks), max(mybreaks))) +
      #guides(fill = "none") +
      # labs(subtitle = mylabel,
      #      x = "Relative Resistance",
      #      y = mylabeller(grep("relative[^IR]+", colnames(data[[i]]),
      #                          value = TRUE))[[1]]) +
      theme_bw() +
      NULL
  )
  dev.off()
  
  ##Now points only (no landscape)
  #Make file
  # png(
  #   paste(sep = "", "./Model_plots/", vars_pattern, mylabel, 
  #         "_pointsonly.png"),
  #   width = 6, height = 5, units = "in", res = 150)
  # 
  # #Make plot
  # print(
  #   ggplot(data[[i]],
  #          aes_string(x = "relativeR",
  #                     y = grep("relative[^IR]+", colnames(data[[i]]),
  #                              value = TRUE))) +
  #     # geom_contour_filled(aes(z = log10(Cell2_population/Cell_population)),
  #     #                     breaks = mybreaks) +
  #     geom_point(aes(color = log10(Cell2_population/Cell_population)),
  #                size = 3) +
  #     scale_x_continuous(trans = "log10") +
  #     scale_y_continuous(trans = "log10") +
  #     scale_fill_viridis_d(drop = FALSE) +
  #     scale_color_continuous(name = "Fitness", type = "viridis",
  #                            limits = c(min(mybreaks), max(mybreaks))) +
  #     guides(fill = "none") +
  #     labs(subtitle = mylabel,
  #          x = "Relative Resistance",
  #          y = mylabeller(grep("relative[^IR]+", colnames(data[[i]]),
  #                              value = TRUE))[[1]]) +
  #     theme_bw() +
  #     NULL
  # )
  # dev.off()
  
  #Plot change in fitness landscape of local relative to global
  
  #First calculate change in frequency bt two landscapes
  # glob_idx <- grep("uniform", names(data)[these_vars_idx])
  # loc_idx <- grep("gPhage", names(data)[these_vars_idx])
  # 
  # glob_temp <- data[[these_vars_idx[glob_idx]]]
  # glob_temp$inv_fitness_glob <- 
  #   glob_temp$Cell2_population/glob_temp$Cell_population
  # 
  # loc_temp <- data[[these_vars_idx[loc_idx]]]
  # loc_temp$inv_fitness_loc <-
  #   loc_temp$Cell2_population/loc_temp$Cell_population
  # 
  # #Drop extra columns
  # glob_temp  <- 
  #   glob_temp[, 
  #             c(match(c("relativeR", "inv_fitness_glob"),
  #                     colnames(glob_temp)),
  #               grep("relative[^IR]+", colnames(glob_temp)))]
  # loc_temp  <- 
  #   loc_temp[, 
  #             c(match(c("relativeR", "inv_fitness_loc"),
  #                     colnames(loc_temp)),
  #               grep("relative[^IR]+", colnames(loc_temp)))]
  # 
  # temp <- dplyr::full_join(glob_temp, loc_temp)
  # temp$delta_freq <- temp$inv_fitness_glob - temp$inv_fitness_loc
  # 
  # #Make file
  # png(
  #   paste(sep = "",
  #         "./Model_plots/",
  #         paste(collapse = "_", strsplit(names(data)[i], split = "_")[[1]][1:3]), 
  #         "_delta.png"),
  #   width = 6, height = 5, units = "in", res = 150)
  # #Make plot
  # print(
  #   ggplot(temp,
  #          aes_string(x = "relativeR",
  #                     y = grep("relative[^IR]+", colnames(temp),
  #                              value = TRUE))) +
  #     geom_contour_filled(aes(z = delta_freq)) +
  #     #geom_point(aes(color = delta_freq), size = 3) +
  #     scale_x_continuous(trans = "log10") +
  #     scale_y_continuous(trans = "log10") +
  #     scale_fill_brewer(type = "div") +
  #     #scale_color_gradient2(name = "Change in Fitness") +
  #     #guides(fill = FALSE) +
  #     labs(x = "Relative Resistance",
  #          y = mylabeller(grep("relative[^IR]+", colnames(temp),
  #                              value = TRUE))[[1]]) +
  #     theme_bw() +
  #     NULL
  # )
  # dev.off()
}

