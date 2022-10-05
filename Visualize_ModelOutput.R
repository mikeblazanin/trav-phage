library(ggplot2)

#Read in modeling results data
data <- lapply(
  X = grep(pattern = "Analysis.*/.*csv", value = TRUE,
           x = list.files("./Modeling/", recursive = TRUE, full.names = TRUE)),
  FUN = read.csv)
names(data) = 
  gsub(pattern = "\\./Modeling//", replacement = "",
       x = grep(pattern = "Modeling//(Analysis.*/.*csv)", value = TRUE,
       x = list.files("./Modeling/", recursive = TRUE, full.names = TRUE)))

data <- lapply(X = data,
               FUN = function(x) {
                 cbind(x, data.frame("relativeR" = 1/x[, "relativeI"]))})

#Make plots
for (i in 1:length(data)) {
  #Get initial distribution information
  if(strsplit(names(data)[i], split = "/")[[1]][1] == "Analysis") {
    mylabel = "Global Parasite Distribution"
  } else if (
    strsplit(names(data)[i], split = "/")[[1]][1] == "Analysis_GaussPhage") {
    mylabel = "Local Parasite Distribution"
  } else if (
    strsplit(names(data)[i], split = "/")[[1]][1] == "Analysis_GaussPhage_Wide") {
    mylabel = "Global Gaussian Parasite Distribution"
  }
  
  #Get limits for log10(Cell2/Cell) across the two dift initial distributions
  vars_pattern <- 
    paste(strsplit(strsplit(names(data)[i], split = "[\\/\\.]")[[1]][2],
                   "_")[[1]][1:3], collapse = "_")
  these_vars_idx <-  grep(pattern = vars_pattern, x = names(data))
  
  obs_vals <- 
    unlist(lapply(data[these_vars_idx],
                  function(x) {x[, "Cell2_population"]/x[, "Cell_population"]}))
  mybreaks <- pretty(x = log10(obs_vals), n = 6)
  
  #axis labels labeller for the non-resistance axis
  mylabeller <-
    as_labeller(c(relativecA = "Relative Attractant Consumption",
                  relativeChi = "Relative Chemotactic Sensitivity",
                  relativecR = "Relative Growth Rate",
                  relativeY = "Relative Growth Yield"))
  
  #Make file
  png(
    paste(sep = "", "./Model_plots/", vars_pattern, mylabel, ".png"),
    width = 6, height = 5, units = "in", res = 150)
  #Make plot
  print(
    ggplot(data[[i]],
           aes_string(x = "relativeR",
                      y = grep("relative[^IR]+", colnames(data[[i]]),
                               value = TRUE))) +
      geom_contour_filled(aes(z = log10(Cell2_population/Cell_population)),
                          breaks = mybreaks) +
      geom_point(aes(color = log10(Cell2_population/Cell_population)),
                 size = 3) +
      scale_x_continuous(trans = "log10") +
      scale_y_continuous(trans = "log2") +
      scale_fill_viridis_d(drop = FALSE) +
      scale_color_continuous(name = "Fitness", type = "viridis",
                             limits = c(min(mybreaks), max(mybreaks))) +
      #guides(fill = "none") +
      labs(subtitle = mylabel,
           x = "Relative Resistance",
           y = mylabeller(grep("relative[^IR]+", colnames(data[[i]]),
                               value = TRUE))[[1]]) +
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

