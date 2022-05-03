data <- lapply(X = paste("./Model_data/", list.files("./Model_data/"), sep = ""),
               FUN = read.csv)
names(data) = list.files("./Model_data/")
data <- lapply(X = data,
               FUN = function(x) {
                 cbind(x, data.frame("relativeR" = 1/x[, "relativeI"]))})

library(ggplot2)

for (i in 1:length(data)) {
  #Get initial distribution information
  if(strsplit(names(data)[i], split = "_")[[1]][4] == "uniform.csv") {
    mylabel = "Global Parasite Distribution"
  } else {mylabel = "Local Parasite Distribution"}
  
  #Get limits for log10(Cell2/Cell) across the two dift initial distributions
  these_vars_idx <- 
    grep(paste(collapse = "_", strsplit(names(data)[i], split = "_")[[1]][1:3]),
         names(data))
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
    paste(sep = "",
          "./Model_plots/",
          paste(collapse = "_", strsplit(names(data)[i], split = "_")[[1]][1:3]), 
          "_", mylabel, ".png"),
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
      scale_y_continuous(trans = "log10") +
      scale_fill_viridis_d(drop = FALSE) +
      scale_color_continuous(name = "Fitness", type = "viridis",
                             limits = c(min(mybreaks), max(mybreaks))) +
      guides(fill = FALSE) +
      labs(subtitle = mylabel,
           x = "Relative Resistance",
           y = mylabeller(grep("relative[^IR]+", colnames(data[[i]]),
                               value = TRUE))[[1]]) +
      theme_bw() +
      NULL
  )
  dev.off()
}