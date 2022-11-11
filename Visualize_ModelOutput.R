library(ggplot2)
library(ggtext)

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
data_mrg$distrib <- 
  factor(data_mrg$distrib,
         levels = c("no_paras", "local", "global", "global_gauss"))

#Make plots
for (vars_manip in unique(data_mrg$vars_manip)) {
  myrows <- which(data_mrg$vars_manip == vars_manip)
  my_data <- data_mrg[myrows, ]
  
  #Get limits for log10(Cell2/Cell) across the dift initial distributions
  mybreaks <- pretty(x = log10(my_data$Cell2_Cell1), n = 6)
  
  #axis labels labeller for the non-resistance axis
  axislabeller <-
    as_labeller(c(relativecA = "Relative Attractant Consumption",
                  relativeChi = "Relative Chemotactic Sensitivity",
                  relativecR = "Relative Growth Rate",
                  relativeY = "Relative Growth Yield"))
  
  #facet labels labeller
  facetlabeller <-
    as_labeller(c(global = "Global Parasites",
                  local = "Local Parasites",
                  global_gauss = "Global Parasites (Gaussian)",
                  no_paras = "No Parasites"))
  
  #Make file
  png(
    paste(sep = "", "./Model_plots/", vars_manip, ".png"),
    width = 6, height = 5, units = "in", res = 150)
  #Make plot
  p <-
    ggplot(my_data,
           aes_string(x = "relativeR", 
                      y = paste("relative", my_data$vars_manip_2[1], sep = ""))) +
      geom_contour_filled(aes(z = log10(Cell2_population/Cell_population)),
                          breaks = mybreaks) +
      geom_point(aes(color = log10(Cell2_population/Cell_population)),
                 size = 3) +
      facet_wrap(~distrib, labeller = facetlabeller) +
      scale_x_continuous(trans = "log10",
                         breaks = 10**c(-2, -1, 0, 1, 2),
                         labels = c("0.01", "0.1", "1", "10", "100")) +
      scale_fill_viridis_d(drop = FALSE) +
      scale_color_continuous(name = "Fitness", type = "viridis",
                             limits = c(min(mybreaks), max(mybreaks))) +
      guides(fill = "none") +
      labs(x = "Relative Resistance") +
      theme_bw() +
      NULL
  if(vars_manip == "I_vs_Y") {
    p <- p + 
      scale_y_continuous(trans = "log2", breaks = 2**c(-.5, 0, .5),
                         labels = c("-0.5", "0", "0.5"),
                         name = "log<sub>2</sub>(Relative Growth Yield)") +
      theme(axis.title.y = element_markdown())
  } else {p <- p + 
    scale_y_continuous(
      trans = "log2",
      name = axislabeller(paste("relative", my_data$vars_manip_2[1], sep = "")))
  }
  print(p)
  dev.off()
  
  #Make main-text figure
  if(vars_manip == "I_vs_Chi") {
    myrows <- which(data_mrg$vars_manip == vars_manip &
                      data_mrg$distrib %in% c("local", "global", "no_paras"))
    my_data <- data_mrg[myrows, ]
    
    #Get limits for log10(Cell2/Cell) across the dift initial distributions
    mybreaks <- pretty(x = log10(my_data$Cell2_Cell1), n = 6)
    
    #Make plot
    png("./Model_plots/I_vs_Chi_maintext.png", width = 8, height = 3, units = "in", res = 150)
    print(
      ggplot(my_data, aes_string(x = "relativeR", y = "relativeChi")) +
        geom_contour_filled(aes(z = log10(Cell2_population/Cell_population)),
                            breaks = mybreaks) +
        geom_point(aes(color = log10(Cell2_population/Cell_population)),
                   size = 3) +
        facet_grid(~distrib, labeller = facetlabeller) +
        scale_x_continuous(trans = "log10",
                           breaks = 10**c(-2, -1, 0, 1, 2),
                           labels = c("0.01", "0.1", "1", "10", "100")) +
        scale_y_continuous(trans = "log2") +
        scale_fill_viridis_d(drop = FALSE) +
        scale_color_continuous(name = "Fitness", type = "viridis",
                               limits = c(min(mybreaks), max(mybreaks))) +
        guides(fill = "none") +
        labs(x = "Relative Resistance", y = "Relative Dispersal") +
        theme_bw() +
        NULL
    )
    dev.off()
  }
}

