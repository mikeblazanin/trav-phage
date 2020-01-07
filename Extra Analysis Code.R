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
  
#Plot curves & peaks
view_peaks <- function(time_data, dens_data, well_contents, plt_point = TRUE,
                       peak_time = NULL, peak_dens = NULL, peak_well_contents,
                       numplots = 9, lwd = 1) {
  #Inputs:  time_data - vector of time values
  #         dens_data - vector of density values, matched with time_data
  #         well_contents - vector of well contents, matched with time_data
  #         peak_time - vector of peaks' time values
  #         peak_dens - vector of peaks' density values, matched with peak_time
  #         peak_well_contents - vectore of peaks' well contents, matched with peak_time
  require(ggplot2)
  for (start_group in seq(from = 1, to = length(unique(well_contents)),
                          by = numplots)) {
    my_wells <- unique(well_contents)[start_group:(start_group+numplots-1)]
    myplot <- ggplot(data = data.frame(
      "dens" = dens_data[well_contents %in% my_wells],
      "time" = time_data[well_contents %in% my_wells],
      "contents" = well_contents[well_contents %in% my_wells]),
      aes(x = time, y = dens)) +
      geom_line(lwd = lwd) +
      facet_wrap(~contents)
    if (plt_point) {
      myplot <- myplot + 
        geom_point(data = data.frame(
          "dens" = peak_dens[peak_well_contents %in% my_wells],
          "time" = peak_time[peak_well_contents %in% my_wells],
          "contents" = peak_well_contents[peak_well_contents %in% my_wells]),
          aes(x = time, y = dens),
          size = 3, pch = 13)
    }
    print(myplot)
  }
}

view_peaks(tidymerge$Timepoint, tidymerge$OD, tidymerge$Well,
           peak_time = outmerge$maxtime, peak_dens = outmerge$max, 
           peak_well_contents = outmerge$Well)


#Run PCA on gc_summarized (no summarization of replicates) ----

#Looks like some variables are skewed:
# max_percap_gr_rate_dens
# pseudo_K
#So normalize by log10
gc_summarized$max_percap_gr_rate_dens_log10 <- 
  log10(gc_summarized$max_percap_gr_rate_dens)
gc_summarized$pseudo_K_log10 <- log10(gc_summarized$pseudo_K)

#Then view for univariate normality of each variable
if (F) {
  for (var in c("first_min", "first_min_time", 
                "max_percap_gr_rate", "max_percap_gr_rate_time", 
                "max_percap_gr_rate_dens_log10", "max_percap_gr_rate_timesincemin",
                "pseudo_K_log10", "pseudo_K_time", 
                "pseudo_K_timesincemin", "pseudo_K_timesince_maxpercap")) {
    hist(as.numeric(gc_summarized[, var]), main = var)
    qqnorm(as.numeric(gc_summarized[, var]), main = var)
  }
}

#Now check for multivariate normality:

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
CSQPlot(gc_summarized[, c("first_min", 
                          #"first_min_time", 
                          "max_percap_gr_rate", 
                          #"max_percap_gr_rate_time", 
                          "max_percap_gr_rate_dens_log10", 
                          "max_percap_gr_rate_timesincemin",
                          "pseudo_K_log10", 
                          #"pseudo_K_time", 
                          #"pseudo_K_timesincemin", 
                          "pseudo_K_timesince_maxpercap")])
#Our data is highly non-normal

#Get principal components
gc_prin_comp <- princomp(gc_summarized[, c("first_min", 
                                           #"first_min_time", 
                                           "max_percap_gr_rate", 
                                           #"max_percap_gr_rate_time", 
                                           "max_percap_gr_rate_dens_log10", 
                                           "max_percap_gr_rate_timesincemin",
                                           "pseudo_K_log10", 
                                           #"pseudo_K_time", 
                                           "pseudo_K_timesincemin"
                                           #"pseudo_K_timesince_maxpercap"
)],
cor = T,
scores = T)

#Print summary
print(summary(gc_prin_comp), digits = 2)
#Looks like PC1 & 2 are ~72% of the variance, with each contributing
# nearly equally

#Not obvious from screeplot what the cutoff, if any, would be
screeplot(gc_prin_comp, type = "lines", main = "Movie PCA Scree Plot")

#Check loadings of Principal components
print(gc_prin_comp$loadings, digits = 2, cutoff = 0)

#PC1 - the 3 densities (first min, dens as max percap, pseudo K)
#PC2 - everything except first min
#PC3 - max percap & pseudo K time since max percap

gc_summarized <- cbind(gc_summarized, gc_prin_comp$scores)

ggplot(data = gc_summarized,
       aes(x = Comp.1, y = Comp.2, color = Treat, shape = Proj)) +
  geom_point() +
  facet_wrap(~Media)

#Definitely not skewed:
# 125 first min Orig
# 125 pseudo K Orig
# 125 first min Rich
#
#Not skewed but outliers:
#7x pseudo K Orig
#7x pseudo K timesincemin Orig
#7x pseudo K Rich
#125 max percap dens Rich
# 
#Maybe skewed:
#7x max percap rate Orig
#7x max percap dens Orig
#7x max percap timesincemin Orig
#7x first min Rich
#7x max percap rate Rich
#7x max percap timesincemin Rich
#7x pseudo K timesincemin Rich
#125 max percap rate Orig
#125 max percap dens Orig
#125 max percap timesincemin Orig
#125 pseudo K timesincemin Orig
#125 max percap rate Rich
#125 max percap timesincemin Rich
#125 pseudo K Rich
#
#Definitely skewed:
#7x first min Orig
#7x max percap dens Rich
#125 pseudo K timesincemin Rich

##See GC Data Normality.rtf

#Split data frame into two projects
gc_noreps_7x <- gc_sum_isols_wide[gc_sum_isols_wide$Proj == "7x", ]
gc_noreps_125 <- gc_sum_isols_wide[gc_sum_isols_wide$Proj == "125", ]

#Take log10 of everything except 125 first min
for (var in c("first_min_avg_Orig", 
              #"first_min_time_avg_Orig",
              "max_percap_gr_rate_avg_Orig",
              #"max_percap_gr_time_avg_Orig",
              "max_percap_gr_dens_avg_Orig",
              "max_percap_gr_timesincemin_avg_Orig",
              "pseudo_K_avg_Orig",
              #"pseudo_K_time_avg_Orig",
              "pseudo_K_timesincemin_avg_Orig",
              #"pseudo_K_timesince_maxpercap_avg_Orig",
              "first_min_avg_Rich",
              #"first_min_time_avg_Rich",
              "max_percap_gr_rate_avg_Rich",
              #"max_percap_gr_time_avg_Rich",
              "max_percap_gr_dens_avg_Rich",
              "max_percap_gr_timesincemin_avg_Rich",
              "pseudo_K_avg_Rich",
              #"pseudo_K_time_avg_Rich",
              "pseudo_K_timesincemin_avg_Rich"
              #"pseudo_K_timesince_maxpercap_avg_Rich"
)) {
  newname <- paste(var, "_log10", sep = "")
  gc_noreps_7x[, newname] <- log10(gc_noreps_7x[, var])
  if (!var %in% c("first_min_avg_Orig", "first_min_avg_Rich")) {
    gc_noreps_125[, newname] <- log10(gc_noreps_125[, var])
  }
}

#Re-check normality
#7x:
for (var_root in c("first_min_avg_", 
                   #"first_min_time_avg_",
                   "max_percap_gr_rate_avg_",
                   #"max_percap_gr_time_avg_",
                   "max_percap_gr_dens_avg_",
                   "max_percap_gr_timesincemin_avg_",
                   "pseudo_K_avg_",
                   #"pseudo_K_time_avg_",
                   "pseudo_K_timesincemin_avg_",
                   #"pseudo_K_timesince_maxpercap_avg_",
                   "first_min_avg_",
                   #"first_min_time_avg_",
                   "max_percap_gr_rate_avg_",
                   #"max_percap_gr_time_avg_",
                   "max_percap_gr_dens_avg_",
                   "max_percap_gr_timesincemin_avg_",
                   "pseudo_K_avg_",
                   #"pseudo_K_time_avg_",
                   "pseudo_K_timesincemin_avg_"
                   #"pseudo_K_timesince_maxpercap_avg_"
)) {
  for (media in c("Orig", "Rich")) {
    var <- paste(var_root, media, "_log10", sep = "")
    qqnorm(as.numeric(gc_noreps_7x[, var]), 
           main = paste("7x", var))
    qqline(as.numeric(gc_noreps_7x[, var]))
    if (var_root != "first_min_avg_") {
      qqnorm(as.numeric(gc_noreps_125[, var]), 
             main = paste("125", var))
      qqline(as.numeric(gc_noreps_125[, var]))
    }
  }
}

#After looking at the outcomes of log10 normalization,
# first min did not improve, use original values
# All three max percap improved in some cases, use log10 values
# Pseudo K did not change, use original values

#Check for multivariate normality



#Looks like our data is pretty non-normal, although the majority of
# it falls within or near the 95% confidence intervals                    

#Get principal components
gc_princomp_7x <- princomp(gc_noreps_7x[, c("first_min_avg_Orig", 
                                            "first_min_avg_Rich",
                                            "max_percap_gr_rate_avg_Orig_log10",
                                            "max_percap_gr_rate_avg_Rich_log10",
                                            "max_percap_gr_dens_avg_Orig_log10",
                                            "max_percap_gr_dens_avg_Rich_log10",
                                            "max_percap_gr_timesincemin_avg_Orig_log10",
                                            "max_percap_gr_timesincemin_avg_Rich_log10",
                                            "pseudo_K_avg_Orig",
                                            "pseudo_K_avg_Rich",
                                            "pseudo_K_timesincemin_avg_Orig",
                                            "pseudo_K_timesincemin_avg_Rich")],
                           cor = T,
                           scores = T)

gc_princomp_125 <- princomp(gc_noreps_125[, c("first_min_avg_Orig", 
                                              "first_min_avg_Rich",
                                              "max_percap_gr_rate_avg_Orig_log10",
                                              "max_percap_gr_rate_avg_Rich_log10",
                                              "max_percap_gr_dens_avg_Orig_log10",
                                              "max_percap_gr_dens_avg_Rich_log10",
                                              "max_percap_gr_timesincemin_avg_Orig_log10",
                                              "max_percap_gr_timesincemin_avg_Rich_log10",
                                              "pseudo_K_avg_Orig",
                                              "pseudo_K_avg_Rich",
                                              "pseudo_K_timesincemin_avg_Orig",
                                              "pseudo_K_timesincemin_avg_Rich")],
                            cor = T,
                            scores = T)

#Print summary
print(summary(gc_princomp_7x), digits = 2)
#PC1 - 36% of variance, PC2 22%, PC3 16%
print(summary(gc_princomp_125), digits = 2)
#PC1 43%, PC2 18% PC3 14%

#Make screeplots
screeplot(gc_princomp_7x, type = "lines", main = "7x PCA Scree Plot")
#Keep 2 or 4?
screeplot(gc_princomp_125, type = "lines", main = "7x PCA Scree Plot")
#Keep 2

#Check loadings of Principal components for 7x
print(gc_princomp_7x$loadings, digits = 2, cutoff = 0)
#PC1 - first min, percap dens, percap time vs percap rate
#PC2 - pseudo K time & dens

#Check loadings of Principal components for 7x
print(gc_princomp_125$loadings, digits = 2, cutoff = 0)
#PC1 - max percap rate vs max percap time & dens, pseudo K time, first min
#PC2 - first min, pseudo K dens vs max percap time & pseudo K time

#Add scores to dataframes
gc_noreps_7x <- cbind(gc_noreps_7x, gc_princomp_7x$scores)
gc_noreps_125 <- cbind(gc_noreps_125, gc_princomp_125$scores)

#Make plots of all isols
ggplot(data = gc_noreps_7x,
       aes(x = Comp.1, y = Comp.2, color = Treat)) +
  geom_point()

ggplot(data = gc_noreps_125,
       aes(x = Comp.1, y = Comp.2, color = Treat)) +
  geom_point()

#Summarize each pop
gc_noreps_7x <- group_by(gc_noreps_7x,
                         "Proj",
                         "Pop",
                         "Treat")
gc_noisols_7x <- summarise(gc_noreps_7x,
                           
)


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

#Old resistance plots ----
# #For local viewing
# ggplot(resis_data[resis_data$Treat != "A", ], aes(x = Treat, y = 1-EOP)) +
#   scale_x_discrete(name = "Treatment") +
#   ylab("Resistance to Phage") +
#   facet_grid(Proj~Pop, labeller = labeller(Proj = my_facet_labels)) +
#   geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.01, 
#                dotsize = 3) +
#   # geom_jitter(width = 0.2, height = 0, size = 2) +
#   theme_bw() + ggtitle("Population") +
#   theme(plot.title = element_text(size = 12, hjust = 0.5), 
#         axis.title = element_text(size = 12),
#         axis.text.x = element_text(color = "black", size = 12)) +
#   geom_hline(yintercept = 0, linetype = "dotted", size = 1.5)
# 
# #For poster
# png("resis_isols.png", width = 14, height = 9, units = "in", res = 300)
# ggplot(resis_data[resis_data$Treat != "A", ], aes(x = Treat, y = 1-EOP)) +
#   scale_x_discrete(name = "Treatment") +
#   ylab("Resistance to Phage") +
#   facet_grid(Proj~Pop, labeller = labeller(Proj = my_facet_labels)) +
#   geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.02, 
#                dotsize = 2) +
#   theme_bw() + ggtitle("Population") +
#   theme(plot.title = element_text(size = 24, hjust = 0.5), 
#         axis.title = element_text(size = 24),
#         axis.text = element_text(color = "black", size = 20),
#         strip.text = element_text(size = 24)) +
#   geom_hline(yintercept = 0, linetype = "dotted", size = 1.5)
# dev.off()
# 
# #resistance vs growth
# gc_resis_data <- merge(resis_data, gc_mppti)
# gc_resis_mppt <- group_by(gc_resis_data, Media, Proj, Pop, Treat)
# gc_resis_mppt <- summarize(gc_resis_mppt, avg_eop = mean(EOP),
#                            avg_gr = mean(gr_max_avg))
# 
# #Make plot of all isolates
# ggplot(gc_resis_data, aes(x = 1-EOP, y = gr_max_avg)) +
#   geom_point() + 
#   facet_grid(Media~., labeller = labeller(Media = my_facet_labels)) +
#   geom_smooth(method = "lm") +
#   labs(x = "Resistance", y = "Per Capita Growth Rate (/hour)")
# 
# summary(lm(gr_max_avg~Media*EOP, data = gc_resis_data))
# 
# #Make plot of all pops
# #For local viewing
# ggplot(gc_resis_mppt, aes(x = 1-avg_eop, y = avg_gr)) +
#   geom_point(size = 2, aes(pch = Treat)) + 
#   scale_shape_manual(values = c(3, 15, 16, 17), name = "Treatment",
#                      breaks = c("A", "C", "G", "L"),
#                      labels = c("WT", "Control", "Global", "Local")) +
#   facet_grid(Media~Proj, labeller = labeller(Media = my_facet_labels,
#                                              Proj = my_facet_labels)) +
#   geom_smooth(method = "lm") +
#   labs(x = "Resistance", y = "Average Per Capita Growth Rate (/hour)") +
#   theme_bw()
# 
# #For poster
# png("resis_gc_tradeoff.png", width = 10, height = 7, units = "in",
#     res = 300)
# ggplot(gc_resis_mppt, aes(x = 1-avg_eop, y = avg_gr)) +
#   geom_point(size = 4, aes(pch = Treat)) + 
#   scale_shape_manual(values = c(3, 15, 16, 17), name = "Treatment",
#                      breaks = c("A", "C", "G", "L"),
#                      labels = c("WT", "Control", "Global", "Local")) +
#   facet_grid(Proj~Media, labeller = labeller(Media = my_facet_labels,
#                                              Proj = my_facet_labels)) +
#   geom_smooth(method = "lm", se = F) +
#   labs(x = "Resistance", y = "Average Maximum Per Capita Growth Rate (/hr)") +
#   theme_bw() +
#   theme(axis.text = element_text(size = 16),
#         strip.text = element_text(size = 20),
#         axis.title = element_text(size = 20),
#         legend.title = element_text(size = 20),
#         legend.text = element_text(size = 16))
# dev.off()
# 
# summary(lm(avg_gr~avg_eop*Media, data = gc_resis_mppt))
# 
