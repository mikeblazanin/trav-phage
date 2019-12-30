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