library(deSolve)
library(ggplot2)

#Define function to calculate derivatives
derivs <- function(time, y, parms, nx, r_mid, r_end, dx, disp_dx) {
  with(as.list(parms), {
    N <- y[1:nx]
    P <- y[(nx+1):(2*nx)]
    R <- y[(2*nx+1):(3*nx)]
    A <- y[(3*nx+1):(4*nx)]
    
    #Adding this incredibly slows down the integration
    # N[N < 0] <- 0
    # P[P < 0] <- 0
    # R[R < 0] <- 0
    # A[A < 0] <- 0
    
    #Calculate diffusion & migration (zero gradient at boundaries)
    flux_N <- -diff(r_end * -D_N * diff(c(N[1], N, N[nx]))/disp_dx)/r_mid/dx -
      chi*N*diff(r_end * diff(c(A[1], A, A[nx]))/disp_dx)/r_mid/dx
    flux_P <- -diff(r_end * -D_P * diff(c(P[1], P, P[nx]))/disp_dx)/r_mid/dx
    flux_R <- -diff(r_end * -D_R * diff(c(R[1], R, R[nx]))/disp_dx)/r_mid/dx
    flux_A <- -diff(r_end * -D_A * diff(c(A[1], A, A[nx]))/disp_dx)/r_mid/dx
    
    #for troubleshooting:
    if (F) {
      par(mfrow = c(2, 1))
      plot(1:length(N), N)
      plot(1:length(N), 
           -diff(r_end * -D_N * diff(c(N[1], N, N[nx]))/disp_dx)/r_mid/dx,
           main = "diffusion")
      lines(1:length(N),
            -chi*N*diff(r_end * diff(c(A[1], A, A[nx]))/disp_dx)/r_mid/dx)
      par(mfrow = c(1, 1))
      plot(-diff(r_end * -D_N * diff(c(N[1], N, N[nx]))/disp_dx)/r_mid/dx,
           -chi*N*diff(r_end * diff(c(A[1], A, A[nx]))/disp_dx)/r_mid/dx,
           xlab = "Diffusion", ylab = "Migration")
    }
    
    #Calculate growth and death
    grow_N <- c_R*yield*R/(R+k_R)*N - i*N*P
    grow_P <- i*(b-1)*N*P
    
    #Calculate consumption
    cons_R <- c_R*R/(R+k_R)*N
    cons_A <- c_A*A/(A+k_A)*N
    
    #Change = diffusion&migration + growth&death + consumption
    dN <- flux_N + grow_N
    dP <- flux_P + grow_P
    dR <- flux_R - cons_R
    dA <- flux_A - cons_A
    
    return(list(c(dN, dP, dR, dA)))
  })
}

#Simple model (for future reference)
if (F) {
  #Calculate physical dimensions
  radius <- 45000                      #um
  nx <- 100                            #concentric circles
  dx <- radius/nx                      #thickness of each circle
  r_mid <- seq(dx/2,by = dx,len = nx)  #distance from center to ea mid-layer
  r_end <- seq(0,by = dx,len = nx+1)   #distance from center to end of each layer
  disp_dx <- dx                        #dispersion distances
  
  inoc_r_n <- 1425
  inoc_r_p <- 5*1425
  
  #Define parameters
  parms <- c(D_N = 50, D_P = 0, D_R = 800, D_A = 800,
             chi = 300,
             c_A = 4*10**-13, c_R = 2*10**-11,
             yield = 10**7,
             k_R = .05, k_A = .001,
             i = 10**-12, b = 50)
  
  #Define init conditions (NPRA)
  y_init = matrix(c(
    rep(250000/round(inoc_r_n/dx), round(inoc_r_n/dx)), rep(0, nx-round(inoc_r_n/dx)),
    rep(25000/round(inoc_r_p/dx), round(inoc_r_p/dx)), rep(0, nx-round(inoc_r_p/dx)),
    rep(2000/nx, nx),
    rep(25/nx, nx)),
    ncol = 4)
  
  #Define times (in seconds)
  times = seq(from = 0, to = 24*60*60, by = 1*60)
  
  #Run model
  yout <- ode.1D(y = y_init, times = times,
                 func = derivs, parms = parms,
                 nspec = 4, names = c("N", "P", "R", "A"),
                 nx = nx, dx=dx, r_mid = r_mid, r_end = r_end, 
                 disp_dx = disp_dx, maxsteps = 25000)
  yout_df <- as.data.frame(yout)
  yout_df <- yout_df[as.numeric(yout_df$time) %% (15*60) == 0, ]
  
  #Reorganize
  colnames(yout_df)[2:ncol(yout_df)] <- 
    c(paste(rep("N", nx), 0:(nx-1), sep = "_"), 
      paste(rep("P", nx), 0:(nx-1), sep = "_"),
      paste(rep("R", nx), 0:(nx-1), sep = "_"),
      paste(rep("A", nx), 0:(nx-1), sep = "_"))
  
  yout_lng <- as.data.frame(tidyr::pivot_longer(yout_df, cols = -time,
                                                names_to = c("pop", "x"),
                                                names_sep = "_",
                                                values_to = "density"))
  
  if (F) {
    ggplot(data = yout_lng[yout_lng$pop == "R", ], 
           aes(x = as.numeric(time), 
               y = as.numeric(x), 
               color = as.numeric(density))) +
      #color = log10(as.numeric(ifelse(density>0, density, 0)+10)))) +
      #geom_contour_filled() +
      geom_point() +
      #geom_line(aes(group = x)) +
      #facet_grid(~pop) +
      #scale_color_continuous(name = "dens") +
      NULL
  }
  
  #Make contour plots
  if (F) {
    n <- ggplot(data = yout_lng[yout_lng$pop == "N", ], 
                aes(x = as.numeric(time)/3600, 
                    y = as.numeric(x), 
                    z = log10(density+1))) +
      geom_contour_filled() +
      facet_grid(~pop)
    p <- ggplot(data = yout_lng[yout_lng$pop == "P", ], 
                aes(x = as.numeric(time)/3600, 
                    y = as.numeric(x), 
                    z = log10(density+1))) +
      geom_contour_filled() +
      facet_grid(~pop)
    r <- ggplot(data = yout_lng[yout_lng$pop == "R", ], 
                aes(x = as.numeric(time)/3600, 
                    y = as.numeric(x), 
                    z = log10(density+1))) +
      geom_contour_filled() +
      facet_grid(~pop) +
      NULL
    a <- ggplot(data = yout_lng[yout_lng$pop == "A", ], 
                aes(x = as.numeric(time)/3600, 
                    y = as.numeric(x), 
                    z = log10(density+1))) +
      geom_contour_filled() +
      facet_grid(~pop) +
      NULL
    
    tiff("model_plot.tiff", width = 10, height = 10, units = "in", res = 300)
    cowplot::plot_grid(n, p, r, a, nrow = 2)
    dev.off()
  }
  
  my_times <- 60*60*c(0, 1, 3, 6, 12, 18, 24)
  ggplot(data = yout_lng[yout_lng$time %in% my_times &
                           yout_lng$pop %in% c("N", "P", "R", "A"), ],
         aes(x = as.numeric(x)*dx/10000, y = as.numeric(density)+1)) +
    facet_wrap(pop ~ ., scales = "free") +
    geom_line(lwd = 1.25, aes(color = as.factor(time/3600))) +
    scale_y_continuous(trans = "log10") +
    geom_hline(yintercept = 1, lty = 2) +
    xlim(NA, 2) +
    scale_color_manual(values = 
                         scales::seq_gradient_pal("red", "blue")(
                           seq(0,1,length.out = length(my_times)))) +
    NULL
  
  if (F) {
    ggplot(data = yout_lng[yout_lng$x %in% c(0:3) &
                             yout_lng$pop %in% c("R", "N", "A"), ],
           aes(x = as.numeric(time), y = as.numeric(density))) +
      geom_line() +
      facet_grid(pop~x, scales = "free")
  }
}

#Define function to run simulations across grid of values
run_sims <- function(inoc_r_n = 1425, inoc_r_p = 1425,
                     init_N_dens = 250000/1425, 
                     init_P_dens = 25000/1425,
                     chi = 300, c_A = 4*10**-13, c_R = 2*10**-11, 
                     yield = 10**7, i = 10**-12, nx = 100, 
                     combinatorial = TRUE, print_info = TRUE) {
  #nx is the number of concentric circles the plate is broken up into
  #inoc_r is the radius, in um, n and p are respectively inoculated into
  #init_N and P are the total inoculated bacterial & phage populations
  # (which will be distributed over inoc_r_n area)
  
  if(any(length(nx) > 1)) {stop("Varying nx is not supported")}
  
  #Save parameter value combinations into dataframe
  if (combinatorial == TRUE) {
  param_combos <- expand.grid(list("inoc_r_n" = inoc_r_n, "inoc_r_p" = inoc_r_p,
                     "init_N_dens" = init_N_dens, "init_P_dens" = init_P_dens,
                     "chi" = chi, "c_A" = c_A, "c_R" = c_R, 
                     "yield" = yield, "i" = i),
                     stringsAsFactors = FALSE)
  num_sims <- nrow(param_combos)
  } else {
    num_sims <- max(sapply(X = list(inoc_r_n, inoc_r_p, init_N_dens, init_P_dens,
                                    chi, c_A, c_R, yield, i), FUN = length))
    
    #Check for parameter lengths being non-divisible with the
    # number of simulations inferred from the longest parameter length
    if (!all(num_sims %% sapply(X = list(inoc_r_n, inoc_r_p, 
                                         init_N_dens, init_P_dens,
                                         chi, c_A, c_R, yield, i), 
                                FUN = length) == 0)) {
      warning("Combinatorial=TRUE but longest param vals length is not a multiple of all other param vals lengths")
    }
    
    #Save parameters into dataframe, replicating shorter parameter
    # vectors as needed to reach # of simulations
    param_combos <- data.frame("inoc_r_n" = rep_len(inoc_r_n, num_sims), 
                               "inoc_r_p" = rep_len(inoc_r_p, num_sims),
                               "init_N_dens" = rep_len(init_N_dens, num_sims), 
                               "init_P_dens" = rep_len(init_P_dens, num_sims),
                               "chi" = rep_len(chi, num_sims), 
                               "c_A" = rep_len(c_A, num_sims), 
                               "c_R" = rep_len(c_R, num_sims), 
                               "yield" = rep_len(yield, num_sims), 
                               "i" = rep_len(i, num_sims),
                               stringsAsFactors = FALSE)
  }
  
  #Define physical parameters of plate (true for all runs)
  radius <- 45000                      #um
  dx <- radius/nx                      #thickness of each circle
  r_mid <- seq(dx/2,by = dx,len = nx)  #distance from center to ea mid-layer
  r_end <- seq(0,by = dx,len = nx+1)   #distance from center to end of each layer
  disp_dx <- dx                        #dispersion distances
  
  #Define times (in seconds) (true for all runs)
  times = seq(from = 0, to = 24*60*60, by = 1*60)
  times_keep <- seq(from = 0, to = 24*60*60, by = 15*60)
  
  #Create output dataframe (yout)
  bigout <- 
    cbind(
      data.frame(
        "uniq_run" = rep(1:nrow(param_combos), each = length(times_keep)),
        "inoc_r_n" = rep(param_combos$inoc_r_n, each = length(times_keep)), 
        "inoc_r_p" = rep(param_combos$inoc_r_p, each = length(times_keep)),
        "init_N_dens" = rep(param_combos$init_N_dens, each = length(times_keep)), 
        "init_P_dens" = rep(param_combos$init_P_dens, each = length(times_keep)),
        "chi" = rep(param_combos$chi, each = length(times_keep)), 
        "c_A" = rep(param_combos$c_A, each = length(times_keep)), 
        "c_R" = rep(param_combos$c_R, each = length(times_keep)), 
        "yield" = rep(param_combos$yield, each = length(times_keep)), 
        "i" = rep(param_combos$i, each = length(times_keep)),
        stringsAsFactors = FALSE),
      as.data.frame(matrix(NA, nrow = length(times_keep)*nrow(param_combos), 
                           ncol = 1+4*nx,
                           dimnames = list(NULL, c("time", 1:(4*nx))))))
  
  #Print number of simulations that will be run
  if(print_info) {
    print(paste(num_sims, "simulations will be run"))
    
    #Save sequence of 10% cutoff points for later reference
    progress_seq <- round(seq(from = 0, to = num_sims, by = num_sims/10))
  }
  
  #Run simulations
  for (myrun in 1:nrow(param_combos)) {
    #Define parameters vector
    parms <- c(D_N = 50, D_P = 0, D_R = 800, D_A = 800,
             chi = param_combos$chi[myrun],
             c_A = param_combos$c_A[myrun], 
             c_R = param_combos$c_R[myrun],
             yield = param_combos$yield[myrun],
             k_R = .05, k_A = .001,
             i = param_combos$i[myrun], b = 50)
    
    #Define initial conditions
    y_init = matrix(c(
      #N
      rep(param_combos$init_N_dens[myrun], 
          round(param_combos$inoc_r_n[myrun]/dx)),
      rep(0, nx - round(param_combos$inoc_r_n[myrun]/dx)),
      #P
      rep(param_combos$init_P_dens[myrun], 
          round(param_combos$inoc_r_p[myrun]/dx)),
      rep(0, nx - round(param_combos$inoc_r_p[myrun]/dx)),
      #R
      rep(2000/nx, nx),
      #A
      rep(25/nx, nx)),
      ncol = 4)
    
    #Run model
    yout <- as.data.frame(ode.1D(y = y_init, times = times,
                                 func = derivs, parms = parms,
                                 nspec = 4, names = c("N", "P", "R", "A"),
                                 nx = nx, dx=dx, r_mid = r_mid, r_end = r_end, 
                                 disp_dx = disp_dx, maxsteps = 25000))
    yout <- yout[yout$time %in% times_keep, ]
    
    #Save results
    bigout[((myrun-1)*length(times_keep)+1):(myrun*length(times_keep)), 
           11:ncol(bigout)] <- yout
    
    #Print progress update
    if (print_info & myrun %in% progress_seq) {
      print(paste((which(progress_seq == myrun)-1)*10,
                  "% completed", sep = ""))
    }
  }
  
  #Rename columns
  #Reorganize
  colnames(bigout)[12:ncol(bigout)] <- 
    c(paste(rep("N", nx), 0:(nx-1), sep = "_"), 
      paste(rep("P", nx), 0:(nx-1), sep = "_"),
      paste(rep("R", nx), 0:(nx-1), sep = "_"),
      paste(rep("A", nx), 0:(nx-1), sep = "_"))
  
  return(list(bigout, param_combos))
}

run1_nx <- 100
if (F) {
  run1 <- run_sims(inoc_r_p = c(0, 1425, 45000),
                   chi = c(300, 400, 500, 600),
                   c_R = c(2, 2.5, 3, 3.5, 4)*10**-11,
                   c_A = c(4, 6, 8)*10**-13,
                   i = 10**c(-12, -14, -16, -18), 
                   nx = run1_nx)
  write.csv(run1[[1]], "./run1_1.csv", row.names = FALSE, quote = FALSE)
  write.csv(run1[[2]], "./run1_2.csv", row.names = FALSE, quote = FALSE)
} else {run1 <- list(read.csv("./run1_1.csv", header = TRUE),
                     read.csv("./run1_2.csv", header = TRUE))}

#Pivot to tidy
run1_lng <- as.data.frame(
  tidyr::pivot_longer(run1[[1]],
                      cols = -c("uniq_run", "inoc_r_n", "inoc_r_p", 
                                "init_N_dens", "init_P_dens", "chi", "c_A", 
                                "c_R", "yield", "i", "time"),
                      names_to = c("pop", "x"),
                      names_sep = "_",
                      values_to = "density"))

#Summarize
run1_sum <- 
  dplyr::summarise(
    dplyr::group_by(run1_lng[run1_lng$time == max(run1_lng$time) &
                               run1_lng$pop == "N", ],
                    uniq_run, inoc_r_n, inoc_r_p, init_N_dens, init_P_dens,
                    chi, c_A, c_R, yield, i),
    tot_dens = sum(density),
    dx = 45000/length(unique(x)),
    percentile_95 = dx*max(which(cumsum(density) < 0.95*tot_dens)),
    percentile_90 = dx*max(which(cumsum(density) < 0.90*tot_dens)),
    first_100cfuml = dx*max(which(density >= 10**2)),
    first_1000cfuml = dx*max(which(density >= 10**3)),
    first_10000cfuml = dx*max(which(density >= 10**4)),
    first_1000cfuml_dens = density[max(which(density >= 10**3))])

#Make profile plots for each run
dir.create("./run1_profile_plots", showWarnings = FALSE)
my_times <- 60*60*c(0, 1, 2, 3, 6, 12, 18, 24)
for (uniq_run in unique(run1_lng$uniq_run)) {
  tiff(paste("./run1_profile_plots/", uniq_run, ".tiff", sep = ""),
       width = 6, height = 4, units = "in", res = 200)
  print(ggplot(data = run1_lng[run1_lng$time %in% my_times &
                                 run1_lng$uniq_run == uniq_run, ],
               aes(x = as.numeric(x)*(45000/run1_nx)/10000, y = as.numeric(density)+1)) +
          facet_wrap(pop ~ ., scales = "free",
                     labeller = labeller(pop = c("A" = "Attractant", "N" = "Bacteria",
                                                 "P" = "Phage", "R" = "Resources"))) +
          geom_line(lwd = 1.25, aes(color = as.factor(time/3600))) +
          geom_point(data = cbind(run1_sum[run1_sum$uniq_run == uniq_run, ],
                                  "pop" = "N"),
                     aes(x = first_1000cfuml/10000,
                         y = first_1000cfuml_dens+1)) +
          scale_y_continuous(trans = "log10") +
          geom_hline(yintercept = 1, lty = 2) +
          scale_color_manual(name = "Time (hrs)",
                             values = 
                               scales::seq_gradient_pal("red", "blue")(
                                 seq(0,1,length.out = length(my_times)))) +
          labs(x = "Position (cm)", y = "Density + 1 (cfu, pfu, or \U003BCmol)") +
          NULL)
  dev.off()
}

#In run 1, when phage are global and infec rate is high (e-12), bact form
# traveling peaks. Otherwise they form expanding fronts.
# (uniq runs 3, 6, 9 after 24 hrs, 12, 15, 18, 21, 24, 27 at multiple high hrs)
                   
tiff("./Modeling_plots/run1_contour1.tiff", width = 10, height = 10,
     units = "in", res = 300)
ggplot(data = run1_sum, aes(x = chi, y = i, z = percentile_95)) +
  geom_contour_filled() +
  facet_grid(inoc_r_p*c_A ~ c_R,
             labeller = labeller(inoc_r_p = c("0" = "Control", 
                                       "1425" = "Local", "45000" = "Global"))) +
  scale_y_continuous(trans = "log10") +
  labs(x = "chemotaxis sensitivity", y = "infection rate",
       subtitle = "Resource Consumption Rate", title = "95th Percentile")
dev.off()

tiff("./Modeling_plots/run1_contour2.tiff", width = 10, height = 10,
     units = "in", res = 300)
ggplot(data = run1_sum, aes(x = chi, y = i, z = percentile_90)) +
  geom_contour_filled() +
  facet_grid(inoc_r_p*c_A ~ c_R,
             labeller = labeller(inoc_r_p = c("0" = "Control", 
                                              "1425" = "Local", "45000" = "Global"))) +
  scale_y_continuous(trans = "log10") +
  labs(x = "chemotaxis sensitivity", y = "infection rate",
       subtitle = "Resource Consumption Rate", title = "90th Percentile")
dev.off()

tiff("./Modeling_plots/run1_contour3.tiff", width = 10, height = 10,
     units = "in", res = 300)
ggplot(data = run1_sum, aes(x = chi, y = i, z = first_100cfuml)) +
  geom_contour_filled() +
  facet_grid(inoc_r_p*c_A ~ c_R,
             labeller = labeller(inoc_r_p = c("0" = "Control", 
                                              "1425" = "Local", "45000" = "Global"))) +
  scale_y_continuous(trans = "log10") +
  labs(x = "chemotaxis sensitivity", y = "infection rate",
       subtitle = "Resource Consumption Rate", 
       title = "Threshold 100 cfu/mL Percentile")
dev.off()

tiff("./Modeling_plots/run1_contour4.tiff", width = 10, height = 10,
     units = "in", res = 300)
ggplot(data = run1_sum, aes(x = chi, y = i, z = first_1000cfuml)) +
  geom_contour_filled() +
  facet_grid(inoc_r_p*c_A ~ c_R,
             labeller = labeller(inoc_r_p = c("0" = "Control", 
                                              "1425" = "Local", "45000" = "Global"))) +
  scale_y_continuous(trans = "log10") +
  labs(x = "chemotaxis sensitivity", y = "infection rate",
       subtitle = "Resource Consumption Rate", 
       title = "Threshold 1000 cfu/mL Percentile")
dev.off()

tiff("./Modeling_plots/run1_contour5.tiff", width = 10, height = 10,
     units = "in", res = 300)
ggplot(data = run1_sum, aes(x = chi, y = i, z = first_10000cfuml)) +
  geom_contour_filled() +
  facet_grid(inoc_r_p*c_A ~ c_R,
             labeller = labeller(inoc_r_p = c("0" = "Control", 
                                              "1425" = "Local", "45000" = "Global"))) +
  scale_y_continuous(trans = "log10") +
  labs(x = "chemotaxis sensitivity", y = "infection rate",
       subtitle = "Resource Consumption Rate", 
       title = "Threshold 10000 cfu/mL Percentile")
dev.off()

for (my_c_R in unique(run1_sum$c_R)) {
  tiff(paste("./Modeling_plots/run1_contour4_c_r=", my_c_R, ".tiff", sep = ""),
       width = 10, height = 10, units = "in", res = 300)
  print(ggplot(data = run1_sum[run1_sum$c_R == my_c_R, ], 
         aes(x = chi, y = i, z = first_1000cfuml)) +
    geom_contour_filled() +
    facet_grid(inoc_r_p ~ c_A,
               labeller = labeller(inoc_r_p = c("0" = "Control", 
                                                "1425" = "Local", "45000" = "Global"))) +
    scale_y_continuous(trans = "log10") +
    labs(x = "chemotaxis sensitivity", y = "infection rate",
         subtitle = "Resource Consumption Rate", 
         title = "Threshold 1000 cfu/mL Percentile"))
  dev.off()
}

#The changes included in run1 (specifically those for c_A and chi)
# aren't enough to actually affect the migration rate
# so changes in migration rate aren't proportional to changes in c_A and chi

run2_nx <- 100
if (F) {
  run2 <- run_sims(inoc_r_p = 0,
                   chi = 300*10**c(0, 2, 4, 6, 8),
                   c_R = 2*10**-11,
                   c_A = 4*10**c(-13, -11, -9, -7),
                   i = 10**-12, 
                   nx = run2_nx)
  write.csv(run2[[1]], "./run2_1.csv", row.names = FALSE, quote = FALSE)
  write.csv(run2[[2]], "./run2_2.csv", row.names = FALSE, quote = FALSE)
} else {run1 <- list(read.csv("./run2_1.csv", header = TRUE),
                     read.csv("./run2_2.csv", header = TRUE))}

#Pivot to tidy
run2_lng <- as.data.frame(
  tidyr::pivot_longer(run2[[1]],
                      cols = -c("uniq_run", "inoc_r_n", "inoc_r_p", 
                                "init_N_dens", "init_P_dens", "chi", "c_A", 
                                "c_R", "yield", "i", "time"),
                      names_to = c("pop", "x"),
                      names_sep = "_",
                      values_to = "density"))

run2_sum <- 
  dplyr::summarise(
    dplyr::group_by(run2_lng[run2_lng$time == max(run2_lng$time) &
                               run2_lng$pop == "N", ],
                    uniq_run, inoc_r_n, inoc_r_p, init_N_dens, init_P_dens,
                    chi, c_A, c_R, yield, i),
    tot_dens = sum(density),
    dx = 45000/length(unique(x)),
    percentile_95 = dx*max(which(cumsum(density) < 0.95*tot_dens)),
    percentile_90 = dx*max(which(cumsum(density) < 0.90*tot_dens)),
    first_1000cfuml = dx*max(which(density >= 10**3)),
    first_10000cfuml = dx*max(which(density >= 10**4)))

tiff("./Modeling_plots/run2_contour4.tiff", width = 10, height = 10,
     units = "in", res = 300)
ggplot(data = run2_sum, aes(x = chi, y = c_A, z = first_1000cfuml)) +
  geom_contour_filled() +
  # facet_grid(inoc_r_p*c_A ~ c_R,
  #            labeller = labeller(inoc_r_p = c("0" = "Control", 
  #                                             "1425" = "Local", "45000" = "Global"))) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +
  # labs(x = "chemotaxis sensitivity", y = "infection rate",
  #      subtitle = "Resource Consumption Rate", 
  #      title = "Threshold 1000 cfu/mL Percentile")
  NULL
dev.off()
