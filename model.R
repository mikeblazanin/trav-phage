#Import & set global options ----
library(deSolve)
library(ggplot2)

make_statplots <- TRUE
make_curveplots <- FALSE

#Define function to calculate derivatives ----
derivs <- function(time, y, parms, nx, r_mid, r_end, dx, disp_dx) {
  #Parms contains:
  #              D_N, D_P, D_R, D_A,
  #              chi1, chi2,
  #              c1_A, c2_A,
  #              c1_R, c2_R,
  #              yield,
  #              k_R, k_A,
  #              i1, i2, b
  
  cen_diff <- function(x) {
    return(diff(c(x[2], x, x[length(x)-1]),
                lag = 2))
  }

  with(as.list(parms), {
    N1 <- y[1:nx]
    N2 <- y[(nx+1):(2*nx)]
    P <- y[(2*nx+1):(3*nx)]
    R <- y[(3*nx+1):(4*nx)]
    A <- y[(4*nx+1):(5*nx)]
    
    #Calculate diffusion & migration (zero gradient at boundaries)
    flux_N1 <- -cen_diff(r_end * -D_N * cen_diff(N1)/disp_dx)/r_mid/dx -
      chi1*cen_diff(r_end * N1 * cen_diff(A)/disp_dx)/r_mid/dx
    flux_N2 <- -cen_diff(r_end * -D_N * cen_diff(N2)/disp_dx)/r_mid/dx -
      chi2*cen_diff(r_end * N2 * cen_diff(A)/disp_dx)/r_mid/dx
    flux_P <- -cen_diff(r_end * -D_P * cen_diff(P)/disp_dx)/r_mid/dx
    flux_R <- -cen_diff(r_end * -D_R * cen_diff(R)/disp_dx)/r_mid/dx
    flux_A <- -cen_diff(r_end * -D_A * cen_diff(A)/disp_dx)/r_mid/dx
    
    #Calculate growth and death
    grow_N1 <- c1_R*yield*R/(R+k_R)*N1 - i1*N1*P
    grow_N2 <- c2_R*yield*R/(R+k_R)*N2 - i2*N2*P
    grow_P <- i1*(b-1)*N1*P + i2*(b-1)*N2*P
    
    #Calculate consumption
    cons_R <- c1_R*R/(R+k_R)*N1 + c2_R*R/(R+k_R)*N2
    cons_A <- c1_A*A/(A+k_A)*N1 + c2_A*A/(A+k_A)*N2
    
    #Change = diffusion&migration + growth&death + consumption
    dN1 <- flux_N1 + grow_N1
    dN2 <- flux_N2 + grow_N2
    dP <- flux_P + grow_P
    dR <- flux_R - cons_R
    dA <- flux_A - cons_A
    
    return(list(c(dN1, dN2, dP, dR, dA)))
  })
}

#Simple model (for future reference) ----
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

#Define function to run simulations across grid of values ----
run_sims <- function(inoc_r_n1 = 1425, inoc_r_n2 = 0, inoc_r_p = 1425,
                     init_N1_dens = 250000/1425, init_N2_dens = 250000/1425,
                     init_P_dens = 25000/1425,
                     chi1 = 300, chi2 = 300,
                     c1_A = 4*10**-13, c2_A = 4*10**-13,
                     c1_R = 2*10**-11, c2_R = 2*10**-11,
                     yield = 10**7, 
                     i1 = 10**-12, i2 = 10**-12, 
                     radius = 45000, nx = 100,
                     times = seq(from = 0, to = 24*60*60, by = 1*60),
                     times_keep = seq(from = 0, to = 24*60*60, by = 15*60),
                     D_N = 50, D_P = 0, D_R = 800, D_A = 800,
                     k_R = 0.05, k_A = 0.001, b = 50,
                     init_R_dens = 2000/nx, init_A_dens = 25/nx,
                     combinatorial = TRUE, print_info = TRUE) {
  
  #nx is the number of concentric circles the plate is broken up into
  #inoc_r is the radius, in um, n and p are respectively inoculated into
  #init_N and P are the total inoculated bacterial & phage populations
  # (which will be distributed over inoc_r_n area)
  
  #Input checks
  stopifnot(
    length(radius) == 1, length(nx) == 1, 
    length(D_N) == 1, length(D_P) == 1, length(D_R) == 1, length(D_A) == 1, 
    length(k_R) == 1, length(k_A) == 1, length(b) == 1,
    length(init_R_dens) == 1, length(init_A_dens) == 1,
    length(times) > 1, all(times_keep %in% times))
  
  #Save parameter value combinations into dataframe
  if (combinatorial == TRUE) {
    param_combos <- expand.grid(list(
      "inoc_r_n1" = inoc_r_n1, "inoc_r_n2" = inoc_r_n2, 
      "inoc_r_p" = inoc_r_p,
      "init_N1_dens" = init_N1_dens, "init_N2_dens" = init_N2_dens,
      "init_P_dens" = init_P_dens,
      "chi1" = chi1, "chi2" = chi2, 
      "c1_A" = c1_A, "c2_A" = c2_A, "c1_R" = c1_R, "c2_R" = c2_R,
      "yield" = yield, "i1" = i1, "i2" = i2),
      stringsAsFactors = FALSE)
    num_sims <- nrow(param_combos)
  } else {
    num_sims <- max(sapply(X = list(inoc_r_n1, inoc_r_n2, inoc_r_p, 
                                    init_N1_dens, init_N2_dens, init_P_dens,
                                    chi1, chi2, c1_A, c2_A, c1_R, c2_R, 
                                    yield, i1, i2), FUN = length))
    
    #Check for parameter lengths being non-divisible with the
    # number of simulations inferred from the longest parameter length
    if (!all(num_sims %% sapply(X = list(inoc_r_n1, inoc_r_n2, inoc_r_p, 
                                         init_N1_dens, init_N2_dens, init_P_dens,
                                         chi1, chi2, c1_A, c2_A, c1_R, c2_R, 
                                         yield, i1, i2), 
                                FUN = length) == 0)) {
      warning("Combinatorial=TRUE but longest param vals length is not a multiple of all other param vals lengths")
    }
    
    #Save parameters into dataframe, replicating shorter parameter
    # vectors as needed to reach # of simulations
    param_combos <- data.frame("inoc_r_n1" = rep_len(inoc_r_n1, num_sims), 
                               "inoc_r_n2" = rep_len(inoc_r_n2, num_sims), 
                               "inoc_r_p" = rep_len(inoc_r_p, num_sims),
                               "init_N1_dens" = rep_len(init_N1_dens, num_sims), 
                               "init_N2_dens" = rep_len(init_N2_dens, num_sims), 
                               "init_P_dens" = rep_len(init_P_dens, num_sims),
                               "chi1" = rep_len(chi1, num_sims), 
                               "chi2" = rep_len(chi2, num_sims), 
                               "c1_A" = rep_len(c1_A, num_sims), 
                               "c2_A" = rep_len(c2_A, num_sims), 
                               "c1_R" = rep_len(c1_R, num_sims), 
                               "c2_R" = rep_len(c2_R, num_sims), 
                               "yield" = rep_len(yield, num_sims), 
                               "i1" = rep_len(i1, num_sims),
                               "i2" = rep_len(i2, num_sims),
                               stringsAsFactors = FALSE)
  }
  
  #Define physical parameters of plate (true for all runs)
  # (radius is default defined in um)
  # (time is default defines in seconds)
  dx <- radius/nx                      #thickness of each circle
  r_mid <- seq(dx/2,by = dx,len = nx)  #distance from center to ea mid-layer
  r_end <- seq(0,by = dx,len = nx+1)   #distance from center to end of each layer
  disp_dx <- dx                        #dispersion distances
  
  #Create output dataframe (yout)
  bigout <- 
    cbind(
      data.frame(
        "uniq_run" = rep(1:nrow(param_combos), each = length(times_keep)),
        "inoc_r_n1" = rep(param_combos$inoc_r_n1, each = length(times_keep)), 
        "inoc_r_n2" = rep(param_combos$inoc_r_n2, each = length(times_keep)), 
        "inoc_r_p" = rep(param_combos$inoc_r_p, each = length(times_keep)),
        "init_N1_dens" = rep(param_combos$init_N1_dens, each = length(times_keep)), 
        "init_N2_dens" = rep(param_combos$init_N2_dens, each = length(times_keep)), 
        "init_P_dens" = rep(param_combos$init_P_dens, each = length(times_keep)),
        "chi1" = rep(param_combos$chi1, each = length(times_keep)), 
        "chi2" = rep(param_combos$chi2, each = length(times_keep)), 
        "c1_A" = rep(param_combos$c1_A, each = length(times_keep)), 
        "c2_A" = rep(param_combos$c2_A, each = length(times_keep)), 
        "c1_R" = rep(param_combos$c1_R, each = length(times_keep)), 
        "c2_R" = rep(param_combos$c2_R, each = length(times_keep)), 
        "yield" = rep(param_combos$yield, each = length(times_keep)), 
        "i1" = rep(param_combos$i1, each = length(times_keep)),
        "i2" = rep(param_combos$i2, each = length(times_keep)),
        stringsAsFactors = FALSE),
      as.data.frame(matrix(NA, nrow = length(times_keep)*nrow(param_combos), 
                           ncol = 1+5*nx,
                           dimnames = list(NULL, c("time", 1:(5*nx))))))
  
  #Print number of simulations that will be run
  if(print_info) {
    print(paste(num_sims, "simulations will be run"))
    
    #Save sequence of 10% cutoff points for later reference
    progress_seq <- round(seq(from = 0, to = num_sims, by = num_sims/10))
  }
  
  #Run simulations
  for (myrun in 1:nrow(param_combos)) {
    #Define parameters vector
    parms <- c(D_N = D_N, D_P = D_P, D_R = D_R, D_A = D_A,
               chi1 = param_combos$chi1[myrun],
               chi2 = param_combos$chi2[myrun],
               c1_A = param_combos$c1_A[myrun], 
               c2_A = param_combos$c2_A[myrun],
               c1_R = param_combos$c1_R[myrun],
               c2_R = param_combos$c2_R[myrun],
               yield = param_combos$yield[myrun],
               k_R = k_R, k_A = k_A,
               i1 = param_combos$i1[myrun], i2 = param_combos$i2[myrun], b = b)
    
    #Define initial conditions
    y_init = matrix(c(
      #N1
      rep(param_combos$init_N1_dens[myrun], 
          round(param_combos$inoc_r_n1[myrun]/dx)),
      rep(0, nx - round(param_combos$inoc_r_n1[myrun]/dx)),
      #N2
      rep(param_combos$init_N2_dens[myrun], 
          round(param_combos$inoc_r_n2[myrun]/dx)),
      rep(0, nx - round(param_combos$inoc_r_n2[myrun]/dx)),
      #P
      rep(param_combos$init_P_dens[myrun], 
          round(param_combos$inoc_r_p[myrun]/dx)),
      rep(0, nx - round(param_combos$inoc_r_p[myrun]/dx)),
      #R
      rep(init_R_dens, nx),
      #A
      rep(init_A_dens, nx)),
      ncol = 4)
    
    #Run model
    yout <- as.data.frame(ode.1D(y = y_init, times = times,
                                 func = derivs, parms = parms,
                                 nspec = 5, names = c("N1", "N2", "P", "R", "A"),
                                 nx = nx, dx=dx, r_mid = r_mid, r_end = r_end, 
                                 disp_dx = disp_dx, maxsteps = 25000))
    yout <- yout[yout$time %in% times_keep, ]
    
    #Save results
    bigout[((myrun-1)*length(times_keep)+1):(myrun*length(times_keep)), 
           17:ncol(bigout)] <- yout
    
    #Print progress update
    if (print_info & myrun %in% progress_seq) {
      print(paste((which(progress_seq == myrun)-1)*10,
                  "% completed", sep = ""))
    }
  }
  
  #Rename columns
  #Reorganize
  colnames(bigout)[18:ncol(bigout)] <- 
    c(paste(rep("N1", nx), 0:(nx-1), sep = "_"), 
      paste(rep("N2", nx), 0:(nx-1), sep = "_"), 
      paste(rep("P", nx), 0:(nx-1), sep = "_"),
      paste(rep("R", nx), 0:(nx-1), sep = "_"),
      paste(rep("A", nx), 0:(nx-1), sep = "_"))
  
  return(list(bigout, param_combos))
}

#Run 1 ----
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
  tidyr::pivot_longer(run1[[1]][run1[[1]]$time %in% 
                                  c(60*60*c(0, 1, 2, 3, 6, 12, 18, 24)), ],
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
    first_1000cfuml_dens = density[max(which(density >= 10**3))],
    # y - y1 = m (x - x1); m = (y2 - y1)/(x2 - x1) = (y2 - y1)/dx
    # Rearrange: x = dx * (y - y1)/(y2 - y1) + x1
    #   (where x, y is any point on that line. Set y = 1000, calc x)
    first_1000cfuml_interp = dx*(1000 - first_1000cfuml_dens)/
      (density[(max(which(density >= 10**3)) + 1)] - first_1000cfuml_dens) +
      first_1000cfuml)

#Make profile plots for each run
dir.create("./run1_profile_plots", showWarnings = FALSE)
my_times <- 60*60*c(0, 1, 2, 3, 6, 12, 18, 24)
if (make_curveplots) {
  for (uniq_run in unique(run1_lng$uniq_run)) {
    tiff(paste("./run1_profile_plots/", uniq_run, ".tiff", sep = ""),
         width = 6, height = 4, units = "in", res = 100)
    print(ggplot(data = run1_lng[run1_lng$time %in% my_times &
                                   run1_lng$uniq_run == uniq_run, ],
                 aes(x = as.numeric(x)*(45000/run1_nx)/10000, y = as.numeric(density)+1)) +
            facet_wrap(pop ~ ., scales = "free",
                       labeller = labeller(pop = c("A" = "Attractant", "N" = "Bacteria",
                                                   "P" = "Phage", "R" = "Resources"))) +
            geom_line(lwd = 1.25, aes(color = as.factor(time/3600))) +
            geom_point(data = cbind(run1_sum[run1_sum$uniq_run == uniq_run, ],
                                    "pop" = "N"),
                       aes(x = first_1000cfuml_interp/10000,
                           y = 1000+1)) +
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
}
  
#In run 1, when phage are global and infec rate is high (e-12), bact form
# traveling peaks. Otherwise they form expanding fronts.
# (uniq runs 3, 6, 9 after 24 hrs, 12, 15, 18, 21, 24, 27 at multiple high hrs)

if (make_statplots) {         
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
  ggplot(data = run1_sum, 
         aes(x = chi, y = i, z = first_1000cfuml_interp/(24*1000))) +
    geom_contour_filled() +
    facet_grid(inoc_r_p*c_A ~ c_R,
               labeller = labeller(inoc_r_p = c("0" = "Control", 
                                                "1425" = "Local", "45000" = "Global"))) +
    scale_y_continuous(trans = "log10") +
    labs(x = "chemotaxis sensitivity", y = "infection rate",
         subtitle = "Resource Consumption Rate", 
         title = "Soft agar growth (mm/hr) [1000 cfu threshold]")
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
                 aes(x = as.numeric(chi), y = as.numeric(i), 
                     z = as.numeric(first_1000cfuml_interp))) +
            geom_contour_filled() +
            facet_grid(inoc_r_p ~ c_A,
                       labeller = labeller(inoc_r_p = c("0" = "Control", 
                                                        "1425" = "Local", "45000" = "Global"))) +
            scale_y_continuous(trans = "log10") +
            labs(x = "chemotaxis sensitivity", y = "infection rate",
                 subtitle = "Attractant Consumption Rate", 
                 title = "Threshold 1000 cfu/mL Percentile"))
    dev.off()
  }
}

#The changes included in run1 (specifically those for c_A and chi)
# aren't enough to actually affect the migration rate
# so changes in migration rate aren't proportional to changes in c_A and chi

#Run 2 ----
#
#Tested chi up to 3300, c_A up to 33e-13, rate only gets to about 0.5 cm/hr
run2_nx <- 100
if (F) {
  run2 <- run_sims(inoc_r_n1 = 1425, inoc_r_n2 = 1425,
                   inoc_r_p = c(0, 1425, 45000),
                   i2 = c(10**-12, 10**-18),
                   chi2 = c(300, 600),
                   nx = run2_nx)
                   
    
    
    # inoc_r_p = c(0, 1425, 45000),
    #                chi = 300,
    #                c_R = 3*10**-11,
    #                c_A = 4*10**c(-13),
    #                i = 10**c(-8, -9, -10, -11, -12),
    #                D_P = 1,
    #                nx = run2_nx)
  write.csv(run2[[1]], "./run2_1.csv", row.names = FALSE, quote = FALSE)
  write.csv(run2[[2]], "./run2_2.csv", row.names = FALSE, quote = FALSE)
} else {run2 <- list(read.csv("./run2_1.csv", header = TRUE),
                     read.csv("./run2_2.csv", header = TRUE))}

#Pivot to tidy
run2_lng <- as.data.frame(
  tidyr::pivot_longer(run2[[1]],
                      cols = -c("uniq_run", "inoc_r_n1", "inoc_r_n2",
                                "inoc_r_p", 
                                "init_N1_dens", "init_N2_dens",
                                "init_P_dens", "chi1", "chi2", "c1_A", "c2_A",
                                "c1_R", "c2_R", "yield", "i1", "i2", "time"),
                      names_to = c("pop", "x"),
                      names_sep = "_",
                      values_to = "density"))

run2_sum <- 
  dplyr::summarise(
    dplyr::group_by(run2_lng[run2_lng$time == max(run2_lng$time) &
                               run2_lng$pop %in% c("N1", "N2"), ],
                    uniq_run, inoc_r_n1, inoc_r_n2, inoc_r_p,
                    init_N1_dens, init_N2_dens, init_P_dens,
                    chi1, chi2, c1_A, c2_A, c1_R, c2_R, 
                    yield, i1, i2, time, pop),
    tot_dens = sum(density),
    dx = 45000/length(unique(x)),
    percentile_95 = dx*max(which(cumsum(density) < 0.95*tot_dens)),
    percentile_90 = dx*max(which(cumsum(density) < 0.90*tot_dens)),
    first_1000cfuml = dx*max(which(density >= 10**3)),
    first_10000cfuml = dx*max(which(density >= 10**4)),
    first_1000cfuml_dens = density[max(which(density >= 10**3))],
    # y - y1 = m (x - x1); m = (y2 - y1)/(x2 - x1) = (y2 - y1)/dx
    # Rearrange: x = dx * (y - y1)/(y2 - y1) + x1
    #   (where x, y is any point on that line. Set y = 1000, calc x)
    first_1000cfuml_interp = dx*(1000 - first_1000cfuml_dens)/
      (density[(max(which(density >= 10**3)) + 1)] - first_1000cfuml_dens) +
      first_1000cfuml)

tiff("./Modeling_plots/run2_contour4.tiff", width = 10, height = 10,
     units = "in", res = 300)
ggplot(data = run2_sum, aes(x = chi, y = c_A, 
                            z = first_1000cfuml_interp/(24*1000))) +
  geom_contour_filled() +
  facet_grid(~ c_R) +
  # facet_grid(inoc_r_p*c_A ~ c_R,
  #            labeller = labeller(inoc_r_p = c("0" = "Control", 
  #                                             "1425" = "Local", "45000" = "Global"))) +
  # scale_y_continuous(trans = "log10") +
  # scale_x_continuous(trans = "log10") +
  labs(x = "chemotaxis sensitivity", y = "A uptake",
       subtitle = "Resource Consumption Rate", 
       title = "Soft agar growth (mm/hr) [1000 cfu threshold]") +
  NULL
dev.off()

#Make profile plots
#Make profile plots for each run
dir.create("./run2_profile_plots", showWarnings = FALSE)
my_times <- 60*60*c(0, 1, 2, 3, 6, 12, 18, 24)
if (make_curveplots) {
  for (uniq_run in unique(run2_lng$uniq_run)) {
    tiff(paste("./run2_profile_plots/", uniq_run, ".tiff", sep = ""),
         width = 6, height = 4, units = "in", res = 200)
    print(ggplot(data = run2_lng[run2_lng$time %in% my_times &
                                   run2_lng$uniq_run == uniq_run, ],
                 aes(x = as.numeric(x)*(45000/run2_nx)/10000, y = as.numeric(density)+1)) +
            facet_wrap(pop ~ ., scales = "free",
                       labeller = labeller(pop = c("A" = "Attractant", "N1" = "Bacteria1",
                                                   "N2" = "Bacteria2",
                                                   "P" = "Phage", "R" = "Resources"))) +
            geom_line(lwd = 1.25, aes(color = as.factor(time/3600))) +
            # geom_point(data = cbind(run2_sum[run2_sum$uniq_run == uniq_run, ],
            #                         "pop" = "N"),
            #            aes(x = first_1000cfuml/10000,
            #                y = first_1000cfuml_dens+1)) +
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
}

##Run 3 ----
# run3_nx <- 100
# if (F) {
#   run3 <- rbind(run_sims(inoc_r_p = c(0, 1425, 45000),
#                          chi = c(300, 600, 900, 1200),
#                          c_R = 2*10**-11,
#                          c_A = 4*10**c(-13),
#                          i = 10**c(-8, -9, -10, -11, -12),
#                          nx = run2_nx,
#                          D_P = 10,
#     
#     run_sims(inoc_r_p = c(0, 1425, 45000),
#                    chi = 300,
#                    c_R = 3*10**-11,
#                    c_A = 4*10**c(-13),
#                    i = 10**c(-8, -9, -10, -11, -12),
#                    nx = run2_nx,
#                    D_N = 0)
#   write.csv(run3[[1]], "./run3_1.csv", row.names = FALSE, quote = FALSE)
#   write.csv(run3[[2]], "./run3_2.csv", row.names = FALSE, quote = FALSE)
# } else {run3 <- list(read.csv("./run3_1.csv", header = TRUE),
#                      read.csv("./run3_2.csv", header = TRUE))}









