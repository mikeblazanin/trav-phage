library(deSolve)
library(ggplot2)

#Define derivs func
derivs <- function(t, y, parms, n_x, dx) {
  #The derivs function must return the derivative of all the variables at a
  # given time, in a list
  
  #Parms contains:  D_N, D_P, D_R, D_A (diffusion coefficients) (um/s)
  #                 chi - chemotaxis sensitivity coeff (um/s)
  #                 c_R, c_A - consumption rate of resources or attractant/cell 
  #                       (mmol/cfu/s)
  #                 yield - cell yield per resource (cfu/mmol)
  #                 k_R, k_A - Michaelis-Menten rate params (mmol/mL)
  #                 N_0, P_0 - cfu or pfu/mL
  #                 i - infection rate (infec/cell/ml/pfu/ml/s)
  #                 b - burst size (pfu/infec)
  #                 N_0, P_0 - initial densities (cfu or pfu/mL)
  #                 R_0, A_0 - initial densities (mmol/mL)
  
  #n_x is the number of grid spaces
  N <- y[1:n_x]
  P <- y[(n_x+1):(2*n_x)]
  R <- y[(2*n_x+1):(3*n_x)]
  A <- y[(3*n_x+1):(4*n_x)]
  
  N[N < 0] <- 0
  P[P < 0] <- 0
  R[R < 0] <- 0
  A[A < 0] <- 0
  
  #dx is the size of each grid space (in um)
  calc_gradient <- function(x, dx) {
    # return(c(diff(c(x[1:2])),
    #           diff(c(x[1:length(x)]), lag = 2)/2,
    #           diff(c(x[(length(x)-1):length(x)]))))
    return(diff(c(x[2], x, x[length(x)-1]), lag = 2)/dx)
    #return(c(diff(x), 0))
  }
  
  #Calculate changes
  with(as.list(parms), {
    grad_N <- calc_gradient(N, dx)
    grad_P <- calc_gradient(P, dx)
    grad_R <- calc_gradient(R, dx)
    grad_A <- calc_gradient(A, dx)
    
    dN = calc_gradient(D_N*grad_N - chi*N*grad_A, dx) +
          c_R*yield*(R/(R+k_R))*N - i*N*P
    dP = D_P*calc_gradient(grad_P, dx) + i*(b-1)*N*P
    dR = D_R*calc_gradient(grad_R, dx) - c_R*(R/(R+k_R))*N
    dA = D_A*calc_gradient(grad_A, dx) - c_A*(A/(A+k_A))*N
    
    return(list(c(dN, dP, dR, dA)))
  })
  #y <- c(as.matrix(yout[9, -1]))
}

#Scale all parameters in units of mols accordingly
# res_scale > 1 scales the unit down such that measures are larger
# i.e. res_scale = 10**3 scales all measures from mmol to umol
res_scale <- 10**12
parms <- c(D_N = 0*50, D_P = 0, D_R = 800, D_A = 800,
           chi = 0*300,
           c_A = 0*res_scale*0*4*10**-16, c_R = 0*res_scale*0*8*10**-11,
           yield = 0*1*10**10/res_scale,
           k_R = res_scale*5*10**-5, k_A = res_scale*1*10**-6,
           i = 0*10**-12, b = 0*50)
n_x <- 4500
dx <- 45000/n_x #in um
#In order of N, P, R, A
# (setting each dx to a cube with volume (in cm^3=mL))
dx_vol <- (dx/10000)**3
y_init = matrix(#c(rep(2450000*dx_vol, round(1425/dx)), rep(0, n_x-round(1425/dx)),
                c(rep(0, n_x),
                  rep(0*dx_vol, n_x),
                  #res_scale*rep(.05*dx_vol, n_x),
                  res_scale*c(rep(.05, 10), rep(0, n_x-10)),
                  res_scale*rep(dx_vol*1*10**-3, n_x)),
                ncol = 4)
#time in seconds
times = seq(from = 0, to = 60, by = 1)
            #to = 24*60*60, by = 5*60)

#Run derivs
yout <- ode.1D(y = y_init, times = times,
               func = derivs, parms = parms, n_x = n_x, dx=dx,
               nspec = 4, dimens = n_x, names = c("N", "P", "R", "A"),
               hmax = 0.25)
yout_df <- as.data.frame(yout)


#Reorganize
colnames(yout_df)[2:ncol(yout_df)] <- 
  c(paste(rep("N", n_x), 1:n_x, sep = "_"), 
    paste(rep("P", n_x), 1:n_x, sep = "_"),
    paste(rep("R", n_x), 1:n_x, sep = "_"),
    paste(rep("A", n_x), 1:n_x, sep = "_"))

yout_lng <- as.data.frame(tidyr::pivot_longer(yout_df, cols = -time,
                               names_to = c("pop", "x"),
                               names_sep = "_",
                               values_to = "density"))

#Plot
# ggplot(data = yout_lng[yout_lng$x %in% 499:502, ],
#        aes(x = time/3600, y = density+10, color = pop)) +
#   geom_line() +
#   facet_grid(pop~x, scales = "free") +
#   scale_y_continuous(trans = "log10")

n <- ggplot(data = yout_lng[yout_lng$pop == "N", ], 
            aes(x = as.numeric(time)/3600, 
                            y = as.numeric(x), 
                z = log10(density+10))) +
  geom_contour_filled() +
  facet_grid(~pop)
p <- ggplot(data = yout_lng[yout_lng$pop == "P", ], 
            aes(x = as.numeric(time)/3600, 
                y = as.numeric(x), 
                z = log10(density+10))) +
  geom_contour_filled() +
  facet_grid(~pop)
r <- ggplot(data = yout_lng[yout_lng$pop == "R", ], 
            aes(x = as.numeric(time)/3600, 
                y = as.numeric(x), 
                color = log10(as.numeric(ifelse(density>0, density, 0)+10)))) +
  #geom_contour_filled() +
  geom_line(aes(group = x)) +
  facet_grid(~pop) +
  scale_color_continuous(name = "dens") +
  NULL
a <- ggplot(data = yout_lng[yout_lng$pop == "A", ], 
            aes(x = as.numeric(time)/3600, 
                y = as.numeric(x), 
                color = log10(as.numeric(ifelse(density>0, density, 0)+10)))) +
  #geom_contour_filled() +
  geom_line(aes(group = x)) +
  facet_grid(~pop) +
  scale_color_continuous(name = "dens") +
  NULL

tiff("model_plot.tiff", width = 10, height = 10, units = "in", res = 300)
cowplot::plot_grid(n, p, r, a, nrow = 2)
dev.off()

tiff("model_plot2.tiff", width = 10, height = 10, units = "in", res = 300)
ggplot(data = yout_lng[yout_lng$time %in% 0:3 &
                         yout_lng$pop == "R", ],
       aes(x = as.numeric(x), y = as.numeric(density)+0)) +
  facet_grid(pop ~ time, scales = "free") +
  geom_line() +
  #scale_y_continuous(trans = "log10") +
  xlim(NA, 60) +
  #geom_hline(yintercept = 10, lty = 2) +
  NULL
dev.off()

ggplot(data = yout_lng[yout_lng$x %in% 7:15, ],
       aes(x = as.numeric(time), y = as.numeric(density)+0)) +
  facet_grid(pop ~ x, scales = "free") +
  geom_line() +
  scale_y_continuous(trans = "log10") +
  #xlim(NA, 500) +
  #geom_hline(yintercept = 10, lty = 2) +
  NULL

yout_lng[yout_lng$pop == "R" & as.numeric(yout_lng$time) %in% c(9, 10) &
           yout_lng$x %in% c(1:100), ]
