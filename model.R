library(deSolve)
library(ggplot2)

#Define derivs func
derivs <- function(t, yvals, parms, n_x, dx) {
  #The derivs function must return the derivative of all the variables at a
  # given time, in a list
  
  #Parms contains:  D_N, D_P, D_R, D_A (diffusion coefficients) (um/s)
  #                 chi - chemotaxis sensitivity coeff (um/s)
  #                 c_R, c_A - consumption rate of resources or attractant/cell 
  #                       (mmol/cfu/s)
  #                 y - cell yield per resource (cfu/mmol)
  #                 k_R, k_A - Michaelis-Menten rate params (mmol/mL)
  #                 N_0, P_0 - cfu or pfu/mL
  #                 i - infection rate (infec/cell/ml/pfu/ml/s)
  #                 b - burst size (pfu/infec)
  #                 N_0, P_0 - initial densities (cfu or pfu/mL)
  #                 R_0, A_0 - initial densities (mmol/mL)

  #n_x is the number of grid spaces
  N <- yvals[1:n_x]
  P <- yvals[(n_x+1):(2*n_x)]
  R <- yvals[(2*n_x+1):(3*n_x)]
  A <- yvals[(3*n_x+1):(4*n_x)]
  
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
          c_R*y*(R/(R+k_R))*N - i*N*P
    dP = D_P*calc_gradient(grad_P, dx) + i*(b-1)*N*P
    dR = D_R*calc_gradient(grad_R, dx) - c_R*(R/(R+k_R))*N
    dA = D_A*calc_gradient(grad_A, dx) - c_A*(A/(A+k_A))*N
    
    return(list(c(dN, dP, dR, dA)))
  })
}

parms <- c(D_N = 50, D_P = 0, D_R = 800, D_A = 800,
           chi = 300,
           c_A = 4*10**-16, c_R = 8*10**-11,
           y = 1*10**10,
           k_R = 5*10**-5, k_A = 1*10**-6,
           i = 10**-12, b = 50)
n_x <- 1000
dx <- 45000/n_x #in um
#In order of N, P, R, A
# (and arbitrarily setting the volume of each dx to 1 mL)
y_init = matrix(c(rep(2450000, round(1425/dx)), rep(0, n_x-round(1425/dx)),
                  rep(0, n_x),
                  rep(.05, n_x),
                  rep(1*10**-3, n_x)),
                ncol = 4)
#time in seconds
times = seq(from = 0, to = 24*60*60, by = 5*60)

#Run derivs
yout <- as.data.frame(ode.1D(y = y_init, times = times,
               func = derivs, parms = parms, n_x = n_x, dx=dx,
               nspec = 4, names = c("N", "P", "R", "A")))

#1st col: time, cols 2:n_x+1 - density of N at spot x
#         then start w/ P cols, then start w/ R cols

#Reorganize
colnames(yout)[2:ncol(yout)] <- 
  c(paste(rep("N", n_x), 1:1000, sep = "_"), 
    paste(rep("P", n_x), 1:1000, sep = "_"),
    paste(rep("R", n_x), 1:1000, sep = "_"),
    paste(rep("A", n_x), 1:1000, sep = "_"))

yout_lng <- as.data.frame(tidyr::pivot_longer(yout, cols = -time,
                               names_to = c("pop", "x"),
                               names_sep = "_",
                               values_to = "density"))

#Plot
ggplot(data = yout_lng[yout_lng$x %in% 499:502, ],
       aes(x = time/3600, y = density+10, color = pop)) +
  geom_line() +
  facet_grid(pop~x, scales = "free") +
  scale_y_continuous(trans = "log10")

n <- ggplot(data = yout_lng[yout_lng$pop == "N", ], 
            aes(x = as.numeric(time)/3600, 
                            y = as.numeric(x), z = log10(density+10))) +
  geom_contour_filled() +
  facet_grid(~pop)
p <- ggplot(data = yout_lng[yout_lng$pop == "P", ], 
            aes(x = as.numeric(time)/3600, 
                y = as.numeric(x), z = log10(density+10))) +
  geom_contour_filled() +
  facet_grid(~pop)
r <- ggplot(data = yout_lng[yout_lng$pop == "R", ], 
            aes(x = as.numeric(time)/3600, 
                y = as.numeric(x), 
                z = as.numeric(log10(density+10)))) +
  geom_contour_filled() +
  facet_grid(~pop)
a <- ggplot(data = yout_lng[yout_lng$pop == "A", ], 
            aes(x = as.numeric(time)/3600, 
                y = as.numeric(x), 
                z = as.numeric(log10(density+10)))) +
  geom_contour_filled() +
  facet_grid(~pop)

tiff("model_plot.tiff", width = 10, height = 10, units = "in", res = 300)
cowplot::plot_grid(n, p, r, a, nrow = 2)
dev.off()

#ggplot(data = yout_lng[yout_lng$time == 