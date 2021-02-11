library(deSolve)
library(ggplot2)

#Define function to calculate derivatives
derivs <- function(time, y, parms, nx, r_mid, r_end, dx, disp_dx) {
  with(as.list(parms), {
    N <- y[1:nx]
    P <- y[(nx+1):(2*nx)]
    R <- y[(2*nx+1):(3*nx)]
    A <- y[(3*nx+1):(4*nx)]
    
    #Calculate diffusion & migration (zero gradient at boundaries)
    #TODO fix by r_end, rr and dr into flux
    flux_N <- -diff(r_end * -D_N * diff(c(N[1], N, N[nx]))/disp_dx)/r_mid/dx -
      chi*N*diff(diff(c(A[1], A, A[nx]))/disp_dx)/r_mid/dx
    flux_P <- -diff(r_end * -D_P * diff(c(P[1], P, P[nx]))/disp_dx)/r_mid/dx
    flux_R <- -diff(r_end * -D_R * diff(c(R[1], R, R[nx]))/disp_dx)/r_mid/dx
    flux_A <- -diff(r_end * -D_A * diff(c(A[1], A, A[nx]))/disp_dx)/r_mid/dx
    
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
    
#Calculate physical dimensions
radius <- 45000                      #um
nx <- 10                           #concentric circles
dx <- radius/nx                     #thickness of each circle
r_mid <- seq(dx/2,by = dx,len = nx) #distance from center to ea mid-layer
r_end <- seq(0,by = dx,len = nx+1)  #distance from center to end of each layer
disp_dx <- dx                        #dispersion distances

res_scale <- 1

#Define parameters
parms <- c(D_N = 50, D_P = 0, D_R = 800, D_A = 800,
           chi = 0*300,
           c_A = res_scale*4*10**-16, c_R = res_scale*8*10**-11,
           yield = 1*10**10/res_scale,
           k_R = res_scale*5*10**-5, k_A = res_scale*1*10**-6,
           i = 0*10**-12, b = 0*50)

#Define init conditions (NPRA)
y_init <- matrix(c(rep(.5, 3), rep(0, nx-3),
                   rep(0, nx),
                   rep(0.05, nx),
                   #c(rep(0.5, 5), rep(0, nx-5)),
                   rep(1*10**-3, nx)),
                 ncol = 4)

#Define times (in seconds)
times = seq(from = 0, to = 24*60*60, by = 60)
#to = 24*60*60, by = 5*60)

#Run model
yout <- ode.1D(y = y_init, times = times,
               func = derivs, parms = parms,
               nspec = 4, names = c("N", "P", "R", "A"),
               nx = nx, dx=dx, r_mid = r_mid, r_end = r_end, 
               disp_dx = disp_dx)
yout_df <- as.data.frame(yout)

#Reorganize
colnames(yout_df)[2:ncol(yout_df)] <- 
  c(paste(rep("N", nx), 1:nx, sep = "_"), 
    paste(rep("P", nx), 1:nx, sep = "_"),
    paste(rep("R", nx), 1:nx, sep = "_"),
    paste(rep("A", nx), 1:nx, sep = "_"))

yout_lng <- as.data.frame(tidyr::pivot_longer(yout_df, cols = -time,
                                              names_to = c("pop", "x"),
                                              names_sep = "_",
                                              values_to = "density"))


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

ggplot(data = yout_lng[yout_lng$pop %in% c("R", "N", "A"), ], 
       aes(x = as.numeric(time), 
           y = as.numeric(x), 
           z = as.numeric(density))) +
  #color = log10(as.numeric(ifelse(density>0, density, 0)+10)))) +
  geom_contour_filled() +
  #geom_line(aes(group = x)) +
  facet_grid(~pop) +
  #scale_color_continuous(name = "dens") +
  NULL

ggplot(data = yout_lng[yout_lng$time %in% summary(yout_lng$time) &
                         yout_lng$pop %in% c("R", "N", "A"), ],
       aes(x = as.numeric(x), y = as.numeric(density)+0)) +
  facet_grid(pop ~ ., scales = "free") +
  geom_line(aes(color = as.factor(time))) +
  #scale_y_continuous(trans = "log10") +
  #xlim(NA, 60) +
  #geom_hline(yintercept = 10, lty = 2) +
  NULL
