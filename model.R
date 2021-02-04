library(deSolve)
library(ggplot2)

#Define derivs func
derivs <- function(t, y, parms, n_x) {
  #The derivs function must return the derivative of all the variables at a
  # given time, in a list
  
  #Parms contains:  D_N, D_P, D_R (diffusion coefficients)
  #                 chi - chemotaxis sensitivity coeff
  #                 a - infection rate
  #                 b - burst size
  #                 c - consumption rate of resources/cell
  #                 m - constant metabolic maintenance cost
  
  #n_x is the number of grid spaces
  N <- y[1:n_x]
  P <- y[(n_x+1):(2*n_x)]
  R <- y[(2*n_x+1):(3*n_x)]
  
  N[N < 0] <- 0
  P[P < 0] <- 0
  R[R < 0] <- 0
  
  calc_gradient <- function(x) {
    return(c(diff(c(x[1:2])),
              diff(c(x[1:length(x)]), lag = 2)/2,
              diff(c(x[(length(x)-1):length(x)]))))
    #return(diff(c(x[1], x, x[length(x)])))
    #return(c(diff(x), 0))
  }
  
  #Calculate changes
  with(as.list(parms), {
    grad_N <- calc_gradient(N)
    grad_P <- calc_gradient(P)
    grad_R <- calc_gradient(R)
    
    #Have a problem where N has to be diff-d twice while the others
    # only once
    # and diff always reduces size of vector by 1
    # so have to decide how to handle edge cases
    
    dN = calc_gradient(D_N*grad_N - chi*N*grad_R) +
          N*(c*R - m) - a*N*P
    dP = D_P*calc_gradient(grad_P) + a*(b-1)*N*P
    dR = D_R*calc_gradient(grad_R) - c*N*R
    
    return(list(c(dN, dP, dR)))
  })
}

#time in hours
parms <- c(D_N = 0, D_P = 0, D_R = 0,
           chi = 0,
           a = 10**-10,
           b = 50,
           c = 10**-9,
           m = 1)
n_x <- 1001
y_init = matrix(c(rep(0, (n_x-1)/2), 10**5, rep(0, (n_x-1)/2),
                  rep(10**2, n_x),
                  rep(10**9, n_x)),
                ncol = 3)
times = seq(from = 0, to = 24, by = 0.1)

#Run derivs
yout <- as.data.frame(ode.1D(y = y_init, times = times,
               func = derivs, parms = parms, n_x = n_x,
               nspec = 3, names = c("N", "P", "R")))

#1st col: time, cols 2:n_x+1 - density of N at spot x
#         then start w/ P cols, then start w/ R cols

#Reorganize
yout_lng <- as.data.frame(tidyr::pivot_longer(yout, cols = -time,
                               names_to = "col",
                               values_to = "density"))
yout_lng$pop <- NA
yout_lng$x <- NA
temp_rows <- which(as.numeric(yout_lng$col) <= (ncol(yout)-1)/3)
yout_lng$pop[temp_rows] <- "N"
yout_lng$x[temp_rows] <- yout_lng$col[temp_rows]
temp_rows <- which(as.numeric(yout_lng$col) > (ncol(yout)-1)/3 & 
                     as.numeric(yout_lng$col) <= 2*(ncol(yout)-1)/3)
yout_lng$pop[temp_rows] <- "P"
yout_lng$x[temp_rows] <- as.numeric(yout_lng$col[temp_rows])-(ncol(yout)-1)/3
temp_rows <- which(as.numeric(yout_lng$col) > 2*(ncol(yout)-1)/3 & 
                     as.numeric(yout_lng$col) <= (ncol(yout)-1))
yout_lng$pop[temp_rows] <- "R"
yout_lng$x[temp_rows] <- as.numeric(yout_lng$col[temp_rows])-2*(ncol(yout)-1)/3
if(any(is.na(yout_lng$pop))) {stop("some cols uncategorized")}

ggplot(data = yout_lng[yout_lng$x == "501", ],
       aes(x = time, y = density+1, color = pop)) +
  geom_line() +
  facet_grid(pop~., scales = "free") +
  scale_y_continuous(trans = "log10")

ggplot(data = yout_lng, aes(x = as.numeric(time), 
                            y = as.numeric(x), z = log10(density+10))) +
  geom_contour_filled() +
  facet_grid(~pop)
