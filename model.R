library(deSolve)

#Define derivs func
derivs <- function(t, y, parms, n_x) {
  #The derivs function must return the derivative of all the variables at a
  # given time, in a list
  
  #Parms contains:  D_N, D_P, D_R (diffusion coefficients)
  #                 chi - chemotaxis sensitivity coeff
  #                 u_N - bact growth rate
  #                 k - bact carrying capacity
  #                 a - infection rate
  #                 b - burst size
  #                 c - consumption rate/cell
  
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
          u_N*N*(1-(N/k)) - a*N*P
    dP = D_P*calc_gradient(grad_P) + a*(b-1)*N*P
    dR = D_R*calc_gradient(grad_R) - c*N
    
    return(list(c(dN, dP, dR)))
  })
}

#time in hours
parms <- c(D_N = 1, D_P = 5, D_R = 10,
           chi = 1,
           u_N = 1,
           k = 10**9,
           a = 10**-12,
           b = 50,
           c = 1)
n_x <- 1001
y_init = matrix(c(rep(0, (n_x-1)/2), 10**6, rep(0, (n_x-1)/2),
                  rep(0, n_x),
                  rep(10, n_x)),
                ncol = 3)
times = seq(from = 0, to = 24, by = 0.1)

#Run derivs
yout <- as.data.frame(ode.1D(y = y_init, times = times,
               func = derivs, parms = parms, n_x = n_x,
               nspec = 3, names = c("N", "P", "R")))

#1st col: time, cols 2:n_x+1 - density of N at spot x
#         then start w/ P cols, then start w/ R cols
plot(y = yout$`501`, x = yout$time)