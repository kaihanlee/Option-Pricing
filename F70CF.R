################################################################
###  Methods for F70CF  ########################################
###    - implicit Finite Difference Method                ######
###       for European options                            ######
################################################################

"implicit_FD_one_step" <- function(pt, p0, r, sigma, deltaT){
  # one step from t+deltaT to t in the implicit FD scheme
  # pt are the prices at time t+deltaTr (end of subperiod)
  # p0 are prices at time t (start of subperiod), 
  #    only p0[1] and p0[M+1] are used: 
  #       - p0[1] is the price for S_t = 0
  #       - p0[M+1] is the option price for the highest possible 
  #                 share price, S_t = S_max 
  # r is the interest rate p.a.
  # sigma is the volatility p.a.
  # deltaT is the length of the time step, length of a subperiod
  # the function returns p0, the option prices at time t, the 
  #    start of subperiod [t, t+deltaT]
  M = length(pt)-1  
  j = 1:(M-1)   # M-1 inner points
  a = 0.5*(r*j - sigma^2 * j^2)*deltaT
  b = 1+(sigma^2 * j^2 + r)*deltaT
  c = -0.5*(sigma^2 * j^2 + r*j)*deltaT
  # construct linear system of equations
  pMat = matrix(0, nrow = M-1, ncol = M-1)
  for (k in j){
    if (k>1) pMat[k,k-1] = a[k]
    pMat[k,k] = b[k]
    if (k<(M-1)) pMat[k,k+1] = c[k]
  }
  p_iplus1 = pt[2:M]    # inner points
  p_iplus1[1]   = p_iplus1[1]   - a[1]*p0[1]
  p_iplus1[M-1] = p_iplus1[M-1] - c[M-1]*p0[M+1]
  p_i = solve(pMat, p_iplus1) 
  p0[2:M] = p_i
  return(p0)
}

"FD_Euro" <- function(p, r, sigma, deltaT){
  # Finite Difference scheme for European derivatives
  # p is the matrix of option prices, only p]1,], p[M+1,] and 
  #    p[,N+1] are used, see comments in implicit_FD_one_step
  # the function returns the matrix p with the derivative prices
  #     at all grid points
  N = ncol(p)-1   # number of time steps
  for (t in (N:1)){ # start with last column (time = N * deltaT)
      p[,t] = implicit_FD_one_step(p[,t+1], p[,t], r, sigma, deltaT)
  }
  return(p)
}

"FD_Amer" <- function(p, r, sigma, deltaT){
  N = ncol(p)-1
  
  for (t in (N:1)){ # start with last column (time = N * deltaT)
    p[,t] = implicit_FD_one_step(p[,t+1], p[,t], r, sigma, deltaT)
    
    # compare all elements in a FD column to the early exercise payoff
    for (j in (2:M)){
      if (p[j,t] <= p[j,t+1]){
        p[j,t] = p[j,t+1]     # take early exercise payoff if it is greater
      }else{
        p[j,t] = p[j,t]     # take put option price if it is greater
      }
    }
    
  }
  return(p)
}

