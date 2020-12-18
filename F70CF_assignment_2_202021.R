# F70CF 2020/21
# assignment 2 - Option pricing

##################################################
####  REPLACE 1234 IN THE FOLLOWING LINE WITH ####
####  THE LAST FOUR DIGITS OF YOUR STUDENT ID ####
set.seed(3031)
##################################################
rm(list = ls())    # delete everything from the work space
source("F70CF.R")  # load functions for implicit FD
####  some constants  ############################
K = 1             # strike price
T = 2.5           # time (in years) to maturity
r = 0.05         # interest rate p.a.
sigma = 0.4 + rexp(1, 10)      # volatility p.a.
S_max = 4         # max share price for FD
deltaS = 0.1
M = S_max / deltaS  # number of share price intervals for FD
deltaT = 1 / 12     # length of a time step in FD method
N = T / deltaT      # number of time intervals for FD
##################################################

#### YOUR CODE SHOULD BEGIN AFTER THIS LINE   ####
### Q1 ###########################################
S0 = deltaS * (0:M)     # create a vector for S0 in range 0 to Smax
d1 = (log(S0/K) + (r + (sigma^2)/2) * (T-0)) / (sigma * sqrt(T-0))    # calculate d1 of BS formula
d2 = d1 - sigma * sqrt(T)    # calculate d2 of BS formula

# calculate put option price using BS formula and pnorm function for d1 and d2
put_value = (K * exp(-r * T) * pnorm(d2, lower.tail=FALSE)) - (S0 * pnorm(d1, lower.tail=FALSE))
plot(S0, put_value, type="l", ylim=0:1, xlab="Risky Asset Price", ylab="Put Option Value / Payoff", 
     main="European Put Option Value / Payoff Diagram")   # plot graph

put_payoff1 = (K - S0) * ((K - S0) > 0)   # calculate put option payoff
lines(S0, put_payoff1, col="red")   # add plot of put option payoff
legend(2.5, 0.85, legend=c("Put Option Value", "Put Option Payoff"), col=c("black", "red"), lwd=c(1:2), cex=1.2)

### Q2 ##########################################
###  Boundary condition for FD for Put option  ####################
St = deltaS * (0:M)     # create a vector for St in range 0 to Smax
t = deltaT * (0:N)     # create a vector for time points in range 0 to 30

P = matrix(NA, nrow = M+1, ncol = N+1)    # create matrix for finite difference method
put_payoff_2a = (K - St) * ((K - St) > 0)   # calculate put option payoff
P[,N+1] = put_payoff_2a     # Boundary condition 1: p(N,j) is the put option payoff 
P[1,] = K * exp(-r * (T-t))   # Boundary condition 2: p(i,0) is when St=0
# European options take discounting into account

P[M+1,] = 0     # Boundary condition 3: p(i,M) is when St=Smax

####### FD prices for Put option ##################################
put_fd_euro = FD_Euro(P,r,sigma,deltaT)   # implicit finite difference method for European option
FD_Euro_grid = cbind(St[(M+1):1], round(put_fd_euro[(M+1):1,],3))   # concatenate St column with FD mesh
print(FD_Euro_grid)   # display FD mesh

plot(FD_Euro_grid[,1], FD_Euro_grid[,2], type="l", ylim=0:1, xlab="Risky Asset Price", ylab="Put Option Value", 
     main="European Put Option Value")    # plot graph
points(St, put_value)   # add plot of the put option price from Q1
legend(2.5, 0.85, legend=c("Implicit FD", "Black-Scholes"), lty=c(1,NA),pch=c(NA,1), cex=1.2)

### Q3 ############################################################
P[1,] = K   # Update boundary condition for American option
put_fd_amer = FD_Amer(P,r,sigma,deltaT)   # implicit finite difference method for American option
FD_Amer_grid = cbind(St[(M+1):1], round(put_fd_amer[(M+1):1,],3))   # concatenate St column with FD mesh
print(FD_Amer_grid)   # display FD mesh
plot(FD_Amer_grid[,1], FD_Amer_grid[,2], type="l", xlab="Risky Asset Price", ylab="Put Option Value", 
      main="Put Option Value (Implicit FD)", ylim=0:1)    # plot graph
points(FD_Euro_grid[,1], FD_Euro_grid[,2], col="blue")   # add plot of the European put option price from Q2
lines(FD_Amer_grid[,1], FD_Amer_grid[,32], col="red", lty=2)   # add plot of the intrinsic value of American put option
legend(2.5, 0.85, legend=c("American", "European", "American (Intrinsic)"), 
       col=c("black","blue","red"), lty=c(1,NA,2), pch=c(NA,1,NA), cex=1.2)



