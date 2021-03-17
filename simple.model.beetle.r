# Overwintering Beetle Model 
# Andre Szejner Sigal
# March 2021

library('plot.matrix')
library('ggplot2')
library('plot3D')

#---- Beetle Parameters ----

x_max  <- 25    # Maximum fat reserves
x_d    <- seq(from = 0, to = x_max, length.out = 100)   #Discretized values of x

#---- Decisions ----

# Beetles either saves energy or invests in cold (# actions = 2)
# Each decision has the same global risk of mortality denoted by mu. (Just a small general risk)
# Each decision has a different risk of cold mortality denoted by mod, depending if it was a cold day or not
# Each decision also has different energy investment denoted by e 
# Decision 1 = save energy, Decision 2 = invest in cold

mu <- c(0.01, 0.01)  # Daily global risk
mod <- c(0.5, 0)     # Mortality modifier for each decision added only if its a bad day
e  <- c(0, 1)        # Cold tolerance investment per day


#---- Winter Parameters ---- 

# Winter consists of 10 Days and each day is divided into 2 time steps (Time). 
# Probability of a cold day is p_c
# Any night has a cost (c_b) and (c_g) 
# The probability of a cold night is given by p_b, and adds mortality if beetles didnt invest in cold

Time <- 2    # Number of time periods per day
Days <- 10   # Number of days in winter
c_g  <- 1    # Cost of a winter day 
c_b  <- 1    # Cost of a winter day (these are legacy code for cost of a bad day and a good day)
p_c  <- 0.1  # Probability of a cold night


# NEW STEP: Assign an ID to Day type, Bad_Day = 1, Good_Day =0
# The goal was to apply high mortality only to the group who chose not to invest, AND if it was a cold night
set.seed(2)
d.type = rbinom(Days, 1, prob = p_c)


##########################################################
# Model starts

# Terminal fitness is the probability of surviving the winter
# First, the empty fitness array
Fit<- array(data = NA, dim = c(length(x_d), Time+1, Days), 
            dimnames = list(x_d, 1:(Time+1), 1:(Days)) )

# Calculate terminal fitness F(x, t+1, d) where x=state 
for (j in 1:length(x_d)) {
  # Note: I use j as the index of an x-value in x_d and x 
  # as the value of x_d[j].
  x <- x_d[j]
  # If your state does not allow you to survive a good night fitness is zero. 
  if (x < c_g) {
    Fit[j, Time+1, Days] <- 0
    # If your state allows you to survive a good night, but not a cold one,
    # you survive with probability p_g = 1 - p_c (the probability of a good night). 
  } else if (x < c_b) {
    Fit[j, Time+1, Days] <- (1-p_c)
    # If your state allows you to survive a cold night, you will always survive. 
  } else {
    Fit[j, Time+1, Days] <- 1
  } # end if-else loop
} # end for loop

as.matrix(Fit[, Time+1, Days]) # Test if it works
# The left column shows state, the right column shows the terminal fitness we just calculated.

#---- Discrete state variable ----
# Interpolate the discrete state variable 
interpolate <- function (x, t, d) {
  
  # Returns the index of the closest discrete x value 
  closest_discrete_x <- function(x) {
    # Note: only the first value is returned if there are multiple 
    # equidistant values.
    return(which(abs(x_d - x) == min(abs(x_d - x)))[1])
  }
  
  # Interpolate between two values a and b. dx is a value between 0 and 1.
  linear_interpolation <- function (a, b, dx) {
    return((1-dx)*a+ b*dx)
  }
  
  # No point in doing interpolation if energy reserves are negative.
  # Bird is dead.
  if (x < 0) { return(0) }
  
  # Figure out between which two discretized values of x our x lies.
  closest <- closest_discrete_x(x)
  if (x < x_d[closest]) {
    j1 <- closest -1
    j2 <- closest
  } else if (x > x_d[closest]){
    j1 <- closest
    j2 <- closest +1
  } 
  # Fitness value for x is already present in the matrix. No need to interpolate.
  else { return( Fit[closest, t, d]) }
  
  # Calculate how x is positioned in relation to x_d[j1] and x_d[j2].
  # 0: Closer to x_d[j1]
  # 1: Closer to x_d[j2]
  delta_x <- (x-x_d[j1])/(x_d[j2]-x_d[j1])
  
  # Interpolate.
  return(linear_interpolation(Fit[j1, t, d], Fit[j2, t, d], delta_x))
}


##########################################################
# ---- Backwards Iteration ---- 

# Array for the decision loop (decisions, H)
H <- array(data = NA, dim = c(length(x_d), Time, Days), dimnames = list(
  x_d, 1:Time, 1:Days))

# Iterate backwards across days at Time+1
for (d in Days:1) {
  # Terminal fitness has already been calculated, 
  # so don't calculate fitness for T+1 on the last day.
  if (d != Days) {
    for (j in 1:length(x_d)) {
      x <- x_d[j]
      
      # Interpolation function used here
      # If a beetle can't survive a good night:
      if (x < c_g) {
        Fit[j, Time+1, d] <- 0
        
        # If a beetle can survive a good night but not a cold night:
        # Fitness = prob. of a good night * (state minus the overnight cost)
      } else if (x < c_b) {
        Fit[j, Time+1, d] <- (1-p_c)*interpolate(x-c_g, 1, d+1)  
        
        # Otherwise is can survive a cold night:
        # Fitness = prob. of a good night * (state minus overnight cost) + prob. of a cold night * (state minus overnight cost)
        # Uses the logic as the HK model
      } else {
        Fit[j, Time+1, d] <- (1-p_c)*interpolate(x-c_g, 1, d+1) + p_c*interpolate(x-c_b, 1, d+1)
      }
    } # end if-else loop
  } # end j loop
  
  # Iterate backwards from max Time for each day
  for (t in Time:1) {
    for (j in 1:length(x_d)) {
      x <- x_d[j]
      # Vector for storing fitness values for each of the decisions
      F_i <- vector(mode = 'numeric', length=2)
      
      # Calculate resulting fitness of choosing each patch
      for (i in 1:2) {
        # Calculating the expected state in the future.
        x_mark <-  x + ((1/Time) * -e[i])   #new state = old state - cold investment
        
        # Calculate the expected fitness given the decision h.
        # BUT FIRST: check what kind of day it was/is/will be
        if (d.type[d] == 1) {  # if it's a cold day, include the mod in mortality mu
          F_i[i] <- (1 - (1/Time)*(mu[i] + mod[i])) * interpolate(x_mark, t+1, d)
        } else {
          F_i[i] <- (1 - (1/Time)*(mu[i])) * interpolate(x_mark, t+1, d)
          
          
        } #end the if day.type loop
      } # end i loop 
      
      # Which decision choice maximizes fitness?
      Fit[j, t, d] <- max(F_i)[1]
      
      # Optimal decision choice is the one that maximizes fitness.
      # In cases where more than one patch shares the same fitness, 
      # the first one (i.e. lower risk) is chosen.
      H[j, t, d] <- which(F_i == max(F_i))[1]
      
    } # end j loop 
  } # end t loop
} # end d loop


##########################################################
#--- Plots ----

# Fitness function at Days-9, Days-5, Days-1
# This figure should be equivalent to 5.1 in Clark and Mangel (1999).

# Create an empty plot for plotting lines.
plot(NA, type="n", 
     xlab="Fat reserves (g)",
     ylab="Fitness, F[x, t]", 
     xlim=c(0, x_max), ylim=c(0, 1))

# Plotting the lines.
lines(x_d, as.vector(Fit[,2,Days-9]), col = "black", lty = 1)
lines(x_d, as.vector(Fit[,2,Days-5]), col = "blue", lty = 1,  lwd = 2)
lines(x_d, as.vector(Fit[,2,Days-1]), col = "red", lty = 1,  lwd = 2)

# Reverse the x-dimension of the array so it is ordered from high to low. 
# Nicer when plotting.
Fit.rev <- Fit[length(x_d):1,,]
H.rev <-   H  [length(x_d):1,,]

# Plotting the optimal decision at any given time.
# black:  Patch 1 (0)  energy saving
# yellow: Patch 2 (1)  cold tolerant
# Should be equivalent to Figure 5.2 in Clark and Mangel (1999).
plot(H.rev[,,Days-9], breaks=c(0.5, 1.5, 2.5), col=c("black", "yellow"),
     xlab = "Time of day",
     ylab = "Fat reserves",
     main = "mix energy/cold Day 1")

plot(H.rev[,,Days-5], breaks=c(0.5, 1.5, 2.5), col=c("black", "yellow"),
     xlab = "Time of day",
     ylab = "Fat reserves",
     main = "mix energy/cold Day 5")

plot(H.rev[,,Days-1], breaks=c(0.5, 1.5, 2.5), col=c("black", "yellow"),
     xlab = "Time of day",
     ylab = "Fat reserves",
     main = "mix energy/cold Day 5")

# Fitness 
plot(Fit.rev[,,Days-9],
     xlab = "Time of day",
     ylab = "Fat reserves")

# Fitness landscape plot.
persp3D(z = Fit.rev[,,Days-9], theta = 135, phi = 45,
        xlab = "State (x)", 
        ylab = "Time (t)",
        zlab = "Fitness (F)")



##########################################################
#---- Forward Simulation ---- 

# Number of days to forward iterate for.
Days  <- 10

# Time/day is given from the model above.
# Initial states for individual (A vector of indexes in x_d).
j_0 <- c(100,75, 55,10)

# Number of individuals to iterate for each x_0. Total number of individuals
# is length(x_0)*N.ind:
N.ind <- 100

# X[j_0, N, D, T]:
# State of individual N with initial state j_0, at day D and time T.
X <- array(data = NA, dim = c(length(j_0), N.ind, Days*Time +1), dimnames = list(
  x_d[j_0],
  1:N.ind,
  1:(Days*Time+1)
))

for (j in 1:length(j_0)) {
  x_0 <- x_d[j_0[j]]
  for (n in 1:N.ind) {
    # Set the initial state.
    X[j, n, 1] <- x_0
    
    # Iterate through the time.
    for (z in 1:(Days*Time) ) {
      t = (z-1) %%  Time  +1 # Time of day
      d = (z-1) %/% Time  +1 # Day
      # The current state.
      x <- X[j, n, z]
      if (x < 0) {
        # Beetle is dead. It will remain dead.
        x_new <- x
      } else {
        # The best decision given the current state (Found by rounding x to the closest item in x_d).
        # DOIT interpolate the decision between the two closest optimal decisions.
        h <- H[which(abs(x_d - x) == min(abs(x_d - x)))[1], t, d]
        
        if (runif(1) <= (1/Time)*(mu[h]+mod[h])) { # Predation risk.
          # Negative x = dead.
          x_new <- -1
        } else {
          investment   <- e[h]
          
          # New state are the previous one minus investment
          x_new <- x - investment
        }
      }
      
      if (t == Time) {
        # If it is the end of the day, nightly costs need to be applied as well.
        if (runif(1) < p_c) {
          # Cold night.
          x_new <- x_new - c_b
          
        } else {
          # Good night.
          x_new <- x_new - c_g
        }
        
        # Set the state for the start of the following day.
      } 
      
      # Set new state.
      X[j, n, z+1] <- x_new
    }
  }
}


## The below plots multiple individuals on the same plot.
# Red vertical lines:     Separates days.
# Light blue solid line:  Fat reserves necessary for surviving a cold night.
# Light blue dotted line: Fat reserves necessary for surviving a good night.

# Create an empty plot, with the necessary xlim and ylim.
plot(NA, type="n", 
     xlab="Time",
     ylab="Fat reserves (g)", xlim=c(1, (Days*Time +1)), ylim=c(0, x_max))

# Cost of good and cold nights
abline(h=c_g, col = "blue", lty = 3)
abline(h=c_g, col = "blue", lty = 1)
# Position of the nights
for (i in 1:Days) {
  abline(v=i*Time, col = "red", lty = 3)
}
# Plot each individual for the chosen x_0.
for (n in 1:N.ind) {
  lines(1:(Days*Time +1), X[1, n, ], col = "black")
}
