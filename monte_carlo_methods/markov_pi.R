# Krauth - Markov pi -- the adult sim on pi calculation on a helipad

# Square size 1 * 1
# Circle inscribed in it radius 1/2 => Area of circle = pi / 4
# Random darts thrown = N ; if N_d in circle then pi / 4 = N_d / N and we can 
# approximate pi

# For this version - the "darts" are replaced by random steps from wherever we are at time t
# If the step takes us outside the square, redo step. If not, move to new location and 
# increment N_d by 1 if qithin circle, else we are in square but not circle

library(ggplot2)
library(tidyr)
library(dplyr)

# start at center (0,0) of the square

# the square coordinates are from top left to bottom right clockwise: (-1/2,1/2), (1/2,1/2)
# (1/2, -1/2) and (-1/2, -1/2)

boundary_x <- 0.5
boundary_y <- 0.5

# check : at any step - abs(x) <= 1/2 and abs(y) <= 1/2 else you are outside square

num_expts <- 10000

N <- 4000
step_size <- 0.3 # random sample both directions in [-step_size, step_size]
approx_pi <- numeric(length = num_expts)

for (expts in 1:num_expts) {
  cat(expts, "\n")
  start_x <- 0.5
  start_y <- 0.5
  
  x <- start_x
  y <- start_y
  
  N_hits <- 0
  
  for (t in 1:N) {
    delta_x <- runif(1, -step_size, step_size)
    delta_y <- runif(1, -step_size, step_size)
    
    if ((abs(x + delta_x) <= boundary_x) && (abs(y + delta_y) <= boundary_y)) {
      x <- x + delta_x
      y <- y + delta_y
      
    }
    
    # inside circle? 
    if (x*x + y*y <= 0.25) { # circle of radius 0.5 centered at (0,0)
      N_hits <- N_hits + 1
    }
    
  }
  value <- N_hits / N
  approx_pi[expts] <- 4 * value
}

# plot the histogram of this

df <- data.frame(
  api = approx_pi
)

# the histo
ggplot(
  data = df, aes(x = api)
) + geom_histogram(alpha = 0.3, fill = "red")


# superimpose a density of the histogram

ggplot(df, aes(x = api)) +
  geom_histogram(aes(y = ..density..), 
                 bins = 30, 
                 fill = "skyblue", 
                 color = "black", 
                 alpha = 0.6) +
  geom_density(color = "red", linewidth = 1.2) +
  theme_minimal() + ggtitle("Markov sims to compute Pi - 10000 experiments and each trial 4000 steps")





