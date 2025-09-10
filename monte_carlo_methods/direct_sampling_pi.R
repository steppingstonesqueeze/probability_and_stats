# direct sampling for Pi

library(ggplot2)

# circle centered at 0, radius 0.5; inscribed inside square of side 1 and sides
# from top left ot bottom right clockwise : 

# the square coordinates are from top left to bottom right clockwise: (-1/2,1/2), (1/2,1/2)
# (1/2, -1/2) and (-1/2, -1/2)

# Algo : sample (x,y) both uar from [-0.5,0,5] ; if distance of (x,y) from (center) <= 0.25
# its a hit

boundary <- 0.5

# check : at any step - abs(x) <= 1/2 and abs(y) <= 1/2 else you are outside square

num_expts <- 10000

N <- 4000

approx_pi <- numeric(length = num_expts)

for (expts in 1:num_expts) {
 
  N_hits <- 0
  
  for (t in 1:N) {
    x <- runif(1, -boundary, boundary)
    y <- runif(1, -boundary, boundary)
    
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
) + geom_histogram(alpha = 0.3, fill = "gray")


# superimpose a density of the histogram

ggplot(df, aes(x = api)) +
  geom_histogram(aes(y = ..density..), 
                 bins = 30, 
                 fill = "green", 
                 color = "black", 
                 alpha = 0.6) +
  geom_density(color = "yellow", linewidth = 1.2) +
  theme_minimal() + ggtitle("Direct sampling to compute Pi - 10000 experiments and each trial 4000 steps")





