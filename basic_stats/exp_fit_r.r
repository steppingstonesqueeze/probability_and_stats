# Read the data
data <- read.csv("/Users/gn/Downloads/histo.csv", header = TRUE)

# Create x values (index positions)
x <- 1:nrow(data)
y <- data[, 1]  # First column values

# Display the data
print("Data points:")
print(data.frame(x = x, y = y))

# Fit exponential model: y = a * exp(b * x)
# We'll use nonlinear least squares (nls)
# Starting values: estimate a from first point, b from decay rate
start_a <- y[1]
start_b <- log(y[2]/y[1])

# Fit the model
exp_model <- nls(y ~ a * exp(b * x), 
                 start = list(a = start_a, b = start_b))

# Display model summary
print("\nExponential Model: y = a * exp(b * x)")
print(summary(exp_model))

# Get coefficients
coef_a <- coef(exp_model)["a"]
coef_b <- coef(exp_model)["b"]
cat("\nFitted equation: y =", round(coef_a, 4), "* exp(", round(coef_b, 4), "* x)\n")

# Calculate R-squared
y_pred <- predict(exp_model)
ss_res <- sum((y - y_pred)^2)
ss_tot <- sum((y - mean(y))^2)
r_squared <- 1 - (ss_res / ss_tot)
cat("R-squared:", round(r_squared, 4), "\n")

# Plot the data and fitted curve
plot(x, y, pch = 19, col = "blue", cex = 1.5,
     main = "Exponential Curve Fitting",
     xlab = "Index", ylab = "Value",
     ylim = c(0, max(y) * 1.1))

# Add fitted curve
x_smooth <- seq(min(x), max(x), length.out = 100)
y_smooth <- coef_a * exp(coef_b * x_smooth)
lines(x_smooth, y_smooth, col = "red", lwd = 2)

# Add legend
legend("topright", 
       legend = c("Data", "Fitted curve"), 
       col = c("blue", "red"), 
       pch = c(19, NA), 
       lty = c(NA, 1),
       lwd = c(NA, 2))

# Add equation to plot
equation_text <- paste0("y = ", round(coef_a, 2), " * exp(", round(coef_b, 3), " * x)")
r2_text <- paste0("R² = ", round(r_squared, 4))
text(x = max(x) * 0.6, y = max(y) * 0.9, 
     labels = paste(equation_text, "\n", r2_text), 
     pos = 4, cex = 0.9)
