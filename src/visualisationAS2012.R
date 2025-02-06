hist(PD_cube$PD, breaks=50, col="lightgray", probability=TRUE, 
     main="Histogram with fitted Probability Curve AST2012",
     xlab = "PD values", ylab="Density",
     xlim = c(min(PD_cube$PD), max(PD_cube$PD)),  # Ensures x-axis covers the range of the data
     ylim = c(0, 0.00025))

density_curve <- density(PD_cube$PD)
lines(density_curve, col = "blue", lwd = 2)

# Calculate the mean and mode
mean_value <- mean(PD_cube$PD)
mode_value <- density_curve$x[which.max(density_curve$y)]

# Add lines for mean and mode
abline(v = mean_value, col = "red", lwd = 2, lty = 2)  # Mean (red, dashed line)
abline(v = mode_value, col = "green", lwd = 2, lty = 2)  # Mode (green, dashed line)

# Add a legend
legend("topright", legend = c("Density Curve", "Mean", "Mode"), 
       col = c("blue", "red", "green"), lwd = 2, lty = 2)

# Manually place the text beneath the legend
# We position the text just below the legend by adjusting the coordinates for the 'text' function
legend_x <- 0.85  # Adjust the x position relative to the plot area
legend_y <- 0.95  # Adjust the y position to be slightly below the legend

# Add text box with mean and mode values under the legend
text(legend_x, legend_y - 0.05, labels = paste("Mean = ", round(mean_value, 2), "\nMode = ", round(mode_value, 2)),
     cex = 1.2, col = "black", adj = 0, pos = 4)  # `adj = 0` aligns text to the left