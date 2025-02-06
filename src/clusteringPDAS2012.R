install.packages("mclust")
library(mclust)

set.seed(42)
gmm_model <- Mclust(PD_cube$PD, G=1:15)
plot(gmm_model, what="BIC")
#Visual inspection of the BIC value indicates that at 7 components we reach a turning point in the BIC

gmm_model <- Mclust(PD_cube$PD, G=7)
partitionedPDvalues <- split(PD_cube$PD, gmm_model$classification)


par(mfrow = c(1, 3), mar = c(4, 4, 1, 1))



hist(PD_cube$PD, breaks = 50, col = "lightgray", probability = TRUE,
     main = "Fitted Gaussian Mixture Model on PD value (AST2012)",
     xlab = "PD values")

for (key in 1:7){
  densityCurve <- density(unlist(partitionedPDvalues[key]))
  densityScaled <- densityCurve$y *gmm_model$parameters$pro[key]
  lines(densityCurve$x, densityScaled, lwd = 2)
}

lines(density(PD_cube$PD), col = "blue", lwd = 2)


