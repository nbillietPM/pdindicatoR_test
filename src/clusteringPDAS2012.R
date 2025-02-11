install.packages("mclust")
install.packages("viridis")
library(mclust)
library(viridis)
set.seed(42)
gmm_model <- Mclust(pdCube$PD, G=1:15)
plot(gmm_model, what="BIC")
#Visual inspection of the BIC value indicates that at 7 components we reach a turning point in the BIC

gmm_model <- Mclust(pdCube$PD, G=7)
partitionedPDvalues <- split(pdCube$PD, gmm_model$classification)
partitionedPDvalues

par(mfrow = c(1, 3), mar = c(4, 4, 1, 1))


pdGeoCube$gmmClass <- gmm_model$classification

mapPD <- ggplot() +
  geom_sf(data = pdGeoCube, mapping=aes(fill = factor(.data$gmmClass))) +  # Ensure gmmClass is treated as a factor
  scale_fill_manual(name = "PD cluster", 
                    values = viridis::viridis(7, option = "magma"),  # Assign distinct colors from the 'magma' palette
                    labels = sapply(1:7, function(i) paste("Cluster", i))) +  # Labels for clusters
  geom_sf(data = belgiumBorder, fill = NA, color = "white", linewidth = 0.5)
mapPD


hist(pdCube$PD, breaks = 50, col = "lightgray", probability = TRUE,
     main = "Fitted Gaussian Mixture Model on PD value (AST2012)",
     xlab = "PD values")

for (key in 1:7){
  densityCurve <- density(unlist(partitionedPDvalues[key]))
  densityScaled <- densityCurve$y *gmm_model$parameters$pro[key]
  lines(densityCurve$x, densityScaled, lwd = 2)
}

lines(density(PD_cube$PD), col = "blue", lwd = 2)


