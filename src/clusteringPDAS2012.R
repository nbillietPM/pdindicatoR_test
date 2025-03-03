install.packages("mclust")
install.packages("viridis")
library(mclust)
library(viridis)
set.seed(42)
gmm_model <- Mclust(pdCube$PD, G=1:15)
plot(gmm_model, what="BIC")
#Visual inspection of the BIC value indicates that at 7 components we reach a turning point in the BIC

numClusters <- 7
gmm_model <- Mclust(pdCube$PD, G=numClusters)

pdGeoCube$gmmClass <- gmm_model$classification

means <- gmm_model$parameters$mean
deviation <- sqrt(gmm_model$parameters$variance$sigma)

cluster_labels <- sapply(1:numClusters, function(i) {
  paste("Cluster", i, "(", round(means[i], 2), ",", round(deviation[i], 2),")")
})

mapPD <- ggplot() +
  geom_sf(data = pdGeoCube, mapping=aes(fill = factor(.data$gmmClass))) +  # Ensure gmmClass is treated as a factor
  scale_fill_manual(name = "PD cluster", 
                    values = viridis::viridis(numClusters, option = "turbo"),  # Assign distinct colors from the 'magma' palette
                    labels = cluster_labels) +
  geom_sf(data = belgiumBorder, fill = NA, color = "white", linewidth = 1)
mapPD

png('img/gmmAST2012Belgie.png',type='cairo', units = 'px', width = 6000, height = 4000, res = 600)
plot(mapPD)
dev.off()

hist(pdCube$PD, breaks = 50, col = "lightgray", probability = TRUE,
     main = "Fitted Gaussian Mixture Model on PD value (AST2012)",
     xlab = "PD values")

for (key in 1:numClusters){
  densityCurve <- density(unlist(partitionedPDvalues[key]))
  densityScaled <- densityCurve$y *gmm_model$parameters$pro[key]
  lines(densityCurve$x, densityScaled, lwd = 2)
}



