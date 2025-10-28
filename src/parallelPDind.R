library(parallel)
library(rgbif)
library(ape)
library(pdindicatoR)
library(dplyr)

ex_data <- retrieve_example_data()
ex_tree <- ex_data$tree
cube <- ex_data$cube
grid <- ex_data$grid

ex_matched <- taxonmatch(ex_tree)
matched_nona <- ex_matched %>%
  dplyr::filter(!is.na(gbif_id))
mcube <- append_ott_id(ex_tree, cube, matched_nona)
mcube <- mcube %>%
  dplyr::filter(!is.na(ott_id))
aggr_cube <- aggregate_cube(mcube)


#Taxonmatch rgbif edition

matchStatistics <- function(matches){
  head(matches)
  exactIdx<-which(matches$matchType=="EXACT")
  fuzzyIdx<-which(matches$matchType=="FUZZY")
  higherrankIdx<-which(matches$matchType=="HIGHERRANK")
  noneIdx<-which(matches$matchType=="NONE")
  
  nmbExact<-length(exactIdx)
  nmbFuzzy<-length(fuzzyIdx)
  nmbHigherRank<-length(higherrankIdx)
  nmbNone<-length(noneIdx)
  nmbTips<-nrow(matches)
  
  print(paste("The number of exact matches:", nmbExact, "(", (nmbExact / nmbTips) * 100, "%)"))
  print(paste("The number of fuzzy matches:", nmbFuzzy, "(", (nmbFuzzy / nmbTips) * 100, "%)"))
  print(paste("The number of higher rank matches:", nmbHigherRank, "(", (nmbHigherRank / nmbTips) * 100, "%)"))
  print(paste("The number of none matches:", nmbNone, "(", (nmbNone / nmbTips) * 100, "%)"))
}

taxonmatchGBIF <- function(phyloTree, nmbCores=1){
  if(nmbCores>detectCores()-1){
        cat("The number of requested cores exceeds the number of available cores. The maximum number of cores is ",detectCores()-1)
        return(invisible(NULL))
  }
  if(nmbCores>1){
    #Generate indices to split the list so that we can process the fragments
    chunkIndices <- cut(seq_along(phyloTree$tip.label), breaks = nmbCores, labels = FALSE)
    #Split up the data based on the indices
    splitList <- split(phyloTree$tip.label, chunkIndices)
    #Depending on the OS different parallelisation is employed
    if (Sys.info()['sysname']=="Windows"){
      #setup computation cluster using the specified number of cores
      computeCluster<-makeCluster(nmbCores, type = "PSOCK")
      #Apply the name_backbone_checklist to each sublist
      tipMatches <- parLapply(computeCluster,
                              splitList,
                              name_backbone_checklist)
      #Combine the separate sublists into a dataframe
      tipMatches<-do.call(bind_rows, tipMatches)
      return(tipMatches)
    } else {
      #setup computation cluster using the specified number of cores
      computeCluster<-makeCluster(nmbCores, type = "FORK")
      #Apply the name_backbone_checklist to each sublist
      tipMatches <- parLapply(computeCluster,
                              splitList,
                              name_backbone_checklist)
      #Combine the separate sublists into a dataframe
      tipMatches<-do.call(bind_rows, tipMatches)
      return(tipMatches)
    }
  } else {
    print("This is the else statement")
    tipMatches <- name_backbone_checklist(phyloTree$tip.label)
    return(tipMatches)
  }
}

#Single core test
startTime <- Sys.time()
results<-taxonmatchGBIF(ex_tree)
endTime <- Sys.time()
endTime-startTime
matchStatistics(results)
results
#Multi core test
startTime <- Sys.time()
results<-taxonmatchGBIF(ex_tree, nmbCores = 4)
endTime <- Sys.time()
endTime-startTime
matchStatistics(results)
results


#Load the tree from the data directory
angiospermTree2020 <- read.nexus("data/phylogeneticAngiosperm/oo_330891/oo_330891.tre")
angiospermTree2020$tip.label<-gsub("WITHDASH", "-", angiospermTree2020$tip.label)
angiospermTree2020$tip.label<-gsub("_", " ", angiospermTree2020$tip.label)

start_time <- Sys.time()
taxonmatch(angiospermTree2020)
end_time <- Sys.time()
end_time-start_time

start_time <- Sys.time()
tipMatches<-taxonmatchGBIF(angiospermTree2020, nmbCores = 8)
end_time <- Sys.time()
end_time-start_time
matchStatistics(tipMatches)
head(tipMatches)
taxonmatch("test",parallel = FALSE)


speciesList <- unique(cube$species)
#check non empty string
speciesList<- speciesList[nzchar(speciesList)]
speciesList

calculateFaithPD <- function(phyloTree, MRCA, speciesList) {
  tipIDs <- match(speciesList, phyloTree$tip.label)
  tipIDs <- tipIDs[!is.na(tipIDs)]  # Remove NAs
  
  # Compute node paths
  nodepath <- lapply(tipIDs, function(id) ape::nodepath(phyloTree, MRCA, id))
  
  # Compute edge IDs
  edge_ids <- lapply(nodepath, function(path) {
    which(phyloTree$edge[, 1] %in% path[-length(path)] & phyloTree$edge[, 2] %in% path[-1])
  })
  
  # unique edge extraction
  edge_ids_unique <- Reduce(union, edge_ids)  # Faster than unique(unlist(...))
  
  # Sum branch lengths
  sum(phyloTree$edge.length[edge_ids_unique], na.rm = TRUE)  # Handle possible NAs efficiently
}


start_time <- Sys.time()
calculateFaithPD(angiospermTree2020, MRCA, speciesList)
end_time <- Sys.time()
end_time-start_time

start_time <- Sys.time()
test<-calculate_faithpd(angiospermTree2020, speciesList, MRCA)
end_time <- Sys.time()
end_time-start_time
test
