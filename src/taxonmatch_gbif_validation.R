library(rgbif)
library(rotl)
library(ape)
library(tidyr)
library(dplyr)

#Load the tree from the data directory
angiospermTree2020 <- read.nexus("data/phylogeneticAngiosperm/oo_330891/oo_330891.tre")

angiospermTree2020$tip.label<-gsub("WITHDASH", "-", angiospermTree2020$tip.label)

tipMatches <- name_backbone_checklist(angiospermTree2020$tip.label)

str(tipMatches[1,])

exactMatches<-dplyr::filter(tipMatches, matchType=="EXACT")
fuzzyMatches<-dplyr::filter(tipMatches, matchType=="FUZZY")
higherrankMatches<-dplyr::filter(tipMatches, matchType=="HIGHERRANK")
noneMatches<-dplyr::filter(tipMatches, matchType=="NONE")

exactIdx<-which(tipMatches$matchType=="EXACT")
fuzzyIdx<-which(tipMatches$matchType=="FUZZY")
higherrankIdx<-which(tipMatches$matchType=="HIGHERRANK")
noneIdx<-which(tipMatches$matchType=="NONE")
  
nmbExact<-nrow(exactIdx)
nmbFuzzy<-nrow(fuzzyIdx)
nmbHigherRank<-nrow(higherrankIdx)
nmbNone<-nrow(noneIdx)
nmbTips<-nrow(tipMatches)

print(paste("The number of exact matches:", nmbExact, "(", (nmbExact / nmbTips) * 100, "%)"))
print(paste("The number of fuzzy matches:", nmbFuzzy, "(", (nmbFuzzy / nmbTips) * 100, "%)"))
print(paste("The number of higher rank matches:", nmbHigherRank, "(", (nmbHigherRank / nmbTips) * 100, "%)"))
print(paste("The number of none matches:", nmbNone, "(", (nmbNone / nmbTips) * 100, "%)"))

fuzzyIdx

fuzzyMatches$species[1]
angiospermTree2020$tip.label[fuzzyIdx[1]]
higherrankMatches$species[1]
angiospermTree2020$tip.label[higherrankIdx[1]]
str(higherrankMatches[1,])


str(name_backbone(name="Tritonia pallida"))
# Retrieve all matches for the name
lookup_results <- name_lookup(query = "Tritonia pallida", rank = "species")

# View results
str(lookup_results$data$kingdom)
str(lookup_results$data[5,])
