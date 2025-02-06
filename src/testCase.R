library(pdindicatoR)
library(ape)

angiospermTree2012 <- read.tree("data/phylogeneticAngiosperm/DaPhnE_01/DaPhnE_01.tre")
angiospermTree2012$tip.label <- lapply(angiospermTree2012$tip.label, function(label) {gsub("_", " ", label)})

angiospermTree2020 <- read.nexus("data/phylogeneticAngiosperm/oo_330891/oo_330891.tre")
angiospermTree2020$tip.label <- lapply(angiospermTree2020$tip.label, function(label) {gsub("_", " ", label)})
#Match all the names in the phylogenetic trees
# AST - AngioSpermTree
startTime <- Sys.time()
matchedAST2012 <- taxonmatch(angiospermTree2012)
endTime <- Sys.time()
timetaken <- endTime - startTime
timetaken
"
> timetaken
Time difference of 15.72772 mins
"
matchedAST2012
"
> matchedAST2012 <- taxonmatch(angiospermTree2012)
Error in .tnrs_match_names(names = names, context_name = context_name,  : 
  Argument ‘names’ must be of class ‘character
"
startTime <- Sys.time()
matchedAST2020 <- taxonmatch(angiospermTree2020)
endTime <- Sys.time()
timetaken <- endTime - startTime
timetaken


testCube <- read.csv("dataScripts/0011344-250127130748423.csv", header=TRUE, sep = "\t")
testCube
