library(pdindicatoR)
library(ape)
library(dplyr)

angiospermTree2020 <- read.nexus("data/phylogeneticAngiosperm/oo_330891/oo_330891.tre")

if (file.exists("data/taxonMatch/AST2020.csv")==FALSE){
  matchedAST2020 <- taxonmatch(angiospermTree2020)
  write.csv(matchedAST2012, "data/taxonMatch/AST2020.csv")
}
else {
  matchedAST2020 <- read.csv("data/taxonMatch/AST2020.csv", header=TRUE, sep=",")
}