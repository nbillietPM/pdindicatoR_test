library(pdindicatoR)
library(ape)
library(dplyr)
library(rotl)

angiospermTree2020 <- read.nexus("data/phylogeneticAngiosperm/oo_330891/oo_330891.tre")
angiospermTree2020$tip.label<-gsub("WITHDASH", "-", angiospermTree2020$tip.label)

if (file.exists("data/taxonMatch/AST2020.csv")==FALSE){
  matchedAST2020 <- taxonmatchLagged(angiospermTree2020,10)
  write.csv(matchedAST2020, "data/taxonMatch/AST2020.csv")
}
else {
  matchedAST2020 <- read.csv("data/taxonMatch/AST2020.csv", header=TRUE, sep=",")
}

matchedAST2020 <- matchedAST2020 %>% rename(gbif_id = gbifID)

matched_nona <- matchedAST2020 %>% dplyr::filter(!is.na(gbif_id))
head(matched_nona)

dataCube <- read.csv("data/cubes/nativeAngioCube.csv", header=TRUE, sep = "\t")

head(dataCube)

colnames(matched_nona)

matched_nona[, c("ott_id", "gbif_id", "unique_name", "orig_tiplabel")]
head(matched_nona)
mcube <- append_ott_id(angiospermTree2020,dataCube, matched_nona)

check_completeness(mcube)

mcube <- mcube %>% dplyr::filter(!is.na(ott_id))