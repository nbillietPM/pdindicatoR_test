library(pdindicatoR)
library(ape)
library(dplyr)

#Read in the older phylogenetic tree
angiospermTree2012 <- read.tree("data/phylogeneticAngiosperm/DaPhnE_01/DaPhnE_01.tre")
#angiospermTree2012$tip.label <- lapply(angiospermTree2012$tip.label, function(label) {gsub("_", " ", label)})

#Perform the taxonmatch function one time and load in the result
if (file.exists("data/taxonMatch/AST2012.csv")==FALSE){
  matchedAST2012 <- taxonmatch(angiospermTree2012)
  write.csv(matchedAST2012, "data/taxonMatch/AST2012.csv")
}
else {
  matchedAST2012 <- read.csv("data/taxonMatch/AST2012.csv", header=TRUE, sep=",")
}

matched_nona <- matchedAST2012 %>%
  dplyr::filter(!is.na(gbif_id))
head(matched_nona)

dataCube <- read.csv("data/cubes/nativeAngioCube.csv", header=TRUE, sep = "\t")

mcube <- append_ott_id(angiospermTree2012,dataCube, matched_nona)
head(mcube)

check_completeness(mcube)

mcube <- mcube %>%
  dplyr::filter(!is.na(ott_id))

head(mcube)

mcube$species %in% angiospermTree2012$tip.label
"
> mcube$species %in% angiospermTree2012$tip.label
   [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
  [20] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
  [39] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
"
head(mcube$species)
"
> head(mcube$species)
[1] "Carlina vulgaris"     "Senecio inaequidens"  "Lepidium didymum"     "Sambucus nigra"       "Senecio inaequidens" 
[6] "Claytonia perfoliata"
"
head(angiospermTree2012$tip.label)
"
> head(angiospermTree2012$tip.label)
[1] "Huperzia_selago"           "Lycopodiella_inundata"     "Lycopodium_clavatum"       "Lycopodium_annotinum"     
[5] "Diphasiastrum_zeilleri"    "Diphasiastrum_tristachyum"
"

#Substitution of the _ character in the tip labels needs to happen
angiospermTree2012$tip.label <- lapply(angiospermTree2012$tip.label, function(label) {gsub(" ", "_", label)})
mcube$species %in% angiospermTree2012$tip.label
"
> mcube$species %in% angiospermTree2012$tip.label
   [1]  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE
  [20]  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
  [39] FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
  [58]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
"
mcube <- mcube[mcube$species %in% angiospermTree2012$tip.label,]

PD_cube <- get_pd_cube(mcube, angiospermTree2012, metric = "faith")


PDcubeToCSV <- function(PDcube, filename){
  # Identify the columns that are stored as lists in the dataframe
  listColumns <- sapply(PDcube, is.list)
  
  # Convert list columns to strings by collapsing the list elements
  PDcube[listColumns] <- lapply(PDcube[listColumns], function(col) {
    sapply(col, function(cell) paste(cell, collapse = ","))
  })
  
  # Write the entire dataframe to CSV, including the modified list columns
  write.csv(PDcube, paste("data/PDcubes/", filename), row.names = FALSE)
}
PDcubeToCSV(PD_cube, "nativePDcubeAST2012.csv")

"
csvToPDcube <- function(filename){
  cubeDF <- read.csv(paste("data/PDcubes/", filename))
  
}

Write a function that converts the csv format of the PD cube back to the useable format
"



"
aggr_cube <- aggregate_cube(mcube, NULL)

aggr_cube
mcube
angiospermTree2012$tip.label
all_matched_sp <- unique(mcube[["orig_tiplabel"]])
all_matched_sp

MRCA <- getMRCA(angiospermTree2012, all_matched_sp)
MRCA

> MRCA <- getMRCA(angiospermTree2012, all_matched_sp)
Error in if (nd == rootnd) break : missing value where TRUE/FALSE needed

missing_species <- setdiff(all_matched_sp, angiospermTree2012$tip.label)
missing_species

Note
  there appears to be some mismatching happening between the tree and the cube
  getMRCA throws an error when comparing the original tiplabel to the tiplabels in the tree
  This gets resolved when the original regular expression substitution of the '_' character 
  in the tip tree labels with a ' ' gets reverted.
"