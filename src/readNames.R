library(dplyr)
acceptedNames <- read.csv("phylogeneticAngiosperm/acceptedNames.csv", header=TRUE, encoding = "UTF-8")

acceptedNames$fullNames <- paste(acceptedNames$Family.to.Genus.Name, 
                                 acceptedNames$Species.name,
                                 sep=" ")

#Certain characters are encoded in such a way (the /xd7 character) provides issues
#substitute the multiple empty spaces with a single space
acceptedNames$fullNames<-lapply(acceptedNames$fullNames, function(text) {
                         text <- iconv(text, from = "windows-1252", to = "UTF-8") 
                         text <- gsub("NA", "", text)
                         gsub("\\s+", " ", text)})

#Filter out the species name that consist of a genus+species
speciesNames <- acceptedNames %>% filter(code.NF6 =="key") %>% pull(fullNames)

invasiveNames <- read.csv("phylogeneticAngiosperm/alienPlantsBelgium.csv", header=TRUE, encoding = "UTF-8")$Name..excl..author.citation.

length(speciesNames[invasiveNames %in% speciesNames])
"
> length(speciesNames[invasiveNames %in% speciesNames])
[1] 610

In this list 610 of the species names that are present are considered invasive according to
the TRIAS list. We will filter out these names and store the rest as considered native species
"

nativeSpecies <- speciesNames[!(invasiveNames %in% speciesNames)]

