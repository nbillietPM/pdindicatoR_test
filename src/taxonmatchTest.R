library(rgbif)
library(ape)
library(rotl)
library(stringr)

#Load the tree from the data directory
angiospermTree2020 <- read.nexus("data/phylogeneticAngiosperm/oo_330891/oo_330891.tre")

#Clean the tip labels that are formatted with 'WITHDASH' with '-' so that they become searchable
angiospermTree2020$tip.label<-gsub("WITHDASH", "-", angiospermTree2020$tip.label)
head(angiospermTree2020$tip.label)
"
> head(angiospermTree2020$tip.label)
[1] "Amborella_trichopoda"  "Trithuria_bibracteata" "Trithuria_submersa"    "Trithuria_austinensis" "Trithuria_australis"  
[6] "Trithuria_filamentosa"
"

nmbSplits=10
#Generate indices to split the list so that we can process the fragments
chunkIndices <- cut(seq_along(angiospermTree2020$tip.label), breaks = nmbSplits, labels = FALSE)
#Split up the data based on the indices
splitList <- split(angiospermTree2020$tip.label, chunkIndices)

#Perform matching on each of the sublists to prevent HTML errors for server overloading
matchedTaxa <- lapply(splitList, function(list){
  result<-tnrs_match_names(list)
  Sys.sleep(5)
  return(result)
})

#Combine the list fragments together in an overarching dataframe
matchedTaxa <- do.call(rbind, matchedTaxa)
head(matchedTaxa)
"
> head(matchedTaxa)
            search_string           unique_name approximate_match     score ott_id is_synonym flags number_matches gbif_id
1.1  amborella_trichopoda  Amborella trichopoda              TRUE 0.9500000 303950      FALSE                    1 3667902
1.2 trithuria_bibracteata Trithuria bibracteata              TRUE 0.9523810  58193      FALSE                    1 2882465
1.3    trithuria_submersa    Trithuria submersa              TRUE 0.9444444 334510      FALSE                    1 2882469
1.4 trithuria_austinensis Trithuria austinensis              TRUE 0.9523810 982723      FALSE                    1 2882466
1.5   trithuria_australis   Trithuria australis              TRUE 0.9473684 344006      FALSE                    3 2882475
1.6 trithuria_filamentosa Trithuria filamentosa              TRUE 0.9523810 856534      FALSE                    1 2882474
"

#Extract the taxonomic information
taxonInformation <- taxonomy_taxon_info(matchedTaxa[!is.na(matchedTaxa$ott_id),]$ott_id)
taxonInformation[[1]]
"
> taxonInformation[[1]]
$flags
list()

$is_suppressed
[1] FALSE

$is_suppressed_from_synth
[1] FALSE

$name
[1] "Amborella trichopoda"

$ott_id
[1] 303950

$rank
[1] "species"

$source
[1] "ott3.7draft2"

$synonyms
list()

$tax_sources
$tax_sources[[1]]
[1] "ncbi:13333"

$tax_sources[[2]]
[1] "gbif:3667902"


$unique_name
[1] "Amborella trichopoda"
"

#Match the instances where the taxon sources are listed and contain a GBIF ID
matchedIdx <- lapply(taxonInformation, function(entry) {
  if (length(entry$tax_sources)>1) {
    return(grep("gbif", entry$tax_sources))
  } else {
    return(NA)  # Handle cases where tax_sources[[2]] is missing
  }
})

#Extract all the indices of instances where the nmb of GBIF id's are more than 1
multiValueIndices <- which(sapply(matchedIdx, length) > 1)

#Extract the corresponding unique names of these instances
multiValueUniqueNames <- lapply(taxonInformation[multiValueIndices], function(entry){return(entry$unique_name)})

multiValueNames[1]
"
> multiValueNames[1]
$`399368`
[1] "Cabomba furcata"
"
#Use the unique names to fetch the GBIF ID's using the GBIF taxonomic backbone
rgbifID_mulVal <- lapply(multiValueNames, function(name){name_backbone(name)$usageKey})
rgbifID_mulVal[1]
"
> rgbifID_mulVal[1]
$`399368`
[1] 4925544
"
#Extract the taxon sources from the multi GBIF instances
multiValueTaxonSources <- lapply(taxonInformation[multiValueIndices], function(entry){return(entry$tax_sources)}
multiValueTaxonSources[1]
"> multiValueTaxonSources[1]
$`399368`
$`399368`[[1]]
[1] "ncbi:296032"

$`399368`[[2]]
[1] "gbif:9341094"

$`399368`[[3]]
[1] "gbif:4925544"

$`399368`[[4]]
[1] "irmng:10892988"
"

#Extract the taxon sources that contain multiple gbif IDs and extract the numerical ID format as int
ott_gbifID<- lapply(taxonInformation[multiValueIndices], function(info){
  lapply(grep("gbif", info$tax_sources, value = TRUE), function(gbifSource){
    as.integer(str_replace(gbifSource, "gbif:", ""))
  })
})

ott_gbifID[[1]]
"
> ott_gbifID[[1]]
[[1]]
[1] 9341094

[[2]]
[1] 4925544
"

#Compare the extracted GBIF ID's with the ID's present in the GBIF backbone


