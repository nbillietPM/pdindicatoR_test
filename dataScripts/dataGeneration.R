#Data Generation
library(rgbif)
library(dplyr)

acceptedNames <- read.csv("data/phylogeneticAngiosperm/acceptedNames.csv", header=TRUE, encoding = "UTF-8")

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

invasiveNames <- read.csv("data/phylogeneticAngiosperm/alienPlantsBelgium.csv", header=TRUE, encoding = "UTF-8")$Name..excl..author.citation.

length(speciesNames[invasiveNames %in% speciesNames])

"
> length(speciesNames[invasiveNames %in% speciesNames])
[1] 610

In this list 610 of the species names that are present are considered invasive according to
the TRIAS list. We will filter out these names and store the rest as considered native species
"
nativeSpecies <- speciesNames[!(invasiveNames %in% speciesNames)]

nativeTaxonBackbone <- lapply(nativeSpecies, name_backbone)
nativeTaxonKeys <- lapply(nativeTaxonBackbone, function(df){df$usageKey})

completeTaxonBackbone <- lapply(speciesNames, name_backbone)
completeTaxonKeys <- lapply(completeTaxonBackbone, function(df){df$usageKey})

#Construct a sql query with the users gbif email account.
nativeSqlQuery <- c("{",
              "  \"sendNotification\": true,",
              "  \"notificationAddresses\": [",
              paste("    \"",GBIF_EMAIL,"\" "),
              "  ],",
              "  \"format\": \"SQL_TSV_ZIP\", ",
              paste("    \"sql\": \"SELECT \\\"year\\\" , GBIF_EEARGCode(1000, decimalLatitude, decimalLongitude, COALESCE(coordinateUncertaintyInMeters, 1000)) AS eeaCellCode, speciesKey, species, establishmentMeans, degreeOfEstablishment, pathway, COUNT(*) AS occurrences, COUNT(DISTINCT recordedBy) AS distinctObservers FROM occurrence WHERE countrycode='BE' AND occurrenceStatus = 'PRESENT' AND hasCoordinate = TRUE AND NOT ARRAY_CONTAINS(issue, 'ZERO_COORDINATE') AND NOT ARRAY_CONTAINS(issue, 'COORDINATE_OUT_OF_RANGE') AND NOT ARRAY_CONTAINS(issue, 'COORDINATE_INVALID') AND NOT ARRAY_CONTAINS(issue, 'COUNTRY_COORDINATE_MISMATCH') AND \\\"year\\\" IS NOT NULL AND taxonKey in (",paste(unlist(nativeTaxonKeys), collapse = ","),") GROUP BY \\\"year\\\", eeaCellCode, speciesKey, species, establishmentMeans, degreeOfEstablishment, pathway ORDER BY \\\"year\\\" DESC, eeaCellCode ASC, speciesKey ASC\""),
              "}")

writeLines(nativeSqlQuery, "dataScripts/nativeAngioCube.json")

completeSqlQuery <- c("{",
                    "  \"sendNotification\": true,",
                    "  \"notificationAddresses\": [",
                    paste("    \"",GBIF_EMAIL,"\" "),
                    "  ],",
                    "  \"format\": \"SQL_TSV_ZIP\", ",
                    paste("    \"sql\": \"SELECT \\\"year\\\" , GBIF_EEARGCode(1000, decimalLatitude, decimalLongitude, COALESCE(coordinateUncertaintyInMeters, 1000)) AS eeaCellCode, speciesKey, species, establishmentMeans, degreeOfEstablishment, pathway, COUNT(*) AS occurrences, COUNT(DISTINCT recordedBy) AS distinctObservers FROM occurrence WHERE countrycode='BE' AND occurrenceStatus = 'PRESENT' AND hasCoordinate = TRUE AND NOT ARRAY_CONTAINS(issue, 'ZERO_COORDINATE') AND NOT ARRAY_CONTAINS(issue, 'COORDINATE_OUT_OF_RANGE') AND NOT ARRAY_CONTAINS(issue, 'COORDINATE_INVALID') AND NOT ARRAY_CONTAINS(issue, 'COUNTRY_COORDINATE_MISMATCH') AND \\\"year\\\" IS NOT NULL AND taxonKey in (",paste(unlist(completeTaxonKeys), collapse = ","),") GROUP BY \\\"year\\\", eeaCellCode, speciesKey, species, establishmentMeans, degreeOfEstablishment, pathway ORDER BY \\\"year\\\" DESC, eeaCellCode ASC, speciesKey ASC\""),
                   "}")

writeLines(completeSqlQuery, "dataScripts/completeAngioCube.json")
