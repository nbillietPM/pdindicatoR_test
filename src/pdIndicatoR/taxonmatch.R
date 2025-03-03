library(parallel)
library(foreach)
library(doParallel)
library(ape)
library(rotl)

#replace the tiplabels that contain the substring 'WITHDASH' with an actual dash '-'
angiospermTree2020$tip.label<-gsub("WITHDASH", "-", angiospermTree2020$tip.label)

Sys.info()["sysname"]

computeCluster <- makeCluster(8)

taxonmatchLagged <- function(tree, nmbSplits){
  #Generate indices to split the list so that we can process the fragments
  chunkIndices <- cut(seq_along(tree$tip.label), breaks = nmbSplits, labels = FALSE)
  #Split up the data based on the indices
  splitList <- split(tree$tip.label, chunkIndices)
  matchedTaxa <- lapply(splitList, function(list){
    result<-tnrs_match_names(list)
    Sys.sleep(5)
    return(result)
    })
  #Combine the list fragments together in an overarching dataframe
  matchedTaxa <- do.call(rbind, matchedTaxa)
  #Extract the taxonomic information
  taxonInformation <- taxonomy_taxon_info(matchedTaxa[!is.na(matchedTaxa$ott_id),]$ott_id)
  #Extract GBIF information from the taxonInformation and store it in the matchedTaxa
  matchedTaxa$gbif_id <- lapply(taxonInformation, function(entry) {
    if (length(entry$tax_sources) >= 2) {
      gbif_value <- entry$tax_sources[[2]]  # Extract GBIF value
      modified_value <- gsub("gbif:", "", gbif_value)  # Apply gsub() to remove "gbif:"
      return(modified_value)
    } else {
      return(NA)  # Handle cases where tax_sources[[2]] is missing
    }
  })
  #Assign the GBIF ID's to the rows that have a OTT ID
  matchedTaxa$gbif_id[!is.na(matchedTaxa$ott_id)] <- as.integer(matchedTaxa$gbif_id)
  original_df <- data.frame(
    #Extract all unique tree labels
    orig_tiplabel = unique(tree$tip.label),
    #Convert all labels to the lowercase format
    search_string = tolower(unique(tree$tip.label)))
  
  #Analogously to the join operator in SQL. Do an inner join by search string
  matchedTaxa <- merge(matchedTaxa, original_df, by = "search_string", all.x = TRUE)
  return(matchedTaxa)
}


#Generate indices to split the list so that we can process the fragments
chunkIndices <- cut(seq_along(tree$tip.label), breaks = nmbSplits, labels = FALSE)
#Split up the data based on the indices
splitList <- split(tree$tip.label, chunkIndices)
matchedTaxa <- lapply(splitList, function(list){
  result<-tnrs_match_names(list)
  Sys.sleep(5)
  return(result)
})

head(matchedTaxa)

taxonInformation <- taxonomy_taxon_info(matchedTaxa[!is.na(matchedTaxa$ott_id),]$ott_id)

gbifID <- lapply(taxonInformation, function(entry) {
  if (length(entry$tax_sources) >= 2) {
    gbif_value <- entry$tax_sources[[2]]  # Extract GBIF value
    modified_value <- gsub("gbif:", "", gbif_value)  # Apply gsub() to remove "gbif:"
    return(modified_value)
  } else {
    return(NA)  # Handle cases where tax_sources[[2]] is missing
  }
})

gbifID <- lapply(taxonInformation, function(taxonInfo){gsub("gbif:", "", taxonInfo$tax_sources[[2]])})

matchedTaxa$gbifID <- unlist(lapply(taxonomy_taxon_info(matchedTaxa[!is.na(matchedTaxa$ott_id),]$ott_id), 
                                    function(taxonInfo){
                                      gsub("gbif:", "", taxonInfo$tax_sources[[2]])
}))








gbifID <- lapply(gbifID, function(id){
  if(is.na(id)==TRUE){
    return(id)  
  }
  else{
    return(as.integer(id))
  }
  })

sum(is.na(gbifID))

matchedTaxa[!is.na(matchedTaxa$ott_id),'gbifID']<-unlist(gbifID)

original_df <- data.frame(
  #Extract all unique tree labels
  orig_tiplabel = unique(angiospermTree2020$tip.label),
  #Convert all labels to the lowercase format
  search_string = tolower(unique(angiospermTree2020$tip.label)))
original_df
#Analogously to the join operator in SQL. Do an inner join by search string
matched_result <- merge(matchedTaxa, original_df, by = "search_string", all.x = TRUE)

matched_result

