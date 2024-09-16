# download GBIF

# install.packages("taxize")
# usethis::edit_r_environ()

library(rgbif)
library(taxize)
library(purrr)
library(dplyr)
library(magrittr) # for %T>% pipe

gbif_taxon_keys <- GUTI_only_sp[[1]] %>% # get the first column of the data.frame
  taxize::get_gbifid_(method="backbone") %>%  # match names to the GBIF backbone to get taxonkeys
  imap(~ .x %>% mutate(original_sciname = .y)) %>%  # add original name back into data.frame
  bind_rows() %T>% # combine all data.frames into one
  readr::write_tsv(path = "all_matches.tsv") %>% # save as side effect for you to inspect if you want
  filter(matchtype == "EXACT" & status == "ACCEPTED") %>% # get only accepted and matched names
  filter(kingdom == "Plantae") %>% # remove anything that might have matched to a non-plant
  pull(usagekey) # get the gbif taxonkeys

occ_download(
  pred_in("taxonKey", gbif_taxon_keys),
  pred_not(pred_in("basisOfRecord",c("FOSSIL_SPECIMEN","LIVING_SPECIMEN"))),
  # pred_gt("elevation", 5000),
  # pred("country", "US"),
  pred("occurrenceStatus","PRESENT"), 
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  pred("hasCoordinate", TRUE),
  pred_gte("year", 1980),
  pred_or(  
    pred_lt("coordinateUncertaintyInMeters",10000),
    pred_isnull("coordinateUncertaintyInMeters")
  ),
  format = "SIMPLE_CSV")


