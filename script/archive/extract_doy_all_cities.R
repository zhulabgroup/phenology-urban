library(dplyr)
library(stringr)

extract_doy_all_cities <- function(taxa, x, threshold) {
  
  #get all doy attribute
  folder_path <- "~/phenology-batchplanet/alldata/PSdata/doy"
  file_names <- list.files(folder_path, full.names = TRUE)
  taxa_file <- file_names[grep(taxa, file_names)] 
  
  taxa_df <- data.frame()
  for (file in taxa_file) {
    city_code <- gsub(str_c(".*/doy_(.*?)_",taxa,"\\.rds"), "\\1", file)
    doy_taxa <- readRDS(file)
    doy_taxa$city <- city_code
    taxa_df <- rbind(taxa_df, doy_taxa)
  } 
  
  #get all position
  folder_path_ts <- "~/phenology-batchplanet/alldata/PSdata/ts"
  file_names_ts <- list.files(folder_path_ts, full.names = TRUE)
  taxa_file_ts <- file_names_ts[grep(taxa, file_names_ts)] # change the taxa 
  
  location_df <- data.frame()
  for (file in taxa_file_ts) {
    city_code <- gsub(str_c(".*/ps_(.*?)_",taxa,"\\.rds"), "\\1", file)
    ts_taxa <- readRDS(file) %>%
      distinct(id, .keep_all = TRUE)
    ts_taxa$city <- city_code
    location_df <- rbind(location_df, ts_taxa)
  } 
  location_df <- subset(location_df, select = c(id, lon, lat, city))
  
  
  down_doy <- subset(taxa_df, year == x & direction == "down" & thres == threshold, 
                     select = c(id, doy, city)) %>% rename(down = doy)
  up_doy <- subset(taxa_df, year == x & direction == "up" & thres == threshold, 
                   select = c(id, doy, city)) %>% rename(up = doy)
  taxa_doy_df <- down_doy %>%
    left_join(up_doy, by = c("id", "city")) %>%
    left_join(location_df, by = c("id", "city"))
  #read in the doy data, left join the ts (which has the lon/lat)
  
  return(taxa_doy_df)
  
}