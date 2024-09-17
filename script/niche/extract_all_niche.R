library(occTest)
library(xfun)

xfun::pkg_attach(c("sf","sp", "dplyr", "terra", "ggplot2", "tidyr", "readxl", "stringr"), install = TRUE) 
library(rlang)

check_exist <- function(x){ 
  arg <- quo_name(enquo(x))
  if(exists(arg, where = .GlobalEnv)){
    return(T)
  } else {
    return(F)
  }
}

###############################                      
### Potential species #########
###############################

little_species <- read.csv("~/lab-data/datasets/vegetation/USTreeAtlas/Little_datatable.csv")
add_city <- read.csv("~/lab-data/proj-phenology-urban/CONUS/climate-niche/guti/add_city.csv")


# GUTI_v1_0 <- read_excel("data/raw/CONUS/climate-niche/guti/GUTI_v1.0.xlsx", 
#                         sheet = "Appendix 3.5", col_names = FALSE, 
#                         skip = 2)
# GUTI_v1_0[1,1] = "Species"
# selected_columns <- GUTI_v1_0 %>%
#   select(where(~ .[1] == "United States of America")) %>%
#   bind_cols(GUTI_v1_0[,1]) %>%  
#   mutate(zero_count = rowSums(. == 0)) %>%
#   subset(zero_count != 200)
# saveRDS(selected_columns, "data/raw/CONUS/climate-niche/guti/us_occ_species.rds")

GUTI_species <- readRDS("~/lab-data/proj-phenology-urban/CONUS/climate-niche/guti/us_occ_species.rds") %>%
  select(...201) %>% slice(-1, -2) %>%
  rename(Latin.Name = ...201)

potential_sp <- inner_join(little_species, GUTI_species, by = "Latin.Name") 
# both recorded in US natural and cities' inventory

little_only_sp <- anti_join(little_species, GUTI_species, by = "Latin.Name") 
# just in US natural env? check with GBIF, if there are records in cities

GUTI_only_sp <- anti_join(GUTI_species, little_species, by = "Latin.Name") %>%
  sample_n(10)
# just in US cities? check with GBIF, if there are original natural habitat, check where it comes from using global GBIF


###############################                      
### Differentiate the urban/natural points ####
############################### 
library(data.table)
df <- fread("~/lab-data/proj-phenology-urban/CONUS/climate-niche/GBIF/us_potential.csv")
df <- df %>%
  dplyr::select("gbifID", "decimalLatitude", "decimalLongitude", "year","genus", "species")

source("script/niche/ext_species_gisa.R")
source("script/niche/run_occtest.R")

# output_env <- rast("~/lab-data/proj-phenology-urban/CONUS/chelsa/chelsa_bio_1km.tif")

path <- "~/lab-data/datasets/climate/WorldClim/wc2.1_10m"
tif_files <- list.files(path, pattern = "bio.*\\.tif$", full.names = TRUE)
output_env <- rast(tif_files)


final_summary_stats <- data.frame()

for (i in 1:nrow(potential_sp)){
  species <- potential_sp$Latin.Name[i]
  us_species_gbif <- ext_species_gisa(dataset = df, species = species)
  
  #natural
  occ_species_natural <- us_species_gbif  %>%
    mutate(datasetName = 'gbif') %>%
    filter(group == "natural")
  
  if (nrow(occ_species_natural) > 0) {
    natural_result <- run_occtest(occ_species = occ_species_natural, output_env = output_env, species = species)
    n_clean_occ <- natural_result$clean_occ
    if (is.null(n_clean_occ)) {
      n_clim_values <- data.frame()
    } else {
      natural_point <- st_as_sf(n_clean_occ[1]$filteredDataset, coords = c("decimalLongitude", "decimalLatitude"), crs = crs(output_env))
      
      n_clim_values <- terra::extract(output_env, natural_point) %>%
        bind_cols(as.data.frame(natural_point)) %>%
        select(1:20) %>%
        mutate(source = "natural")
    }
  } else {
    n_clim_values <- data.frame() 
  }
  
  # urban
  occ_species_urban <- us_species_gbif  %>%
    mutate(datasetName = 'gbif') %>%
    filter(group == "urban") 
  
  if (nrow(occ_species_urban) > 0) {
    urban_result <- run_occtest(occ_species = occ_species_urban, output_env = output_env, species = species)
    u_clean_occ <- urban_result$clean_occ
    if (is.null(u_clean_occ)) {
      u_clim_values <- data.frame()
    } else {
      urban_point <- st_as_sf(u_clean_occ[1]$filteredDataset, coords = c("decimalLongitude", "decimalLatitude"), crs = crs(output_env))
      
      u_clim_values <- terra::extract(output_env, urban_point) %>%
        bind_cols(as.data.frame(urban_point)) %>%
        select(1:20) %>%
        mutate(source = "urban")
    }
  } else {
    u_clim_values <- data.frame() 
  }
  
  # GUTI
  selected_columns <- readRDS("~/lab-data/proj-phenology-urban/CONUS/climate-niche/guti/us_occ_species.rds")
  
  select_rows <- selected_columns %>%
    filter(selected_columns[, 201] == species)  %>%             # define the species
    bind_rows(selected_columns[c(1,2),])
  
  df_without_last_col <- select_rows %>% select(-ncol(select_rows))
  guti <- t(df_without_last_col) %>% as.data.frame()
  colnames(guti) <- c("Occurance", "Country", "Cities")
  us_db <- guti %>%
    filter(Occurance == 1) %>%
    filter(Country == "United States of America") %>%
    left_join(add_city, by = c("Cities" = "City")) %>%
    filter(str_detect(Country.y, "USA")) %>%
    st_as_sf(coords = c("Longitude","Latitude"), crs = 4326)
  
  guti_clim_values <- terra::extract(output_env, us_db) %>%
    bind_cols(as.data.frame(us_db)) %>%
    select(1:20) %>%
    mutate(source = "urban_guti")
  
  # little
  
  shp_name <- little_species %>%
    filter(Latin.Name == species) %>%
    pull("SHP..")
  
  shp_path <- paste0("~/lab-data/datasets/vegetation/USTreeAtlas/shp/", shp_name, "/", shp_name, ".shp")
  
  little_map <- st_read(shp_path) %>%
    st_set_crs(crs(output_env)) # Acer saccharinum
  
  if (check_exist(natural_point)){
    natural_point_little <- st_intersection(natural_point, little_map)
    
    n_l_clim_values <- terra::extract(output_env, natural_point_little)%>%
      bind_cols(as.data.frame(natural_point_little)) %>%
      select(1:20) %>%
      mutate(source = "natural_little")
  } else {
    n_l_clim_values <- data.frame()}
  
  all_clim_values <- bind_rows(n_clim_values, u_clim_values, guti_clim_values, n_l_clim_values)
  
  summary_stats <- all_clim_values %>%
    group_by(source) %>%
    summarise(across(starts_with("wc"), list(
      min = ~ min(., na.rm = TRUE),
      max = ~ max(., na.rm = TRUE),
      mean = ~ mean(., na.rm = TRUE),
      median = ~ median(., na.rm = TRUE),
      quantile_95 = ~ quantile(., 0.95, na.rm = TRUE),
      quantile_5 = ~ quantile(., 0.05, na.rm = TRUE)
    ))) %>%
    mutate(species = species)
  
  final_summary_stats <- bind_rows(final_summary_stats, summary_stats)
  saveRDS(final_summary_stats, "data/proc/conus/climate_niche/worldclim/final_summary_stats.rds")
}


######################
###### Plot ##########
######################

final_summary_stats <- readRDS("~/phenology-urban/data/proc/conus/climate_niche/final_summary_stats.rds") 

filter_4_source <-  final_summary_stats %>%
  group_by(species) %>%
  count(species) %>%
  ungroup() %>%
  filter(n == 4) %>%
  right_join(final_summary_stats, by = "species") %>%
  na.omit() %>%
  subset(source %in% c("natural", "urban"))

source("script/niche/plot_seg.R")
plot_seg(data = filter_4_source, biom = "bio5_")