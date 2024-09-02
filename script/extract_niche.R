library(occTest)
library(sp)
library(sf)
library(terra)
library(ggplot2)
library(tidyr)
library(dplyr)

###############################                      
### Differentiate the urban/natural points
############################### 

us_gisa <- rast("~/phenology-urban/data/raw/CONUS/GISA/us_gisa.tif")
us_saccharum_gbif1 <- us_saccharum_gbif %>%
  rename(lon = decimalLongitude, lat = decimalLatitude)
us_saccharum_gbif1 <- vect(us_saccharum_gbif1[, c("lon", "lat")], crs = crs(us_gisa))

# Extract the GISA values
gbif_gisa_values2 <- terra::extract(us_gisa, us_saccharum_gbif1)
us_saccharum_gbif <- us_saccharum_gbif %>%
  mutate(gisa = gbif_gisa_values[, 2]) %>%
  mutate(gisa = recode(gisa,
                       `1` = 1972, `2` = 1978, `3` = 1985, `4` = 1986, `5` = 1987, 
                       `6` = 1988, `7` = 1989, `8` = 1990, `9` = 1991, `10` = 1992, 
                       `11` = 1993, `12` = 1994, `13` = 1995, `14` = 1996, `15` = 1997, 
                       `16` = 1998, `17` = 1999, `18` = 2000, `19` = 2001, `20` = 2002, 
                       `21` = 2003, `22` = 2004, `23` = 2005, `24` = 2006, `25` = 2007, 
                       `26` = 2008, `27` = 2009, `28` = 2010, `29` = 2011, `30` = 2012, 
                       `31` = 2013, `32` = 2014, `33` = 2015, `34` = 2016, `35` = 2017, 
                       `36` = 2018, `37` = 2019))
# Natrual: record before urbanization and always natural
# Urban: record after urbanization
us_saccharum_gbif <- us_saccharum_gbif %>%
  mutate(group = case_when(
    gisa == 0 ~ "natural",
    gisa != 0 & gisa <= year ~ "urban",
    gisa != 0 & gisa > year ~ "natural"
  )) %>%
  filter(!(gisa == 0 & year > 2019)) # remove the natrual records after 2019, because we cannot determine whether it is urban or not

# Visualization
us_boundary <- st_read("~/phenology-urban/data/raw/CONUS/boundary/us_boundarydissolve.shp")
us_saccharum_gbif1 <- st_as_sf(us_saccharum_gbif, coords = c("decimalLongitude", "decimalLatitude"), crs = crs(us_gisa))
ggplot(us_saccharum_gbif1) +
  geom_sf(aes(color = group), size = 0.02) +
  geom_sf(data=us_boundary,color="grey60",fill=NA) +
  labs(x = "Longitude", y = "Latitude", color = "Group") +
  theme_minimal()


###############################
### Prepare th envRaster
###############################

# #load the raster of the environmental variables
# tif_dir <- "/Volumes/seas-zhukai/datasets/climate/CHELSA/climatology/"
# tif_files <- list.files(tif_dir, pattern = "bio\\d{1,2}_.*\\.tif$", full.names = TRUE)
# ras_env <- rast(tif_files)
# original_names <- names(ras_env)
# new_names <- sub("CHELSA_(bio\\d{1,2})_.*", "\\1", original_names)
# names(ras_env) <- new_names
# 
# # mask by the us boundary
# us_boundary <- st_read("~/phenology-urban/data/raw/CONUS/boundary/us_boundarydissolve.shp") %>%
#   st_transform(crs(ras_env))
# ras_env <- crop(ras_env, us_boundary)
output_env <- rast("~/phenology-urban/data/raw/CONUS/chelsa/merged_rasters.tif")


###############################                      
### Example: Acer saccharum: natural points
###############################  

# Test
# do coordinates filtering
# do geo test (outlier, accuracy), no country status
# do institution locality test (botanic garden, biodiversity institutions or GBIF headquarters)
# do human detection to double check


#load data
occ_species_natural <- us_saccharum_gbif  %>%
  mutate(datasetName = 'gbif') %>%
  filter(group == "natural")

names(occ_species_natural)[4] <- 'eventDate'

default_settings <- defaultSettings()
my_analysis_settings <- default_settings$analysisSettings

#modify settings
# Decide whether conduct landuse and human influence tests - only do institution locality test
# Because I rely on GISA data to filter urban distributed points (I could do it for double check)
my_analysis_settings$humanDetection$doHumanDetection <- F
my_analysis_settings$landUse$doLandUse <- F
my_analysis_settings$institutionLocality$doInstitutionLocality <- T

# Decide whether conduct country and range analysis
my_analysis_settings$centroidDetection$doCentroidDetection <- F
my_analysis_settings$countryStatusRange$doCountryRecordAnalysis <- F
my_analysis_settings$countryStatusRange$doRangeAnalysis <-  F


# Decide whether conduct time accuracy test
my_analysis_settings$doTimeAccuracy <- F

# Decide whether conduct environment outliers test
my_analysis_settings$geoOutliers$doGeoOutliers <- F 
my_analysis_settings$envOutliers$doEnvOutliers <- F      # environment outliers: identified in niche calculation

myGeoOutliers = list (
  doGeoOutliers=TRUE,
  methodGeoOutliers=c('alphaHull', 'distance', 'median', 'grubbs','dixon','rosner','mcp'),  # 'quantilesamplingcorr': Error: Cannot parse CN-TW into a known Country
  alpha.parameter = 2,
  mcp_percSample =95
)

my_analysis_settings$geoOutliers <- myGeoOutliers     # geographical outliers

my_writeout_settings <- set_writeout(output.dir = "~/phenology-urban/data/proc/conus/climate_niche/acer_saccharum_natrual" ,writeAllOutput = T)
#run workflow
occ_species_test <- occTest(sp.name = 'Acer saccharum',
                            sp.table = occ_species_natural,
                            r.env =output_env,
                            analysisSettings = my_analysis_settings,
                            writeoutSettings = my_writeout_settings,
                            verbose = T)
table(occ_species_test$geoOutliers_score,exclude = F)

selMM = occFilter(df = occ_species_test,errorAcceptance = 'strict')
# if (errorAcceptance == 'strict')    {colMetaData$errorThreshold =0.2}
# if (errorAcceptance == 'majority')  {colMetaData$errorThreshold =0.5}
# if (errorAcceptance == 'relaxed')   {colMetaData$errorThreshold =0.7}
names(selMM)
selMM[['filteredDataset']]
head(selMM[['summaryStats']])


# change the rules
# custom_rules<-readRDS(system.file('ext/fieldMetadata.rds',package='occTest')
# custom_rules$errorThreshold<-rep(0.5,nrow(custom_rules))
# custom_rules[custom_rules$testBlock=="lu","errorThreshold"]<-0.7
# custom_rules[custom_rules$testBlock=="geo","errorThreshold"]<-0.2
# 
# selMM = occFilter(df = outMartesMartes, custom = custom_rules)
list_of_plots <- plot(occ_species_test, 
                      occFilter_list = selMM, 
                      show_plot = T)

list_of_plots




###############################                      
### Example: Acer saccharum: urban
############################### 

# Test
# do coordinates filtering
# do geo test (outlier, accuracy), no country status
# do institution locality test and human detection to double check


# Test
# do coordinates filtering
# do geo test (outlier, accuracy), no country status
# do institution locality test (botanic garden, biodiversity institutions or GBIF headquarters)
# do human detection to double check


#load data
occ_species_urban <- us_saccharum_gbif  %>%
  mutate(datasetName = 'gbif') %>%
  filter(group == "urban")

names(occ_species_urban)[4] <- 'eventDate'

default_settings <- defaultSettings()
my_analysis_settings <- default_settings$analysisSettings

#modify settings
# Decide whether conduct landuse and human influence tests - only do institution locality test
# Because I rely on GISA data to filter urban distributed points (I could do it for double check)
my_analysis_settings$humanDetection$doHumanDetection <- F
my_analysis_settings$landUse$doLandUse <- F

# Not decide whether conduct institution locality test (keep it, because it still located in urban env)
my_analysis_settings$institutionLocality$doInstitutionLocality <- F

# Decide whether conduct country and range analysis
my_analysis_settings$centroidDetection$doCentroidDetection <- F
my_analysis_settings$countryStatusRange$doCountryRecordAnalysis <- F
my_analysis_settings$countryStatusRange$doRangeAnalysis <- F

# Decide whether conduct time accuracy test
my_analysis_settings$doTimeAccuracy <- F

# Decide whether conduct environment outliers test
my_analysis_settings$envOutliers$doEnvOutliers <- F      # environment outliers: identified in niche calculation
my_analysis_settings$geoOutliers$doGeoOutliers <- F

myGeoOutliers = list (
  doGeoOutliers=TRUE,
  methodGeoOutliers=c('alphaHull', 'distance', 'median', 'grubbs','dixon','rosner','mcp'),  # 'quantilesamplingcorr': Error: Cannot parse CN-TW into a known Country
  alpha.parameter = 2,
  mcp_percSample =95
)

my_analysis_settings$geoOutliers <- myGeoOutliers     # geographical outliers

my_writeout_settings <- set_writeout(output.dir = "~/phenology-urban/data/proc/conus/climate_niche/acer_saccharum_urban" ,writeAllOutput = T)
#run workflow
occ_species_test_urban <- occTest(sp.name = 'Acer saccharum',
                                  sp.table = occ_species_urban,
                                  r.env =output_env,
                                  analysisSettings = my_analysis_settings,
                                  writeoutSettings = my_writeout_settings,
                                  verbose = T)
table(occ_species_test_urban$Reason,exclude = F)

urban = occFilter(df = occ_species_test_urban,errorAcceptance = 'strict')
# if (errorAcceptance == 'strict')    {colMetaData$errorThreshold =0.2}
# if (errorAcceptance == 'majority')  {colMetaData$errorThreshold =0.5}
# if (errorAcceptance == 'relaxed')   {colMetaData$errorThreshold =0.7}
names(urban)
urban[['filteredDataset']]
head(urban[['summaryStats']])


# change the rules
# custom_rules<-readRDS(system.file('ext/fieldMetadata.rds',package='occTest')
# custom_rules$errorThreshold<-rep(0.5,nrow(custom_rules))
# custom_rules[custom_rules$testBlock=="lu","errorThreshold"]<-0.7
# custom_rules[custom_rules$testBlock=="geo","errorThreshold"]<-0.2
# 
# selMM = occFilter(df = outMartesMartes, custom = custom_rules)
list_of_plots <- plot(occ_species_test_urban, 
                      occFilter_list = urban, 
                      show_plot = T)

list_of_plots



##################################
###### Niche difference ##########
##################################

# check the niche difference after the filter test

# natural/urban before and after
original_point <- st_as_sf(occ_species_urban, coords = c("decimalLongitude", "decimalLatitude"), crs = crs(output_env))
extracted_values <- terra::extract(output_env, original_point)

filter_point <- st_as_sf(urban[1]$filteredDataset, coords = c("decimalLongitude", "decimalLatitude"), crs = crs(output_env))
filter_values <- terra::extract(output_env, filter_point)

##BIO1 = Annual Mean Temperature
##BIO5 = Max Temperature of Warmest Month
##BIO6 = Min Temperature of Coldest Month
##BIO12 = Annual Precipitation
##BIO17 = Precipitation of Driest Quarter

# calculate the range of the environmental variables
variables <- c("bio1", "bio5", "bio6")

original_ranges <- data.frame(
  variable = variables,
  min_value = sapply(variables, function(var) quantile(extracted_values[[var]], 0.01, na.rm = TRUE)),
  max_value = sapply(variables, function(var) quantile(extracted_values[[var]], 0.99, na.rm = TRUE))
)
filter_ranges <- data.frame(
  variable = variables,
  min_value = sapply(variables, function(var) quantile(filter_values[[var]], 0.01, na.rm = TRUE)),
  max_value = sapply(variables, function(var) quantile(filter_values[[var]], 0.99, na.rm = TRUE))
)


comparison_df <- original_ranges %>%
  rename(min_original = min_value, max_original = max_value) %>%
  inner_join(filter_ranges %>% rename(min_filter = min_value, max_filter = max_value), by = "variable")

# visualization
comparison_df_long <- comparison_df %>%
  pivot_longer(cols = c(min_original, max_original, min_filter, max_filter),
               names_to = c("set", "type"),
               names_sep = "_") %>%
  pivot_wider(names_from = set, values_from = value)

ggplot(comparison_df_long) +
  geom_segment(aes(x = as.numeric(as.factor(variable)) + ifelse(type == "original", -0.1, 0.1),
                   xend = as.numeric(as.factor(variable)) + ifelse(type == "original", -0.1, 0.1),
                   y = min, yend = max, color = type),
               linewidth = 1, lineend = "round") +
  scale_x_continuous(breaks = 1:length(unique(comparison_df_long$variable)),
                     labels = unique(comparison_df_long$variable)) +
  labs(x = "Variable", y = "Range", title = "Urban before and after testing") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("original" = "blue", "filter" = "red")) +
  coord_flip()




# natural/urban comparison
# natural
original_point <- st_as_sf(selMM[1]$filteredDataset, coords = c("decimalLongitude", "decimalLatitude"), crs = crs(output_env))
  # occ_species_test %>%
  # filter(Reason == "Duplicated--GridCell" | is.na(Reason) ) %>%
  # st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = crs(output_env))
extracted_values <- terra::extract(output_env, original_point)

# urban
filter_point <- st_as_sf(urban[1]$filteredDataset, coords = c("decimalLongitude", "decimalLatitude"), crs = crs(output_env))
  # occ_species_test_urban %>%
  # filter(Reason == "Duplicated--GridCell" | is.na(Reason) ) %>%
  # st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = crs(output_env))
filter_values <- terra::extract(output_env, filter_point)

# st_as_sf(urban[1]$filteredDataset, coords = c("decimalLongitude", "decimalLatitude"), crs = crs(output_env))

# calculate the range of the environmental variables
variables <- c("bio12", "bio17")

extracted_values_df <- extracted_values %>% mutate(source = "natural")
filter_values_df <- filter_values %>% mutate(source = "urban")
combined_values <- bind_rows(extracted_values_df, filter_values_df)

combined_values_long <- combined_values %>%
  select(all_of(variables), source) %>%
  pivot_longer(cols = -source, names_to = "variable", values_to = "value")

# boxplot
ggplot(combined_values_long, aes(x = variable, y = value, fill = source)) +
  geom_boxplot() +
  labs(title = "Natural vs. Urban",
       x = "Environmental Variable", y = "Value") +
  theme_minimal() +
  scale_fill_manual(values = c("natural" = "darkslategray4", "urban" = "chocolate"))




original_ranges <- data.frame(
  variable = variables,
  min_value = sapply(variables, function(var) quantile(extracted_values[[var]], 0, na.rm = TRUE)),
  max_value = sapply(variables, function(var) quantile(extracted_values[[var]], 1, na.rm = TRUE))
)
filter_ranges <- data.frame(
  variable = variables,
  min_value = sapply(variables, function(var) quantile(filter_values[[var]], 0, na.rm = TRUE)),
  max_value = sapply(variables, function(var) quantile(filter_values[[var]], 1, na.rm = TRUE))
)

# visualization
comparison_df <- original_ranges %>%
  rename(min_natural = min_value, max_natural = max_value) %>%
  inner_join(filter_ranges %>% rename(min_urban = min_value, max_urban = max_value), by = "variable")

# visualization
comparison_df_long <- comparison_df %>%
  pivot_longer(cols = c(min_natural, max_natural, min_urban, max_urban),
               names_to = c("set", "type"),
               names_sep = "_") %>%
  pivot_wider(names_from = set, values_from = value)

ggplot(comparison_df_long) +
  geom_segment(aes(x = as.numeric(as.factor(variable)) + ifelse(type == "natural", -0.1, 0.1),
                   xend = as.numeric(as.factor(variable)) + ifelse(type == "natural", -0.1, 0.1),
                   y = min, yend = max, color = type),
               linewidth = 1, lineend = "round") +
  scale_x_continuous(breaks = 1:length(unique(comparison_df_long$variable)),
                     labels = unique(comparison_df_long$variable)) +
  labs(x = "Variable", y = "Range", title = "Urban vs. Natural") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("natural" = "darkslategray4", "urban" = "chocolate")) +
  coord_flip()
