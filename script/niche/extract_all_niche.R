library(occTest)
library(sp)
library(sf)
library(terra)
library(ggplot2)
library(tidyr)
library(dplyr)
library(readxl)
library(stringr)

###############################                      
### Potential species #########
###############################

little_species <- read.csv("/Volumes/seas-zhukai/datasets/vegetation/USTreeAtlas/Little_datatable.csv")
add_city <- read_csv("data/raw/CONUS/climate-niche/guti/add_city.csv")


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

GUTI_species <- readRDS("data/raw/CONUS/climate-niche/guti/us_occ_species.rds") %>%
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
df <- fread("~/phenology-urban/data/raw/CONUS/climate-niche/GBIF/us_potential.csv")
df <- df %>%
  dplyr::select("gbifID", "decimalLatitude", "decimalLongitude", "year","genus", "species")

source("script/niche/ext_species_gisa.R")
source("script/niche/run_occtest.R")

output_env <- rast("~/phenology-urban/data/raw/CONUS/chelsa/chelsa_bio_1km.tif")

final_summary_stats <- data.frame()
for (i in 1:2){
  species <- potential_sp$Latin.Name[i]
  us_species_gbif <- ext_species_gisa(dataset = df, species = species)
  
  #natural
  occ_species_natural <- us_species_gbif  %>%
    mutate(datasetName = 'gbif') %>%
    filter(group == "natural")
  
  natural_result <- run_occtest(occ_species = occ_species_natural, output_env = output_env, species = species)
  
  n_occ_species_test <- natural_result$occ_species_test
  n_clean_occ <- natural_result$clean_occ
  
  natural_point <- 
    st_as_sf(n_clean_occ[1]$filteredDataset, coords = c("decimalLongitude", "decimalLatitude"), crs = crs(output_env))
  
  n_clim_values <- terra::extract(output_env, natural_point)%>%
    bind_cols(as.data.frame(natural_point)) %>%
    select(1:20) %>%
    mutate(source = "natural")
  
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
  selected_columns <- readRDS("data/raw/CONUS/climate-niche/guti/us_occ_species.rds")
  
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
  
  shp_path <- paste0("/Volumes/seas-zhukai/datasets/vegetation/USTreeAtlas/shp/", shp_name, "/", shp_name, ".shp")
  
  little_map <- st_read(shp_path) %>%
    st_set_crs(crs(output_env)) # Acer saccharinum
  
  natural_point_little <- st_intersection(natural_point, little_map)
  
  n_l_clim_values <- terra::extract(output_env, natural_point_little)%>%
    bind_cols(as.data.frame(natural_point_little)) %>%
    select(1:20) %>%
    mutate(source = "natural_little")
  
  
  all_clim_values <- bind_rows(n_clim_values, u_clim_values, guti_clim_values, n_l_clim_values)

  summary_stats <- all_clim_values %>%
    group_by(source) %>%
    summarise(across(starts_with("bio"), list(
      min = ~ min(., na.rm = TRUE),
      max = ~ max(., na.rm = TRUE),
      mean = ~ mean(., na.rm = TRUE),
      median = ~ median(., na.rm = TRUE),
      quantile_95 = ~ quantile(., 0.95, na.rm = TRUE),
      quantile_5 = ~ quantile(., 0.05, na.rm = TRUE)
    ))) %>%
    mutate(species = species)

  final_summary_stats <- bind_rows(final_summary_stats, summary_stats)
  saveRDS(final_summary_stats, "data/proc/conus/climate_niche/final_summary_stats.rds")
}
    

######################
###### Plot ##########
######################


library(tibble)
#chelsa
data <- tribble(
  ~variable, ~source, ~mean_value, ~median_value, ~p99_value, ~p1_value, ~max_value, ~min_value,
  "bio1", "gbif - natural", 11.2, 11.1, 19.1, 4.85, 23.4, 3.25,
  "bio1", "gbif - urban", 11.7, 12.1, 18.6, 6.05, 20.5, 3.65,
  "bio5", "gbif - natural", 28.7, 28.4, 34.3, 24.0, 37.7, 20.9,
  "bio5", "gbif - urban", 28.6, 28.1, 34.3, 24.9, 38.0, 20.0,
  "bio6", "gbif - natural", -6.10, -5.95, 5.99, -15.4, 14.1, -18.5,
  "bio6", "gbif - urban", -4.87, -3.95, 4.45, -14.3, 8.95, -17.5,
  "bio1", "guti", 10.39951, 10.25, 17.644, 6.252, 18.75, 6.25,
  "bio5", "guti", 27.53835, 27.55, 32.944, 23.456, 34.35, 21.35,
  "bio6", "guti", -6.470388, -6.45, 6.648, -14.848, 8.35, -14.85,
  "bio1", "gbif - little", 11.01146, 10.75, 17.83, 5.75, 19.45, 4.55,
  "bio5", "gbif - little", 28.55017, 28.25, 33.24, 24.16, 33.65, 22.55,
  "bio6", "gbif - little", -6.587493, -6.45, 3.22, -14.95, 6.55, -16.45
)



# 定义 dodge 的宽度
dodge <- position_dodge(width = 0.5)

ggplot(data, aes(x = interaction(source, variable), y = mean_value, color = source)) +
  geom_linerange(aes(ymin = min_value, ymax = max_value), linewidth = 1, show.legend = TRUE) +  # 绘制 min/max 范围
  geom_point(aes(y = median_value, shape = "median"), size = 3, color = "black", show.legend = TRUE) +  # 中值点
  geom_point(aes(y = mean_value, shape = "mean"), size = 3, color = "red", show.legend = TRUE) +  # 平均值点
  geom_point(aes(y = p1_value, shape = "p1"), size = 2, color = "blue", show.legend = TRUE) +  # p1 点
  geom_point(aes(y = p99_value, shape = "p99"), size = 2, color = "blue", show.legend = TRUE) +  # p99 点
  labs(title = "Comparison of Environmental Variables Across Sources",
       x = "Environmental Variable and Source", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # x轴标签倾斜
  scale_color_manual(values = c("gbif - natural" = "darkgreen", "gbif - urban" = "brown", 
                                "guti" = "purple", "gbif - little" = "orange")) +
  scale_shape_manual(name = "Statistics",
                     values = c("mean" = 16, "median" = 15, "p1" = 17, "p99" = 17),
                     labels = c("Mean", "Median", "P1", "P99")) +  # 添加 mean/median/p1/p99 的图例
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +  # 确保图例显示，并调整顺序
  coord_flip()  # 翻转坐标轴


