library(occTest)
library(sp)
library(sf)
library(terra)
library(ggplot2)
library(tidyr)
library(dplyr)
library(readxl)

###############################                      
### Potential species #########
###############################

little_species <- read.csv("/Volumes/seas-zhukai/datasets/vegetation/USTreeAtlas/Little_datatable.csv")


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

us_species_gbif <- ext_species_gisa(dataset = df, species = "Acer saccharinum")

# us_species_gbif <- result$us_species_gbif
# distribution <- result$plot
# 
# print(distribution)


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
#   st_transform(crs(output_env))
# ras_env <- crop(ras_env, us_boundary)
output_env <- rast("~/phenology-urban/data/raw/CONUS/chelsa/chelsa_bio_1km.tif")


###############################                      
### Example: Acer saccharum ####
###############################  

# Test
# do coordinates filtering
# do geo test (outlier, accuracy), no country status
# do institution locality test (botanic garden, biodiversity institutions or GBIF headquarters)
# do human detection to double check

source("script/niche/run_occtest.R")


###############################                      
### Example: Acer saccharinum: natural ####
############################### 

occ_species_natural <- us_species_gbif  %>%
  mutate(datasetName = 'gbif') %>%
  filter(group == "natural")

natural_result <- run_occtest(occ_species = occ_species_natural, output_env = output_env, species = "Acer saccharinum")

# result
n_occ_species_test <- natural_result$occ_species_test
n_clean_occ <- natural_result$clean_occ
table(n_occ_species_test$Reason)

# change the rules
# custom_rules<-readRDS(system.file('ext/fieldMetadata.rds',package='occTest')
# custom_rules$errorThreshold<-rep(0.5,nrow(custom_rules))
# custom_rules[custom_rules$testBlock=="lu","errorThreshold"]<-0.7
# custom_rules[custom_rules$testBlock=="geo","errorThreshold"]<-0.2
# 
# selMM = occFilter(df = outMartesMartes, custom = custom_rules)
list_of_plots <- plot(n_occ_species_test, 
                      occFilter_list = n_clean_occ, 
                      show_plot = T)

list_of_plots

###############################                      
### Example: Acer saccharinum: urban #####
############################### 

# Test
# do coordinates filtering
# do geo test (outlier, accuracy), no country status
# do institution locality test and human detection to double check

occ_species_urban <- us_species_gbif  %>%
  mutate(datasetName = 'gbif') %>%
  filter(group == "urban")

urban_result <- run_occurrence_analysis(occ_species = occ_species_urban, output_env = output_env, species = "Acer saccharinum")

# result
u_occ_species_test <- urban_result$occ_species_test
u_clean_occ <- urban_result$clean_occ
table(u_occ_species_test$Reason)
# change the rules
# custom_rules<-readRDS(system.file('ext/fieldMetadata.rds',package='occTest')
# custom_rules$errorThreshold<-rep(0.5,nrow(custom_rules))
# custom_rules[custom_rules$testBlock=="lu","errorThreshold"]<-0.7
# custom_rules[custom_rules$testBlock=="geo","errorThreshold"]<-0.2
# 
# selMM = occFilter(df = outMartesMartes, custom = custom_rules)
list_of_plots <- plot(u_occ_species_test, 
                      occFilter_list = u_clean_occ, 
                      show_plot = T)

list_of_plots


##################################
###### Niche difference ##########
##################################

# natural/urban comparison
# natural
natural_point <- 
  st_as_sf(n_clean_occ[1]$filteredDataset, coords = c("decimalLongitude", "decimalLatitude"), crs = crs(output_env))
  # n_occ_species_test %>%
  # filter(Reason == "Duplicated--GridCell" | is.na(Reason) ) %>%
  # st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = crs(output_env))

# only remove not te/exact dup coor+year
natural_point <- n_occ_species_test %>%
  filter(Reason != "not terrestrial habitat" | is.na(Reason)) %>%
  distinct(eventDate, decimalLongitude, decimalLatitude, .keep_all = TRUE) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = crs(output_env))

n_clim_values <- terra::extract(output_env, natural_point)%>%
  bind_cols(as.data.frame(natural_point))





# urban
urban_point <- 
  st_as_sf(u_clean_occ[1]$filteredDataset, coords = c("decimalLongitude", "decimalLatitude"), crs = crs(output_env))
  # u_occ_species_test %>%
  # filter(Reason == "Duplicated--GridCell" | is.na(Reason) ) %>%
  # st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = crs(output_env))
urban_point <- u_occ_species_test %>%
  filter(Reason != "not terrestrial habitat" | is.na(Reason)) %>%
  distinct(eventDate, decimalLongitude, decimalLatitude, .keep_all = TRUE)  %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = crs(output_env))

u_clim_values <- terra::extract(output_env, urban_point) %>%
  bind_cols(as.data.frame(urban_point))

# st_as_sf(urban[1]$filteredDataset, coords = c("decimalLongitude", "decimalLatitude"), crs = crs(output_env))

# calculate the range of the environmental variables
variables <- c("bio1", "bio5", "bio6")

n_clim_values_df <- n_clim_values %>% mutate(source = "natural")
u_clim_values_df <- u_clim_values %>% 
  mutate(source = "urban")
  filter(eventDate %in% c(2015, 2005, 1995)) %>%
  mutate(source = case_when(
    eventDate == 1995 ~ "urban-1995",
    eventDate == 2005 ~ "urban-2005",
    eventDate == 2015 ~ "urban-2015"
  ))

combined_values <- bind_rows(n_clim_values_df, u_clim_values_df)

combined_values_long <- combined_values %>%
  select(all_of(variables), source) %>%
  pivot_longer(cols = -source, names_to = "variable", values_to = "value")

# boxplot
ggplot(combined_values_long, aes(x = variable, y = value, fill = source)) +
  geom_boxplot() +
  labs(title = "Natural vs. Urban",
       x = "Environmental Variable", y = "Value") +
  scale_y_continuous(limits = c(12, 12.5)) +
  theme_minimal() +
  scale_fill_manual(values = c("natural" = "darkslategray4", 
                               "urban-1995" = "chocolate4",
                               "urban-2005" = "chocolate3",
                               "urban-2015" = "chocolate1"))

# plot t

summary_stats <- combined_values_long %>%
  group_by(variable, source) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),      # 计算均值
    median_value = median(value, na.rm = TRUE),  # 计算中位数
    p99_value = quantile(value, 0.99, na.rm = TRUE),  # 计算99%分位数
    p1_value = quantile(value, 0.01, na.rm = TRUE),    # 计算1%分位数
    max_value = max(value, na.rm = TRUE),        # 计算最大值
    min_value = min(value, na.rm = TRUE)         # 计算最小值
  )

# 查看结果
print(summary_stats)


natural_ranges <- data.frame(
  variable = variables,
  min_value = sapply(variables, function(var) quantile(n_clim_values[[var]], 0, na.rm = TRUE)),
  max_value = sapply(variables, function(var) quantile(n_clim_values[[var]], 1, na.rm = TRUE))
)
urban_ranges <- data.frame(
  variable = variables,
  min_value = sapply(variables, function(var) quantile(u_clim_values[[var]], 0, na.rm = TRUE)),
  max_value = sapply(variables, function(var) quantile(u_clim_values[[var]], 1, na.rm = TRUE))
)

# visualization
comparison_df <- natural_ranges %>%
  rename(min_natural = min_value, max_natural = max_value) %>%
  inner_join(urban_ranges %>% rename(min_urban = min_value, max_urban = max_value), by = "variable")

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

##################################
###### compare to GUTI: Acer saccharum ####
##################################



selected_columns <- ("data/raw/CONUS/climate-niche/guti/us_occ_species.rds")
select_rows <- selected_columns %>%
  filter(selected_columns[, 201] == "Acer saccharinum")  %>%             # define the species
  bind_rows(selected_columns[c(1,2),])

df_without_last_col <- select_rows %>% select(-ncol(select_rows))
guti <- t(df_without_last_col) %>% as.data.frame()
colnames(guti) <- c("Occurance", "Country", "Cities")


add_city <- read_csv("data/raw/CONUS/climate-niche/guti/add_city.csv")

# instead here should use city boundary - get all climate grids
us_db <- guti%>%
  filter(Occurance == 1) %>%
  filter(Country == "United States of America") %>%
  left_join(add_city, by = c("Cities" = "City")) %>%
  filter(str_detect(Country.y, "USA")) %>%
  st_as_sf(coords = c("Longitude","Latitude"), crs = 4326)

plot <- ggplot() +
  geom_sf(data = us_boundary, fill = "white", color = "black") +
  geom_sf(data = us_db, size = 2, color = "#00B2EE") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Acer saccharinum in GUTI", x = "Longitude", y = "Latitude") +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 24),
        plot.background = element_rect(fill = "transparent", color = NA))
plot

guti_clim_values <- terra::extract(output_env, us_db) %>%
  bind_cols(as.data.frame(us_db))

summary_stats <- guti_clim_values[, c(2,3,4)] %>%
  summarise(across(everything(), list(
    mean = ~ mean(.x, na.rm = TRUE),
    median = ~ median(.x, na.rm = TRUE),
    p99 = ~ quantile(.x, 0.99, na.rm = TRUE),
    p1 = ~ quantile(.x, 0.01, na.rm = TRUE),
    max = ~ max(.x, na.rm = TRUE),
    min = ~ min(.x, na.rm = TRUE)
  )))
print(summary_stats)


##################################
###### Filter use little's map ####### 
##################################

# load the little's map
little_map <- st_read("/Volumes/seas-zhukai/datasets/vegetation/USTreeAtlas/shp/acersacc/acersacc.shp") %>%
  st_set_crs(crs(output_env)) # Acer saccharinum

map <- ggplot() +
  geom_sf(data = little_map, aes(fill = "Little's map"), color = "grey60") +
  scale_fill_manual(values = c("Little's map" = "grey90")) +
  geom_sf(data = us_species_gbif1, size = 0.02, aes(color  = group)) +
  scale_color_manual(values = c("natural" = "darkslategray4", "urban" = "chocolate"), 
                     labels = c("natural" = "GBIF - Natural", "urban" = "GBIF - Urban")) +
  geom_sf(data = us_db, size = 2, aes(shape = "GUTI - City"), color = "magenta3") +
  scale_shape_manual(values = c("GUTI - City" = 3)) +
  geom_sf(data=us_boundary,color="grey60",fill=NA) +
  labs(x = "Longitude", y = "Latitude", 
       color = "Group", fill = "Little's map", shape = "Legend") +
  theme_minimal()

natural_point_little <- st_intersection(natural_point, little_map)

n_l_clim_values <- terra::extract(output_env, natural_point_little)%>%
  bind_cols(as.data.frame(natural_point_little))

summary_stats <- n_l_clim_values[, c(2,3,4)] %>%
  summarise(across(everything(), list(
    mean = ~ mean(.x, na.rm = TRUE),
    median = ~ median(.x, na.rm = TRUE),
    p99 = ~ quantile(.x, 0.99, na.rm = TRUE),
    p1 = ~ quantile(.x, 0.01, na.rm = TRUE),
    max = ~ max(.x, na.rm = TRUE),
    min = ~ min(.x, na.rm = TRUE)
  )))
summary_stats



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
# worldclim
data <- tribble(
  ~variable, ~source, ~mean_value, ~median_value, ~p99_value, ~p1_value, ~max_value, ~min_value,
  "bio1", "gbif - natural", 10.5, 10.2, 19.0, 3.55, 23.6, 2.63,
  "bio1", "gbif - urban", 11.0, 10.7, 19.0, 4.25, 20.9, 2.71,
  "bio5", "gbif - natural", 29.5, 29.2, 34.7, 24.7, 37.2, 23.1,
  "bio5", "gbif - urban", 29.6, 29.3, 35.5, 24.9, 37.2, 21.6,
  "bio6", "gbif - natural", -8.40, -8.20, 3.94, -21.2, 12.1, -23.4,
  "bio6", "gbif - urban", -7.55, -7.52, 3.99, -19.4, 6.23, -22.9,
  "bio1", "guti", 9.772536, 9.89575, 18.4901, 4.779036, 18.74046, 4.722469,
  "bio5", "guti", 28.27661, 28.50175, 32.57762, 24.00812, 34.15425, 22.78857,
  "bio6", "guti", -8.081136, -7.962865, 6.334125, -15.88733, 7.830075, -16.3735,
  "bio1", "gbif - little", 10.36991, 10.047, 17.74517, 4.660585, 19.37501, 3.783552,
  "bio5", "gbif - little", 29.31268, 28.98875, 33.8081, 25.16708, 33.98575, 23.73625,
  "bio6", "gbif - little", -8.88696, -8.4755, 1.148455, -18.91378, 4.2395, -21.4185
)


# 定义 dodge 的宽度
dodge <- position_dodge(width = 0.5)

# 绘制图形
ggplot(data, aes(x = interaction(source, variable), y = mean_value, color = source)) +
  geom_linerange(aes(ymin = min_value, ymax = max_value), size = 1) +  # 绘制 min/max 范围
  geom_point(aes(y = median_value), shape = 15, size = 3, color = "black") +  # 中值点
  geom_point(aes(y = mean_value), shape = 16, size = 3, color = "red") +  # 平均值点
  geom_point(aes(y = p1_value), shape = 17, size = 2, color = "blue") +  # p1 点
  geom_point(aes(y = p99_value), shape = 17, size = 2, color = "blue") +  # p99 点
  labs(title = "Comparison of Environmental Variables Across Sources",
       x = "Environmental Variable and Source", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # x轴标签倾斜
  scale_color_manual(values = c("gbif - natural" = "darkgreen", "gbif - urban" = "brown", 
                                "guti" = "purple", "gbif - little" = "orange")) +
  coord_flip() 


library(ggplot2)

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


