library(dplyr)

# data3 <- read.delim("~/phenology-urban/data/raw/CONUS/climate-niche/GBIF/Acer_20240617_1980_2024/occurrence.txt",
#                    sep = "\t",header = TRUE, quote = "", fill = T, fileEncoding = "UTF-8",stringsAsFactors = FALSE)
# saveRDS(data3, "~/phenology-urban/data/raw/CONUS/climate-niche/GBIF/us_acer_gbif.rds")
data <- readRDS("~/phenology-urban/data/raw/CONUS/climate-niche/GBIF/us_acer_gbif.rds") %>% 
  distinct() %>%
  filter(year <= 2023) %>%
  filter(species = "Acer saccharum")%>%
  dplyr::select("gbifID", "decimalLatitude", "decimalLongitude", "year","genus", "specificEpithet", "species")  

## check the distribution of laughlin's data


library(readr)
laughlin <- read_csv("data/raw/CONUS/climate-niche/Laughlin24/occurrences.csv") %>%
  filter(species == "Alnus rubra")

max_bio1 <- laughlin %>%
  group_by(species) %>%
  filter(bio1 == max(bio1)) %>%
  slice(1) %>%
  mutate(type = "max") %>%
  ungroup()

min_bio1 <- laughlin %>%
  group_by(species) %>%
  filter(bio1 == min(bio1)) %>%
  slice(1) %>%
  mutate(type = "min") %>%
  ungroup()

check <- bind_rows(max_bio1, min_bio1) %>%
  group_by(species, source) %>%
  summarise(count = n()) %>%
  ungroup()


variables <- c("bio1", "bio5", "bio6")

p_df <- laughlin %>% mutate(group = "potential")
r_df <- laughlin %>% 
  filter(source == "bien") %>%
  mutate(group = "realized")
combined_values <- bind_rows(p_df, r_df)

combined_values_long <- combined_values %>%
  select(all_of(variables), group) %>%
  pivot_longer(cols = -group, names_to = "variable", values_to = "value")

# boxplot
ggplot(combined_values_long, aes(x = variable, y = value, fill = group)) +
  geom_boxplot() +
  labs(title = "Natural vs. Urban",
       x = "Environmental Variable", y = "Value") +
  theme_minimal() +
  scale_fill_manual(values = c("potential" = "darkslategray4", "realized" = "chocolate"))


p_ranges <- data.frame(
  variable = variables,
  min_value = sapply(variables, function(var) quantile(p_df[[var]], 0.01, na.rm = TRUE)),
  max_value = sapply(variables, function(var) quantile(p_df[[var]], 0.99, na.rm = TRUE))
)
r_ranges <- data.frame(
  variable = variables,
  min_value = sapply(variables, function(var) quantile(r_df[[var]], 0.01, na.rm = TRUE)),
  max_value = sapply(variables, function(var) quantile(r_df[[var]], 0.99, na.rm = TRUE))
)

# visualization
comparison_df <- p_ranges %>%
  rename(min_p = min_value, max_p = max_value) %>%
  inner_join(r_ranges %>% rename(min_r = min_value, max_r = max_value), by = "variable")

# visualization
comparison_df_long <- comparison_df %>%
  pivot_longer(cols = c(min_p, max_p, min_r, max_r),
               names_to = c("set", "type"),
               names_sep = "_") %>%
  pivot_wider(names_from = set, values_from = value)

ggplot(comparison_df_long) +
  geom_segment(aes(x = as.numeric(as.factor(variable)) + ifelse(type == "p", -0.1, 0.1),
                   xend = as.numeric(as.factor(variable)) + ifelse(type == "p", -0.1, 0.1),
                   y = min, yend = max, color = type),
               linewidth = 1, lineend = "round") +
  scale_x_continuous(breaks = 1:length(unique(comparison_df_long$variable)),
                     labels = unique(comparison_df_long$variable)) +
  labs(x = "Variable", y = "Range", title = "Potential vs. Realized") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("p" = "darkslategray4", "r" = "chocolate")) +
  coord_flip()


tif_dir <- "/Volumes/seas-zhukai/datasets/climate/WorldClim/wc2.1_10m/"
tif_files <- list.files(tif_dir, pattern = "bio_(1|5|6)\\.tif$", full.names = TRUE)

ras_env <- rast(tif_files)
original_names <- names(ras_env)
new_names <- sub("CHELSA_(bio\\d{1,2})_.*", "\\1", original_names)
names(ras_env) <- c("bio1", "bio5", "bio6")

# mask by the us boundary
us_boundary <- st_read("~/phenology-urban/data/raw/CONUS/boundary/us_boundarydissolve.shp") %>%
  st_transform(crs(ras_env))
output_env <- crop(ras_env, us_boundary)


raster_df <- as.data.frame(ny_gisa, xy = TRUE, na.rm = TRUE)
head(raster_df)

ggplot(raster_df, aes(x = x, y = y, fill = GISA_v02_1972_2019_45N_70W)) +
  geom_tile() +  # 使用 tile 绘制栅格
  scale_fill_gradientn(
    colors = c("darkslategray4", "chocolate4", "orange"),
    values = c(0, 1/35, 1),  # 设置颜色渐变比例
    breaks = c(0, 35),  # 设置图例显示的断点
    labels = c("0", "1-35")) +
  labs(title = "GISA Raster", fill = "Value") +
  coord_fixed() +  # 保持x和y的比例
  theme_minimal()
