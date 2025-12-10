city = "NewYork"

input_sos <- readRDS("~/phenology-urban/data/proc/model_output/nooutlier_SOS_input.rds") |>
  filter(city == !!city)
uhat_daily_day <- readRDS(paste0("~/phenology-urban/data/proc/tree_climate/Grid_daily/",city ,"_day.rds"))
uhat_daily_night <- readRDS(paste0("~/phenology-urban/data/proc/tree_climate/Grid_daily/",city ,"_night.rds"))

# night temp

test <- (input_sos) |>
  filter(canonicalname == "Acer rubrum") |>
  mutate(
    date_start = as.Date(start - 1, origin = paste0(year, "-01-01")),
    date_end   = as.Date(end   - 1, origin = paste0(year, "-01-01")),
    date_doy   = as.Date(doy   - 1, origin = paste0(year, "-01-01"))
  ) |>
  left_join(
    uhat_daily_night, by = "id"
  )

date_cols <- grep("^20\\d{2}_\\d{2}_\\d{2}$", names(test), value = TRUE)
df_long <- test %>%
  pivot_longer(
    cols = all_of(date_cols),
    names_to  = "date_str",
    values_to = "temp"
  ) %>%
  filter(!is.na(temp)) %>%
  mutate(date = ymd(gsub("_", "-", date_str)),
         GDD = pmax(temp - 5, 0))

df_gdd_3month <- df_long %>%
  group_by(id, year, doy, temp_winter, prcp_winter, srad_winter) %>%
  filter(date < date_doy,
         date >= date_doy - days(90)) |>
  summarise(
    mean_temp = mean(temp, na.rm = TRUE),
    GDD_sum = sum(GDD, na.rm = TRUE),
    .groups = "drop"
  )
  
df_gdd_greenup <- df_long %>%
  group_by(id, year, doy, temp_winter, prcp_winter, srad_winter) %>%
  filter(date < date_doy,
         date >= date_start) |>
  summarise(
    mean_temp = mean(temp, na.rm = TRUE),
    GDD_sum = sum(GDD, na.rm = TRUE),
    .groups = "drop"
  )


ggplot(df_gdd_greenup, aes(x = mean_temp, y = temp_winter)) +
  geom_point() +
  labs(title = paste0("Acer rubrum SOS vs Nighttime GDD (3 months before SOS)"))

