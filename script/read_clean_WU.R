library(zhulabtools)
load_packages(c("dplyr","lubridate", "ggplot2", "sf", "readr", "purrr", "stringr", "tidyr"))


################### Combine all the weather data from different stations ###################

read_clean_WU <- function(city = "NY", run_checks = TRUE) {
  file_path <- paste0("~/urban-cooling/data/raw/WU/", city, "/data")
  file_list <- list.files(path = file_path, pattern = "*.csv", full.names = TRUE)
  all_data <- data.frame()
  
  for (file in file_list) {
    temp_data <- read.csv(file)
    temp_data$name <- tools::file_path_sans_ext(basename(file))
    all_data <- bind_rows(all_data, temp_data)
  }
  
  all_data <- all_data %>%
    mutate(
      Date = ifelse(
        grepl("/\\d{2}$", Date),  #
        sub("/(\\d{2})$", "/20\\1", Date),  #
        Date
      ),
      Date = as.Date(Date, format = "%m/%d/%Y") 
    )
  
  
  all_temp <- all_data %>% select(Date, HighTemp, AvgTemp, LowTemp, Sum_mm, name)
  
  columns_to_clean <- c("Sum_mm", "HighTemp", "AvgTemp", "LowTemp")
  
  for (col in columns_to_clean) {
    if (col %in% names(all_temp)) {
      if (col == "Sum_mm") {
        all_temp[[col]] <- ifelse(
          grepl("in$", all_temp[[col]]),                      
          as.numeric(sub(" in", "", all_temp[[col]])) * 25.4, 
          NA
        )
      } else {
        all_temp[[col]] <- ifelse(
          grepl("F$", all_temp[[col]]),  
          (as.numeric(sub(" F", "", all_temp[[col]])) - 32) * 5 / 9, 
          NA
        )
      }
    }
  }
  all_temp <- all_temp %>% 
    filter(!is.na(Sum_mm)) %>%
    filter(!is.na(HighTemp)) %>%
    filter(!is.na(AvgTemp)) %>%
    filter(!is.na(LowTemp))
  
  if (run_checks) {
    ################### Clean the data (identify outliers) ###################
    
    ### 1. Identical Value Sequences ###
    
    check_identical <- function(col, is_precip = FALSE) {
      rle_result <- rle(col)
      if (is_precip) {
        indices <- which(rle_result$lengths >= 20 & rle_result$values > 0)
      } else {
        indices <- which(rle_result$lengths >= 20)
      }
      rle_result$values <- rep(FALSE, length(rle_result$values))
      if (length(indices) > 0) {
        rle_result$values[indices] <- TRUE
      }
      flag <- inverse.rle(rle_result)
      return(flag)
    }
    
    all_data_iden <- all_temp %>%
      group_by(name) %>%
      arrange(Date) %>%
      mutate(
        HighTemp_Identical = check_identical(HighTemp),
        AvgTemp_Identical = check_identical(AvgTemp),
        LowTemp_Identical = check_identical(LowTemp),
        Precip_Identical = check_identical(Sum_mm, is_precip = TRUE)
      ) %>%
      ungroup()
    
    ### 2. Extreme Precipitation ###
    
    check_extreme_prcp <- function(data, station_name, targetDate, row_index) {
      if (row_index %% 1000 == 0) {
        message(glue::glue("Processing row {row_index}"))
      }
      data %>%
        filter(name == station_name) %>%
        mutate(
          month = month(Date)
        ) %>%
        filter(
          (month == month(targetDate))
        ) %>%
        mutate(
          Sum_mm_Extreme = (Sum_mm >= sum(Sum_mm[Date != targetDate], na.rm = TRUE) + 300)) %>%
        filter(Date == targetDate)
    }
    
    all_data_iden_ext <- all_data_iden %>%
      select(Date, name) %>%
      # filter(name == "KNYNEWYO1585") %>%
      mutate(row_index = row_number()) %>%
      pmap_dfr(~ check_extreme_prcp(all_data_iden, ..2, ..1, ..3))
      
    ### 3. Extreme Temperature ###
    
    check_extreme_temp <- function(data, station_name, targetDate, row_index) {
      if (row_index %% 1000 == 0) {
        message(glue::glue("Processing row {row_index}"))
      }
      data %>%
        filter(name == station_name) %>%
        mutate(
          month_day = format(Date, "%m-%d"),
          start_day = format(targetDate - 15, "%m-%d"),
          end_day = format(targetDate + 15, "%m-%d")
        ) %>%
        filter(
          (start_day > end_day & (month_day >= start_day | month_day <= end_day)) |
            (start_day <= end_day & month_day >= start_day & month_day <= end_day)
        ) %>%
        mutate(
          HighTemp_Extreme = abs((HighTemp - mean(HighTemp)) / sd(HighTemp)) >= 6,
          AvgTemp_Extreme = abs((AvgTemp - mean(AvgTemp)) / sd(AvgTemp)) >= 6,
          LowTemp_Extreme = abs((LowTemp - mean(LowTemp)) / sd(LowTemp)) >= 6
          ) %>%
        filter(Date == targetDate)
    }
    
    all_data_iden_ext_tp <- all_data_iden_ext %>%
      select(Date, name) %>%
      mutate(row_index = row_number()) %>%
      pmap_dfr(~ check_extreme_temp(all_data_iden_ext, ..2, ..1, ..3))%>%
      select(-c(month_day, start_day, end_day))
    
    ### 4. Time Consistancy ###
    
    all_data_iden_ext_tcon <- all_data_iden_ext_tp %>%
      group_by(name) %>%
      arrange(Date) %>%
      mutate(
        HighTemp_SpikeDip = abs(HighTemp - lag(HighTemp)) >= 25 & abs(HighTemp - lead(HighTemp)) >= 25,
        AvgTemp_SpikeDip = abs(AvgTemp - lag(AvgTemp)) >= 25 & abs(AvgTemp - lead(AvgTemp)) >= 25,
        LowTemp_SpikeDip = abs(LowTemp - lag(LowTemp)) >= 25 & abs(LowTemp - lead(LowTemp)) >= 25
        )
    
    ### 5. Spatial Consistancy ###
    
    all_data_iden_ext_con <- all_data_iden_ext_tcon %>%
      group_by(Date) %>%
      mutate(
        HighTemp_Mean_Others = map_dbl(name, ~ mean(HighTemp[name != .], na.rm = TRUE)),
        AvgTemp_Mean_Others = map_dbl(name, ~ mean(AvgTemp[name != .], na.rm = TRUE)),
        LowTemp_Mean_Others = map_dbl(name, ~ mean(LowTemp[name != .], na.rm = TRUE)),
        HighTemp_Anomaly = abs(HighTemp - HighTemp_Mean_Others) > 20,
        AvgTemp_Anomaly = abs(AvgTemp - AvgTemp_Mean_Others) > 20,
        LowTemp_Anomaly = abs(LowTemp - LowTemp_Mean_Others) > 20,
        Precip_Mean_Others = map_dbl(name, ~ mean(Sum_mm[name != . & Sum_mm != 0], na.rm = TRUE)),
        Precip_Magnitude_Anomaly = if_else(
          is.na(Precip_Mean_Others), 
          Sum_mm > 100, 
          (Sum_mm - Precip_Mean_Others) >= 100
        )
      ) %>%
      ungroup() %>%
     select(-c(HighTemp_Mean_Others, LowTemp_Mean_Others,
               AvgTemp_Mean_Others, Precip_Mean_Others, month))
      
    
    ### 6. Summary ###
    
    final <- all_data_iden_ext_con %>%
      mutate(across(
        c(HighTemp_Identical, LowTemp_Identical, AvgTemp_Identical, Precip_Identical, 
          HighTemp_Extreme, AvgTemp_Extreme, LowTemp_Extreme, Sum_mm_Extreme, 
          HighTemp_SpikeDip, AvgTemp_SpikeDip, LowTemp_SpikeDip, 
          HighTemp_Anomaly, AvgTemp_Anomaly, LowTemp_Anomaly),
        ~ replace_na(., FALSE)
      )) %>%
      mutate(
        HighTemp_Flag = (HighTemp_Identical | HighTemp_Extreme | HighTemp_SpikeDip | HighTemp_Anomaly),
        AvgTemp_Flag = (AvgTemp_Identical | AvgTemp_Extreme | AvgTemp_SpikeDip | AvgTemp_Anomaly),
        LowTemp_Flag = (LowTemp_Identical | LowTemp_Extreme | LowTemp_SpikeDip | LowTemp_Anomaly),
        Prcp_Flag = (Precip_Identical | Sum_mm_Extreme | Precip_Magnitude_Anomaly)
      ) %>%
      select(-c(HighTemp_Identical, LowTemp_Identical, AvgTemp_Identical, Precip_Identical, 
                HighTemp_Extreme, AvgTemp_Extreme, LowTemp_Extreme, Sum_mm_Extreme, 
                HighTemp_SpikeDip, AvgTemp_SpikeDip, LowTemp_SpikeDip, 
                HighTemp_Anomaly, AvgTemp_Anomaly, LowTemp_Anomaly, Precip_Magnitude_Anomaly)) %>%
      mutate(
        LowTemp = case_when(LowTemp_Flag ~ NA_real_, TRUE ~ LowTemp),
        AvgTemp = case_when(AvgTemp_Flag ~ NA_real_, TRUE ~ AvgTemp),
        HighTemp = case_when(HighTemp_Flag ~ NA_real_, TRUE ~ HighTemp),
        Sum_mm = case_when(Prcp_Flag ~ NA_real_, TRUE ~ Sum_mm)
      )
    
    saveRDS(final, "~/urban-cooling/data/raw/WU/NY/NY_wu.rds")
    return(final)
  }
  return(all_temp)
}



# nyboundary <- st_read("~/urban-cooling/data/raw/NYC/boundary/nybb_dissolved.shp") %>%
#   st_transform(4326)
# 
# # transfer the df to sf object
# wu_location <- read_csv("data/raw/NYC/WU/location.csv") %>%
#   st_as_sf(coords = c("Lon", "Lat"), crs = 4326)
# 
# wu_buffer <- st_buffer(wu_location, dist = 100)
# 
# intersects <- st_intersects(exp_all, wu_buffer)
# points_in_buffer <- exp_all[lengths(intersects) > 0, ]
# 
# table(points_in_buffer$change, points_in_buffer$Year)
# print(nrow(points_in_buffer))
# 
# ggplot() +
#   geom_sf(data = nyboundary, fill = "grey") +
#   geom_sf(data = wu_location, color = "blue", size = 1) +
#   geom_sf(data = points_in_buffer, color = "red", size = 0.1) +
#   theme_minimal()
# 
# 
# # plot the location of the weather stations
# table(points_in_buffer$change, points_in_buffer$Year)
# print(nrow(points_in_buffer))
# 
# ggplot() +
#   geom_sf(data = nyboundary, fill = "grey") +
#   geom_sf(data = wu_location, color = "blue", size = 1) +
#   geom_sf(data = points_in_buffer, color = "red", size = 0.1) +
#   theme_minimal()
# 
#   
# # show the trend
# 
# ggplot(KNJBAYON28, aes(x = Date, y = AvgTemp)) +
#   geom_line() +
#   labs(title = "Daily average Temp. (Station KNJBAYON28)",
#        x = "Date",
#        y = "Temperature (C)") +
#   theme_minimal()
