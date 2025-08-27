library(zhulabtools)
load_packages(c("tidyverse","lubridate", "ggplot2", "sf", "vroom"))


################### Combine all the weather data from different stations ###################
parse_mm <- function(x) {
  case_when(
    str_detect(x, "in$") ~ parse_number(x) * 25.4,
    str_detect(x, "mm$") ~ parse_number(x),
    !str_detect(x, "[A-Za-z]$") ~ suppressWarnings(as.numeric(x)),
    TRUE ~ NA_real_
  )
}
parse_C <- function(x) {
  case_when(
    str_detect(x, "F$") ~ (parse_number(x) - 32) * 5/9,
    str_detect(x, "C$") ~ parse_number(x),
    !str_detect(x, "[A-Za-z]$") ~ suppressWarnings(as.numeric(x)),
    TRUE ~ NA_real_
  )
}

flag_identical_runs <- function(x, min_len = 20L, require_positive = FALSE) {
  r <- rle(x)
  keep <- r$lengths >= min_len
  if (require_positive) keep <- keep & (r$values > 0 & !is.na(r$values))
  rep(keep, r$lengths)
}

read_clean_WU <- function(city, run_checks = TRUE) {
  file_path <- "~/phenology-urban/data/raw/CONUS/WU"
  dir_path  <- file.path(file_path, city)
  files     <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE)
  files     <- files[basename(files) != "location.csv"]
  if (length(files) == 0) stop("No station CSVs found for city: ", city)
  print(paste0("Clean ", length(files), " station files for city: ", city))
  all_data <- vroom(files,
                    id = "src",                                 
                    col_types = readr::cols(.default = readr::col_character()),
                    .name_repair = "unique",
                    show_col_types = FALSE,
                    progress = FALSE
  ) %>%
    mutate(name = tools::file_path_sans_ext(basename(src))) %>%
    dplyr::select(-src) %>% 
    dplyr::select(name, Date, HighTemp, AvgTemp, LowTemp, Sum_mm)
  
  all_num <- all_data %>%
    mutate(
      Sum_mm  = parse_mm(Sum_mm),
      HighTemp= parse_C(HighTemp),
      AvgTemp = parse_C(AvgTemp),
      LowTemp = parse_C(LowTemp)
    ) %>%
    filter(!is.na(Date))
  
  if (!run_checks) {
    saveRDS(all_temp, paste0("~/phenology-urban/data/raw/CONUS/WU_clean/", city, "_wu.rds"))
    return(all_temp)
  }
  
  if (run_checks) {
    ################### Clean the data (identify outliers) ###################
    
    ### 1. Identical Value Sequences ###
    print("Starting identical value sequence check...")
    all_data_iden <- all_num %>%
      arrange(name, Date) %>%
      group_by(name) %>%
      mutate(
        HighTemp_Identical = flag_identical_runs(HighTemp),
        AvgTemp_Identical  = flag_identical_runs(AvgTemp),
        LowTemp_Identical  = flag_identical_runs(LowTemp),
        Precip_Identical   = flag_identical_runs(Sum_mm, require_positive = TRUE)
      ) %>%
      ungroup()
    
    ### 2. Extreme Precipitation ###
    print("Starting extreme precipitation check...")
    all_data_iden_ext <- all_data_iden %>%
      mutate(M = month(Date)) %>%
      group_by(name, M) %>%
      mutate(
        .tot_month = sum(Sum_mm, na.rm = TRUE),
        .oth_month = pmax(.tot_month - Sum_mm, 0),
        Sum_mm_Extreme = Sum_mm >= (.oth_month + 300)
      ) %>%
      ungroup() %>%
      select(-.tot_month, -.oth_month, -M)
    
    ### 3. Extreme Temperature ###
    print("Starting extreme temperature check...")
    rollsum_circ <- function(x, k = 31L) {
      stopifnot(length(x) > 0, k %% 2 == 1)
      n <- length(x); h <- (k - 1L) / 2L
      xx <- c(x, x, x)
      cs <- c(0, cumsum(xx))
      starts <- (n + 1L - h):(n + n - h)
      ends   <- (n + 1L + h):(n + n + h)
      cs[ends + 1L] - cs[starts]
    }
    
    build_mu_sd <- function(df, col) {
      df %>%
        mutate(doy = yday(Date)) %>%
        group_by(name, doy) %>%
        summarise(
          s   = sum(.data[[col]], na.rm = TRUE),
          ss  = sum((.data[[col]])^2, na.rm = TRUE),
          nn  = sum(!is.na(.data[[col]])),
          .groups = "drop"
        ) %>%
        group_by(name) %>%
        tidyr::complete(doy = 1:366, fill = list(s = 0, ss = 0, nn = 0)) %>%
        arrange(doy, .by_group = TRUE) %>%
        mutate(
          s_w   = rollsum_circ(s, 31L),
          ss_w  = rollsum_circ(ss, 31L),
          n_w   = rollsum_circ(nn, 31L),
          mu    = ifelse(n_w > 0, s_w / n_w, NA_real_),
          sd    = ifelse(n_w > 1, sqrt(pmax(ss_w / n_w - mu^2, 0)), NA_real_)
        ) %>%
        select(name, doy, mu, sd) %>%
        ungroup()
    }
    
    mu_sd_high <- build_mu_sd(all_data_iden_ext, "HighTemp") %>% rename(mu_H = mu, sd_H = sd)
    mu_sd_avg  <- build_mu_sd(all_data_iden_ext, "AvgTemp")  %>% rename(mu_A = mu, sd_A = sd)
    mu_sd_low  <- build_mu_sd(all_data_iden_ext, "LowTemp")  %>% rename(mu_L = mu, sd_L = sd)
    
    all_data_iden_ext_tp <- all_data_iden_ext %>%
      mutate(doy = yday(Date)) %>%
      left_join(mu_sd_high, by = c("name", "doy")) %>%
      left_join(mu_sd_avg,  by = c("name", "doy")) %>%
      left_join(mu_sd_low,  by = c("name", "doy")) %>%
      mutate(
        HighTemp_Extreme = if_else(!is.na(sd_H) & sd_H > 0, abs((HighTemp - mu_H) / sd_H) >= 6, FALSE),
        AvgTemp_Extreme  = if_else(!is.na(sd_A) & sd_A > 0, abs((AvgTemp  - mu_A) / sd_A) >= 6, FALSE),
        LowTemp_Extreme  = if_else(!is.na(sd_L) & sd_L > 0, abs((LowTemp  - mu_L) / sd_L) >= 6, FALSE)
      ) %>%
      select(-doy, -mu_H, -sd_H, -mu_A, -sd_A, -mu_L, -sd_L)
    
    
    # check_extreme_temp <- function(data, station_name, targetDate, row_index) {
    #   if (row_index %% 1000 == 0) message(glue::glue("Processing row {row_index}"))
    #   targetDate <- as.Date(
    #     targetDate,
    #     tryFormats = c("%m/%d/%Y", "%m/%d/%y")
    #   )
    #   data %>%
    #     filter(name == station_name) %>%
    #     mutate(
    #       Date = as.Date(Date),
    #       month_day = format(Date, format = "%m-%d"),
    #       start_day = format(targetDate - 15, format = "%m-%d"),
    #       end_day   = format(targetDate + 15, format = "%m-%d")
    #     ) %>%
    #     filter(
    #       (start_day > end_day & (month_day >= start_day | month_day <= end_day)) |
    #         (start_day <= end_day & month_day >= start_day & month_day <= end_day)
    #     ) %>%
    #     {
    #       muH <- mean(.$HighTemp, na.rm = TRUE); sdH <- sd(.$HighTemp, na.rm = TRUE)
    #       muA <- mean(.$AvgTemp , na.rm = TRUE); sdA <- sd(.$AvgTemp , na.rm = TRUE)
    #       muL <- mean(.$LowTemp , na.rm = TRUE); sdL <- sd(.$LowTemp , na.rm = TRUE)
    #       mutate(.,
    #              HighTemp_Extreme = sdH > 0 & abs((HighTemp - muH)/sdH) >= 6,
    #              AvgTemp_Extreme  = sdA > 0 & abs((AvgTemp  - muA)/sdA) >= 6,
    #              LowTemp_Extreme  = sdL > 0 & abs((LowTemp  - muL)/sdL) >= 6
    #       )
    #     } %>%
    #     filter(Date == targetDate) %>%
    #     select(-month_day, -start_day, -end_day)
    # }
    # 
    # all_data_iden_ext_tp <- all_data_iden_ext %>%
    #   select(Date, name) %>%
    #   mutate(row_index = row_number(),
    #          Date = as.Date(Date, tryFormats = c("%Y-%m-%d", "%m/%d/%Y", "%m/%d/%y"))) %>%
    #   pmap_dfr(~ check_extreme_temp(all_data_iden_ext, ..2, ..1, ..3))%>%
    #   select(-c(month_day, start_day, end_day))
    # 
    ### 4. Time Consistancy ###
    print("Starting time consistency check...")
    all_data_iden_ext_tcon <- all_data_iden_ext_tp %>%
      arrange(name, Date) %>%
      group_by(name) %>%
      mutate(
        HighTemp_SpikeDip = abs(HighTemp - lag(HighTemp)) >= 25 & abs(HighTemp - lead(HighTemp)) >= 25,
        AvgTemp_SpikeDip  = abs(AvgTemp  - lag(AvgTemp))  >= 25 & abs(AvgTemp  - lead(AvgTemp))  >= 25,
        LowTemp_SpikeDip  = abs(LowTemp  - lag(LowTemp))  >= 25 & abs(LowTemp  - lead(LowTemp))  >= 25
      ) %>%
      ungroup()
    
    ### 5. Spatial Consistancy ###
    print("Starting spatial consistency check...")
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
               AvgTemp_Mean_Others, Precip_Mean_Others))
  
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
    
    saveRDS(final, paste0("~/phenology-urban/data/raw/CONUS/WU_clean/",city,"_wu.rds"))
    all_temp = final
  }
  return(all_temp)
}

# example of running the function
read_clean_WU(city = "Albuquerque", run_checks = TRUE)

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
