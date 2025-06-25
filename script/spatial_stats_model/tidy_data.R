valid_climate <- function(points_in_buffer_joined_d, phe_type){
  if (phe_type == "SOS") {
    phe_metric <- "doy"
    direction <- "up"
    thres <- 0.5
    model_input <- points_in_buffer_joined_d %>%
      filter(direction == !!direction & thres == !!thres) %>%
      select(id, taxa,genus, lon, lat, year, doy, geometry, ends_with("winter")) %>%
      drop_na()
  } else if (phe_type == "EOS") {
    phe_metric <- "doy"
    direction <- "down"
    thres <- 0.5
    model_input <- points_in_buffer_joined_d %>%
      filter(direction == !!direction & thres == !!thres) %>%
      select(id, taxa, genus, lon, lat, year, doy, geometry, ends_with("summer")) %>%
      drop_na()
  } else if (phe_type == "GreenUp") {
    phe_metric <- "doy_diff"
    direction <- "up"
    diff <- points_in_buffer_joined_d %>%
      filter(direction == !!direction & thres %in% c(0.2, 0.8)) %>%
      group_by(id, year) %>%
      summarise(
        doy_diff = abs(diff(doy[order(thres)])),
        .groups = "drop"
      )
    model_input = points_in_buffer_joined_d %>%
      filter(direction == !!direction & thres == 0.2) %>%
      left_join(diff, by = c("year", "id")) %>%
      select(id, taxa, genus, lon, lat, year, doy_diff, geometry, ends_with("spring")) %>%
      drop_na()
    
  } else {
    direction <- "down"
    diff <- points_in_buffer_joined_d %>%
      filter(direction == !!direction & thres %in% c(0.2, 0.8)) %>%
      group_by(id, year) %>%
      summarise(
        doy_diff = abs(diff(doy[order(thres)])),
        .groups = "drop"
      )
    model_input = points_in_buffer_joined_d %>%
      filter(direction == !!direction & thres == 0.2) %>%
      left_join(diff, by = c("year", "id")) %>%
      select(id, taxa, genus, lon, lat, year, doy_diff, geometry, ends_with("fall")) %>%
      drop_na()
  }
  
  return(model_input)
}