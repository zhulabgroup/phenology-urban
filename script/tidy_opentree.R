library(jsonlite)
library(purrr)
library(stringr)
library(tidyverse)
clean <- F

tidy_opentree <- function() {
  if (clean == T){  
    sources_bounds <- fromJSON(here("data/raw/CONUS/StreetTrees/OpenTree/source-stats.json"))
    city_info <- fromJSON(here("data/raw/CONUS/StreetTrees/OpenTree/sources-out.json"))
    
    city_points <- imap_dfr(sources_bounds, function(city_data, city_id) {
      if (!is.null(city_data$bounds)) {
        bounds <- unlist(city_data$bounds)
        # change the order of the coordinates for Cupertino
        if (city_id == "cupertino") {
          bounds <- c(-122.0902, 37.29152, -121.9958, 37.33817)  # xmin, ymin, xmax, ymax
        }
        xmin <- bounds[1]
        ymin <- bounds[2]
        xmax <- bounds[3]
        ymax <- bounds[4]
        
        x_center <- if (xmin == 0) xmax else if (xmax == 0) xmin else (xmin + xmax) / 2
        y_center <- if (ymin == 0) ymax else if (ymax == 0) ymin else (ymin + ymax) / 2
        
        point <- st_point(c(x_center, y_center))
        polygon <- st_polygon(list(matrix(c(
          xmin, ymin,
          xmin, ymax,
          xmax, ymax,
          xmax, ymin,
          xmin, ymin
        ), ncol = 2, byrow = TRUE)))
        
        
      } else {
        point <- st_geometrycollection()
        polygon <- st_geometrycollection()
      }
      
      keep_count <- if (!is.null(city_data$keepCount)) city_data$keepCount else NA_integer_
      
      match_row <- city_info %>%
        filter(tolower(id) == tolower(city_id) | tolower(short) == tolower(city_id)) %>%
        slice(1)
      
      country_name <- if (nrow(match_row) > 0) match_row$country else NA_character_
      city_name <- if (nrow(match_row) > 0) match_row$short else NA_character_
      
      tibble(
        id = city_id,
        TreeCount = keep_count,
        country = country_name,
        city_name = city_name,
        center_geom = st_sfc(point),
        boundary_geom = st_sfc(polygon)
      )
    })
    st_geometry(city_points) <- "center_geom"
    # Bind center_geom back explicitly, then summarize
    city_points_sf <- city_points |>
      filter(country == "USA") |>
      filter(id != "washington_me") |>  # Remove invalid/duplicate city
      mutate(
        city_name = str_replace_all(city_name, " ", ""),
        city_name = if_else(city_name == "Portland,Oregon", "Portland", city_name),
        city_name = if_else(city_name == "UNT", "Denton", city_name),
        city_name = if_else(city_name == "StAugustine", "St.Augustine", city_name)
      ) |>
      group_by(city_name) |>
      summarise(
        TreeCount = sum(TreeCount),
        geometry = st_sfc(first(center_geom)),  # Use st_sfc() to keep geometry shape consistent
        boundary_geom = st_union(boundary_geom),
        .groups = "drop"
      )
    # set crs
    st_crs(city_points_sf) <- 4326
    
    city_points_sf$city_name[city_points_sf$city_name == "Colombus"] <- "Columbus"
    # Change 2 points from the projected coordinates to WGS84
    city_points_sf$geometry[city_points_sf$city_name == "St.Augustine"] <- st_sfc(st_point(c(-81.318, 29.89611)))
    city_points_sf$geometry[city_points_sf$city_name == "CapeCoral"] <- st_sfc(st_point(c(-81.97942, 26.61865)))
    # Change the coordinates of Washington, DC, San Francisco, and Provindence
    city_points_sf$geometry[city_points_sf$city_name == "WashingtonDC"] <- st_sfc(st_point(c(-77.0369, 38.9072)))
    city_points_sf$geometry[city_points_sf$city_name == "Providence"] <- st_sfc(st_point(c(-71.4128, 41.8240)))
    city_points_sf$geometry[city_points_sf$city_name == "SanFrancisco"] <- st_sfc(st_point(c(-122.435, 37.770)))
    
    city_polygons <- st_read("data/raw/CONUS/boundary/USCities/tl_2024_us_uac20.shp") |>
      separate(NAME20, into = c("city_raw", "state"), sep = ", ", remove = FALSE) |>
      # Split rows if city contains '--'
      mutate(city_list = str_split(city_raw, "--")) |>
      unnest(city_list) |>
      # Clean up column names
      rename(city = city_list) |>
      # remove space in the name
      mutate(city = str_replace_all(city, "\\s+", "")) |>
      select(NAME20, state, city, ALAND20,AWATER20, INTPTLAT20, INTPTLON20) |>
      mutate(city = if_else(NAME20 == "Washington--Arlington, DC--VA--MD" & city == "Washington",
                            "WashingtonDC",
                            city))
    
    states_polygon <- st_read(here("data/raw/CONUS/boundary/USStates/s_18mr25.shp")) |>
      select(STATE) |>
      st_transform(crs = 4326)
    
    # Now cross and filter
    sf_use_s2(FALSE)
    matched_polygons <- city_points_sf |>
      rename(geometry_point = geometry) |>
      st_join(states_polygon, join = st_within) |>
      tidyr::crossing(
        city_polygons |> rename(geometry_poly = geometry)
      ) |>
      filter(
        str_detect(city, fixed(city_name, ignore_case = TRUE)),
        str_detect(state, fixed(STATE, ignore_case = TRUE))
      ) |>
      filter(city != "DentonSouthwest") |>
      mutate(geometry = geometry_poly) |>  # use matched city polygon geometry
      select(city_name, TreeCount, geometry) |>
      st_as_sf()
    
    boundary_overlap <- city_points_sf %>%
      mutate(boundary_geom = st_make_valid(boundary_geom)) %>%
      group_split(city_name) %>%
      map_dfr(function(city_row) {
        city_id <- city_row$city_name[1]
        city_boundary <- st_transform(st_sfc(city_row$boundary_geom[[1]], crs = st_crs(city_points_sf)), 3857)
        poly_match <- matched_polygons %>% 
          filter(city_name == city_id) %>%
          st_transform(3857)
        if (nrow(poly_match) == 0) return(NULL)
        if (city_id %in% c("St.Augustine", "CapeCoral")) {
          # Use the full matched polygon, no intersection
          geom <- st_union(poly_match$geometry)
        } else {
          # Proper intersection
          city_boundary <- st_transform(
            st_sfc(city_row$boundary_geom[[1]], crs = st_crs(city_points_sf)), 3857
          )
          geom <- st_intersection(city_boundary, st_union(poly_match$geometry))
        }
        
        tibble(
          city_name = city_id,
          TreeCount = city_row$TreeCount[1],
          geometry = st_sfc(geom)
        )
      }) %>%
      st_as_sf(crs = 3857) %>%
      st_transform(4326)  # optional: back-transform to lat/lon
    
    
    unmatched_fallback <- city_points_sf |>
      filter(!city_name %in% boundary_overlap$city_name) |>
      mutate(geometry = boundary_geom) |>
      select(city_name, TreeCount, geometry) |>
      st_set_crs(4326)
    city_open = bind_rows(boundary_overlap, unmatched_fallback)
    
    saveRDS(city_open, here("~/phenology-urban/data/raw/CONUS/StreetTrees/OpenTree/city_open.rds"))}
  
  city_open = readRDS(here("~/phenology-urban/data/raw/CONUS/StreetTrees/OpenTree/city_open.rds"))
  return(city_open)
}


