fix_map_ratio <- function(xlim, ylim, ratio = 1.5) {
  xmid <- mean(xlim)
  ymid <- mean(ylim)
  
  xr <- diff(xlim)
  yr <- diff(ylim)
  
  current_ratio <- xr / yr
  
  if (current_ratio < ratio) {
    xr_new <- yr * ratio
    xlim <- c(xmid - xr_new / 2, xmid + xr_new / 2)
  } else {
    yr_new <- xr / ratio
    ylim <- c(ymid - yr_new / 2, ymid + yr_new / 2)
  }
  
  list(xlim = xlim, ylim = ylim)
}

search_boundary <- function(city) {
  boundary <- switch(
    city,
    "NewYork" = st_read(here("data/input/boundary/USCities/Example_city/NY/boundary/nybb_dissolved.shp"), quiet = TRUE),
    "DT" = st_read(here("data/input/boundary/USCities/Example_city/DT/City_of_Detroit_Boundary/City_of_Detroit_Boundary.shp"), quiet = TRUE),
    "Seattle" = st_read(here("data/input/boundary/USCities/Example_city/ST/seattle_boundary/CityUGA.shp"), quiet = TRUE) |>
      filter(CITY_NM == "Seattle"),
    "Houston" = st_read(here("data/input/boundary/USCities/Example_city/HT/HGAC_City_Boundaries/HGAC_City_Boundaries.shp"), quiet = TRUE) |>
      filter(NAME == "Houston") |>
      st_make_valid(),  # change to some distance
    "TP" = st_read(here("data/input/boundary/USCities/Example_city/TP/Tampa_Planning_Districts/Tampa_Planning_Districts.shp"), quiet = TRUE) |>
      st_make_valid() |> st_union(),
    "Denver" = st_read(here("data/input/boundary/USCities/Example_city/DV/DRCOG_Municipalities/DRCOG_Municipalities.shp"), quiet = TRUE) |>
      filter(city == "Denver"),
    "Austin" = st_read(here("data/input/boundary/USCities/Example_city/AT/TxDOT_City_Boundaries/Cities.shp"), quiet = TRUE) |>
      filter(CITY_NM == "Austin"),
    "Rochester" = st_read(here("data/input/boundary/USCities/Example_city/Rochester/City_of_Rochester_Border/City_of_Rochester_Border.shp"), quiet = TRUE),
    "Stockton" = st_read(here("data/input/boundary/USCities/Example_city/Stockton/City_Limits/geo_export_95d77263-4f44-4046-b5b9-23c00947fd1f.shp"), quiet = TRUE) |>
      filter(name == "STOCKTON"),
    "LosAngeles" = st_read(here("data/input/boundary/USCities/Example_city/LosAngeles/City_Boundary.shp"), quiet = TRUE),
    "SanFrancisco" = st_read(here("data/input/boundary/USCities/Example_city/SF/Planning_Neighborhood_Groups_Map/geo_export_a4e97979-9d5a-4474-8af3-073be4dffae4.shp"), quiet = TRUE) |>
      st_make_valid() |> st_union(),
    "Knoxville" = st_read(here("data/input/boundary/USCities/Example_city/Knoxville/TN_City_Boundaries/TN_City_Boundaries.shp"), quiet = TRUE) |>
      filter(NAME == "Knoxville"),
    "WashingtonDC" = st_read(here("data/input/boundary/USCities/Example_city/WashingtonDC/DC_Boundary/DC_Boundary.shp"), quiet = TRUE),
    "Columbus" = st_read(here("data/input/boundary/USCities/Example_city/Columbus/Corporate_Boundary/Corporate_Boundary.shp"), quiet = TRUE),
    "Arlington" = st_read(here("data/input/boundary/USCities/Example_city/Arlington/Municipal_Boundary/Municipal_Boundary.shp"), quiet = TRUE)|> 
      filter(MUNICIPALN == "ARLINGTON"),
    "Sacramento" = st_read(here("data/input/boundary/USCities/Example_city/Sacramento/City_and_County_Boundary/City_and_County_Boundary.shp"), quiet = TRUE) |>
      filter(CITY_NAME == "SACRAMENTO") |>
      st_make_valid(),
    "Atlanta" = st_read(here("data/input/boundary/USCities/Example_city/Atlanta/Quadrants/Quadrants.shp"), quiet = TRUE) |>
      st_make_valid() |> st_union()  |>
      nngeo::st_remove_holes(),
    "NewOrleans" = st_read(here("data/input/boundary/USCities/Example_city/NewOrleans/Orleans Parish Boundary_20260416/geo_export_2201ea0f-3a1c-4066-aec3-a8e995d07a7a.shp"), quiet = TRUE),
    "Portland" = st_read(here("data/input/boundary/USCities/Example_city/Portland/City_Boundaries/City_Boundaries.shp"), quiet = TRUE) |>
      filter(CITYNAME == "Portland") |>
      st_make_valid(),
    "Bloomington" = st_read(here("data/input/boundary/USCities/Example_city/Bloomington/Bloomington Municipal Boundary_20260421/geo_export_951c8648-2ab5-41e2-820c-1a4d268e2537.shp"), quiet = TRUE),
    "Norfolk" = st_read(here("data/input/boundary/USCities/Example_city/Norfolk/City_Limit/City_Limit_-_City_of_Norfolk.shp"), quiet = TRUE),
    stop("Invalid city abbreviation.")
  )
  boundary <- st_transform(boundary, crs = 4326)  # ensure WGS84
  return(boundary)
}
