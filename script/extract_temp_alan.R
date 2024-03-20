library(sp)
library(sf)
library(raster)

extract_temp_alan <- function(point_tree) {
  
  ## insert boundary
  roi_sf <- sf::read_sf('data/raw/CONUS/boundary/us_boundarydissolve.shp')
  ny_boundary <- st_transform(roi_sf, crs = st_crs(4326))
  
  
  ## insert alan
  alan_2018 <- raster("~/phenology-urban/data/raw/CONUS/ALAN/VNP46A4_t2018.tif")
  ## Remove very low values of NTL; can be considered noise
  alan_2018[alan_2018 <= 2] <- 0
  ## Distribution is skewed, so log
  alan_2018 <- log(alan_2018 + 1)
  
  
  ##insert lst
  # lstimg <- raster("~/phenology-urban/data/raw/NY/LST_L8/year_mean/2018.tif")
  # lst <- lstimg |> mask(ny_boundary)
  
  
  #insert daymet
  tavg_ann <- raster("data/proc/daymet/conus/tavg_ann_2018.tif")
  tavg_spring <- raster("data/proc/daymet/conus/tavg_spring_2018.tif")
  tavg_summer <- raster("data/proc/daymet/conus/tavg_summer_2018.tif")
  
  
  #extract
  # extract_lst <- extract(lst, point_tree)
  # point_tree$lst_2018 <- extract_lst
  
  extract_alan <- extract(alan_2018, point_tree)
  point_tree$alan_2018 <- extract_alan
  
  extract_tavg_ann <- extract(tavg_ann, point_tree)
  point_tree$tavg_ann_2018 <- extract_tavg_ann
  
  extract_tavg_spring <- extract(tavg_spring, point_tree)
  point_tree$tavg_spr_2018 <- extract_tavg_spring
  
  extract_tavg_summer <- extract(tavg_summer, point_tree)
  point_tree$tavg_sum_2018 <- extract_tavg_summer
  
  
  
  # doy_taxa <- readRDS(str_c("~/phenology-batchplanet/alldata/PSdata/urban/doy/doy_NY_",taxa,".rds"))
  # down_doy <- doy_taxa %>%
  #   filter(year == 2018, direction == "down", thres == 0.5) %>%
  #   select(id, down_05_2018 = doy)
  # up_doy <- doy_taxa %>%
  #   filter(year == 2018, direction == "up", thres == 0.5) %>%
  #   select(id, up_05_2018 = doy)
  # doy_taxa <- down_doy %>%
  #   left_join(up_doy, by = "id")
  # 
  # point_tree <- left_join(point_tree, doy_taxa, by = "id")
  
  point_tree <- point_tree %>%
    filter(!is.na(tavg_ann_2018) & is.finite(tavg_ann_2018)) %>%
    filter(!is.na(alan_2018))
    # filter(!is.na(down_05_2018))
  
  return(point_tree)
  
}
