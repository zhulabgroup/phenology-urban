library(sp)
library(sf)
library(raster)

extract_temp_alan <- function(point_tree, str_taxa) {
  taxa_point <- point_tree %>% filter(str_detect(taxa, str_taxa))
  
  
  ## insert boundary
  roi_sf <- sf::read_sf('data/raw/NY/boundary/nybb_dissolved.shp')
  ny_boundary <- st_transform(roi_sf, crs = st_crs(4326))
  
  
  
  ## insert alan
  alan_2018 <- raster("~/phenology-urban/data/raw/NY/ALAN/VNP46A4_t2018.tif")
  ## Remove very low values of NTL; can be considered noise
  alan_2018[alan_2018 <= 2] <- 0
  ## Distribution is skewed, so log
  alan_2018 <- log(alan_2018 + 1)
  
  
  ##insert lst
  lstimg <- raster("~/phenology-urban/data/raw/NY/LST_L8/year_mean/2018.tif")
  lst <- lstimg |> mask(ny_boundary)
  
  
  #insert daymet
  tavg_ann <- raster("data/proc/daymet/tavg_ann_2018.tif")
  tavg_spring <- raster("data/proc/daymet/tavg_spring_2018.tif")
  tavg_summer <- raster("data/proc/daymet/tavg_summer_2018.tif")
  
  
  #extract
  extract_lst <- extract(lst, taxa_point)
  taxa_point$lst_2018 <- extract_lst
  
  extract_alan <- extract(alan_2018, taxa_point)
  taxa_point$alan_2018 <- extract_alan
  
  extract_tavg_ann <- extract(tavg_ann, taxa_point)
  taxa_point$tavg_ann_2018 <- extract_tavg_ann
  
  extract_tavg_spring <- extract(tavg_spring, taxa_point)
  taxa_point$tavg_spr_2018 <- extract_tavg_spring
  
  extract_tavg_summer <- extract(tavg_summer, taxa_point)
  taxa_point$tavg_sum_2018 <- extract_tavg_summer
  
  
  
  doy_taxa <- readRDS(str_c("~/phenology-batchplanet/alldata/PSdata/urban/doy/doy_NY_",taxa,".rds"))
  down_doy <- doy_taxa %>%
    filter(year == 2018, direction == "down", thres == 0.5) %>%
    select(id, down_05_2018 = doy)
  up_doy <- doy_taxa %>%
    filter(year == 2018, direction == "up", thres == 0.5) %>%
    select(id, up_05_2018 = doy)
  doy_taxa <- down_doy %>%
    left_join(up_doy, by = "id")
  
  taxa_point <- left_join(taxa_point, doy_taxa, by = "id")
  
  taxa_point <- taxa_point %>%
    filter(!is.na(tavg_ann_2018) & is.finite(tavg_ann_2018)) %>%
    filter(!is.na(alan_2018)) %>%
    filter(!is.na(down_05_2018))
  
  return(taxa_point)
  
}
