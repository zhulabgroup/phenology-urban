inventory_dir <- here("data/input/tree_list/")
output_dir    <- here("data/input/tree_list/summary/")

file_list <- list.files(inventory_dir, pattern = "\\.rds$", full.names = TRUE)

all_centroids <- list()
all_bboxes    <- list()

for (f in file_list) {
  message("Processing: ", f)
  obj <- tryCatch(readRDS(f), error = function(e) { message("  ! read error: ", e$message); return(NULL) })
  if (is.null(obj)) next
  
  sf_pts <- NULL
  if (inherits(obj, "sf")) {
    sf_pts <- obj
  } else {
    cn <- names(obj)
    lon_col <- cn[str_detect(tolower(cn), "^(lon|longitude|x)$")][1]
    lat_col <- cn[str_detect(tolower(cn), "^(lat|latitude|y)$")][1]
    
    if (!is.na(lon_col) && !is.na(lat_col)) {
      df <- obj %>% filter(!is.na(.data[[lon_col]]), !is.na(.data[[lat_col]]))
      if (nrow(df) == 0) next
      sf_pts <- st_as_sf(df, coords = c(lon_col, lat_col), crs = 4326, remove = FALSE)
    } else {
      next
    }
  }
  
  if (nrow(sf_pts) == 0) next

  
  center_pt <- sf_pts %>%
    st_geometry() %>%
    st_union() %>%
    st_centroid() %>%
    st_as_sf()
  
  st_crs(center_pt) <- st_crs(sf_pts)
  
  
  bbox_poly <- st_as_sfc(st_bbox(sf_pts))
  st_crs(bbox_poly) <- st_crs(sf_pts)
  bbox_poly <- st_as_sf(bbox_poly)
  base <- sub(".rds", "", basename(f))
  center_pt$site <- base
  center_pt$num_tree <- nrow(obj)
  bbox_poly$site  <- base
  all_centroids[[base]] <- center_pt
  all_bboxes[[base]]    <- bbox_poly
}


all_centroids_sf <- do.call(rbind, all_centroids)
all_bboxes_sf    <- do.call(rbind, all_bboxes)
write_rds(all_centroids_sf, file.path(output_dir, "inventory_centroids.rds"))
write_rds(all_bboxes_sf,    file.path(output_dir, "inventory_bboxes.rds"))

