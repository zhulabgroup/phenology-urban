library(dplyr)
library(sp)
library(sf)
library(foreach)
library(progress)
library(doSNOW)
library(readr)

faster_intersect <- function(x, cliplayer){
  
  x_geo <- sf::st_geometry(x)
  number_p <- nrow(cliplayer)
  cliplayer_list <- split(cliplayer, seq(nrow(cliplayer)))
  
  cl <- makeCluster(23, outfile = "")
  registerDoSNOW(cl)
  intersection <- foreach(cliplayer_row = cliplayer_list) %dopar% {
    test <- sf::st_intersects(x_geo, sf::st_geometry(cliplayer_row)) # test if polygons overlap clipping layer
    print(cliplayer_row$id)
    # if result is empty (o: FALSE), out of study area, we delete
    x_inter <- x[lengths(test) > 0,]  # exclude polygons that don't even intersect
    if (nrow(x_inter) > 0) { # only if there are polygons left:
      # now we check which polygons are fully within as we don't need to clip them at all
      test2 <- sf::st_covered_by(x_inter, sf::st_geometry(cliplayer_row))
      # and we intersect the rest
      x_cut <- x_inter[lengths(test2) == 0,]  # extract objects to cut
      x_within <- x_inter[lengths(test2) != 0,]
      x_within <- dplyr::mutate(x_within,
                                id = cliplayer_row$id,
                                taxa = cliplayer_row$taxa,
                                site = cliplayer_row$site)
      # object within, but need to give attribute
      # perform the intersection
      cut2 <- suppressWarnings(
        sf::st_intersection(x_cut, cliplayer_row))
      cut_result <- sf::st_make_valid(cut2)
      x_intersect <- rbind(cut_result, x_within)
      x_intersect$shape_area <- sf::st_area(x_intersect)
      }
    
    as.data.frame(x_intersect)
    x_intersect
    
  }
  intersection_results <- bind_rows(intersection)
  stopCluster(cl)
  
  write_rds(intersection_results, "result/temp/Morus_intersect_200.rds")
}


