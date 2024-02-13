checkgeometry <- function(x, target = "POLYGON"){
  
  if (!target %in% c("POLYGON", "POINT")) stop("Target geometry must be POLYGON or POINT")
  
  # if we're dealing with a list, only perform on the elements that need it
  if (inherits(x, "list")){
    
    x <- lapply(x,
                function(y){
                  if (sf::st_geometry_type(y, by_geometry = FALSE) == target){
                    return(sf::st_make_valid(y))  # if tile already has the right geometry, do nothing but validate
                  } else {  # otherwise fix geometry and cast to target (multi to single part)
                    y <- sf::st_make_valid(y) %>%
                      # remove geometries that are valid but incorrect type
                      # (e.g line and points from polygon clipping)
                      # Note: this WILL lead to small holes in some places
                      dplyr::filter(sf::st_is(., if(target == "POLYGON"){
                        c("POLYGON","MULTIPOLYGON")} else {c("POINT","MULTIPOINT")})) %>%
                      sf::st_cast(to = ifelse(target == "POLYGON", "MULTIPOLYGON", "MULTIPOINT"), warn = FALSE) %>%
                      sf::st_cast(to = target, warn = FALSE)
                    return(y)
                  }
                })
    
  } else if(inherits(x, "sf")){
    if (sf::st_geometry_type(x, by_geometry = FALSE) != target){
      x <- sf::st_make_valid(x) %>%
        # remove geometries that are valid but incorrect type (e.g line and points from polygon clipping)
        dplyr::filter(sf::st_is(., if(target == "POLYGON"){
          c("POLYGON","MULTIPOLYGON")} else {c("POINT","MULTIPOINT")})) %>%
        sf::st_cast(to = ifelse(target == "POLYGON", "MULTIPOLYGON", "MULTIPOINT"), warn =FALSE) %>%
        sf::st_cast(to = target, warn = FALSE)} else {
          x <- sf::st_make_valid(x)
        }
  }
  return(x)
}