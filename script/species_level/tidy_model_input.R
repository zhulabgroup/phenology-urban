prepare_data <- function(phe_type, points_in_buffer, all_sites_temp, srad_daily) {
  points_in_buffer_joined <- points_in_buffer %>%
    left_join((srad_daily %>%
                 select(id, year, season, Avg_srad) %>%
                 pivot_wider(names_from = season, values_from = Avg_srad, names_prefix = "srad_")), 
              by = c("id", "year")) %>%
    left_join(all_sites_temp %>%
                select(name, year, season, Avg_AvgTemp) %>%
                pivot_wider(names_from = season, values_from = Avg_AvgTemp, names_prefix = "temp_"), 
              by = c("Site" = "name", "year")) %>%
    left_join(all_sites_temp %>%
                select(name, year, season, Sum_Sum_mm) %>%
                pivot_wider(names_from = season, values_from = Sum_Sum_mm, names_prefix = "prcp_"), 
              by = c("Site" = "name", "year"))
  
  if (phe_type == "SOS") {
    phe_metric <- "doy"
    direction <- "up"
    thres <- 0.5
    model_input <- points_in_buffer_joined %>%
      filter(direction == !!direction & thres == !!thres) %>%
      select(id, taxa, lon, lat, year, doy, Site, geometry, ends_with("winter")) %>%
      drop_na()
  } else if (phe_type == "EOS") {
    phe_metric <- "doy"
    direction <- "down"
    thres <- 0.5
    model_input <- points_in_buffer_joined %>%
      filter(direction == !!direction & thres == !!thres) %>%
      select(id, taxa, lon, lat, year, doy, Site, geometry, ends_with("summer")) %>%
      drop_na()
  } else if (phe_type == "GreenUp") {
    phe_metric <- "doy_diff"
    direction <- "up"
    diff <- points_in_buffer_joined %>%
      filter(direction == !!direction & thres %in% c(0.2, 0.8)) %>%
      group_by(id, year) %>%
      summarise(
        doy_diff = abs(diff(doy[order(thres)])),
        .groups = "drop"
      )
    model_input = points_in_buffer_joined %>%
      filter(direction == !!direction & thres == 0.2) %>%
      left_join(diff, by = c("year", "id")) %>%
      select(id, taxa, lon, lat, year, Site, doy_diff, geometry, ends_with("spring")) %>%
      drop_na()
    
  } else {
    direction <- "down"
    diff <- points_in_buffer_joined %>%
      filter(direction == !!direction & thres %in% c(0.2, 0.8)) %>%
      group_by(id, year) %>%
      summarise(
        doy_diff = abs(diff(doy[order(thres)])),
        .groups = "drop"
      )
    model_input = points_in_buffer_joined %>%
      filter(direction == !!direction & thres == 0.2) %>%
      left_join(diff, by = c("year", "id")) %>%
      select(id, taxa, lon, lat, year, Site, doy_diff, geometry, ends_with("fall")) %>%
      drop_na()
  }
  
  combine_exact_fuzzy <- function(df, threshold) {
    exact <- df |> filter(matchtype == "EXACT") |>
      drop_na(specieskey) # should have the exact species information
    fuzzy <- df |>
      filter(matchtype == "FUZZY", !(original_sciname %in% exact$original_sciname)) |>
      drop_na(specieskey) |>
      group_by(original_sciname) |> 
      arrange(desc(status == "ACCEPTED"), desc(confidence)) |> # sort by status, confidence and similarity
      slice(1) |>
      ungroup() |>
      add_count(original_sciname, name = "n_dup") |> # check remaining duplicates
      filter(n_dup == 1) |> # keep only unique ones, since we can not make final desicion
      select(-n_dup) |>
      filter(confidence >= threshold * 100)
    bind_rows(exact, fuzzy) |> 
      group_by(original_sciname) |>
      arrange(desc(status == "ACCEPTED")) |> 
      slice(1)
  }
  
  gbif_taxon_keys <- model_input$taxa %>% # get the first column of the data.frame
    unique()  %>% # remove duplicates
    get_gbifid_(method="backbone") %>%  # match names to the GBIF backbone to get taxonkeys
    imap(~ .x %>% mutate(original_sciname = .y))  %>%   # add original name back into data.frame
    bind_rows()  %>% 
    filter(status %in% c("ACCEPTED", "SYNONYM") & kingdom == "Plantae" & matchtype != "HIGHERRANK") %>% # get only accepted names and remove anything that might have matched to a non-plant
    combine_exact_fuzzy(threshold = 0.9) %>%
    select(canonicalname,original_sciname) # check for missing values
  
  top_10_species <- model_input %>%
    left_join(gbif_taxon_keys, by = c("taxa" = "original_sciname")) %>% # join the GBIF taxon keys
    filter(canonicalname %in% (model_input %>%
                                 left_join(gbif_taxon_keys, by = c("taxa" = "original_sciname")) %>%
                                 drop_na(canonicalname) %>% 
                                 group_by(canonicalname) %>%
                                 summarise(n = n(), .groups = "drop") %>%
                                 arrange(desc(n)) %>%
                                 # head(10) %>%
                                 pull(canonicalname))) # filter to the taxa in the original data)
  
  return(top_10_species)
}