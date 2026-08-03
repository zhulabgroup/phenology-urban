###################################### Clean TRY #################################

trait_data_raw <- rtry_import(here("data/input/tree_list/summary/TRY_street_tree.txt"))|>
  rtry_remove_dup() |>
  rtry_select_row(TraitID == 37) |>
  rtry_select_col(AccSpeciesName, OriglName, OrigValueStr, OrigUnitStr, Comment)

phe_type <- trait_data_raw |>
  filter(!str_detect(OriglName, regex("shoot phenology", ignore_case = TRUE))) |>
  filter(OrigUnitStr != "days")

phe_type$OrigValueStr <- tolower(phe_type$OrigValueStr) # make lowercase

comment_subset <- phe_type |>
  filter(!is.na(Comment) & Comment != "") |>
  mutate(type = case_when(
    str_to_lower(OrigValueStr) == "semi" ~ "semi",
    str_to_lower(OrigValueStr) == "sd" ~ "semi deciduous",
    str_to_lower(OrigValueStr) == "d" ~ "deciduous",
    str_to_lower(OrigValueStr) == "e" ~ "evergreen",
    str_to_lower(OrigValueStr) == "dc" ~ "deciduous",
    str_to_lower(OrigValueStr) == "ev" ~ "evergreen",
    str_to_lower(OrigValueStr) == "bd" ~ "deciduous broadleaved",
    
    str_to_lower(OrigValueStr) == "w" &  str_detect(Comment, "W=wintergreen") ~ "wintergreen",
    str_to_lower(OrigValueStr) == "w" &  str_detect(Comment, "W=winter deciduous") ~ "deciduous",
    
    # keep origlvalue
    str_detect(Comment, "native to site") |
      str_detect(Comment, "3\\.22") | 
      str_detect(Comment, "based on TRY") |
      str_detect(Comment, "Deciduous D / evergreen E / semi-deciduous SD") |
      str_detect(Comment, "Deciduousness as adult") |
      (Comment == "Evergreen or Deciduous") |
      (Comment == "evergreen, winter deciduous, drought deciduous, semideciduous") |
      str_detect(Comment, "Indicator of deciduousness of individual") |
      (Comment == "leaf phenology") |
      (Comment == "Leaf phenology") |
      str_detect(Comment, "PFT") |
      str_detect(Comment, "from the TRY database") ~ OrigValueStr,
    
    # 0/1 evergreen vs deciduous
    str_detect(Comment, "0 = no, 1 = yes") & OriglName == "evergreen" ~ if_else(OrigValueStr == "1", "evergreen", "deciduous"),
    str_detect(Comment, "0 = no, 1 = yes") & OriglName == "deciduous" ~ if_else(OrigValueStr == "1", "deciduous", "evergreen"),
    
    # Phenological type: evergreen: Y deciduous: N
    str_detect(Comment, "Phenological type: evergreen: Y") ~ if_else(
      str_to_lower(OrigValueStr) == "y", "evergreen", "deciduous"
    ),
    
    # 0= evergreen coniferous, 1= deciduous broadleaved, 2= deciduous coniferous
    str_detect(Comment, "0= evergreen coniferous") ~ case_when(
      OrigValueStr == 0 ~ "evergreen",
      OrigValueStr == 1 ~ "deciduous",
      OrigValueStr == 2 ~ "deciduous",
      TRUE ~ NA_character_
    ),
    
    # 1= aestival green; 2= partial evergreen; 3= evergreen
    str_detect(Comment, "1= aestival green") ~ case_when(
      OrigValueStr == 1 ~ "deciduous",
      OrigValueStr == 2 ~ "semi evergreen",
      OrigValueStr == 3 ~ "evergreen",
      TRUE ~ NA_character_
    ),
    TRUE ~ NA_character_
  )
  ) |>
  filter(!type %in% c("variable deciduous", "deciduous/evergreen", "nonevergreen")) |>
  filter(str_detect(type, "deciduous") | str_detect(type, "evergreen") |
           str_detect(type, "aestival") | str_detect(type, "vernal") | str_detect(type, "hibernal")) |>
  mutate(type = if_else(type == "deciduous broadleaved", "deciduous", type),
         type = if_else(type == "partly evergreen", "semi-evergreen", type),
         type = if_else(type == "semi deciduous", "semi-deciduous", type))

non_comment_subset <- phe_type |>
  filter(Comment == "") |>
  mutate(OriglName = str_to_lower(OriglName)) |>
  mutate(type = case_when(
    str_to_lower(OrigValueStr) == "semi" ~ "semi",
    str_to_lower(OrigValueStr) == "s" ~ "semi",
    str_to_lower(OrigValueStr) == "sd" ~ "semi deciduous",
    str_to_lower(OrigValueStr) == "sg" ~ "semi evergreen",
    str_to_lower(OrigValueStr) == "d" ~ "deciduous",
    str_to_lower(OrigValueStr) == "e" ~ "evergreen",
    str_to_lower(OrigValueStr) == "dc" ~ "deciduous",
    str_to_lower(OrigValueStr) == "ev" ~ "evergreen",
    
    # Phenological type: evergreen: Y deciduous: N
    (OriglName == "deciduous") ~ if_else(
      str_to_lower(OrigValueStr) %in% c("y","yes"),"deciduous", "-"
    ),
    (OriglName == "evergreen") ~ case_when(
      str_to_lower(OrigValueStr) %in% c("y", "yes") ~ "evergreen",
      str_to_lower(OrigValueStr) %in% c("n", "no")  ~ "-",
      TRUE ~ str_to_lower(OrigValueStr)
    ),
    (OriglName == "leaf phenology: deciduous") ~ if_else(
      str_to_lower(OrigValueStr) %in% c("y","yes"),"deciduous", "-"
    ),
    (OriglName == "leaf phenology: evergreen") ~ if_else(
      str_to_lower(OrigValueStr) %in% c("y","yes"),"evergreen", "-" 
    ),
    (OriglName == "leaf phenology: leaf exchanger") ~ if_else(
      str_to_lower(OrigValueStr) %in% c("y","yes"),"exchanger", "-"
    ),
    (OriglName == "leaf phenology: semi-deciduous") ~ if_else(
      str_to_lower(OrigValueStr) %in% c("y","yes"),"semi-deciduous", "-" 
    ),
    TRUE ~ OrigValueStr
  )
  ) |>
  filter(!str_detect(type, "http")) |>
  filter(!str_detect(type, "deciduous or evergreen")) |>
  filter(!type %in% c("y/n", "", "-", "aphyllous", "d_ev","cladophylls", "deciduous/evergreen", "evergreen/deciduous", "leafless","other", 
                      "phyllodium", "seasonal biomass peak", "stem succulent")) |>
  mutate(type = if_else(type == "always overwintering green", "hibernal", type),
         type = if_else(type == "always persistent green", "evergreen", type),
         type = if_else(type == "always spring green", "vernal", type),
         type = if_else(type == "always summer green", "aestival", type),
         type = if_else(type == "caducifolio", "deciduous", type),
         type = if_else(type %in% c("deciduous broad-leaved", "deciduous needle-leaved", "cold deciduous", "deciduous scale-like", "deciduous/needlelike",
                                    "decioduous", "decuduous", "drought deciduous", "drought-deciduous", "winter deciduous"), "deciduous", type),
         type = if_else(type %in% c("evergeen", "everggreen", "evergreen broad-leaved", "evergreen grass", "evergreen needle-leaved", "evergreen scale-like",
                                    "evergren", "perennifolio"), "evergreen", type),
         type = if_else(type ==  "winter semi-deciduous", "semi-deciduous", type),
         type = if_else(type == "semi-evergreen/drought-deciduous", "semi-evergreen/deciduous", type),
         type = if_else(type == "deciduous or semi-evergreen broad-leaved", "deciduous/semi-evergreen", type),
         type = if_else(type == "drought semi-deciduous", "semi-deciduous", type),
         type = if_else(type == "deciduous-semideciduous", "deciduous/semi-deciduous", type),
         type = if_else(type == "semidecicuous", "semi-deciduous", type)) |>
  mutate(
    type = str_replace_all(type, "semi ", "semi-"),
    type = str_replace_all(type, "evergren", "evergreen"),
    type = str_replace_all(type, "semievergreen", "semi-evergreen"), 
    type = str_replace_all(type, "semideciduous", "semi-deciduous"),
    type = str_replace_all(type, " to ", "/")) |>
  mutate(type = if_else(type == "semi-evergreen/deciduous", "deciduous/semi-evergreen", type),
         type = if_else(type == "semi-evergreen/evergreen", "evergreen/semi-evergreen", type),
         type = if_else(type == "semi-deciduous/evergreen", "evergreen/semi-deciduous", type))

phe_type_clean <- bind_rows(comment_subset, non_comment_subset)

gbif_taxon_keys <- phe_type_clean$AccSpeciesName|> # get the first column of the data.frame
  unique() |> # remove duplicates
  name_backbone_checklist(kingdom = "Plantae")|>
  filter(kingdom == "Plantae") |>
  filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
  mutate(canonicalName = if_else(status == "SYNONYM", species, canonicalName)) |>
  select(canonicalName, verbatim_name) # check for missing values

phe_type_clean_gbif <- phe_type_clean |>
  left_join(gbif_taxon_keys, by = c("AccSpeciesName" = "verbatim_name")) |>
  filter(!is.na(canonicalName)) |>
  group_by(canonicalName, type)|>
  summarise(num_obs = n(), .groups = "drop")

write_rds(phe_type_clean_gbif, here("data/input/tree_list/summary/species_phe_types.rds"))
