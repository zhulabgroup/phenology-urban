library(tidyverse)
library(readr)
library(rgbif)
library(rtry)
library(here)
library(sf)

################################# Rules #################################
# Clean species name
# Keep deciduous species
# Reduce samples to 2000/species
# Unique location'

# Format: id | taxa | lat | lon | site(city_short) | genus

save_dir <- "data/raw/CONUS/StreetTrees/cleaned_ps_input/"

###################################### Clean TRY #################################

check <- read_rds(here(paste0(save_dir, "summary/species_phe_types.rds"))) |>
  filter(!type %in% c("hibernal", "vernal", "aestival")) |>
  group_by(canonicalName) |>
  arrange(desc(num_obs), type, .by_group = TRUE) |>
  summarise(
    types = paste(unique(type), collapse = ", "),
    types_with_counts = if (n_distinct(type) > 1) {
      paste0(type, " (", num_obs, ")", collapse = ", ")
    } else {
      NA_character_
    },
    single_type = first(type),
    .groups = "drop"
  )

################################# NY DT DV HT ST TP AT #################################
# Remove NY DV HT AT ST since they come from OpenTree  

metadata <- read_csv("~/lab-data/datasets/vegetation/PS/urban/metadata.csv") |>
  mutate(genus = str_extract(taxa, "^[^ ]+")) |>
  distinct(taxa, genus, lat, lon, site)

gbif_taxon_keys <- metadata$taxa|> # get the first column of the data.frame
  unique() |> # remove duplicates
  name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
  filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
  filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
  mutate(canonicalName = if_else(status == "SYNONYM", species, canonicalName)) |>
  select(canonicalName, verbatim_name) # check for missing values

# Join to TRY data and find the deciduous species

set.seed(123) # for reproducibility

metadata_clean_species <- metadata|>
  left_join(gbif_taxon_keys, by = c("taxa" = "verbatim_name"))|> # join the GBIF taxon keys
  drop_na(canonicalName)|> 
  rename(canonicalname = canonicalName) |>
  group_by(canonicalname, site) |>
  filter(n() >= 30) |>
  group_modify(~ {
    if (nrow(.x) > 2000) {
      .x |> slice_sample(n = 2000)
    } else {
      .x                          
    }
  }) |>
  ungroup() |>
  left_join(check, by = c("canonicalname" = "canonicalName")) |>
  filter(str_detect(single_type, "deciduous")) |>
  select(-taxa, -types_with_counts, -types, -single_type) |>
  rename(taxa = canonicalname)

metadata_clean_species |>
  group_split(site) |>  
  purrr::walk(function(df) {
    site_name <- unique(df$site)
    df <- df |>
      mutate(id = row_number())
    write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
  })

################################# McCoy's #################################
# Remove duplicate: Detroit, Tampa

elife_trees_clean <- readRDS(here("data/raw/CONUS/StreetTrees/elife/elife_trees_clean.rds")) |>
  filter(!is.na(NAME20)) |>
  st_transform(4326) |>
  filter(!condition %in% c("dead/dying", "dead")) |>
  select(scientific_name, condition, site = city_combine.x, state = state_abb.x) |>
  mutate(lon = st_coordinates(geometry)[, 1],
         lat = st_coordinates(geometry)[, 2]) |>
  st_drop_geometry()

gbif_taxon_keys <- elife_trees_clean$scientific_name|> # get the first column of the data.frame
  unique() |> # remove duplicates
  name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
  filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
  filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
  mutate(canonicalName = if_else(status == "SYNONYM", species, canonicalName)) |>
  select(canonicalName, verbatim_name) # check for missing values

# Join to TRY data and find the deciduous species

set.seed(123) # for reproducibility

elife_trees_clean_phe <- elife_trees_clean|>
  left_join(gbif_taxon_keys, by = c("scientific_name" = "verbatim_name"))|> # join the GBIF taxon keys
  drop_na(canonicalName)|> 
  rename(canonicalname = canonicalName) |>
  group_by(canonicalname, site) |>
  filter(n() >= 30) |>
  group_modify(~ {
    if (nrow(.x) > 2000) {
      .x |> slice_sample(n = 2000)
    } else {
      .x                          
    }
  }) |>
  ungroup() |>
  left_join(check, by = c("canonicalname" = "canonicalName")) |>
  filter(str_detect(single_type, "deciduous")) |>
  select(-scientific_name, -types_with_counts, -types, -single_type) |>
  mutate(genus = str_extract(canonicalname, "^[^ ]+")) |>
  rename(taxa = canonicalname)

elife_trees_clean_phe |>
  group_split(site) |>  
  purrr::walk(function(df) {
    site_name <- gsub("\\.", "_", unique(df$site))
    df <- df |>
      mutate(id = row_number())
    write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
  })

################################# OpenTree.org #################################

# Allentown 
Allentown <- FALSE
if(Allentown){
  raw_Allentown <- st_read(here("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Allentown/City_Trees.shp")) |>
    st_transform(4326) |>
    mutate(lon = st_coordinates(geometry)[, 1],
           lat = st_coordinates(geometry)[, 2],
           site = "Allentown") |>
    select(TR_GENUS, TR_COMMON, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_Allentown$TR_GENUS|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  # gbif_taxon_keys_common <- raw_Allentown |> # get the first column of the data.frame
  #   left_join(gbif_taxon_keys, by = c("TR_GENUS" = "verbatim_name"))|> # join the GBIF taxon keys
  #   distinct(TR_COMMON, canonicalName) |>
  #   filter(is.na(canonicalName)) |>
  #   pull(TR_COMMON) |>
  #   comm2sci()
  # 
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility

  clean_Allentown <- raw_Allentown|>
    left_join(gbif_taxon_keys, by = c("TR_GENUS" = "verbatim_name"))|> # join the GBIF taxon keys
    mutate(canonicalName = if_else((is.na(canonicalName) & TR_COMMON == "Apple"), "Malus domestica", canonicalName),
           canonicalName = if_else((is.na(canonicalName) & TR_COMMON == "Sugar Maple"), "Acer saccharum", canonicalName),
           canonicalName = if_else((is.na(canonicalName) & TR_COMMON == "Red Maple"), "Acer rubrum", canonicalName)) |>
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    st_drop_geometry() |>
    select(-geometry, -TR_GENUS, -TR_COMMON, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Allentown |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Amherst
Amherst <- FALSE
if(Amherst){
  raw_Amherst <- st_read(here("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Amherst/Street_Trees.shp")) |>
    st_transform(4326) |>
    mutate(lon = st_coordinates(geometry)[, 1],
           lat = st_coordinates(geometry)[, 2],
           site = "Amherst",
           Species = str_replace(tolower(Species), "^(.)", toupper)) |>
    filter(!is.na(Species)) |>
    filter(Species != "_") |>
    st_drop_geometry() |>
    select(Species, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_Amherst$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  # If there is common name
  gbif_taxon_keys_common <- raw_Amherst |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(Species, canonicalName) |>
    filter(is.na(canonicalName)) |>
    pull(Species)
  comm2sci_list <- map(
    gbif_taxon_keys_common,
    function(x) {
      Sys.sleep(1) 
      tryCatch(
        comm2sci(x), 
        error = function(e) NULL  
      )
    }
  )
  comm2sci <- imap_dfr(comm2sci_list, function(inner_list, outer_name) {
    if (length(inner_list) == 0) {
      tibble(common_group = outer_name,
             common_name  = NA_character_,
             scientific_name = NA_character_)
    } else {
      map_dfr(names(inner_list), function(comm_name) {
        sci_names <- inner_list[[comm_name]]
        if (length(sci_names) == 0) sci_names <- NA_character_
        tibble(common_group = outer_name,
               common_name = comm_name,
               scientific_name = sci_names)
      })
    }
  }) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = common_name)
  
  gbif_taxon_keys <- gbif_taxon_keys |>
    left_join(comm2sci, by = "verbatim_name") |>
    mutate(
      canonicalName = if_else(!is.na(scientific_name), scientific_name, canonicalName)
    ) |>
    select(-scientific_name)
  
  
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility
  
  clean_Amherst <- raw_Amherst|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Species, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Amherst |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Boulder
Boulder <- FALSE
if(Boulder){
  raw_Boulder <- st_read(here("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Boulder/Tree_Inventory_Open_Data.shp")) |>
    st_transform(4326) |>
    mutate(lon = st_coordinates(geometry)[, 1],
           lat = st_coordinates(geometry)[, 2],
           site = "Boulder",
           Species = LATINNAME) |>
    filter(!is.na(Species)) |>
    st_drop_geometry() |>
    select(Species, COMMONNAME, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_Boulder$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  # If there is common name
  gbif_taxon_keys_common <- raw_Boulder |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(Species,COMMONNAME, canonicalName) |>
    filter(is.na(canonicalName)) |>
    pull(COMMONNAME)
  gbif_taxon_keys_common_swapped <- str_replace(
    gbif_taxon_keys_common,
    "^([^,]+),\\s*(.+)$",   # 逗号前后分组
    "\\2 \\1"              # 交换顺序
  )
  comm2sci_list <- map(
    gbif_taxon_keys_common_swapped,
    function(x) {
      Sys.sleep(1) 
      tryCatch(
        comm2sci(x), 
        error = function(e) NULL  
      )
    }
  )
  comm2sci <- imap_dfr(comm2sci_list, function(inner_list, outer_name) {
    if (length(inner_list) == 0) {
      tibble(common_group = outer_name,
             common_name  = NA_character_,
             scientific_name = NA_character_)
    } else {
      map_dfr(names(inner_list), function(comm_name) {
        sci_names <- inner_list[[comm_name]]
        if (length(sci_names) == 0) sci_names <- NA_character_
        tibble(common_group = outer_name,
               common_name = comm_name,
               scientific_name = sci_names)
      })
    }
  }) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = common_name)
  
  common_name <-  raw_Boulder |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(Species,COMMONNAME, canonicalName) |>
    filter(is.na(canonicalName)) |>
    mutate(new_common = str_replace(
      COMMONNAME,
      "^([^,]+),\\s*(.+)$",   # 逗号前后分组
      "\\2 \\1"              # 交换顺序
    )) |>
    left_join(comm2sci, by = c("new_common" = "verbatim_name")) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = COMMONNAME)
  
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility
  
  clean_Boulder <- raw_Boulder|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    left_join(common_name, by = c("COMMONNAME" = "verbatim_name")) |>
    mutate(canonicalName = if_else(is.na(canonicalName), scientific_name, canonicalName)) |>
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-COMMONNAME, -scientific_name, -Species, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Boulder |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Bozeman
Bozeman <- FALSE
if(Bozeman){
  raw_Bozeman <- read_csv(here("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Bozeman.csv")) |>
    select(lon = X,
           lat = Y,
           genus = Genus,
           taxa = Botanical_Name,
           Common_Name) |>
    filter(!(is.na(genus) & is.na(taxa) & is.na(Common_Name))) |>
    distinct()
  
  gbif_taxon_keys <- raw_Bozeman$taxa|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  # If there is common name
  gbif_taxon_keys_common <- raw_Bozeman |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("taxa" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(taxa,Common_Name, canonicalName) |>
    filter(is.na(canonicalName)) |>
    pull(Common_Name)
  gbif_taxon_keys_common_swapped <- str_replace(
    gbif_taxon_keys_common,
    "^([^,]+),\\s*(.+)$",   # 逗号前后分组
    "\\2 \\1"              # 交换顺序
  )
  comm2sci_list <- map(
    gbif_taxon_keys_common_swapped,
    function(x) {
      Sys.sleep(1) 
      tryCatch(
        comm2sci(x), 
        error = function(e) NULL  
      )
    }
  )
  comm2sci <- imap_dfr(comm2sci_list, function(inner_list, outer_name) {
    if (length(inner_list) == 0) {
      tibble(common_group = outer_name,
             common_name  = NA_character_,
             scientific_name = NA_character_)
    } else {
      map_dfr(names(inner_list), function(comm_name) {
        sci_names <- inner_list[[comm_name]]
        if (length(sci_names) == 0) sci_names <- NA_character_
        tibble(common_group = outer_name,
               common_name = comm_name,
               scientific_name = sci_names)
      })
    }
  }) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = common_name) |>
    distinct()
  
  common_name <-  raw_Bozeman |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("taxa" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(taxa,Common_Name, canonicalName) |>
    filter(is.na(canonicalName)) |>
    mutate(new_common = str_replace(
      Common_Name,
      "^([^,]+),\\s*(.+)$",   # 逗号前后分组
      "\\2 \\1"              # 交换顺序
    )) |>
    left_join(comm2sci, by = c("new_common" = "verbatim_name")) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = Common_Name) |>
    distinct()
  
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility
  
  clean_Bozeman <- raw_Bozeman|>
    left_join(gbif_taxon_keys, by = c("taxa" = "verbatim_name"))|> # join the GBIF taxon keys
    left_join(common_name, by = c("Common_Name" = "verbatim_name")) |>
    mutate(canonicalName = if_else(is.na(canonicalName), scientific_name, canonicalName),
           site = "Bozeman") |>
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Common_Name, -scientific_name, -taxa, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Bozeman |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Champaign
Champaign <- FALSE
if(Champaign){
  raw_Champaign <- st_read(here("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Champaign/City_Owned_Trees.shp")) |>
    st_transform(4326) |>
    mutate(lon = st_coordinates(geometry)[, 1],
           lat = st_coordinates(geometry)[, 2],
           site = "Champaign",
           Species = SPP) |>
    st_drop_geometry() |>
    filter(!(is.na(SPP) & is.na(COMMON) )) |>
    select(Species, COMMON, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_Champaign$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  # If there is common name
  gbif_taxon_keys_common <- raw_Champaign |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(Species,COMMON, canonicalName) |>
    filter(is.na(canonicalName)) |>
    pull(COMMON)
  gbif_taxon_keys_common_swapped <- str_replace(
    gbif_taxon_keys_common,
    "^([^,]+),\\s*(.+)$",   # 逗号前后分组
    "\\2 \\1"              # 交换顺序
  )
  comm2sci_list <- map(
    gbif_taxon_keys_common_swapped,
    function(x) {
      Sys.sleep(1) 
      tryCatch(
        comm2sci(x), 
        error = function(e) NULL  
      )
    }
  )
  comm2sci <- imap_dfr(comm2sci_list, function(inner_list, outer_name) {
    if (length(inner_list) == 0) {
      tibble(common_group = outer_name,
             common_name  = NA_character_,
             scientific_name = NA_character_)
    } else {
      map_dfr(names(inner_list), function(comm_name) {
        sci_names <- inner_list[[comm_name]]
        if (length(sci_names) == 0) sci_names <- NA_character_
        tibble(common_group = outer_name,
               common_name = comm_name,
               scientific_name = sci_names)
      })
    }
  }) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = common_name) |>
    distinct()
  
  common_name <-  raw_Champaign |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(Species,COMMON, canonicalName) |>
    filter(is.na(canonicalName)) |>
    mutate(new_common = str_replace(
      COMMON,
      "^([^,]+),\\s*(.+)$",   # 逗号前后分组
      "\\2 \\1"              # 交换顺序
    )) |>
    left_join(comm2sci, by = c("new_common" = "verbatim_name")) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = COMMON) |>
    distinct()
  
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility
  
  clean_Champaign <- raw_Champaign|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    left_join(common_name, by = c("COMMON" = "verbatim_name")) |>
    mutate(canonicalName = if_else(is.na(canonicalName), scientific_name, canonicalName)) |>
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-COMMON, -Species, -scientific_name, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Champaign |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Charlottesville
Charlottesville <- FALSE
if(Charlottesville){
  raw_Charlottesville <- st_read(here("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Charlottesville/Tree_Inventory_Point.shp")) |>
    st_transform(4326) |>
    mutate(lon = st_coordinates(geometry)[, 1],
           lat = st_coordinates(geometry)[, 2],
           site = "Charlottesville",
           Species = str_c(Genus, Species, sep = " ")) |>
    st_drop_geometry() |>
    filter(!(is.na(Species) & is.na(Common_Nam) )) |>
    select(Species, Common_Nam, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_Charlottesville$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  # If there is common name
  gbif_taxon_keys_common <- raw_Charlottesville |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(Species,Common_Nam, canonicalName) |>
    filter(is.na(canonicalName)) |>
    pull(Common_Nam)
  
  comm2sci_list <- map(
    gbif_taxon_keys_common,
    function(x) {
      Sys.sleep(1) 
      tryCatch(
        comm2sci(x), 
        error = function(e) NULL  
      )
    }
  )
  comm2sci <- imap_dfr(comm2sci_list, function(inner_list, outer_name) {
    if (length(inner_list) == 0) {
      tibble(common_group = outer_name,
             common_name  = NA_character_,
             scientific_name = NA_character_)
    } else {
      map_dfr(names(inner_list), function(comm_name) {
        sci_names <- inner_list[[comm_name]]
        if (length(sci_names) == 0) sci_names <- NA_character_
        tibble(common_group = outer_name,
               common_name = comm_name,
               scientific_name = sci_names)
      })
    }
  }) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = common_name) |>
    distinct()
  
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility
  
  clean_Charlottesville <- raw_Charlottesville|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    left_join(comm2sci, by = c("Common_Nam" = "verbatim_name")) |>
    mutate(canonicalName = if_else(is.na(canonicalName), scientific_name, canonicalName)) |>
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Common_Nam, -Species, -scientific_name, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Charlottesville |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Philadelphia
Philadelphia <- FALSE
if(Philadelphia){
  raw_Philadelphia <- st_read(here("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Philadelphia/ppr_tree_inventory_2024.shp")) |>
    st_transform(4326) |>
    mutate(lon = st_coordinates(geometry)[, 1],
           lat = st_coordinates(geometry)[, 2],
           site = "Philadelphia") |>
    separate(tree_name, into = c("Species", "Common_Nam"), sep = " - | – \\s*", remove = FALSE) |>
    st_drop_geometry() |>
    filter(!(is.na(Species) & is.na(Common_Nam) )) |>
    select(Species, Common_Nam, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_Philadelphia$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  # If there is common name
  gbif_taxon_keys_common <- raw_Philadelphia |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(Species,Common_Nam, canonicalName) |>
    filter(is.na(canonicalName)) |>
    pull(Common_Nam) |>
    tolower()
  comm2sci_list <- map(
    gbif_taxon_keys_common,
    function(x) {
      Sys.sleep(1) 
      tryCatch(
        comm2sci(x), 
        error = function(e) NULL  
      )
    }
  )
  comm2sci <- imap_dfr(comm2sci_list, function(inner_list, outer_name) {
    if (length(inner_list) == 0) {
      tibble(common_group = outer_name,
             common_name  = NA_character_,
             scientific_name = NA_character_)
    } else {
      map_dfr(names(inner_list), function(comm_name) {
        sci_names <- inner_list[[comm_name]]
        if (length(sci_names) == 0) sci_names <- NA_character_
        tibble(common_group = outer_name,
               common_name = comm_name,
               scientific_name = sci_names)
      })
    }
  }) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = common_name) |>
    distinct()
  comm2sci$verbatim_name <- toupper(comm2sci$verbatim_name)
  
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility
  
  clean_Philadelphia <- raw_Philadelphia|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    left_join(comm2sci, by = c("Common_Nam" = "verbatim_name")) |>
    mutate(canonicalName = if_else(is.na(canonicalName), scientific_name, canonicalName)) |>
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Common_Nam, -Species, -scientific_name, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Philadelphia |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# StAugustine
StAugustine <- FALSE
if(StAugustine){
  raw_StAugustine <- st_read(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_StAugustine/Trees.shp")) |>
    st_transform(4326) |>
    mutate(lon = st_coordinates(geometry)[, 1],
           lat = st_coordinates(geometry)[, 2],
           site = "StAugustine",
           Species = SPP) |>
    st_drop_geometry() |>
    filter(!(is.na(Species))) |>
    select(Species, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_StAugustine$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility
  
  clean_StAugustine <- raw_StAugustine|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Species, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_StAugustine |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Escondido (Opentree)
Escondido <- FALSE
if(Escondido){
  raw_Escondido <- st_read(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Escondido/Tree_Inventory.shp")) |>
    st_transform(4326) |>
    mutate(lon = st_coordinates(geometry)[, 1],
           lat = st_coordinates(geometry)[, 2],
           site = "Escondido",
           Species = BOTANICAL_,
           Common_Nam = COMMON_NAM) |>
    filter(!(is.na(Species) & is.na(Common_Nam))) |>
    select(Species, Common_Nam, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_Escondido$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  # If there is common name
  gbif_taxon_keys_common <- raw_Escondido |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(Species,Common_Nam, canonicalName) |>
    filter(is.na(canonicalName)) |>
    pull(Common_Nam) |>
    tolower()
  
  comm2sci_list <- map(
    gbif_taxon_keys_common,
    function(x) {
      Sys.sleep(1) 
      tryCatch(
        comm2sci(x), 
        error = function(e) NULL  
      )
    }
  )
  comm2sci <- imap_dfr(comm2sci_list, function(inner_list, outer_name) {
    if (length(inner_list) == 0) {
      tibble(common_group = outer_name,
             common_name  = NA_character_,
             scientific_name = NA_character_)
    } else {
      map_dfr(names(inner_list), function(comm_name) {
        sci_names <- inner_list[[comm_name]]
        if (length(sci_names) == 0) sci_names <- NA_character_
        tibble(common_group = outer_name,
               common_name = comm_name,
               scientific_name = sci_names)
      })
    }
  }) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = common_name) |>
    distinct()
  
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility
  
  clean_Escondido <- raw_Escondido|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    left_join(comm2sci, by = c("Common_Nam" = "verbatim_name")) |>
    mutate(canonicalName = if_else(is.na(canonicalName), scientific_name, canonicalName)) |>
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Common_Nam, -Species, -scientific_name, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Escondido |>
    st_drop_geometry() |>
    select(-geometry) |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Cambridge (Opentree)
Cambridge <- FALSE
if(Cambridge){
  raw_Cambridge <- st_read(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Cambridge/ENVIRONMENTAL_StreetTrees.shp")) |>
    st_transform(4326) |>
    filter(SiteType == "Tree") |>
    mutate(lon = st_coordinates(geometry)[, 1],
           lat = st_coordinates(geometry)[, 2],
           site = "Cambridge",
           Species = Scientific,
           Common_Nam = CommonName) |>
    st_drop_geometry() |>
    filter(is.na(RemovalRea))|>
    filter(!(is.na(Species) & is.na(Common_Nam))) |>
    select(Species, Common_Nam, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_Cambridge$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  # If there is common name
  gbif_taxon_keys_common <- raw_Cambridge |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(Species,Common_Nam, canonicalName) |>
    filter(is.na(canonicalName)) |>
    pull(Common_Nam)
  
  comm2sci_list <- map(
    gbif_taxon_keys_common,
    function(x) {
      Sys.sleep(1) 
      tryCatch(
        comm2sci(x), 
        error = function(e) NULL  
      )
    }
  )
  comm2sci <- imap_dfr(comm2sci_list, function(inner_list, outer_name) {
    if (length(inner_list) == 0) {
      tibble(common_group = outer_name,
             common_name  = NA_character_,
             scientific_name = NA_character_)
    } else {
      map_dfr(names(inner_list), function(comm_name) {
        sci_names <- inner_list[[comm_name]]
        if (length(sci_names) == 0) sci_names <- NA_character_
        tibble(common_group = outer_name,
               common_name = comm_name,
               scientific_name = sci_names)
      })
    }
  }) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = common_name) |>
    distinct()
  
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility
  
  clean_Cambridge <- raw_Cambridge|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    left_join(comm2sci, by = c("Common_Nam" = "verbatim_name")) |>
    mutate(canonicalName = if_else(is.na(canonicalName), scientific_name, canonicalName)) |>
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Common_Nam, -Species, -scientific_name, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Cambridge |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

################################# GUTI inferred #################################

# Boise
Boise <- FALSE
if(Boise){
  raw_Boise <- st_read(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Boise/Street_Trees.shp")) |>
    st_transform(4326) |>
    mutate(lon = st_coordinates(geometry)[, 1],
           lat = st_coordinates(geometry)[, 2],
           site = "Boise",
           Species = tolower(str_c(Genus, Species, sep = " "))) |>
    st_drop_geometry() |>
    filter(!(is.na(Species) & is.na(Common_Nam))) |>
    select(Species, Common_Nam, site, lon, lat) |>
    distinct()
  
  
  gbif_taxon_keys <- raw_Boise$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  # If there is common name
  gbif_taxon_keys_common <- raw_Boise |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(Species,Common_Nam, canonicalName) |>
    filter(is.na(canonicalName)) |>
    pull(Common_Nam)
  gbif_taxon_keys_common_swapped <- str_replace(
    gbif_taxon_keys_common,
    "^([^,]+),\\s*(.+)$",   # 逗号前后分组
    "\\2 \\1"              # 交换顺序
  )
  comm2sci_list <- map(
    gbif_taxon_keys_common_swapped,
    function(x) {
      Sys.sleep(1) 
      tryCatch(
        comm2sci(x), 
        error = function(e) NULL  
      )
    }
  )
  comm2sci <- imap_dfr(comm2sci_list, function(inner_list, outer_name) {
    if (length(inner_list) == 0) {
      tibble(common_group = outer_name,
             common_name  = NA_character_,
             scientific_name = NA_character_)
    } else {
      map_dfr(names(inner_list), function(comm_name) {
        sci_names <- inner_list[[comm_name]]
        if (length(sci_names) == 0) sci_names <- NA_character_
        tibble(common_group = outer_name,
               common_name = comm_name,
               scientific_name = sci_names)
      })
    }
  }) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = common_name) |>
    distinct()
  
  common_name <-  raw_Boise |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(Species,Common_Nam, canonicalName) |>
    filter(is.na(canonicalName)) |>
    mutate(new_common = str_replace(
      Common_Nam,
      "^([^,]+),\\s*(.+)$",   # 逗号前后分组
      "\\2 \\1"              # 交换顺序
    )) |>
    left_join(comm2sci, by = c("new_common" = "verbatim_name")) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = Common_Nam) |>
    distinct()
  
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility
  
  clean_Boise <- raw_Boise|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    left_join(common_name, by = c("Common_Nam" = "verbatim_name")) |>
    mutate(canonicalName = if_else(is.na(canonicalName), scientific_name, canonicalName)) |>
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Common_Nam, -Species, -scientific_name, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Boise |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Boston
Boston <- FALSE
if(Boston){
  raw_Boston <- st_read(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Boston/BPRD_Trees.shp")) |>
    st_transform(4326) |>
    mutate(lon = st_coordinates(geometry)[, 1],
           lat = st_coordinates(geometry)[, 2],
           site = "Boston",
           Species = spp_bot,
           Common_Nam = spp_com) |>
    st_drop_geometry() |>
    filter(!(is.na(Species) & is.na(Common_Nam))) |>
    select(Species, Common_Nam, site, lon, lat) |>
    distinct()
  
  
  gbif_taxon_keys <- raw_Boston$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  # If there is common name
  gbif_taxon_keys_common <- raw_Boston |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(Species,Common_Nam, canonicalName) |>
    filter(is.na(canonicalName)) |>
    pull(Common_Nam)
  gbif_taxon_keys_common_swapped <- str_replace(
    gbif_taxon_keys_common,
    "^([^,]+),\\s*(.+)$",   # 逗号前后分组
    "\\2 \\1"              # 交换顺序
  )
  comm2sci_list <- map(
    gbif_taxon_keys_common_swapped,
    function(x) {
      Sys.sleep(1) 
      tryCatch(
        comm2sci(x), 
        error = function(e) NULL  
      )
    }
  )
  comm2sci <- imap_dfr(comm2sci_list, function(inner_list, outer_name) {
    if (length(inner_list) == 0) {
      tibble(common_group = outer_name,
             common_name  = NA_character_,
             scientific_name = NA_character_)
    } else {
      map_dfr(names(inner_list), function(comm_name) {
        sci_names <- inner_list[[comm_name]]
        if (length(sci_names) == 0) sci_names <- NA_character_
        tibble(common_group = outer_name,
               common_name = comm_name,
               scientific_name = sci_names)
      })
    }
  }) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = common_name) |>
    distinct()
  
  common_name <-  raw_Boston |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(Species,Common_Nam, canonicalName) |>
    filter(is.na(canonicalName)) |>
    mutate(new_common = str_replace(
      Common_Nam,
      "^([^,]+),\\s*(.+)$",   # 逗号前后分组
      "\\2 \\1"              # 交换顺序
    )) |>
    left_join(comm2sci, by = c("new_common" = "verbatim_name")) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = Common_Nam) |>
    distinct()
  
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility
  
  clean_Boston <- raw_Boston|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    left_join(common_name, by = c("Common_Nam" = "verbatim_name")) |>
    mutate(canonicalName = if_else(is.na(canonicalName), scientific_name, canonicalName)) |>
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Common_Nam, -Species, -scientific_name, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Boston |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Davenport
Davenport <- FALSE
if(Davenport){
  raw_Davenport <- st_read(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Davenport/Trees.shp")) |>
    st_transform(4326) |>
    mutate(lon = st_coordinates(geometry)[, 1],
           lat = st_coordinates(geometry)[, 2],
           site = "Davenport",
           Species = Scientific) |>
    st_drop_geometry() |>
    filter(!(is.na(Species) & is.na(Common_Nam))) |>
    filter(lon>-97)|>
    select(Species, Common_Nam, site, lon, lat) |>
    distinct()
  
  
  gbif_taxon_keys <- raw_Davenport$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  # If there is common name
  gbif_taxon_keys_common <- raw_Davenport |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(Species,Common_Nam, canonicalName) |>
    filter(is.na(canonicalName)) |>
    pull(Common_Nam)
  
  comm2sci_list <- map(
    gbif_taxon_keys_common,
    function(x) {
      Sys.sleep(1) 
      tryCatch(
        comm2sci(x), 
        error = function(e) NULL  
      )
    }
  )
  comm2sci <- imap_dfr(comm2sci_list, function(inner_list, outer_name) {
    if (length(inner_list) == 0) {
      tibble(common_group = outer_name,
             common_name  = NA_character_,
             scientific_name = NA_character_)
    } else {
      map_dfr(names(inner_list), function(comm_name) {
        sci_names <- inner_list[[comm_name]]
        if (length(sci_names) == 0) sci_names <- NA_character_
        tibble(common_group = outer_name,
               common_name = comm_name,
               scientific_name = sci_names)
      })
    }
  }) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = common_name) |>
    distinct()
  
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility
  
  clean_Davenport <- raw_Davenport|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    left_join(comm2sci, by = c("Common_Nam" = "verbatim_name")) |>
    mutate(canonicalName = if_else(is.na(canonicalName), scientific_name, canonicalName)) |>
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Common_Nam, -Species, -scientific_name, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Davenport |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Ithaca
Ithaca <- FALSE
if(Ithaca){
  raw_Ithaca <- st_read(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Ithaca/City_Managed_Trees.shp")) |>
    st_transform(4326) |>
    mutate(lon = st_coordinates(geometry)[, 1],
           lat = st_coordinates(geometry)[, 2],
           site = "Ithaca",
           Species = tolower(SPP_bot),
           Common_Nam = tolower(SPP_com)) |>
    st_drop_geometry() |>
    filter(!(is.na(Species) & is.na(Common_Nam))) |>
    filter(lon>-97)|>
    select(Species, Common_Nam, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_Ithaca$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  # If there is common name
  gbif_taxon_keys_common <- raw_Ithaca |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(Species,Common_Nam, canonicalName) |>
    filter(is.na(canonicalName)) |>
    pull(Common_Nam)
  gbif_taxon_keys_common_swapped <- str_replace(
    gbif_taxon_keys_common,
    "^([^,]+),\\s*(.+)$",   # 逗号前后分组
    "\\2 \\1"              # 交换顺序
  )
  comm2sci_list <- map(
    gbif_taxon_keys_common_swapped,
    function(x) {
      Sys.sleep(1) 
      tryCatch(
        comm2sci(x), 
        error = function(e) NULL  
      )
    }
  )
  comm2sci <- imap_dfr(comm2sci_list, function(inner_list, outer_name) {
    if (length(inner_list) == 0) {
      tibble(common_group = outer_name,
             common_name  = NA_character_,
             scientific_name = NA_character_)
    } else {
      map_dfr(names(inner_list), function(comm_name) {
        sci_names <- inner_list[[comm_name]]
        if (length(sci_names) == 0) sci_names <- NA_character_
        tibble(common_group = outer_name,
               common_name = comm_name,
               scientific_name = sci_names)
      })
    }
  }) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = common_name) |>
    distinct()
  
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility
  
  clean_Ithaca <- raw_Ithaca|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    left_join(comm2sci, by = c("Common_Nam" = "verbatim_name")) |>
    mutate(canonicalName = if_else(is.na(canonicalName), scientific_name, canonicalName)) |>
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Common_Nam, -Species, -scientific_name, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Ithaca |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Riverside
Riverside <- FALSE
if(Riverside){
  raw_Riverside <- st_read(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Riverside/Street_Trees.shp")) |>
    st_transform(4326) |>
    mutate(lon = st_coordinates(geometry)[, 1],
           lat = st_coordinates(geometry)[, 2],
           site = "Riverside",
           Species = BotanicalN,
           Common_Nam = CommonName) |>
    st_drop_geometry() |>
    filter(!(is.na(Species) & is.na(Common_Nam))) |>
    filter(Species != "Vacant site")|>
    select(Species, Common_Nam, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_Riverside$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  # If there is common name
  gbif_taxon_keys_common <- raw_Riverside |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(Species,Common_Nam, canonicalName) |>
    filter(is.na(canonicalName)) |>
    pull(Common_Nam)
  
  comm2sci_list <- map(
    gbif_taxon_keys_common,
    function(x) {
      Sys.sleep(1) 
      tryCatch(
        comm2sci(x), 
        error = function(e) NULL  
      )
    }
  )
  comm2sci <- imap_dfr(comm2sci_list, function(inner_list, outer_name) {
    if (length(inner_list) == 0) {
      tibble(common_group = outer_name,
             common_name  = NA_character_,
             scientific_name = NA_character_)
    } else {
      map_dfr(names(inner_list), function(comm_name) {
        sci_names <- inner_list[[comm_name]]
        if (length(sci_names) == 0) sci_names <- NA_character_
        tibble(common_group = outer_name,
               common_name = comm_name,
               scientific_name = sci_names)
      })
    }
  }) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = common_name) |>
    distinct()
  
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility
  
  clean_Riverside <- raw_Riverside|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    left_join(comm2sci, by = c("Common_Nam" = "verbatim_name")) |>
    mutate(canonicalName = if_else(is.na(canonicalName), scientific_name, canonicalName)) |>
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Common_Nam, -Species, -scientific_name, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Riverside |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Spokane
Spokane <- FALSE
if(Spokane){
  raw_Spokane <- st_read(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Spokane/Tree_Inventory.shp")) |>
    st_transform(4326) |>
    mutate(lon = st_coordinates(geometry)[, 1],
           lat = st_coordinates(geometry)[, 2],
           site = "Spokane",
           Common_Nam = species) |>
    st_drop_geometry() |>
    filter(!(is.na(Common_Nam))) |>
    select(Common_Nam, site, lon, lat) |>
    distinct()
  
  # If there is common name
  gbif_taxon_keys_common_swapped <- str_replace(
    raw_Spokane |> select(Common_Nam) |> distinct() |> pull(),
    "^([^,]+),\\s*(.+)$",   # 逗号前后分组
    "\\2 \\1"              # 交换顺序
  )
  
  comm2sci_list <- map(
    gbif_taxon_keys_common_swapped,
    function(x) {
      Sys.sleep(1) 
      tryCatch(
        comm2sci(x), 
        error = function(e) NULL  
      )
    }
  )
  comm2sci <- imap_dfr(comm2sci_list, function(inner_list, outer_name) {
    if (length(inner_list) == 0) {
      tibble(common_group = outer_name,
             common_name  = NA_character_,
             scientific_name = NA_character_)
    } else {
      map_dfr(names(inner_list), function(comm_name) {
        sci_names <- inner_list[[comm_name]]
        if (length(sci_names) == 0) sci_names <- NA_character_
        tibble(common_group = outer_name,
               common_name = comm_name,
               scientific_name = sci_names)
      })
    }
  }) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = common_name) |>
    distinct()
  
  common_name <-  raw_Spokane |> # get the first column of the data.frame
    mutate(new_common = str_replace(
      Common_Nam,
      "^([^,]+),\\s*(.+)$",   # 逗号前后分组
      "\\2 \\1"              # 交换顺序
    )) |>
    left_join(comm2sci, by = c("new_common" = "verbatim_name")) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = Common_Nam) |>
    distinct()
  
  gbif_taxon_keys <- common_name$scientific_name|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) |> # check for missing values 
    left_join(common_name, by = c("verbatim_name" = "scientific_name")) |>
    select(canonicalName, verbatim_name = verbatim_name.y)
  
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility
  
  clean_Spokane <- raw_Spokane|>
    left_join(gbif_taxon_keys, by = c("Common_Nam" = "verbatim_name")) |>
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Common_Nam, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Spokane |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# NewOrleans
NewOrleans <- FALSE
if(NewOrleans){
  raw_NewOrleans <- read_csv(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_NewOrleans.csv")) |>
    mutate(lon = as.numeric(str_extract(the_geom, "(?<=\\()[^ ]+")),
           lat = as.numeric(str_extract(the_geom, "(?<= )[0-9\\.-]+(?=\\))")),
           site = "NewOrleans",
           Species = str_remove(`Botanical Name`, "\\s*\\([^)]*\\)"),
           Common_Nam = str_remove(`Common Name`, "\\s*\\([^)]*\\)")) |>
    filter(!(is.na(Species) & is.na(Common_Nam))) |>
    select(Species, Common_Nam, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_NewOrleans$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  # If there is common name
  gbif_taxon_keys_common <- raw_NewOrleans |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(Species,Common_Nam, canonicalName) |>
    filter(is.na(canonicalName)) |>
    pull(Common_Nam)
  
  comm2sci_list <- map(
    gbif_taxon_keys_common,
    function(x) {
      Sys.sleep(1) 
      tryCatch(
        comm2sci(x), 
        error = function(e) NULL  
      )
    }
  )
  comm2sci <- imap_dfr(comm2sci_list, function(inner_list, outer_name) {
    if (length(inner_list) == 0) {
      tibble(common_group = outer_name,
             common_name  = NA_character_,
             scientific_name = NA_character_)
    } else {
      map_dfr(names(inner_list), function(comm_name) {
        sci_names <- inner_list[[comm_name]]
        if (length(sci_names) == 0) sci_names <- NA_character_
        tibble(common_group = outer_name,
               common_name = comm_name,
               scientific_name = sci_names)
      })
    }
  }) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = common_name) |>
    distinct()
  
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility
  
  clean_NewOrleans <- raw_NewOrleans|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    left_join(comm2sci, by = c("Common_Nam" = "verbatim_name")) |>
    mutate(canonicalName = if_else(is.na(canonicalName), scientific_name, canonicalName)) |>
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Common_Nam, -Species, -scientific_name, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_NewOrleans |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Berkeley
Berkeley <- FALSE
if(Berkeley){
  raw_Berkeley <- st_read(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Berkeley/Berkeley_trees_2013.shp")) |>
    st_transform(4326) |>
    mutate(lon = st_coordinates(geometry)[, 1],
           lat = st_coordinates(geometry)[, 2],
           site = "Berkeley",
           Species = SPECIES) |>
    st_drop_geometry() |>
    filter(Species != "Planting Site")|>
    filter(!(is.na(Species) & (Common_Nam == "#N/A"))) |>
    select(Species, Common_Nam, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_Berkeley$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  # If there is common name
  gbif_taxon_keys_common <- raw_Berkeley |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(Species,Common_Nam, canonicalName) |>
    filter(is.na(canonicalName)) |>
    pull(Common_Nam)
  
  comm2sci_list <- map(
    gbif_taxon_keys_common,
    function(x) {
      Sys.sleep(1) 
      tryCatch(
        comm2sci(x), 
        error = function(e) NULL  
      )
    }
  )
  comm2sci <- imap_dfr(comm2sci_list, function(inner_list, outer_name) {
    if (length(inner_list) == 0) {
      tibble(common_group = outer_name,
             common_name  = NA_character_,
             scientific_name = NA_character_)
    } else {
      map_dfr(names(inner_list), function(comm_name) {
        sci_names <- inner_list[[comm_name]]
        if (length(sci_names) == 0) sci_names <- NA_character_
        tibble(common_group = outer_name,
               common_name = comm_name,
               scientific_name = sci_names)
      })
    }
  }) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = common_name) |>
    distinct()
  
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility
  
  clean_Berkeley <- raw_Berkeley|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    left_join(comm2sci, by = c("Common_Nam" = "verbatim_name")) |>
    mutate(canonicalName = if_else(is.na(canonicalName), scientific_name, canonicalName)) |>
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Common_Nam, -Species, -scientific_name, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Berkeley |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Tempe
Tempe <- FALSE
if(Tempe){
  raw_Tempe <- st_read(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Tempe/Tree_Inventory.shp")) |>
    st_transform(4326) |>
    mutate(lon = st_coordinates(geometry)[, 1],
           lat = st_coordinates(geometry)[, 2],
           site = "Tempe",
           Species = Species_Na) |>
    st_drop_geometry() |>
    filter(!(is.na(Species))) |>
    select(Species, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_Tempe$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility
  
  clean_Tempe <- raw_Tempe|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Species, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Tempe |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Livingston
Livingston <- FALSE
if(Livingston){
  raw_Livingston <- st_read(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Livingston_MT/LivingstonTreeSurvey2022.shp")) |>
    st_transform(4326) |>
    mutate(lon = st_coordinates(geometry)[, 1],
           lat = st_coordinates(geometry)[, 2],
           site = "Livingston",
           Species = TreeSpecie,
           Common_Nam = CommonName) |>
    st_drop_geometry() |>
    filter(!(is.na(Species) & ! is.na(Common_Nam))) |>
    select(Species, Common_Nam, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_Livingston$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  # If there is common name
  gbif_taxon_keys_common <- raw_Livingston |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(Species,Common_Nam, canonicalName) |>
    filter(is.na(canonicalName)) |>
    pull(Common_Nam)
  
  comm2sci_list <- map(
    gbif_taxon_keys_common,
    function(x) {
      Sys.sleep(1) 
      tryCatch(
        comm2sci(x), 
        error = function(e) NULL  
      )
    }
  )
  comm2sci <- imap_dfr(comm2sci_list, function(inner_list, outer_name) {
    if (length(inner_list) == 0) {
      tibble(common_group = outer_name,
             common_name  = NA_character_,
             scientific_name = NA_character_)
    } else {
      map_dfr(names(inner_list), function(comm_name) {
        sci_names <- inner_list[[comm_name]]
        if (length(sci_names) == 0) sci_names <- NA_character_
        tibble(common_group = outer_name,
               common_name = comm_name,
               scientific_name = sci_names)
      })
    }
  }) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = common_name) |>
    distinct()
  
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility
  
  clean_Livingston <- raw_Livingston|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    left_join(comm2sci, by = c("Common_Nam" = "verbatim_name")) |>
    mutate(canonicalName = if_else(is.na(canonicalName), scientific_name, canonicalName)) |>
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Common_Nam, -Species, -scientific_name, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Livingston |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Bloomington
Bloomington <- FALSE
if(Bloomington){
  raw_Bloomington <- read_csv(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Bloomington.csv")) |>
    mutate(lon = as.numeric(str_extract(Geometry, "(?<=\\()[^ ]+")),
           lat = as.numeric(str_extract(Geometry, "(?<= )[0-9\\.-]+(?=\\))")),
           site = "Bloomington",
           Common_Nam = `Common Species`) |>
    filter(!(is.na(Common_Nam))) |>
    filter(Condition != "Dead") |>
    select(Common_Nam, site, lon, lat) |>
    distinct()
  
  # If there is common name
  gbif_taxon_keys_common_swapped <- str_replace(
    raw_Bloomington |> select(Common_Nam) |> distinct() |> pull(),
    "^([^,]+),\\s*(.+)$",   # 逗号前后分组
    "\\2 \\1"              # 交换顺序
  )
  
  comm2sci_list <- map(
    gbif_taxon_keys_common_swapped,
    function(x) {
      Sys.sleep(1) 
      tryCatch(
        comm2sci(x), 
        error = function(e) NULL  
      )
    }
  )
  comm2sci <- imap_dfr(comm2sci_list, function(inner_list, outer_name) {
    if (length(inner_list) == 0) {
      tibble(common_group = outer_name,
             common_name  = NA_character_,
             scientific_name = NA_character_)
    } else {
      map_dfr(names(inner_list), function(comm_name) {
        sci_names <- inner_list[[comm_name]]
        if (length(sci_names) == 0) sci_names <- NA_character_
        tibble(common_group = outer_name,
               common_name = comm_name,
               scientific_name = sci_names)
      })
    }
  }) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = common_name) |>
    distinct()
  
  common_name <-  raw_Bloomington |> # get the first column of the data.frame
    mutate(new_common = str_replace(
      Common_Nam,
      "^([^,]+),\\s*(.+)$",   # 逗号前后分组
      "\\2 \\1"              # 交换顺序
    )) |>
    left_join(comm2sci, by = c("new_common" = "verbatim_name")) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = Common_Nam) |>
    distinct()
  
  gbif_taxon_keys <- common_name$scientific_name|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) |> # check for missing values 
    left_join(common_name, by = c("verbatim_name" = "scientific_name")) |>
    select(canonicalName, verbatim_name = verbatim_name.y)
  
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility
  
  clean_Bloomington <- raw_Bloomington|>
    left_join(gbif_taxon_keys, by = c("Common_Nam" = "verbatim_name")) |>
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Common_Nam, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Bloomington |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

############################ Add more ############################

# Norfolk
Norfolk <- FALSE
if(Norfolk){
  raw_Norfolk <- read_csv(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Norfolk.csv")) |>
    mutate(lon = Longitude,,
           lat = Latitude,
           site = "Norfolk",
           Species = str_c(Genus, Species, sep = " "),
           Common_Nam = `Common Name`) |>
    filter(!(is.na(Species) & is.na(Common_Nam))) |>
    select(Species, Common_Nam, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_Norfolk$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  # If there is common name
  gbif_taxon_keys_common <- raw_Norfolk |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(Species,Common_Nam, canonicalName) |>
    filter(is.na(canonicalName)) |>
    pull(Common_Nam)
  gbif_taxon_keys_common_swapped <- str_replace(
    gbif_taxon_keys_common,
    "^([^,]+),\\s*(.+)$",   # 逗号前后分组
    "\\2 \\1"              # 交换顺序
  )
  comm2sci_list <- map(
    gbif_taxon_keys_common_swapped,
    function(x) {
      Sys.sleep(1) 
      tryCatch(
        comm2sci(x), 
        error = function(e) NULL  
      )
    }
  )
  comm2sci <- imap_dfr(comm2sci_list, function(inner_list, outer_name) {
    if (length(inner_list) == 0) {
      tibble(common_group = outer_name,
             common_name  = NA_character_,
             scientific_name = NA_character_)
    } else {
      map_dfr(names(inner_list), function(comm_name) {
        sci_names <- inner_list[[comm_name]]
        if (length(sci_names) == 0) sci_names <- NA_character_
        tibble(common_group = outer_name,
               common_name = comm_name,
               scientific_name = sci_names)
      })
    }
  }) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = common_name) |>
    distinct()
  
  common_name <-  raw_Norfolk |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(Species,Common_Nam, canonicalName) |>
    filter(is.na(canonicalName)) |>
    mutate(new_common = str_replace(
      Common_Nam,
      "^([^,]+),\\s*(.+)$",   # 逗号前后分组
      "\\2 \\1"              # 交换顺序
    )) |>
    left_join(comm2sci, by = c("new_common" = "verbatim_name")) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = Common_Nam) |>
    distinct()
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility
  
  clean_Norfolk <- raw_Norfolk|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    left_join(common_name, by = c("Common_Nam" = "verbatim_name")) |>
    mutate(canonicalName = if_else(is.na(canonicalName), scientific_name, canonicalName)) |>
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Common_Nam, -Species, -scientific_name, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Norfolk |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Mesa
Mesa <- FALSE
if(Mesa){
  raw_Mesa <- st_read(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Mesa/geo_export_ec872989-c496-49b7-a2d9-75fc799a03b8.shp")) |>
    st_transform(4326) |>
    mutate(lon = st_coordinates(geometry)[, 1],
           lat = st_coordinates(geometry)[, 2],
           site = "Mesa",
           Species = botanical_,
           Common_Nam = common_nam) |>
    filter(!(is.na(Species) & is.na(Common_Nam))) |>
    filter(Species != "Vacant site") |>
    st_drop_geometry() |>
    select(Species, Common_Nam, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_Mesa$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  # If there is common name
  gbif_taxon_keys_common <- raw_Mesa |> # get the first column of the data.frame
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    distinct(Species,Common_Nam, canonicalName) |>
    filter(is.na(canonicalName)) |>
    pull(Common_Nam)
  comm2sci_list <- map(
    gbif_taxon_keys_common,
    function(x) {
      Sys.sleep(1) 
      tryCatch(
        comm2sci(x), 
        error = function(e) NULL  
      )
    }
  )
  comm2sci <- imap_dfr(comm2sci_list, function(inner_list, outer_name) {
    if (length(inner_list) == 0) {
      tibble(common_group = outer_name,
             common_name  = NA_character_,
             scientific_name = NA_character_)
    } else {
      map_dfr(names(inner_list), function(comm_name) {
        sci_names <- inner_list[[comm_name]]
        if (length(sci_names) == 0) sci_names <- NA_character_
        tibble(common_group = outer_name,
               common_name = comm_name,
               scientific_name = sci_names)
      })
    }
  }) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = common_name) |>
    distinct()
  
  set.seed(123) # for reproducibility
  
  clean_Mesa <- raw_Mesa|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> # join the GBIF taxon keys
    left_join(comm2sci, by = c("Common_Nam" = "verbatim_name")) |>
    mutate(canonicalName = if_else(is.na(canonicalName), scientific_name, canonicalName)) |>
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Common_Nam, -Species, -scientific_name, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Mesa |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Utah
UT <- FALSE
if(UT){
  raw_UT <- st_read(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_UT/UrbanTreeInventory.shp")) |>
    st_transform(4326) |>
    mutate(lon = st_coordinates(geometry)[, 1],
           lat = st_coordinates(geometry)[, 2],
           site = Community) |>
    filter(!(is.na(Species) | (Species == "Unknown"))) |>
    filter(is.na(removed) |  (removed == "no")) |>
    filter(TreeType == "Deciduous") |>
    st_drop_geometry() |>
    select(Common_Nam = Species, site, lon, lat) |>
    filter(site %in% c( "Moab", "St. George", "Nephi", "Murray")) |>
    mutate(Common_Nam = if_else(Common_Nam == "Tree-of-Heaven", "Tree of Heaven", Common_Nam)) |>
    distinct()
  
  # If there is common name
  gbif_taxon_keys_common <- raw_UT |> # get the first column of the data.frame
    distinct(Common_Nam) |>
    pull(Common_Nam)
  gbif_taxon_keys_common_swapped <- str_replace(
    gbif_taxon_keys_common,
    "^([^,]+)-\\s*(.+)$",   # 逗号前后分组
    "\\2 \\1"              # 交换顺序
  )
  comm2sci_list <- map(
    gbif_taxon_keys_common_swapped,
    function(x) {
      Sys.sleep(1) 
      tryCatch(
        comm2sci(x), 
        error = function(e) NULL  
      )
    }
  )
  comm2sci <- imap_dfr(comm2sci_list, function(inner_list, outer_name) {
    if (length(inner_list) == 0) {
      tibble(common_group = outer_name,
             common_name  = NA_character_,
             scientific_name = NA_character_)
    } else {
      map_dfr(names(inner_list), function(comm_name) {
        sci_names <- inner_list[[comm_name]]
        if (length(sci_names) == 0) sci_names <- NA_character_
        tibble(common_group = outer_name,
               common_name = comm_name,
               scientific_name = sci_names)
      })
    }
  }) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = common_name) |>
    distinct()
  
  
  common_name <-  raw_UT |> # get the first column of the data.frame
    mutate(new_common = str_replace(
      Common_Nam,
      "^([^,]+)-\\s*(.+)$",   # 逗号前后分组
      "\\2 \\1"              # 交换顺序
    )) |>
    left_join(comm2sci, by = c("new_common" = "verbatim_name")) |>
    filter(!is.na(scientific_name)) |>
    select(scientific_name, verbatim_name = Common_Nam) |>
    distinct()
  
  gbif_taxon_keys <- common_name$scientific_name|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) |> # check for missing values 
    left_join(common_name, by = c("verbatim_name" = "scientific_name")) |>
    select(canonicalName, verbatim_name = verbatim_name.y)
  
  # Join to TRY data and find the deciduous species
  
  set.seed(123) # for reproducibility
  
  clean_UT <- raw_UT|>
    left_join(gbif_taxon_keys, by = c("Common_Nam" = "verbatim_name")) |>
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Common_Nam, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"),
           site = if_else(site == "St. George", "St_George", site))
  
  clean_UT |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Newburgh
Newburgh <- FALSE
if(Newburgh){
  raw_Newburgh <- read_csv(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Newburgh.csv")) |>
    mutate(lon = Longitude,
           lat = Latitude,
           site = "Newburgh",
           Species = str_extract(Species, "\\(([^)]+)\\)") |> 
             str_replace_all("[()]", "")) |>
    filter(!(is.na(Species) | Species=="stump")) |>
    select(Species, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_Newburgh$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  set.seed(123) # for reproducibility
  
  clean_Newburgh <- raw_Newburgh|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> 
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Species, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Newburgh |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Reading
Reading <- FALSE
if(Reading){
  raw_Reading <- read_csv(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Reading.csv")) |>
    mutate(lon = Longitude,
           lat = Latitude,
           site = "Reading",
           Species = str_extract(Species, "\\(([^)]+)\\)") |> 
             str_replace_all("[()]", "")) |>
    filter(!(is.na(Species) | Species=="stump")) |>
    filter(`Tree Present` != "No") |>
    select(Species, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_Reading$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  set.seed(123) # for reproducibility
  
  clean_Reading <- raw_Reading|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> 
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Species, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Reading |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Erie
Erie <- FALSE
if(Erie){
  raw_Erie <- read_csv(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Erie.csv")) |>
    mutate(lon = Longitude,
           lat = Latitude,
           site = "Erie",
           Species = str_extract(Species, "\\(([^)]+)\\)") |> 
             str_replace_all("[()]", "")) |>
    filter(!(is.na(Species) | Species=="Unknown Tree")) |>
    filter(`Tree Present` != "No") |>
    select(Species, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_Erie$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  set.seed(123) # for reproducibility
  
  clean_Erie <- raw_Erie|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> 
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Species, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Erie |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Lancaster
Lancaster <- FALSE
if(Lancaster){
  raw_Lancaster <- read_csv(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Lancaster.csv")) |>
    mutate(lon = Longitude,
           lat = Latitude,
           site = "Lancaster",
           Species = str_extract(Species, "\\(([^)]+)\\)") |> 
             str_replace_all("[()]", "")) |>
    filter(!(is.na(Species) | Species=="Unknown Tree")) |>
    filter(`Tree Present` != "No") |>
    select(Species, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_Lancaster$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  set.seed(123) # for reproducibility
  
  clean_Lancaster <- raw_Lancaster|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> 
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Species, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Lancaster |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Harrisburg
Harrisburg <- FALSE
if(Harrisburg){
  raw_Harrisburg <- read_csv(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Harrisburg.csv")) |>
    mutate(lon = Longitude,
           lat = Latitude,
           site = "Harrisburg",
           Species = str_extract(Species, "\\(([^)]+)\\)") |> 
             str_replace_all("[()]", "")) |>
    filter(!(is.na(Species) | Species=="Unknown Tree")) |>
    filter(`Tree Present` != "No") |>
    select(Species, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_Harrisburg$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  set.seed(123) # for reproducibility
  
  clean_Harrisburg <- raw_Harrisburg|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> 
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Species, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Harrisburg |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# FortCollins
FortCollins <- FALSE
if(FortCollins){
  raw_FortCollins <- read_csv(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_FortCollins.csv")) |>
    mutate(lon = Longitude,
           lat = Latitude,
           site = "FortCollins",
           Species = str_extract(Species, "\\(([^)]+)\\)") |> 
             str_replace_all("[()]", "")) |>
    filter(!(is.na(Species) | Species=="Vacant Site")) |>
    select(Species, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_FortCollins$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  set.seed(123) # for reproducibility
  
  clean_FortCollins <- raw_FortCollins|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> 
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Species, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_FortCollins |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

# Binghamton
Binghamton <- FALSE
if(Binghamton){
  raw_Binghamton <- read_csv(("data/raw/CONUS/StreetTrees/raw/Tree_Inventory_Binghamton.csv")) |>
    mutate(lon = Longitude,
           lat = Latitude,
           site = "Binghamton",
           Species = str_extract(Species, "\\(([^)]+)\\)") |> 
             str_replace_all("[()]", "")) |>
    filter(!(is.na(Species) | Species=="vacant site small")) |>
    select(Species, site, lon, lat) |>
    distinct()
  
  gbif_taxon_keys <- raw_Binghamton$Species|> # get the first column of the data.frame
    unique() |> # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")|>
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    # filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES")|> # get only accepted names and remove anything that might have matched to a non-species
    mutate(canonicalName = if_else(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES", species, NA_character_)) |>
    select(canonicalName, verbatim_name) # check for missing values
  
  set.seed(123) # for reproducibility
  
  clean_Binghamton <- raw_Binghamton|>
    left_join(gbif_taxon_keys, by = c("Species" = "verbatim_name"))|> 
    drop_na(canonicalName)|> 
    rename(canonicalname = canonicalName) |>
    group_by(canonicalname, site) |>
    filter(n() >= 30) |>
    group_modify(~ {
      if (nrow(.x) > 2000) {
        .x |> slice_sample(n = 2000)
      } else {
        .x                          
      }
    }) |>
    ungroup() |>
    left_join(check, by = c("canonicalname" = "canonicalName")) |>
    filter(str_detect(single_type, "deciduous")) |>
    select(-Species, -types_with_counts, -types, -single_type) |>
    rename(taxa = canonicalname) |>
    mutate(genus = str_extract(taxa, "^[^ ]+"))
  
  clean_Binghamton |>
    group_split(site) |>  
    purrr::walk(function(df) {
      site_name <- gsub("\\.", "_", unique(df$site))
      df <- df |>
        mutate(id = row_number())
      write_rds(df, here(paste0(save_dir, site_name, ".rds")), compress = "gz")
    })
}

