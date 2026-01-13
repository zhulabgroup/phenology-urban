library(zhulabtools)
# Check and  load required packages
load_packages(c("sf","viridis", "terra","patchwork","tidyverse", "ggnewscale", "plotbiomes","patchwork", "terra",  "taxize","metafor", "V.PhyloMaker2", "ggtree", "ggpubr", "ggtreeExtra", "brms", "ggridges", "tidybayes", "readr", "cowplot", "rgbif", "zhulabgbif", "cli", "tidyverse", "broom"))

inventory_dir <- "data/raw/CONUS/StreetTrees/cleaned_ps_input/"
file_list <- list.files(inventory_dir, pattern = "\\.rds$", full.names = TRUE)

all_trees <- map_dfr(file_list, function(file) {
  city <- read_rds(file) |>
    select(taxa, site) |>
    rename(species = taxa)
})  |>
  rename(city = site) |>
  mutate(city = ifelse(city == "St_George", "StGeorge", city)) |>
  mutate(city = ifelse(city == "St.Louis", "StLouis", city))      

################### Phylogeny tree ####################

summary_num <- all_trees |>
  group_by(species) |>
  summarise(num = n()) |>
  ungroup() |>
  arrange(desc(num))

load("~/phenology-urban/data/proc/model_output/final_models_clean_outlier/inla/spatial_spde/tree_species.RData")

D <- summary_num %>%
  mutate(species = gsub("_", " ", species)) %>%
  filter(species %in% phylo_tree$tip.label)

sp_data <- sort(unique(na.omit(D$species)))
sp_tree <- sort(phylo_tree$tip.label)
keep_sp <- intersect(sp_data, sp_tree)

D2 <- D %>%
  filter(species %in% keep_sp) %>%
  droplevels() %>%
  mutate(
    species = factor(species, levels = keep_sp)
  )

tr2 <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, keep_sp))

all_tree <- ggtree(tr2, layout = "circular")+
  geom_tiplab(offset = 55, size =2)
all_ring <- as.data.frame(summary_num) |>
  mutate(label = species)

p_ring <- all_tree +
  geom_fruit(
    data  = all_ring,
    geom  = geom_col,
    mapping = aes(
      y = label,
      x = 1,  
      fill = num 
    ),
    orientation = "y",
    width  = 1,
    offset = 0.05,
    pwidth = 0.1,
    axis.params = list(axis="none", text=FALSE, line=FALSE, ticks=FALSE)
  ) +
  scale_fill_gradient(
    name = "Number of trees",
    low  = "white",
    high = "darkblue" 
  )
  theme_void()

p_ring


################### Safety margin for 1 species ####################

# chelsa_dir <- "~/lab-data/datasets/climate/CHELSA/bioclim"
# tree_dir <- "~/phenology-urban/data/raw/CONUS/StreetTrees/cleaned_ps_input"
# tree_chelsa_dir <- "~/phenology-urban/data/raw/CONUS/StreetTrees/tree_bioclimate"
# chelsa_ras <- load_chelsa(chelsa_vars = c("bio1", "bio12"), chelsa_dir = chelsa_dir)
# process_gbif_chelsa(chelsa_ras, tree_dir, tree_chelsa_dir, overwrite = FALSE)
# 
# tree_chelsa_files <- list.files(tree_chelsa_dir, pattern = "\\.rds$", full.names = TRUE)
# tree_chelsa_summary_all <- map_dfr(tree_chelsa_files, function(tree_chelsa) {
#   read_rds(tree_chelsa) |>
#     group_by(site) |>
#     summarise(
#       bio1_mean   = mean(bio1, na.rm = TRUE),
#       bio1_median = median(bio1, na.rm = TRUE),
#       bio1_p95    = quantile(bio1, 0.95, na.rm = TRUE),
#       bio1_p5     = quantile(bio1, 0.05, na.rm = TRUE),
#       bio12_mean   = mean(bio12, na.rm = TRUE),
#       bio12_median = median(bio12, na.rm = TRUE),
#       bio12_p95    = quantile(bio12, 0.95, na.rm = TRUE),
#       bio12_p5     = quantile(bio12, 0.05, na.rm = TRUE),
#       n_records   = n(),
#       .groups = "drop"
#     )
# }) |>
#   rename(bio1 = bio1_median,
#          bio12 = bio12_median)
# write_rds(tree_chelsa_summary_all, "~/phenology-urban/data/raw/CONUS/StreetTrees/cleaned_ps_input/summary/city_bioclimate.rds", compress = "gz")

city_feature <- readRDS("~/phenology-urban/data/raw/CONUS/StreetTrees/cleaned_ps_input/summary/city_feature.rds") |>
  select(site) |>
  left_join(read_rds("~/phenology-urban/data/raw/CONUS/StreetTrees/cleaned_ps_input/summary/city_bioclimate.rds") |>
              select(site, bio1, bio12), by = "site") |>
  rename(city = site) |>
  mutate(city = ifelse(city == "St_George", "StGeorge", city)) |>
  mutate(city = ifelse(city == "St.Louis", "StLouis", city)) |>              
  st_drop_geometry()

species_exp <- "Acer rubrum"
usage_key_exp <- species_exp |> 
  name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta") |>
  select(usageKey) |>
  pull()

### Read in occurrence data and join with city features
occ_species_exp <- read_rds(paste0("~/phenology-urban/data/proc/tree_origin/chelsa/specieskey_", usage_key_exp, ".rds"))
city_species_exp <- all_trees |>
  filter(species == species_exp) |>
  select(city) |>
  distinct() |>
  left_join(city_feature, by = "city")


plot_safety_margin <- function(occ_species_exp, city_species_exp, city_list, climate_var, var_label, species_exp) {
  # choose which tail to mark
  q_prob <- if (climate_var == "bio12") 0.05 else 0.95
  q_lab  <- paste0(q_prob*100, "%")
  
  # extract vectors
  occ_vals <- occ_species_exp[[climate_var]]
  city_df  <- city_species_exp |>
    select(city, all_of(climate_var)) |>
    filter(city %in% city_list) |>
    rename(x = all_of(climate_var))
  q_val <- quantile(occ_vals, probs = q_prob, na.rm = TRUE)
  if (is.null(var_label)) var_label <- climate_var
  
  ggplot() +
    # density for occurrence background
    geom_density(
      data = tibble(x = occ_vals) %>% filter(is.finite(x)),
      aes(x = x),
      adjust = 2,
      linewidth = 1.0,
      fill = "grey85",
      color = "grey25",
      alpha = 0.8
    ) +
    # percentile line (95% for bio1, 5% for bio12)
    geom_vline(
      xintercept = q_val,
      linetype = "dashed",
      linewidth = 1.1,
      color = "black"
    ) +
    # city lines
    geom_vline(
      data = city_df,
      aes(xintercept = x),
      color = "dodgerblue4",
      alpha = 0.7,
      linewidth = 0.7
    ) +
    labs(x = var_label, y = "Density") +
    theme_classic(base_size = 13)
}

### Plot density and safety margin
p_safety_margin_temp <- plot_safety_margin(
  occ_species_exp = occ_species_exp,
  city_species_exp = city_species_exp,
  city_list = c("Denver", "TP"),
  climate_var = "bio1",
  var_label = "Mean annual temperature (°C)",
  species_exp = species_exp
)

p_safety_margin_prcp <- plot_safety_margin(
  occ_species_exp = occ_species_exp,
  city_species_exp = city_species_exp,
  city_list = c("Denver", "TP"),
  climate_var = "bio12",
  var_label = "Mean annual precipitation (mm)",
  species_exp = species_exp
)

p_safety_margin_temp/p_safety_margin_prcp 
