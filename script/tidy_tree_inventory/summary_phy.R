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
  ) +
  theme_void() +
  theme(  legend.position = "bottom",
          legend.box = "horizontal",
          legend.key.width = unit(3, "cm"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 22)
  )

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
  biom_columns <- occ_species_exp %>%
    group_by(species) %>%
    summarise(
      mean_col    = mean(.data[[climate_var]], na.rm = TRUE),
      median_col  = median(.data[[climate_var]], na.rm = TRUE),
      p5_col      = quantile(.data[[climate_var]], 0.05, na.rm = TRUE),
      p95_col     = quantile(.data[[climate_var]], 0.95, na.rm = TRUE),
      p1_col      = quantile(.data[[climate_var]], 0.01, na.rm = TRUE),
      p99_col     = quantile(.data[[climate_var]], 0.99, na.rm = TRUE),
      .groups     = "drop"
    )
  
  ggplot() +
    geom_linerange(data = biom_columns, aes(x = species,ymin = p1_col, ymax = p99_col), color = "grey30", linewidth = 5, alpha = 0.2) + 
    geom_point(data = biom_columns, aes(x = species,y = median_col, shape = "median"), size = 4) +  
    geom_point(data = biom_columns, aes(x = species,y = mean_col, shape = "mean"), color = "grey30", size = 4, alpha =0.5) + 
    geom_point(data = biom_columns, aes(x = species,y = p5_col, shape = "p5"), color = "grey30", size = 4) +  
    geom_point(data = biom_columns, aes(x = species,y = p95_col, shape = "p95"), color = "grey30", size = 4) + 
    labs(y = var_label) +
    scale_shape_manual(name = "Climatic niche metrics",
                       values = c("mean" = 4, "median" = 20, "p5" = 124, "p95" = 124),
                       labels = c("Mean", "Median", "0.05 quantile", "0.95 quantile")) +  
    guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) + 
    coord_flip()+
    # percentile line (95% for bio1, 5% for bio12)
    geom_hline(yintercept = q_val, linewidth = 1.1, color = "grey30", ) +
    # city lines
    geom_hline(data = city_df, aes(yintercept = x), linetype = "dashed", color = "dodgerblue4", alpha = 0.7, linewidth = 1.1) +
    theme_classic(base_size = 13) +
    theme(legend.position = "bottom",
      axis.text.x = element_text( hjust = 1),
          axis.ticks.y = element_line(linewidth = 0),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())
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

p_safety_margin_temp/p_safety_margin_prcp +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
