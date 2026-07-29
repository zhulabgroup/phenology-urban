pacman::p_load("V.PhyloMaker2", "ggtree", "ggpubr", "ggtreeExtra")

load(here("release/input/tree_list/summary/tree_species.RData"))

plot_variance_ring <- function(phe_type) {
  plot_df <- read_rds(here(paste0("release/intermediate/variance_decompose/variance_partition_", phe_type, ".rds"))) %>%
    select(species, prop_city2, prop_within2, phe_type, sample_city) %>%
    filter(is.finite(prop_city2), is.finite(prop_within2)) %>%
    mutate(
      prop_city2 = round(prop_city2 * 100, 2),
      prop_within2 = round(prop_within2 * 100, 2),
      group = if_else(prop_city2 < abs(prop_within2),
                      "Within-city > Across-city",
                      "Within-city < Across-city")
    ) %>%
    mutate(
      species = factor(species, levels = unique(species))
    ) |>
    pivot_longer(
      cols = c(prop_within2, prop_city2),
      names_to = "component",
      values_to = "prop"
    ) %>%
    mutate(
      component = recode(component,
                         prop_within2 = "Within-city",
                         prop_city2   = "Across-city")
    )
  
  
  D <- plot_df %>%
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
  
  all_tree <- ggtree(tr2, layout = "fan", open.angle = 10)+
    geom_tiplab(offset = 55, size =0, linesize = 0.3)
  all_ring <- as.data.frame(plot_df) |>
    mutate(label = species)
  
  p_ring <- all_tree +
    geom_fruit(
      data = all_ring,
      geom = geom_col,
      mapping = aes(
        y = label,
        x = prop,
        fill = component
      ),
      orientation = "y",
      width = 1,
      offset = 0.05,
      pwidth = 0.5,
      axis.params = list(
        axis = "x",
        text.size = 4.5,
        hjust = 1,
        vjust = 1.5,
        nbreak = 3
      ),
      grid.params = list()
    ) +
    scale_fill_manual(
      values = c("Across-city" = "grey80", "Within-city" = "grey40"),
      name = "Variance proportion (%)",
      labels = c("Between-city", "Within-city")
    )
  
  return(p_ring)
}
