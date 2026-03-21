plot_heatmap <- function(coef_summary_all, level = "species", city_order_bio = "bio1", phe_type = "SOS") {
  df <- coef_summary_all %>%
    filter(phe_type == !!phe_type) %>%
    group_by(city, .data[[level]]) %>%
    summarise(
      Estimate = mean(Estimate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(city_feature, by = "city") %>%
    mutate(
      city = case_when(
        city == "DT" ~ "Detroit",
        city == "TP" ~ "Tampa",
        TRUE ~ city
      )
    )
    
  if (level == "species") {
    class_tree <- phylo_tree
  } else if (level == "genus") {
    class_tree <- tree_genus
  } else if (level == "family") {
    class_tree <- tree_family
    }
  
  df_filtered <- df %>%
    filter(.data[[level]] %in% class_tree$tip.label)
  class_used   <- unique(df_filtered[[level]])
  tree_sub  <- drop.tip(class_tree,
                        setdiff(class_tree$tip.label, class_used))
  class_order <- tree_sub$tip.label 
  city_order <- df_filtered %>%
    distinct(city, .data[[city_order_bio]]) %>%
    arrange(.data[[city_order_bio]]) %>%
    pull(city)
  
  coef_df_plot <- df_filtered |>
    complete(city, .data[[level]]) |>
    mutate(
      city    = factor(city, levels = city_order),
      !!level := factor(.data[[level]], levels = rev(class_order))
    )
  level_title <- case_when(
    level == "species" ~ "Species",
    level == "genus"   ~ "Genus",
    level == "family"  ~ "Family"
  )
  
  p_heat <- ggplot(coef_df_plot, aes(x = .data[[level]], y = city, fill = Estimate)) +
    geom_tile() +
    scale_fill_gradient2(low = "#FCA636FF", mid = "white", high = "#0D0887FF", midpoint = 0,
                         name = "Estimate", na.value = "grey90") +
    theme_minimal(base_size = 10) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      panel.grid = element_blank(),
      axis.title.x = element_text(size = 18),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 18),
      legend.key.size = unit(20, "pt")
    ) +
    labs(fill = expression(beta[temp]),
         x = level_title)
  
  city_strip <- ggplot(coef_df_plot |>
                         group_by(city) |>
                         summarise(Estimate = mean(Estimate, na.rm = TRUE)),
                       aes(x = 1, y = city, fill = Estimate)) +
    geom_tile() +
    scale_fill_gradient2(low = "#FCA636FF", mid = "white", high = "#0D0887FF", 
                         midpoint = 0,  limits = c(range(coef_df_plot$Estimate, na.rm = TRUE)[[1]], 
                                                   range(coef_df_plot$Estimate, na.rm = TRUE)[[2]]), 
                         name = "Estimate", na.value = "grey90") +
    theme_minimal(base_size = 10) +
    theme(
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      axis.title.y = element_text(size = 18),
      legend.position = "none"
    ) +
    labs(y = "City")
  
  tree_plot <- ggtree(tree_sub) +
    layout_dendrogram()
  
  species_strip <- ggplot(coef_df_plot |>
                            group_by(.data[[level]]) |>
                            summarise(Estimate = mean(Estimate, na.rm = TRUE)),
                          aes(x = .data[[level]], y = 1, fill = Estimate)) +
    geom_tile() +
    scale_fill_gradient2(low = "#FCA636FF", mid = "white", high = "#0D0887FF", 
                         midpoint = 0,  limits = c(range(coef_df_plot$Estimate, na.rm = TRUE)[[1]], 
                                                   range(coef_df_plot$Estimate, na.rm = TRUE)[[2]]), 
                         name = "Estimate", na.value = "grey90") +
    theme_minimal(base_size = 10) +
    theme(
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none"
    )
  
  p_left <- plot_spacer() / plot_spacer() / city_strip +
    plot_layout(heights = c(2, 1, 20))
  p_right <- tree_plot / species_strip / p_heat +
    plot_layout(heights = c(2, 1, 20))
  p_combined <- (p_left | p_right) +
    plot_layout(widths = c(1, 20)) +
    plot_annotation(
      subtitle = paste0(phe_type)
    )
  
  return(p_combined)
}








plot_heatmap_sig <- function(coef_summary_all, level = "species", city_order_bio = "bio1", phe_type = "SOS") {
  df <- coef_summary_all %>%
    filter(phe_type == !!phe_type) %>%
    group_by(city, .data[[level]], SigAlpha, Color) %>%
    summarise(
      Estimate = mean(Estimate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(city_feature, by = "city") %>%
    mutate(
      city = case_when(
        city == "DT" ~ "Detroit",
        city == "TP" ~ "Tampa",
        TRUE ~ city
      )
    )
  
  if (level == "species") {
    class_tree <- phylo_tree
  } else if (level == "genus") {
    class_tree <- tree_genus
  } else if (level == "family") {
    class_tree <- tree_family
  }
  
  df_filtered <- df %>%
    filter(.data[[level]] %in% class_tree$tip.label)
  class_used   <- unique(df_filtered[[level]])
  tree_sub  <- drop.tip(class_tree,
                        setdiff(class_tree$tip.label, class_used))
  class_order <- tree_sub$tip.label 
  
  city_order <- df_filtered %>%
    distinct(city, .data[[city_order_bio]]) %>%
    arrange(desc(.data[[city_order_bio]])) %>%
    pull(city)
  
  coef_df_plot <- df_filtered |>
    complete(city, .data[[level]]) |>
    mutate(
      city    = factor(city, levels = city_order),
      !!level := factor(.data[[level]], levels = rev(class_order))
    ) |>
    mutate(
      Color = case_when(
        is.na(Color) ~ "grey95",
        TRUE ~ Color
      ),
      SigAlpha = case_when(
        is.na(SigAlpha) ~ 1,
        TRUE ~ SigAlpha
      ),
      city    = factor(city, levels = rev(city_order)),
      !!level := factor(.data[[level]], levels = rev(class_order))
    )
  level_title <- case_when(
    level == "species" ~ "Species",
    level == "genus"   ~ "Genus",
    level == "family"  ~ "Family"
  )
  
  p_heat <- ggplot(coef_df_plot, aes(x = .data[[level]], y = city, fill = Color, alpha = SigAlpha)) +
    geom_tile() +
    scale_fill_identity(na.value = "grey95") +
    scale_alpha_identity(na.value = "grey95") +
    # scale_fill_gradient2(low = "#FCA636FF", mid = "white", high = "#0D0887FF", midpoint = 0,
    #                      name = "Estimate", na.value = "grey95") +
    theme_minimal(base_size = 10) +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      # axis.title = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      legend.position = "none",
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 18)
    ) # +
    # labs(y = "City", x = level_title)
  
  tree_plot <- ggtree(tree_sub) +
    layout_dendrogram()
  
  p_right <- tree_plot / p_heat +
    plot_layout(heights = c(2, 30))
  p_combined <- (p_right)
  
  return(p_combined)
}

