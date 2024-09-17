plot_seg <- function(data, biom, compare){
  
  filter_4_source <-  data %>%
    group_by(species) %>%
    count(species) %>%
    ungroup() %>%
    filter(n == 4) %>%                          ####
    right_join(data, by = "species") %>%
    na.omit() %>%
    subset(source %in% compare)
  
  
  dodge <- position_dodge(width = 0.3)
  
  biom_columns <- filter_4_source %>%
    select(starts_with(biom)) %>%
    cbind(filter_4_source %>% select("source", "species"))
  
  min_col <- sym(paste0(biom, "min"))
  max_col <- sym(paste0(biom, "max"))
  median_col <- sym(paste0(biom, "median"))
  mean_col <- sym(paste0(biom, "mean"))
  p5_col <- sym(paste0(biom, "quantile_5"))
  p95_col <- sym(paste0(biom, "quantile_95"))
  
  
  sorted_species <- biom_columns %>%
    filter(source == compare[1]) %>%
    arrange(desc(!!median_col)) %>%
    pull(species)
  
  biom_columns <- biom_columns %>%
    mutate(species_order = factor(species, levels = sorted_species)) %>%
    arrange(species_order, source)
  
  
  plot <- ggplot(biom_columns, aes(x = species_order, color = source, group = source)) +
    geom_linerange(aes(ymin = !!p5_col, ymax = !!p95_col), linewidth = 1.5, alpha = 0.2) +  # 绘制 min/max 范围
    geom_point(aes(y = !!median_col, shape = "median", color = source), size = 1.2) +  # 中值点
    geom_point(aes(y = !!mean_col, shape = "mean", color = source), size = 1.2, alpha =0.5) +  # 平均值点
    geom_point(aes(y = !!p5_col, shape = "p5", color = source), size = 1.5) +  # p5 点
    geom_point(aes(y = !!p95_col, shape = "p95", color = source), size = 1.5) +  # p95 点
    labs(title = paste0("Comparison of ", biom, "value between natural and urban niches"),
         x = "Species", y = "Value") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 5)) +  
    scale_color_manual(name = "Source",
                       
                       values = c("blue", "red2"),
                       labels = c(compare[1], compare[2])) +  
    # ,"guti" = "purple", "gbif - little" = "orange"
    
    scale_shape_manual(name = "Statistics",
                       values = c("mean" = 4, "median" = 20, "p5" = 124, "p95" = 124),
                       labels = c("Mean", "Median", "P5", "P95")) +  
    
    guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) + 
    coord_flip()
  
  return(plot)
}