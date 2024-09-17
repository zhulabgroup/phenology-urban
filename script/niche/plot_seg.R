plot_seg <- function(data, biom){
  dodge <- position_dodge(width = 0.3)
  
  biom_columns <- data %>%
    select(starts_with(biom)) %>%
    cbind(data %>% select("source", "species"))
  
  ymin_col <- sym(paste0(biom, "min"))
  ymax_col <- sym(paste0(biom, "max"))
  median_col <- sym(paste0(biom, "median"))
  mean_col <- sym(paste0(biom, "mean"))
  p5_col <- sym(paste0(biom, "quantile_5"))
  p95_col <- sym(paste0(biom, "quantile_95"))
  
  
  sorted_species <- biom_columns %>%
    filter(source == "natural") %>%
    arrange(desc(!!median_col)) %>%
    pull(species)
  
  biom_columns <- biom_columns %>%
    mutate(species_order = factor(species, levels = sorted_species)) %>%
    arrange(species_order, source)
  
  
  plot <- ggplot(biom_columns, aes(x = species_order, color = source, group = source)) +
    geom_linerange(aes(ymin = !!ymin_col, ymax = !!ymax_col), linewidth = 1.5, alpha = 0.2) +  # 绘制 min/max 范围
    geom_point(aes(y = !!median_col, shape = "median", color = source), size = 1.5) +  # 中值点
    geom_point(aes(y = !!mean_col, shape = "mean", color = source), size = 1.5) +  # 平均值点
    geom_point(aes(y = !!p5_col, shape = "p5", color = source), size = 1.5) +  # p5 点
    geom_point(aes(y = !!p95_col, shape = "p95", color = source), size = 1.5) +  # p95 点
    labs(title = "Comparison of Environmental Variables Across Sources",
         x = "Environmental Variable and Source", y = "Value") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # x轴标签倾斜
    scale_color_manual(name = "Source",
                       
                       values = c("natural" = "blue", "urban" = "red2"),
                       labels = c("Natural", "Urban"),) +
    # ,"guti" = "purple", "gbif - little" = "orange"
    
    scale_shape_manual(name = "Statistics",
                       values = c("mean" = 4, "median" = 5, "p5" = 124, "p95" = 124),
                       labels = c("Mean", "Median", "P5", "P95")) +  # 添加 mean/median/p1/p99 的图例
    
    guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +  # 确保图例显示，并调整顺序
    coord_flip()  # 翻转坐标轴
  
  return(plot)
}