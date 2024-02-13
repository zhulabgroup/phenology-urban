make_group_display <- function(point, trend, temp){
  
  #making groups
  if (temp == "lst_2018") {
    bins <- seq(floor(min(point[["lst_2018"]])),
                ceiling(max(point[["lst_2018"]])), by = 2)
  } else {
    bins <- seq(floor(min(point[[temp]])),
                ceiling(max(point[[temp]])), by = 0.1)
  }
 

  point <- point %>% mutate(bin = cut(point[[temp]], breaks = bins))
  
  point <- point %>%
    mutate(group = ifelse(alan_2018 <= quantile(alan_2018, 0.25, na.rm = TRUE), "Low ALAN",
                          ifelse(alan_2018 >= quantile(alan_2018, 0.75, na.rm = TRUE), "High ALAN", "Medium ALAN"))) # different temp
  
  point$group <- factor(point$group, levels = c("High ALAN", "Medium ALAN", "Low ALAN"))
  
  
  #wilcox_result <- betula_200_bld_forstats_alan %>% group_by(bin) %>%
  #  filter("High ALAN" %in% group & "Low ALAN" %in% group) %>%
  #  summarize(p_value = ifelse(n_distinct(group) == 2, wilcox.test(up_05_2018 ~ group)$p.value, NA))
  #significance_results <- data.frame(bin = wilcox_result$bin, p_value = wilcox_result$p_value)

  p <- ggplot(point, aes(x = bin, y = .data[[trend]], fill = factor(group))) +
    geom_boxplot(position = "dodge", width = 0.7, alpha = 0.5) +
    labs(x = temp, y = trend, fill = "Group") +
    scale_fill_manual(values = c("#FF9999", "#FFCC99", "#66FF66")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(str_c("2018, ", taxa))
  #+
  #  geom_text(data = significance_results, aes(x = bin, y = 1.05 * max(betula_200_bld_forstats_alan$up_05_2018), label = ifelse(p_value < 0.05, "*", "")), size = 4, vjust = -0.5)

  return(p)
  }