plot_diff_scatter <- function(data, biom, compare) {
  
  filter_4_source <-  data %>%
    group_by(species) %>%
    count(species) %>%
    ungroup() %>%
    filter(n == 4) %>%                          ####
    right_join(data, by = "species") %>%
    na.omit() %>%
    subset(source %in% compare)
  
  
  
  wider_data <- filter_4_source %>%
    select(species, source, biom) %>%
    pivot_wider(names_from = source, 
                values_from = biom)
  
  diff_data <- wider_data %>%
    mutate(
      diff = !!sym(compare[1]) - !!sym(compare[2])  # Use !!sym() to evaluate column names
    )
  
  model <- lm(as.formula(paste("diff", "~", compare[1])), diff_data)
  coef <- coef(model)
  equation_text <- sprintf("y = %.2fx + %.2f", coef[2], coef[1])
  r_squared <- summary(model)$r.squared
  p_value <- summary(model)$coefficients[2, 4]
  r2_text <- paste0("R² = ", round(r_squared, 3))
  p_text <- paste0("p-value = ", format(p_value, digits = 3, scientific = TRUE))
  
  
  plot <- ggplot(diff_data, aes(x = !!sym(compare[1]), y = diff)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE, level = 0.95,  fill = "gray", alpha = 0.5) +
    labs(
      title = biom, 
      subtitle = paste(
        equation_text, "; ", r2_text, "; ", p_text
      ),
      x = "Natural", 
      y = "Natural - Urban"
    ) +
    theme_classic()
  
  return(plot)
}



