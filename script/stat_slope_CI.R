analyze_stat_slopeCI <- function(data, temp_col, doy_col, title) {
  
  library(ggplot2)
  
  plot_sca_c <- ggplot(data, aes_string(x = temp_col, y = doy_col, color = "group")) +    
    scale_color_manual(values = c("#EE9572", "#00B2EE")) +
    geom_point(alpha = 0, stroke = NA) + 
    theme_classic() +
    labs(x = "Mean summer temperature (°C)", y = "Colored leaves (Day of year)") +             
    theme(axis.text = element_text(size = 20), 
          axis.title = element_text(size = 22),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
          axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(r = 10)),
          legend.position = "none",
          plot.title = element_text(size = 22)) +
    scale_x_continuous(breaks = pretty(data[[temp_col]])) +
    scale_y_continuous(breaks = pretty(data[[doy_col]])) +
    geom_smooth(method = "lm", se = TRUE, aes(group = group), formula = y ~ x) +  
    ggtitle(title)
  
  model_group1 <- lm(as.formula(paste(doy_col, "~", temp_col)), data = data, subset = group == "No ALAN")
  model_group2 <- lm(as.formula(paste(doy_col, "~", temp_col)), data = data, subset = group == "ALAN")
  
  coef_group1 <- coef(model_group1)
  coef_group2 <- coef(model_group2)
  
  equation_text1 <- sprintf("y = %.2fx + %.2f", coef_group1[2], coef_group1[1])
  equation_text2 <- sprintf("y = %.2fx + %.2f", coef_group2[2], coef_group2[1])
  
  confint_group1 <- confint(model_group1, temp_col, level = 0.95)
  confint_group2 <- confint(model_group2, temp_col, level = 0.95)
  
  x_pos <- max(data[[temp_col]]) - diff(range(data[[temp_col]])) * 0.2
  y_pos_group1 <- predict(model_group1, newdata = setNames(data.frame(x_pos), temp_col)) * 1.07
  y_pos_group2 <- predict(model_group2, newdata = setNames(data.frame(x_pos), temp_col)) *0.93
  
  ci_text1 <- sprintf("CI: [%.2f, %.2f]", confint_group1[1,1], confint_group1[1,2])
  ci_text2 <- sprintf("CI: [%.2f, %.2f]", confint_group2[1,1], confint_group2[1,2])
  
  plot_sca_c <- plot_sca_c + 
    annotate("text", x = x_pos, y = y_pos_group1, label = ci_text1, size = 6, color = "#00B2EE") +
    annotate("text", x = x_pos, y = y_pos_group2, label = ci_text2, size = 6,  color = "#EE9572") +
    annotate("text", x = x_pos, y = y_pos_group1 +10, label = equation_text1, size = 6, color = "#00B2EE") +
    annotate("text", x = x_pos, y = y_pos_group2 -10, label = equation_text2, size = 6, color = "#EE9572")
  
  
  return(plot_sca_c)
}

# 