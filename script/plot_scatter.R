plot_function <- function(dataframe, family, x_ind, y_ind, label_x, label_y) {
  plot <- ggplot(dataframe, aes(x = x_ind, y = y_ind, color = Species)) +    # city
    geom_point(alpha = 0.2,  stroke = NA) + 
    theme_classic() +
    labs(x = label_x, y = label_y) +  #city             
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
    scale_x_continuous(breaks = pretty(x_ind)) +
    scale_y_continuous(breaks = pretty(y_ind)) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 18)
    ) +
    geom_smooth(method = "lm", se = FALSE, aes(group = species), formula = y ~ x) +  #city
    # geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", formula = y ~ x) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~~")),
                 formula = y ~ x,
                 parse = TRUE,
                 label.x = "right", label.y = "top",  
                 size = 2.5,
                 geom = "text",   
                 hjust = 1,       
                 vjust = -0.5     
    ) +
    
    ggtitle(str_c("2018, ", family)) +
    theme(plot.title = element_text(size = 22))
  
 # ggsave(str_c("result/alan/",label_x, "_", label_y, "_",taxa,"_allcity2018.jpg"), plot, width = 8, height = 6, units = "in")
  
  return(plot)
}