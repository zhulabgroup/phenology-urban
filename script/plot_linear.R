plot_function <- function(dataframe, taxa, x_ind, y_ind, label_x, label_y) {
  plot <- ggplot(dataframe, aes(x = x_ind, y = y_ind)) +    # change
    geom_hex(bins = 60, color = NA) +  
    scale_fill_gradient(low = "lightcoral", high = "darkred") +
    theme_classic() +
    labs(x = label_x, y = label_y, fill = "Point Count") +               # change
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
    scale_x_continuous(breaks = pretty(x_ind)) +       # change
    scale_y_continuous(breaks = pretty(y_ind)) +  # change
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 18)
    ) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", formula = y ~ x) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~~")),
                 formula = y ~ x,
                 parse = TRUE,
                 label.x = "right", label.y = 0.9,
                 size = 4, color = "blue") +
    ggtitle(str_c("2018, ", taxa)) +
    theme(plot.title = element_text(size = 22))
  
#  ggsave(str_c("result/bld_index/",label_x, "_", label_y,".jpg"), plot, width = 8, height = 6, units = "in")

    return(plot)
  }